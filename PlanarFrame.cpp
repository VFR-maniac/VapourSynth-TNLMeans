/*
**   My PlanarFrame class
**
**   Copyright (C) 2005-2006 Kevin Stone
**
**   This program is free software; you can redistribute it and/or modify
**   it under the terms of the GNU General Public License as published by
**   the Free Software Foundation; either version 2 of the License, or
**   (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**   but WITHOUT ANY WARRANTY; without even the implied warranty of
**   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**   GNU General Public License for more details.
**
**   You should have received a copy of the GNU General Public License
**   along with this program; if not, write to the Free Software
**   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "VapourSynth.h"
#include "AlignedMemory.h"
#include "PlanarFrame.h"

#include <cstring>

inline void PlanarFrame::BitBlt
(
    void       *dstp,
    int         dst_stride,
    const void *srcp,
    int         src_stride,
    int         row_size,
    int         height
)
{
    if( height <= 0 )
        return;
    if( src_stride == dst_stride
     && src_stride == row_size )
        std::memcpy( dstp, srcp, row_size * height );
    else
    {
        int i;
        unsigned char *srcp8 = static_cast<unsigned char *>(const_cast<void *>(srcp));
        unsigned char *dstp8 = static_cast<unsigned char *>                   (dstp);
        for( i = 0; i < height; ++i )
        {
            std::memcpy( dstp8, srcp8, row_size );
            srcp8 += src_stride;
            dstp8 += dst_stride;
        }
    }
}

void PlanarFrame::copyFrom( const VSFrameRef *frame, const VSAPI *vsapi )
{
    if( !y || !u || !v ) return;
    if( vsapi->getFrameFormat( frame )->id == pfYUV420P8 )
    {
        BitBlt( y, ypitch,
                vsapi->getReadPtr   ( frame, 0 ), vsapi->getStride     ( frame, 0 ),
                vsapi->getFrameWidth( frame, 0 ), vsapi->getFrameHeight( frame, 0 ) );
        BitBlt( u, uvpitch,
                vsapi->getReadPtr   ( frame, 1 ), vsapi->getStride     ( frame, 1 ),
                vsapi->getFrameWidth( frame, 1 ), vsapi->getFrameHeight( frame, 1 ) );
        BitBlt( v, uvpitch,
                vsapi->getReadPtr   ( frame, 2 ), vsapi->getStride     ( frame, 2 ),
                vsapi->getFrameWidth( frame, 2 ), vsapi->getFrameHeight( frame, 2 ) );
    }
}

void PlanarFrame::copyTo( VSFrameRef *frame, const VSAPI *vsapi )
{
    if( !y || !u || !v ) return;
    if( vsapi->getFrameFormat( frame )->id == pfYUV420P8 )
    {
        BitBlt( vsapi->getWritePtr( frame, 0 ), vsapi->getStride( frame, 0 ), y,  ypitch,  ywidth,  yheight );
        BitBlt( vsapi->getWritePtr( frame, 1 ), vsapi->getStride( frame, 1 ), u, uvpitch, uvwidth, uvheight );
        BitBlt( vsapi->getWritePtr( frame, 2 ), vsapi->getStride( frame, 2 ), v, uvpitch, uvwidth, uvheight );
    }
}

unsigned char *PlanarFrame::GetPtr( int plane )
{
    if( plane == 0 ) return y;
    if( plane == 1 ) return u;
    return v;
}

int PlanarFrame::GetWidth( int plane )
{
    if( plane == 0 )
        return ywidth;
    else
        return uvwidth;
}

int PlanarFrame::GetHeight( int plane )
{
    if( plane == 0 )
        return yheight;
    else
        return uvheight;
}

int PlanarFrame::GetPitch( int plane )
{
    if( plane == 0 )
        return ypitch;
    else
        return uvpitch;
}

PlanarFrame::PlanarFrame( const VSVideoInfo &viInfo )
{
    ypitch = uvpitch = 0;
    ywidth = uvwidth = 0;
    yheight = uvheight = 0;
    y = u = v = nullptr;
    allocSpace( viInfo );
}

PlanarFrame::~PlanarFrame()
{
    freeSpace();
}

bool PlanarFrame::allocSpace( const VSVideoInfo &viInfo )
{
    freeSpace();
    int height = viInfo.height;
    int width  = viInfo.width;
    if( viInfo.format->id == pfYUV420P8 )
    {
        static const int min_alignment = 32;
        ypitch  = width + ((width % min_alignment) == 0 ? 0 : min_alignment - (width % min_alignment));
        ywidth  = width;
        yheight = height;
        width  >>= 1;
        height >>= 1;
        uvpitch  = width + ((width % min_alignment) == 0 ? 0 : min_alignment - (width % min_alignment));
        uvwidth  = width;
        uvheight = height;
        y = static_cast<unsigned char *>(AlignedMemory::alloc(  ypitch * yheight,  min_alignment ));
        if( !y ) return false;
        u = static_cast<unsigned char *>(AlignedMemory::alloc( uvpitch * uvheight, min_alignment ));
        if( !u ) return false;
        v = static_cast<unsigned char *>(AlignedMemory::alloc( uvpitch * uvheight, min_alignment ));
        if( !v ) return false;
        return true;
    }
    return false;
}

void PlanarFrame::freeSpace()
{
    if( y ) { AlignedMemory::free( y ); y = nullptr; }
    if( u ) { AlignedMemory::free( u ); u = nullptr; }
    if( v ) { AlignedMemory::free( v ); v = nullptr; }
}
