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
    for( int i = 0; i < numPlanes; ++i )
        if( data[i] )
            BitBlt( data[i], stride[i],
                    vsapi->getReadPtr   ( frame, i ), vsapi->getStride     ( frame, i ),
                    vsapi->getFrameWidth( frame, i ), vsapi->getFrameHeight( frame, i ) );
}

void PlanarFrame::copyTo( VSFrameRef *frame, const VSAPI *vsapi )
{
    for( int i = 0; i < numPlanes; ++i )
        if( data[i] )
            BitBlt( vsapi->getWritePtr( frame, i ), vsapi->getStride( frame, i ),
                    data[i], stride[i], width[i], height[i] );
}

unsigned char *PlanarFrame::GetPtr( int plane )
{
    return plane < numPlanes ? data[plane] : nullptr;
}

int PlanarFrame::GetWidth( int plane )
{
    return plane < numPlanes ? width[plane] : 0;
}

int PlanarFrame::GetHeight( int plane )
{
    return plane < numPlanes ? height[plane] : 0;
}

int PlanarFrame::GetPitch( int plane )
{
    return plane < numPlanes ? stride[plane] : 0;
}

PlanarFrame::PlanarFrame( const VSVideoInfo &viInfo )
{
    stride[0] = stride[1] = stride[2] = 0;
    width [0] = width [1] = width [2] = 0;
    height[0] = height[1] = height[2] = 0;
    data  [0] = data  [1] = data  [2] = nullptr;
    numPlanes = viInfo.format->numPlanes;
    allocSpace( viInfo );
}

PlanarFrame::~PlanarFrame()
{
    freeSpace();
}

bool PlanarFrame::allocSpace( const VSVideoInfo &viInfo )
{
    freeSpace();
    for( int i = 0; i < numPlanes; ++i )
    {
        static const int min_alignment = 32;
        width [i] = viInfo.width  >> (i ? viInfo.format->subSamplingW : 0);
        height[i] = viInfo.height >> (i ? viInfo.format->subSamplingH : 0);
        stride[i] = width[i] + ((width[i] % min_alignment) == 0 ? 0 : min_alignment - (width[i] % min_alignment));
        data[i] = static_cast<unsigned char *>(AlignedMemory::alloc( stride[i] * height[i], min_alignment ));
        if( !data[i] ) return false;
    }
    return true;
}

void PlanarFrame::freeSpace()
{
    if( data[0] ) { AlignedMemory::free( data[0] ); data[0] = nullptr; }
    if( data[1] ) { AlignedMemory::free( data[1] ); data[1] = nullptr; }
    if( data[2] ) { AlignedMemory::free( data[2] ); data[2] = nullptr; }
}
