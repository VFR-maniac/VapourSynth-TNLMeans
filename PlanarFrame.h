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

#ifndef __PlanarFrame_H__
#define __PlanarFrame_H__

class PlanarFrame
{
private:
    int numPlanes;
    int stride[3];
    int width [3];
    int height[3];
    unsigned char *data[3];
    bool allocSpace( const VSVideoInfo &viInfo );
    void freeSpace();
    inline void BitBlt( void *dstp, int dst_stride, const void *srcp, int src_stride, int row_size, int height );

public:
    void copyTo( VSFrameRef *frame, const VSAPI *vsapi );
    unsigned char * GetPtr   ( int plane );
    int             GetWidth ( int plane );
    int             GetHeight( int plane );
    int             GetPitch ( int plane );
    /* Constructors */
    PlanarFrame( const VSVideoInfo &viInfo );
    /* Destructor */
    ~PlanarFrame();
};

#endif