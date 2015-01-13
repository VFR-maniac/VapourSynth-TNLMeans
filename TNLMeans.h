/*
**                 TNLMeans for VapourSynth
**
**   TNLMeans is an implementation of the NL-means denoising algorithm.
**   Aside from the original method, TNLMeans also supports extension
**   into 3D, a faster, block based approach, and a multiscale version.
**
**   Copyright (C) 2006-2007 Kevin Stone
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

#include <cmath>        /* std::exp */

struct SDATA
{
    double *weights;
    double *sums;
    double *wmaxs;
};

class nlFrame
{
public:
    int          fnum;
    PlanarFrame *pf;
    SDATA      **ds;
    int         *dsa;
    nlFrame( bool _useblocks, int _size, const VSVideoInfo &vi );
    ~nlFrame();
    void setFNum( int i );
};

class nlCache
{
public:
    nlFrame **frames;
    int start_pos, size;
    nlCache( int _size, bool _useblocks, const VSVideoInfo &vi );
    ~nlCache();
    void resetCacheStart( int first, int last );
    int  getCachePos    ( int n );
    void clearDS        ( nlFrame *nl );
};

class TNLMeans
{
private:
    int                Ax, Ay, Az;
    int                Sx, Sy;
    int                Bx, By;
    int                Sxd, Syd, Sxa;
    int                Bxd, Byd, Bxa;
    int                Axd, Ayd, Axa, Azdm1;
    double             a, a2;
    double             h, hin, h2in;
    double            *sumsb, *weightsb, *gw;
    bool               ssd;
    nlCache           *fc;
    SDATA             *ds;
    PlanarFrame       *dstPF, *srcPFr;
    int mapn( int n );
    inline double GetSSD( const unsigned char *s1, const unsigned char *s2, const double *gwT, const int k ) { return (s1[k] - s2[k]) * (s1[k] - s2[k]) * gwT[k]; };
    inline double GetSAD( const unsigned char *s1, const unsigned char *s2, const double *gwT, const int k ) { return std::abs( s1[k] - s2[k] ) * gwT[k]; };
    inline double GetSSDWeight( const double diff, const double gweights ) { return std::exp( (diff / gweights) * h2in ); };
    inline double GetSADWeight( const double diff, const double gweights ) { return std::exp( (diff / gweights) * hin ); };
    VSFrameRef *CopyTo          ( int n, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    template < int ssd > VSFrameRef *GetFrameByMethod( int n, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    template < int ssd > VSFrameRef *GetFrameWZ      ( int n, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    template < int ssd > VSFrameRef *GetFrameWZB     ( int n, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    template < int ssd > VSFrameRef *GetFrameWOZ     ( int n, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    template < int ssd > VSFrameRef *GetFrameWOZB    ( int n, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );

public:
    VSVideoInfo vi;
    VSNodeRef  *node;
    void RequestFrame( int n, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    VSFrameRef *GetFrame( int n, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    /* Constructor */
    TNLMeans
    (
        int _Ax, int _Ay, int _Az,
        int _Sx, int _Sy,
        int _Bx, int _By,
        double _a, double _h, bool ssd,
        const VSMap *in,
        VSMap       *out,
        const VSAPI *vsapi
    );
    /* Desctructor */
    ~TNLMeans();
};
