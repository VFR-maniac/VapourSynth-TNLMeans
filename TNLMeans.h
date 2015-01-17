/*
**                 TNLMeans for VapourSynth
**
**   TNLMeans is an implementation of the NL-means denoising algorithm.
**   Aside from the original method, TNLMeans also supports extension
**   into 3D, a faster, block based approach, and a multiscale version.
**
**   Copyright (C) 2006-2007 Kevin Stone
**   Copyright (C) 2015      Yusuke Nakamura
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

#include <cmath>
#include <cstring>
#include <algorithm>
#include <limits>

#ifdef __MINGW32__
#include "mingw.thread.h"
#include "mingw.mutex.h"
#else
#include <thread>
#include <mutex>
#endif

struct SDATA
{
    double *weights;
    double *sums;
    double *wmaxs;
};

class nlFrame
{
public:
    int               fnum;
    const VSAPI      *vsapi;
    const VSFrameRef *pf;
    SDATA           **ds;
    int              *dsa;
    nlFrame( bool _useblocks, int _size, const VSVideoInfo &vi, const VSAPI *_vsapi );
    ~nlFrame();
    void setFNum( int i );
};

class nlCache
{
public:
    nlFrame **frames;
    int start_pos, size;
    nlCache( int _size, bool _useblocks, const VSVideoInfo &vi, const VSAPI *vsapi );
    ~nlCache();
    void resetCacheStart( int first, int last );
    int  getCachePos    ( int n );
    void clearDS        ( nlFrame *nl );
};

class nlThread
{
public:
    int      active;
    double  *sumsb, *weightsb, *gw;
    nlCache *fc;
    SDATA   *ds;
    nlThread();
    ~nlThread();
};

class TNLMeans
{
private:
    int       Ax, Ay, Az;
    int       Sx, Sy;
    int       Bx, By;
    int       Sxd, Syd, Sxa;
    int       Bxd, Byd, Bxa;
    int       Axd, Ayd, Axa, Azdm1;
    double    a, a2;
    double    h, hin, h2in;
    bool      ssd;
    int       numThreads;
    nlThread *threads;
    std::mutex mtx;
    int mapn( int n );
    inline double GetSSD( const unsigned char *s1, const unsigned char *s2, const double *gwT, const int k ) { return (s1[k] - s2[k]) * (s1[k] - s2[k]) * gwT[k]; };
    inline double GetSAD( const unsigned char *s1, const unsigned char *s2, const double *gwT, const int k ) { return std::abs( s1[k] - s2[k] ) * gwT[k]; };
    inline double GetSSDWeight( const double diff, const double gweights ) { return std::exp( (diff / gweights) * h2in ); };
    inline double GetSADWeight( const double diff, const double gweights ) { return std::exp( (diff / gweights) * hin ); };
    VSFrameRef *newVideoFrame( int n, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    template < int ssd > VSFrameRef *GetFrameByMethod( int n, int thread, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    template < int ssd > VSFrameRef *GetFrameWZ      ( int n, int thread, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    template < int ssd > VSFrameRef *GetFrameWZB     ( int n, int thread, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    template < int ssd > VSFrameRef *GetFrameWOZ     ( int n, int thread, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    template < int ssd > VSFrameRef *GetFrameWOZB    ( int n, int thread, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );

public:
    VSVideoInfo vi;
    VSNodeRef  *node;
    void RequestFrame( int n, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    VSFrameRef *GetFrame( int n, VSFrameContext *frame_ctx, VSCore *core, const VSAPI *vsapi );
    typedef class {} bad_param;
    typedef class {} bad_alloc;
    /* Constructor */
    TNLMeans
    (
        int _Ax, int _Ay, int _Az,
        int _Sx, int _Sy,
        int _Bx, int _By,
        double _a, double _h, bool ssd,
        const VSMap *in,
        VSMap       *out,
        VSCore      *core,
        const VSAPI *vsapi
    );
    /* Destructor */
    ~TNLMeans();
};

static inline void fill_zero_d( double *x, size_t n )
{
    if( std::numeric_limits<double>::is_iec559 )
        std::memset( x, 0, n * sizeof(double) );
    else
        std::fill_n( x, n, 0.0 );
}
