/*--------------------------------------------------------------------
 *    Basil / Sybil:   c_funcs.c   1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define ANSI_DECLARATORS 1

void hostname_(name, len)
char *name;
int *len;
{
    char *tmp = 0;
    if ( (tmp = getenv("HOSTNAME")) ||
            (tmp = getenv("COMPUTERNAME")) ) {
        strncpy(name,tmp,*len-1);
        *len=strlen(name);
    }

}

/*
 *   returns yy/mm/dd in day[]
 *   returns seconds after midnight if flag=0
 *   else seconds since last call with flag=1
 */

int datime_(day,secs,flag)
int *day, *secs, *flag;
{
    time_t now;
    struct tm *timedat;
    static long jdays,elapsedsecs;

    time(&now);
    timedat = localtime(&now);
    day[0] = timedat->tm_year;
    day[1] = timedat->tm_mon+1;
    day[2] = timedat->tm_mday;

    *secs = timedat->tm_hour*3600 + timedat->tm_min*60 + timedat->tm_sec;
    if (*flag) {
        *secs -= (elapsedsecs-(timedat->tm_yday - jdays)*86400);
    }
    else {
        elapsedsecs = *secs;
        jdays = timedat->tm_yday;
    }
    return(0);
}

/* find maximum of a and b */
#define MAX(a,b)        (((a)>(b))?(a):(b))
#define MIN(a,b)        (((a)<(b))?(a):(b))
#define ABS(a)          (((a)<0) ? -(a) : (a))

#ifdef  ANSI_DECLARATORS
int PntOnLine(int ipx,int ipy,int iqx,int iqy,int itx,int ity)

#else
int PntOnLine(ipx, ipy, iqx, iqy, itx, ity)
int ipx, ipy, iqx, iqy, itx, ity;
#endif
{    
    long px,py,qx,qy,tx,ty;
long f,max;

    px=(long)ipx;
    py=(long)ipy;
    qx=(long)iqx;
    qy=(long)iqy;
    tx=(long)itx;
    ty=(long)ity;
    f= (qy-py)*(tx-px)-(ty-py)*(qx-px);
    max = MAX(ABS(qx-px), ABS(qy-py));
    if ( ABS((qy-py)*(tx-px)-(ty-py)*(qx-px)) >=
        (MAX(ABS(qx-px), ABS(qy-py)))) return(0);
    if (((qx<px)&&(px<tx)) || ((qy<py)&&(py<ty))) return(1);
    if (((tx<px)&&(px<qx)) || ((ty<py)&&(py<qy))) return(1);
    if (((px<qx)&&(qx<tx)) || ((py<qy)&&(qy<ty))) return(3);
    if (((tx<qx)&&(qx<px)) || ((ty<qy)&&(qy<py))) return(3);
    return(2);
}

#ifdef  ANSI_DECLARATORS
int pntonsegment_(float *x, float *y,float *x1, float *y1,
                   float *x2, float *y2)
#else
int pntonsegment_(x, y,x1, y1, x2, y2)
float *x, *y,*x1, *y1, *x2, *y2;
#endif
{
    int res=0;

#if XY
    int ix1, ix2, iy1, iy2, ix, iy;
    if (*x1<0) ix1 = (int)((*x1-5E-7)*1.0E6);
    else ix1 = (int)((*x1+5E-7)*1.0E6);
    if (*y1<0) iy1 = (int)((*y1-5E-7)*1.0E6);
    else iy1 = (int)((*y1+5E-7)*1.0E6);
    if (*x2<0) ix2 = (int)((*x2-5E-7)*1.0E6);
    else ix2 = (int)((*x2+5E-7)*1.0E6);
    if (*y2<0) iy2 = (int)((*y2-5E-7)*1.0E6);
    else iy2 = (int)((*y2+5E-7)*1.0E6);
    if (*x<0) ix = (int)((*x-5E-7)*1.0E6);
    else  ix = (int)((*x+5E-7)*1.0E6);
    if (*y<0) iy = (int)((*y-5E-7)*1.0E6);
    else  iy = (int)((*y+5E-7)*1.0E6);
#endif
    float eps = 5.0e-6;
    float xmin,xmax,ymin,ymax;
    float val1,val2;

    if (*x1 < *x2) {
        xmin = *x1;
        xmax = *x2;
    }
    else {
        xmin = *x2;
        xmax = *x1;
    }
    if (*y1 < *y2) {
        ymin = *y1;
        ymax = *y2;
    }
    else {
        ymin = *y2;
        ymax = *y1;
    }
    if ((*x > xmin-eps)&&(*x < xmax+eps)&&
        (*y > ymin-eps)&&(*y < ymax+eps)) {
        val1 = ABS((*y-*y1)*(xmax-xmin));
        val2 = ABS((ymax-ymin)*(*x-*x1));
        if (val1<=(val2+eps) && val1>=(val2-eps)) res=2;
    }
    return((res==2)?1:0);
}

/*
 * Graphics Gem
 */
/*
 * conflict with wingdi.h if using mingw (see above)
 * but we always define WINDING so just wrap the line
 * to get rid of compiler warning
 */
#ifndef __WIN32__
#endif
#define WINDING

#define X	0
#define Y	1

#ifndef TRUE
#define TRUE	1
#define FALSE	0
#endif


/* test if a & b are within epsilon.  Favors cases where a < b */
#define Near(a,b,eps)	( ((b)-(eps)<(a)) && ((a)-(eps)<(b)) )


/* ======= Crossings algorithm ============================================ */

/* Shoot a test ray along +X axis.  The strategy, from MacMartin, is to
 * compare vertex Y values to the testing point's Y and quickly discard
 * edges which are entirely to one side of the test ray.
 *
 * Input 2D polygon _pgon_ with _numverts_ number of vertices and test point
 * _point_, returns 1 if inside, 0 if outside.	WINDING and CONVEX can be
 * defined for this test.
 */
int crossingstest_( double  *pgon, int *nverts, double  point[2],
                    int *res )
{
#ifdef	WINDING
register int	crossings ;
#endif
register int	j, yflag0, yflag1, inside_flag, xflag0, numverts ;
register double ty, tx, *vtx0, *vtx1 ;
#ifdef	CONVEX
register int	line_flag ;
#endif

    numverts = *nverts;
    tx = point[X] ;
    ty = point[Y] ;

    vtx0 = pgon+ 2*(numverts-1) ;
    /* get test bit for above/below X axis */
    yflag0 = ( vtx0[Y] >= ty ) ;
    vtx1 = pgon;

#ifdef	WINDING
    crossings = 0 ;
#else
    inside_flag = 0 ;
#endif
#ifdef	CONVEX
    line_flag = 0 ;
#endif
    for ( j = numverts+1 ; --j ; ) {

	yflag1 = ( vtx1[Y] >= ty ) ;
	/* check if endpoints straddle (are on opposite sides) of X axis
	 * (i.e. the Y's differ); if so, +X ray could intersect this edge.
	 */
	if ( yflag0 != yflag1 ) {
	    xflag0 = ( vtx0[X] >= tx ) ;
	    /* check if endpoints are on same side of the Y axis (i.e. X's
	     * are the same); if so, it's easy to test if edge hits or misses.
	     */
	    if ( xflag0 == ( vtx1[X] >= tx ) ) {

		/* if edge's X values both right of the point, must hit */
#ifdef	WINDING
		if ( xflag0 ) crossings += ( yflag0 ? -1 : 1 ) ;
#else
		if ( xflag0 ) inside_flag = !inside_flag ;
#endif
	    } else {
		/* compute intersection of pgon segment with +X ray, note
		 * if >= point's X; if so, the ray hits it.
		 */
		if ( (vtx1[X] - (vtx1[Y]-ty)*
		     ( vtx0[X]-vtx1[X])/(vtx0[Y]-vtx1[Y])) >= tx ) {
#ifdef	WINDING
		    crossings += ( yflag0 ? -1 : 1 ) ;
#else
		    inside_flag = !inside_flag ;
#endif
		}
	    }
#ifdef	CONVEX
	    /* if this is second edge hit, then done testing */
	    if ( line_flag ) goto Exit ;

	    /* note that one edge has been hit by the ray's line */
	    line_flag = TRUE ;
#endif
	}

	/* move to next pair of vertices, retaining info as possible */
	yflag0 = yflag1 ;
	vtx0 = vtx1 ;
	vtx1 += 2 ;
    }
#ifdef	CONVEX
    Exit: ;
#endif
#ifdef	WINDING
    /* test if crossings is not zero */
    inside_flag = (crossings != 0) ;
#endif

    *res = inside_flag;
    return( inside_flag ) ;
}
