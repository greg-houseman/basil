/*--------------------------------------------------------------------
 *    Basil / Sybil:   plotdefs.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#ifndef _plottype_h
/* plot_type */
#define XYPLOT    2
#define PRFL      3
#define ARROWS    4
#define CNTRS     5
#define BLANK     6
#define PS        7

#define PRFL1D   31
#define PRFL2D   32
#define PRFLMRK  33

/* plot_description - use bitmask */
#define LINES	  2
#define SHADE     4
#define DIM1      8
#define DIM2_X   16
#define DIM2_Y   32

/* orientation of colourbar for shaded contour plots */
#define CB_VERT	  42
#define CB_HORIZ  41
#define CB_NONE   40

/* check LINECNTOUR (arrow.f,pref77.c) if this order is changed */
#define CNTR_MIN_U 0
#define CNTR_MAX_U 1
#define CNTR_SCL   2
#define CNTR_LVL_U 3
#define CNTR_STEP_U  4
#define MAX_CNTRS  5
#define NUM_CNTRS  6
#define CNTR_MIN   7
#define CNTR_MAX   8
#define CNTR_STPL  9
#define MRK_PRFL  10
#define CNTR_LVL  11
#define CNTR_STEP 12

#define MAXCNTRVALS  13

#define PRFL_MIN_U 0
#define PRFL_MAX_U 1
#define PRFL_SCL   2
#define PRFL_SHFT  3
#define PRFL_MIN   4
#define PRFL_MAX   5
#define PRFL_X1    6 /* maintain this order */
#define PRFL_Y1    7 /* for these 4 points- */
#define PRFL_X2    8 /* x1,y1,x2,y2         */
#define PRFL_Y2    9
#define NUM_PTS   10
#define LWR_LIM   11
#define UPR_LIM   12

#define MAXPRFLVALS   13

/* flags used by motif app */
#define CNTRVALS   10
#define PRFLVALS   11
#define PRFLPTS    12
#define PRFLDIR    13
#define ZOOMVALS   14

/*
 * correspond to xsyb numbers - sent to PLMESH
 */
#define ELMNT             0
#define ELMNTNUM          1
#define BNDRY             3
#define ELMNTDBL          4
#define VISCMSH           5
#define VISCDBL           6
#define BNDRYDBL          7
#define LGELMNT           8
#define LGBNDRY           9
#define INTBND           10
#define INTBNDDBL        11
#define SEMSH            12
#define SEDBL            13
#define INTLBNDRY        14
#define INTLBNDRYDBL     15
#define _plottype_h
#endif
