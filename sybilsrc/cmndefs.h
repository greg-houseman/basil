/*--------------------------------------------------------------------
 *    Basil / Sybil:   cmndefs.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
#ifndef _S_cmndefs_h
#define _S_cmndefs_h

#define SYB_FILENAME_MAX    255
#define USEFL_SIZE            8
#define FOREGROUND            1
#define BACKGROUND            0
#define PLOT_COLOURS         64
#define MAX_COLOURS   PLOT_COLOURS+USEFL_SIZE

#define GREY_MAP	   100
#define STD_MAP 	   101
#define ABS_MAP 	   102

#define DEFAULT_OPTIONS_FILE "sybil.in"
#define DEFAULT_LOG_FILE "sybil.log"
#define DEFAULT_LOG_EXT ".log"
#define DEFAULT_COLUMNS    0  /* in cells */
#define DEFAULT_ROWS       0  /* in cells */
#define DEFAULT_FONTHGT   12
#define DEFAULT_LGE_FONTHGT   24

/* window defaults to A4 (29.7x21.0cm) */
/* reduce to 28.4x19.7cm to allow for printer clipping */
/* USLetter (21.59x27.94cm) reduce to 20.29x26.64 */
#define DEFAULTWIDTH 197.0
#define DEFAULTHEIGHT 284.0
#define USLETTERWIDTH 202.9
#define USLETTERHEIGHT 266.4
#define A4_PAPER           0
#define US_PAPER           1

#define PORTRAIT           0
#define LANDSCAPE          1
 
#define MINWIDTH          1
#define MAXWIDTH         10    /* max no. cells/row */
#define MINHEIGHT         1
#define MAXHEIGHT        10    /* max no. cells/col */

#define NNX3_SIZE		423
#define MAX_LABEL_LEN 80
#define MAX_NUM_LEN 12
#define MAXNAME 30
#define PEN_UP 3
#define PEN_DN 2
#define J_CENTRE 0
#define J_BASE   1
#define J_TOP    2
#define J_LEFT   1
#define J_RIGHT  2
#define TIC_NONE 3

#define XLMARG    0.15
#define XRMARG    0.2
#define YMARG    0.1
#define PG_XMARG  0.0
#define PG_YMARG  0.0
/*
#define XLMARG    0.05
#define XRMARG    0.05
#define YMARG    0.05
*/
 
/* plot parameters */
#define SYB_VELOCITY      0
#define SYB_STRAIN        1
#define SYB_STRESS        2
#define STRN_MRKR         3
#define THICKNESS         4 /* don't break old files */
#define LAYER             5
#define DEFORM            6
#define MESH              7
#define BBOX              8
#define LGMESH            9
#define DENSITY          10
#define ROTATION         11
#define ELLE             12
#define GRAVITY          13

/* fix these when dashln() coordinated in xapk and lppak */
#define SOLID             8
#define DASH              1
#define DOT               7

#define FILE_OPT          0
#define DATA_IN           1
#define XY_PLOT           2
#define LOCATE            6
#define OPTIONS           7
#define HELP              8

#define TITLE             20
#define OPTS              21
#define COMMENT           22
#define LABEL             23

/* pwindo indices */
#define XINIT   0
#define YINIT   1
#define XMIN    2
#define XMAX    3
#define YMIN    4
#define YMAX    5
#define SCLX    6
#define SCLY    7
#define SIZE    8
#define XCMIN   9
#define XCMAX   10
#define YCMIN   11
#define YCMAX   12
#define ULEFT   13
#define URGHT   14
#define UBASE   15
#define UTOP    16
#define XCMINREF 17
#define XCMAXREF 18
#define YCMINREF 19
#define YCMAXREF 20
#define XMARGE 21
#define YMARGE 22

#define PWINDO_ENTRIES   23

/* mesh indices */
#define NX0    0
#define NX1    1
#define NX2    2
#define NX3    3
#define NY0    4
#define NY1    5
#define NY2    6
#define NY3    7
#define NXY    8
#define NNX3   9
#define NNY3  10
#define NP3	  11
#define NWK	  12 /* not used */
#define ICON  13

#define MESH_ENTRIES   14


#ifndef MIN
#define MIN(x,y)    ((x) < (y) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x,y)    ((x) > (y) ? (x) : (y))
#endif

#endif
