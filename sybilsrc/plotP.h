/*--------------------------------------------------------------------
 *    Basil / Sybil:   plotP.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/* 
 * plotP.h - Private definitions for Plot widget
 */

#ifndef _PlotP_h
#define _PlotP_h

/*
 * Include private header file of superclass.
 * In this case, however, Primitive's header file is
 * thrown into XmP.h.
 */
#include <Xm/Xm.h>
#if XmVersion<1002
#include <X11/IntrinsicP.h>
#include <X11/CoreP.h>
#else
#include <Xm/XmP.h>   
#include <Xm/PrimitiveP.h>
#endif

/*
 * Include public header file for this widget.
 */
#include "plot.h"
#include "types.h"

/* New fields for the Plot widget class record */

typedef struct {
    int make_compiler_happy;    /* keep compiler happy */
} PlotClassPart;

/* Full class record declaration */
typedef struct _PlotClassRec {
    CoreClassPart    core_class;
#if XmVersion>1001
    XmPrimitiveClassPart    primitive_class;
#endif
    PlotClassPart    plot_class;
} PlotClassRec;

extern PlotClassRec plotClassRec;

/* New fields for the Plot widget record */
typedef struct {
    /* resources */
    Pixel    foreground;
    XtCallbackList callback;    /* application installed callback function(s) */
    Dimension    pixmap_width_in_cells;
    Dimension     pixmap_height_in_cells;
    Dimension	 cell_width_in_pixels;
    Dimension	 cell_height_in_pixels;
    float	 cell_width_in_cm;
    float	 cell_height_in_cm;
    int cur_x, cur_y;  /* position of visible corner in big pixmap */
    cell_data *cell;    /* array for data for each cell*/
    save_node *plots; /* list of plots, in display order */
    label_node *labels; /* list of optional labels */
    Boolean show_all;  /* whether bitmap should display entire bitmap */

    /* private state */
    Dimension    pixmap_width_in_pixels;
    Dimension    pixmap_height_in_pixels;
    Pixmap big_picture;
    GC    draw_gc;
    GC    undraw_gc;
    GC    border_gc;
    Boolean user_allocated;  /* whether user allocated cell array */
} PlotPart;

/*
 * Full instance record declaration
 */
typedef struct _PlotRec {
    CorePart            core;
#if XmVersion>1001
    XmPrimitivePart        primitive;
#endif
    PlotPart        plot;
} PlotRec;

#endif /* _PlotP_h */
