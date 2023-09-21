/*--------------------------------------------------------------------
 *    Basil / Sybil:   plot.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#ifndef _Plot_h
#define _Plot_h

/*
 * Plot Widget public include file
 */

/*
 * The public header file for the immediate superclass normally
 * must be included.  However, not in this case because the public
 * header file for Primitive is in Xm.h, which is already included
 * in all Motif applications.
 */

/* #include <Xm/Superclass.h>  */

/* 
 * This public structure is used as call_data to the callback.
 * It passes the x, y position of the cell toggled (in units of
 * cells, not pixels) and a mode flag that indicates whether the
 * cell was turned on (1) or off (0).
 */
typedef struct {
	int mode;
	int newx;
	int newy;
} PlotPointInfo;

#define XtNtoggleCallback "toggleCallback"
#define XtNcellWidthInPixels "cellWidthInPixels"
#define XtNcellHeightInPixels "cellHeightInPixels"
#define XtNcellWidthInCm "cellWidthInCm"
#define XtNcellHeightInCm "cellHeightInCm"
#define XtNpixmapWidthInCells "pixmapWidthInCells"
#define XtNpixmapHeightInCells "pixmapHeightInCells"
#define XtNcurX "curX"
#define XtNcurY "curY"
#define XtNcellArray "cellArray"
#define XtNplotList "plotList"
#define XtNlabelList "labelList"
#define XtNshowEntireBitmap "showEntireBitmap"
 
#define XtCToggleCallback "ToggleCallback"
#define XtCCellWidthInPixels "CellWidthInPixels"
#define XtCCellHeightInPixels "CellHeightInPixels"
#define XtCCellWidthInCm "CellWidthInCm"
#define XtCCellHeightInCm "CellHeightInCm"
#define XtCPixmapWidthInCells "PixmapWidthInCells"
#define XtCPixmapHeightInCells "PixmapHeightInCells"
#define XtCCurX "CurX"
#define XtCCurY "CurY"
#define XtCCellArray "CellArray"
#define XtCPlotList "PlotList"
#define XtCLabelList "LabelList"
#define XtCShowEntireBitmap "ShowEntireBitmap"

/*extern char *BitmapEditGetArray();  w */
    /* Widget w; */

/* Class record constants */

extern WidgetClass plotWidgetClass;

typedef struct _PlotClassRec *PlotWidgetClass;
typedef struct _PlotRec      *PlotWidget;

/* modes for drawing */
#define DRAW 1
#define UNDRAW 0
 
extern int Save_option();
extern int Write_plot_data(Widget w,FILE *fp);
extern void UndrawCellBorder();
extern void DrawCellBorder();
extern void UpdateCell();
extern void DrawLine();

extern int Mesh[MESH_ENTRIES];
extern float Pwindo[PWINDO_ENTRIES];
#endif /* _Plot_h */
/* DON'T ADD AFTER THIS #endif */
