/*--------------------------------------------------------------------
     Basil / Sybil:   plot.c  1.1  1 Octoner 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#include <Xm/Xm.h>
#if XmVersion>1001
#include <Xm/XmP.h>
#endif
#include <X11/StringDefs.h>
#include <X11/Intrinsic.h>
#include <X11/Xatom.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cmndefs.h"
#include "plotP.h"
#include "plot.h"
#include "types.h"
#include "globals.h"
#include "strain.h"
#include "deform.h"
#include "string_utils.h"
#include "log.h"
#include "errnum.h"
#include "error.h"

#define XY 0

#define DEFAULT_CELL_WIDTH   200  /* in pixels */
#define DEFAULT_CELL_HEIGHT  100  /* in pixels */

#define DEFAULT_CELL_WIDTH_CM   148.5  /* in cm */
#define DEFAULT_CELL_HEIGHT_CM  DEFAULTHEIGHT/DEFAULT_PIXMAP_HEIGHT/* in cm */

/* values for instance variable is_drawn */
#define DRAWN 1
#define UNDRAWN 0

#define MAXLINES 10   /* max of horiz or vertical cells */
#define SCROLLBARWIDTH 15

#define offset(field) XtOffsetOf(PlotRec, field)

int  Clear_saved_data();
int  Clear_all_saved_data();
int  FindCell();
int  store_label();
int  write_options();
void CreateBigPixmap();
void cms_to_pixels();
void InitBigPixmap();
extern void init_box();
extern void formatnumber_();
extern void setpencolor_();
extern void plotu_();
extern void xpak_init();

static XtResource resources[] = {
/* not needed in Motif - already defined by Primitive.
 *    {
 *   XtNforeground, 
 *   XtCForeground, 
 *   XtRPixel, 
 *   sizeof(Pixel),
 *   offset(plot.foreground), 
 *   XtRString, 
 *   XtDefaultForeground
 *    },
 */
     {
    XtNtoggleCallback, 
    XtCToggleCallback, 
    XtRCallback, 
    sizeof(XtPointer),
    offset(plot.callback), 
    XtRCallback, 
    NULL
     },
     {
    XtNpixmapWidthInCells, 
    XtCPixmapWidthInCells, 
    XtRDimension, 
    sizeof(Dimension),
    offset(plot.pixmap_width_in_cells), 
    XtRImmediate, 
    (XtPointer)DEFAULT_COLUMNS
     },
     {
    XtNpixmapHeightInCells, 
    XtCPixmapHeightInCells, 
    XtRDimension, 
    sizeof(Dimension),
    offset(plot.pixmap_height_in_cells), 
    XtRImmediate, 
    (XtPointer)DEFAULT_ROWS
     },
     {
    XtNcellWidthInPixels, 
    XtCCellWidthInPixels, 
    XtRDimension, sizeof(Dimension),
    offset(plot.cell_width_in_pixels), 
    XtRImmediate, 
    (XtPointer)DEFAULT_CELL_WIDTH
     },
     {
    XtNcellHeightInPixels, 
    XtCCellHeightInPixels, 
    XtRDimension, sizeof(Dimension),
    offset(plot.cell_height_in_pixels), 
    XtRImmediate, 
    (XtPointer)DEFAULT_CELL_HEIGHT
     },
/*
     {
    XtNcellWidthInCm, 
    XtCCellWidthInCm, 
    XtRFloat, sizeof(float),
    offset(plot.cell_width_in_cm), 
    XtRImmediate, 
    (XtPointer)DEFAULT_CELL_WIDTH_CM
     },
     {
    XtNcellHeightInCm, 
    XtCCellHeightInCm, 
    XtRFloat, sizeof(float),
    offset(plot.cell_height_in_cm), 
    XtRImmediate, 
    (XtPointer)DEFAULT_CELL_HEIGHT_CM
     },
*/
     {
    XtNcurX, 
    XtCCurX, 
    XtRInt, 
    sizeof(int),
    offset(plot.cur_x), 
    XtRImmediate, 
    (XtPointer) 0
     },
     {
    XtNcurY, 
    XtCCurY, 
    XtRInt, 
    sizeof(int),
    offset(plot.cur_y), 
    XtRImmediate, 
    (XtPointer) 0
     },
     {
    XtNplotList, 
    XtCPlotList, 
    XtRString, 
    sizeof(String),
    offset(plot.plots), 
    XtRImmediate, 
    (XtPointer) 0
     },
     {
    XtNlabelList, 
    XtCLabelList, 
    XtRString, 
    sizeof(String),
    offset(plot.labels), 
    XtRImmediate, 
    (XtPointer) 0
     },
     {
    XtNcellArray, 
    XtCCellArray, 
    XtRString, 
    sizeof(String),
    offset(plot.cell), 
    XtRImmediate, 
    (XtPointer) 0
     },
     {
    XtNshowEntireBitmap, 
    XtCShowEntireBitmap, 
    XtRBoolean, 
    sizeof(Boolean),
    offset(plot.show_all), 
    XtRImmediate, 
    (XtPointer) TRUE
     },
};

/* Declaration of methods */

static void Initialize();
static void Redisplay();
static void Destroy();
static void Resize();
static Boolean SetValues();
static XtGeometryResult QueryGeometry();

/* these Core methods not needed by BitmapEdit:
 *
 * static void ClassInitialize();
 * static void Realize();
 */

/* the following are private functions unique to plotArea */
static void DrawBorder(), ChangeCellSize();
/*  static void DoCell();  */

/* the following are actions of plotArea */
static void UndrawCell(), OptionalLabel();

static char defaultTranslations[] =
    "<Btn1Down>:    OptionalLabel()        \n\
    <Btn2Down>:    UndrawCell()            \n\
    <Btn2Motion>:  UndrawCell()";

static XtActionsRec actions[] = {
        {"OptionalLabel", OptionalLabel},
        {"UndrawCell", UndrawCell},
};

/* definition in plot.h */
static PlotPointInfo info;

PlotClassRec plotClassRec = {
    {
    /* core_class fields */
#if XmVersion>1001
    /* superclass        */ (WidgetClass) &xmPrimitiveClassRec,
#else
    /* superclass        */ (WidgetClass) &coreClassRec,
#endif
    /* class_name        */ "Plot",
    /* widget_size       */ sizeof(PlotRec),
    /* class_initialize      */ NULL,
    /* class_part_initialize     */ NULL,
    /* class_inited          */ FALSE,
    /* initialize        */ Initialize,
    /* initialize_hook       */ NULL,
    /* realize           */ XtInheritRealize,
    /* actions           */ actions,
    /* num_actions       */ XtNumber(actions),
    /* resources         */ resources,
    /* num_resources         */ XtNumber(resources),
    /* xrm_class         */ NULLQUARK,
    /* compress_motion       */ TRUE,
    /* compress_exposure     */ XtExposeCompressMultiple,
    /* compress_enterleave   */ TRUE,
    /* visible_interest      */ FALSE,
    /* destroy           */ Destroy,
    /* resize            */ Resize,
    /* expose            */ Redisplay,
    /* set_values        */ SetValues,
    /* set_values_hook       */ NULL,
    /* set_values_almost     */ XtInheritSetValuesAlmost,
    /* get_values_hook       */ NULL,
    /* accept_focus      */ NULL,
    /* version           */ XtVersion,
    /* callback_private      */ NULL,
    /* tm_table          */ defaultTranslations,
    /* query_geometry        */ QueryGeometry,
    /* display_accelerator       */ XtInheritDisplayAccelerator,
    /* extension                 */ NULL
    },
#if XmVersion>1001
    {  /* Primitive class fields */
    /* border_highlight           _XtInherit,       */
    /* border_unhighlight         _XtInherit,    */
    /* border_highlight   */ (XtWidgetProc)   _XtInherit,       
    /* border_unhighlight */ (XtWidgetProc)   _XtInherit,    
    /* translations       */        XtInheritTranslations,      
    /* arm_and_activate   */        NULL,             
    /* syn resources      */        NULL,           
    /* num_syn_resources  */        0, 
    /* extension          */        NULL,             
    },
#endif
    {
    /* extension          */        0,
    },
};

WidgetClass plotWidgetClass = (WidgetClass) & plotClassRec;

float Screen_scale_factor;
extern struct {float psca_x, psca_y, ixo, iyo, iox, ioy;} a00000_;

void UndrawCellBorder();
void DrawCellBorder();
void UpdateCell();
void Change_initial_options();
void LocateOptionalLabel();
void Init_save_node();    /**** GH, 26/1/18 to avoid conflicting types error ****/
void clear_X(), ChangeCellConfig(), ConfigureCells();
extern void set_box_dimensions(), Init_Data();
extern void ReadLabel();

static void GetBorderGC(w)
Widget w;
{
    PlotWidget cw = (PlotWidget) w;
    XGCValues values;
    XtGCMask mask = GCForeground | GCBackground | GCDashOffset | 
            GCDashList | GCLineStyle;
    extern GC Border_gc;

    /* 
     * Set foreground and background to those of widget
     * border is drawn as dotted line
     */
#if XmVersion<1002
    values.foreground = 1;
    values.background = 0;
#else
    values.foreground = cw->primitive.foreground;
    values.background = cw->core.background_pixel;
#endif
    values.dashes = 1;
    values.dash_offset = 0;
    values.line_style = LineOnOffDash;

    Border_gc = cw->plot.border_gc = XCreateGC(XtDisplay(cw), 
             cw->plot.big_picture, mask, &values);
}

static void GetDrawGC(w)
Widget w;
{
    PlotWidget cw = (PlotWidget) w;
    XGCValues values;
    XtGCMask mask = GCForeground | GCBackground;
    extern GC Draw_gc;

    /* 
     * Set foreground and background
     */
#if XmVersion<1002
    values.foreground = 1;
    values.background = 0;
#else
    values.foreground = cw->primitive.foreground;
    values.background = cw->core.background_pixel;
#endif

    Draw_gc = cw->plot.draw_gc = XCreateGC(XtDisplay(cw), 
             cw->plot.big_picture, mask, &values);
}

static void GetUndrawGC(w)
Widget w;
{
    PlotWidget cw = (PlotWidget) w;
    XGCValues values;
    XtGCMask mask = GCForeground | GCBackground;
    extern GC Undraw_gc;

    /* 
     * Set foreground and background to bg and fg of widget
     */
    values.foreground = cw->core.background_pixel;
    values.background = cw->primitive.foreground;

    Undraw_gc = cw->plot.undraw_gc = XCreateGC(XtDisplay(cw), 
              cw->plot.big_picture, mask, &values);
}

/*
static void GetCopyGC(w)
Widget w;
{
    PlotWidget cw = (PlotWidget) w;
    XGCValues values;
    XtGCMask mask = GCForeground | GCBackground;
    extern GC Copy_gc;

#if XmVersion<1002
    values.foreground = 0;
    values.background = 1;
#else
    values.foreground = cw->primitive.foreground;
    values.background = cw->core.background_pixel;
#endif

    Copy_gc = cw->plot.copy_gc = XtGetGC(cw, mask, &values);
}
*/

/* ARGSUSED */
static void Initialize(treq, tnew, args, num_args)
Widget treq, tnew;
ArgList args;
Cardinal *num_args;
{
    char msg[100];

    unsigned char warn=0;
    PlotWidget new = (PlotWidget) tnew;
    new->plot.cur_x = 0;
    new->plot.cur_y = 0;
    strcpy(msg,"");

    /* 
     *  Check instance values set by resources that may be invalid. 
     */

    if (new->plot.pixmap_width_in_cells==0) /* not set by command line */
         new->plot.pixmap_width_in_cells = (Dimension)Settings.columns;
    if (new->plot.pixmap_width_in_cells < 1) {
        warn=1;
        new->plot.pixmap_width_in_cells = Settings.columns = 1;
        Initial_Settings.columns = Settings.columns;
    }
    else 
        Initial_Settings.columns=Settings.columns=
                                 new->plot.pixmap_width_in_cells;

    if (new->plot.pixmap_height_in_cells==0) /* not set by command line */
         new->plot.pixmap_height_in_cells = (Dimension)Settings.rows;
    if (new->plot.pixmap_height_in_cells < 1) {
        warn=1;
        new->plot.pixmap_height_in_cells = Settings.rows = 1;
        Initial_Settings.rows = Settings.rows;
    }
    else 
        Initial_Settings.rows=Settings.rows=
                                 new->plot.pixmap_height_in_cells;

    if (warn) {
        sprintf(msg,
        "sybil: columns and/or rows invalid or not set (using %d x %d).",
        new->plot.pixmap_width_in_cells,new->plot.pixmap_height_in_cells); 
        XtWarning(msg);
    }

    ConfigureCells(new);

    GetDrawGC((Widget)new);
    GetUndrawGC((Widget)new);
    GetBorderGC((Widget)new);

    InitBigPixmap(new);
}

/* ARGSUSED */
static void Redisplay(w, event)
Widget w;
XExposeEvent *event;
{
    PlotWidget cw = (PlotWidget) w;
    register int x, y;
    unsigned int width, height;
    if (!XtIsRealized(w))
    return;

    if (event) {  /* called from btn-event or expose */
        x = event->x;
        y = event->y; 
        width = event->width;
        height =  event->height;
    } 
    else {        /* called because complete redraw */
        x = 0;
        y = 0; 
        width = cw->core.width;
        height = cw->core.height;
        width = cw->plot.pixmap_width_in_pixels;
        height = cw->plot.pixmap_height_in_pixels;
    }

/*
    XCopyArea(XtDisplay(cw), cw->plot.big_picture, 
                XtWindow(cw), cw->plot.draw_gc, x + 
                cw->plot.cur_x, y + cw->plot.cur_y, 
                width, height, x, y);
*/
    XCopyArea(XtDisplay(cw), cw->plot.big_picture, 
                XtWindow(cw), cw->plot.draw_gc, x, y, 
                width, height, x, y);

    if (Settings.mark_cell) DrawBorder(cw->plot.draw_gc, DRAW,
                                            cw, Plot_info.curr_cell);
}

/* ARGSUSED */
static Boolean SetValues(current, request, new, args, num_args)
Widget current, request, new;
ArgList args;
Cardinal *num_args;
{
    PlotWidget curcw = (PlotWidget) current;
    PlotWidget newcw = (PlotWidget) new;
    Boolean do_redisplay = False;

/*
    if (curcw->primitive.foreground != newcw->primitive.foreground) {
        XtReleaseGC(curcw, curcw->plot.copy_gc);
        GetCopyGC(newcw);
        do_redisplay = True;
    }
*/

    if ((curcw->plot.cur_x != newcw->plot.cur_x) || 
            (curcw->plot.cur_y != newcw->plot.cur_y))
        do_redisplay = True;

    if (curcw->plot.cell_height_in_pixels != 
            newcw->plot.cell_height_in_pixels) {
        ChangeCellSize((Widget)curcw, newcw->plot.cell_width_in_pixels,
                        newcw->plot.cell_height_in_pixels);
        do_redisplay = True;
    }

    if (curcw->plot.pixmap_width_in_cells != 
            newcw->plot.pixmap_width_in_cells)  {
        newcw->plot.pixmap_width_in_cells = 
                curcw->plot.pixmap_width_in_cells;
        XtWarning("Sybil: pixmap_width_in_cells cannot be set by XtSetValues.\n");
    }

    if (curcw->plot.pixmap_height_in_cells != 
            newcw->plot.pixmap_height_in_cells) {
        newcw->plot.pixmap_height_in_cells = 
                curcw->plot.pixmap_height_in_cells;
        XtWarning("Sybil: pixmap_height_in_cells cannot be set by XtSetValues.\n");
    }

    return do_redisplay;
}

void clear_X(w)
Widget w;
{
    Destroy(w);
}

static void Destroy(w)
Widget w;
{
    PlotWidget cw = (PlotWidget) w;
    if (cw->plot.big_picture)
        XFreePixmap(XtDisplay(cw), cw->plot.big_picture);

    if (cw->plot.draw_gc)
        XFreeGC(XtDisplay(cw), cw->plot.draw_gc);

    if (cw->plot.undraw_gc)
        XFreeGC(XtDisplay(cw), cw->plot.undraw_gc);

    if (cw->plot.border_gc)
        XFreeGC(XtDisplay(cw), cw->plot.border_gc);

    /*
     * Free memory allocated to save plots drawn to each cell.
     */
    Clear_all_saved_data(cw);
    /*
     * Free memory allocated with Calloc.  This was done
     * only if application didn't supply cell array.
     */
    if (!cw->plot.user_allocated)
        XtFree((XtPointer)cw->plot.cell);
}

static void Reinitialize( tnew )
Widget tnew;
{
    char msg[100];

    unsigned char warn=0;
    PlotWidget new = (PlotWidget) tnew;
    new->plot.cur_x = 0;
    new->plot.cur_y = 0;

    /* 
     *  Check instance values set by resources that may be invalid. 
     */

    new->plot.pixmap_width_in_cells = (Dimension)Settings.columns;
    new->plot.pixmap_height_in_cells = (Dimension)Settings.rows;

    if (warn) {
        sprintf(msg,
        "sybil: columns and/or rows invalid or not set (using %d x %d).",
        new->plot.pixmap_width_in_cells,new->plot.pixmap_height_in_cells); 
        XtWarning(msg);
    }
    ConfigureCells( new );

    InitBigPixmap(new);
}

void ConfigureCells(new)
PlotWidget new;
{
    int i,j,indx;
    int dumiro,dumlpl,dumils;
    float pixel_hgt,pixel_wdth,width,height;
    float viewable = 0.85;
    float x_offset,y_offset,tmp;
    float pix_per_mm_vert, pix_per_mm_horiz;
    Screen *scr;

    Screen_scale_factor = 1.0;
    scr = XtScreen( new );
    if (Settings.paper_size==US_PAPER) {
        height = USLETTERHEIGHT;
        width = USLETTERWIDTH;
    }
    else {
        height = DEFAULTHEIGHT;
        width = DEFAULTWIDTH;
    }
    pix_per_mm_vert = (float)HeightOfScreen(scr)/HeightMMOfScreen(scr);
    pix_per_mm_horiz = (float)WidthOfScreen(scr)/WidthMMOfScreen(scr);
    if (Settings.orient==LANDSCAPE) {
        tmp = height;
        height = width;
        width = tmp;
    }
    pixel_hgt = pix_per_mm_vert *height;
    if ((int)pixel_hgt>(HeightOfScreen(scr)*viewable)) {
        pixel_hgt = (float)HeightOfScreen(scr)*viewable;
        Screen_scale_factor= (pixel_hgt/pix_per_mm_vert)/height;
    }
    pixel_wdth = pix_per_mm_horiz* width * Screen_scale_factor;
    x_offset = pixel_wdth * Settings.page_xmargin;
    y_offset = pixel_hgt * Settings.page_ymargin;
    /*
     * allow for main title at top of page
     */
    Plot_info.title_offset += y_offset;
    new->plot.cell_height_in_pixels = (Dimension)(
                                    (pixel_hgt-(float)Plot_info.title_offset
                                       - y_offset)/
                                    (float)new->plot.pixmap_height_in_cells);
    new->plot.cell_width_in_pixels = (Dimension)((pixel_wdth - x_offset*2)/
                                    (float)new->plot.pixmap_width_in_cells);
    /*
     * cell dimensions in cm, before scaling to fit on screen
     */
    new->plot.cell_height_in_cm = (float)new->plot.cell_height_in_pixels/
                     ((float)HeightOfScreen(scr)/HeightMMOfScreen(scr)*10)/
                                   Screen_scale_factor;
    new->plot.cell_width_in_cm = (float)new->plot.cell_width_in_pixels/
                     ((float)WidthOfScreen(scr)/WidthMMOfScreen(scr)*10)/
                                   Screen_scale_factor;

    if ((new->plot.cur_x < 0) ||  (new->plot.cur_y < 0)) {
        XtWarning("sybil: cur_x and cur_y must be non-negative (using 0, 0).");
        new->plot.cur_x = 0;
        new->plot.cur_y = 0;
    }

    if (new->plot.cell == NULL) {
        new->plot.cell = (cell_data *)XtCalloc( 
                new->plot.pixmap_width_in_cells * 
                new->plot.pixmap_height_in_cells, sizeof(cell_data));
        new->plot.user_allocated = False;
    }
    else
        new->plot.user_allocated = True;  /* user supplied cell array */

    new->plot.pixmap_width_in_pixels = (Dimension)pixel_wdth;

    new->plot.pixmap_height_in_pixels = (Dimension)pixel_hgt;

    /* 
     * Motif Primitive sets width and height to provide enough room for
     * the highlight and shadow around a widget.  Newsyb
     * doesn't use these features.  A widget that did use these
     * features would *add* its desired dimensions to those set
     * by Primitive.  To use this widget with another widget set, remove
     * the following two lines. 
     */
    new->core.width = 0;
    new->core.height = 0;

    if (new->core.width == 0) {
        if (new->plot.show_all == False)
            new->core.width = (new->plot.pixmap_width_in_pixels 
                    > (Dimension)pixel_wdth) ? pixel_wdth : 
                    (new->plot.pixmap_width_in_pixels);
        else
            new->core.width = new->plot.pixmap_width_in_pixels;
    }

    if (new->core.height == 0) {
        if (new->plot.show_all == False)
            new->core.height = 
                    (new->plot.pixmap_height_in_pixels > 
                    (Dimension)pixel_hgt) ? pixel_hgt : 
                    (new->plot.pixmap_height_in_pixels);
        else
            new->core.height = new->plot.pixmap_height_in_pixels;
    }

    /* index cells 0 -> columns(rows) -1 */
    Plot_info.max_row=new->plot.pixmap_height_in_cells-1;
    Plot_info.max_col=new->plot.pixmap_width_in_cells-1;
    Plot_info.pixels_per_cm=(float)HeightOfScreen(scr)/
                                            HeightMMOfScreen(scr)*10;

dumiro=(int)Settings.orient; dumlpl=11; dumils=2;
xpak_init( &dumiro,&dumlpl,&dumils,XtDisplay(new) );

a00000_.psca_x=(float)WidthOfScreen(scr)/WidthMMOfScreen(scr)*10/*
*Screen_scale_factor*/;
a00000_.psca_y=Plot_info.pixels_per_cm/* *Screen_scale_factor*/;
a00000_.iyo=new->plot.pixmap_height_in_pixels;

/***FIX this when data read from saved file *****/
        /* set values for each cell */
    for (i=0,indx=0;i<new->plot.pixmap_height_in_cells;i++) {
        for (j=0;j<new->plot.pixmap_width_in_cells;j++,indx++) {
            new->plot.cell[indx].drawn=0;
            new->plot.cell[indx].row=i;
            new->plot.cell[indx].col=j;
            new->plot.cell[indx].rect.x=j*new->plot.cell_width_in_pixels
                                            + x_offset;
            new->plot.cell[indx].rect.y=i*new->plot.cell_height_in_pixels +
                                        Plot_info.title_offset;
            new->plot.cell[indx].rect_cm_x = new->plot.cell[indx].rect.x/
                                              a00000_.psca_x;
            new->plot.cell[indx].rect_cm_y = 
               (new->plot.pixmap_height_in_pixels-new->plot.cell[indx].rect.y)/
                                              a00000_.psca_y;
            new->plot.cell[indx].rect.width=new->plot.cell_width_in_pixels-1;
            new->plot.cell[indx].rect.height=new->plot.cell_height_in_pixels-1;
        }
    }


    new->plot.cell_height_in_cm = (float)new->plot.cell_height_in_pixels/
                                   a00000_.psca_y;
    new->plot.cell_width_in_cm = (float)new->plot.cell_width_in_pixels/
                                   a00000_.psca_x;
    set_box_dimensions( Pwindo,
                        new->plot.cell[Plot_info.curr_cell].rect_cm_x,
                        new->plot.cell[Plot_info.curr_cell].rect_cm_y,
                        new->plot.cell_width_in_cm,
                        new->plot.cell_height_in_cm );

    CreateBigPixmap(new);
}

void ChangeCellConfig(cw)
PlotWidget cw;
{

    if (cw->plot.big_picture)
        XFreePixmap(XtDisplay(cw), cw->plot.big_picture);

    /*
     * Free memory allocated to save plots drawn to each cell.
     */
    Clear_all_saved_data(cw);
    /*
     * Free memory allocated with Calloc.  This was done
     * only if application didn't supply cell array.
     */
    if (!cw->plot.user_allocated) {
        XtFree((XtPointer)cw->plot.cell);
        cw->plot.cell = NULL;
    }

    if (Plot_info.title) {
        free(Plot_info.title);
        Plot_info.title = NULL;
    }
    Plot_info.title_offset = DEFAULT_LGE_FONTHGT * 1.5;
    Plot_info.curr_cell = 0;
    Reinitialize((Widget)cw);

    InitBigPixmap(cw);
}

/*
 * clears outline of previous cell, outlines current cell
 */
void ChangeCell( prev,new )
int prev,new;
{
    extern PlotWidget plotArea;
    PlotWidget cw = plotArea;

    if (Settings.mark_cell && XtIsRealized((Widget)cw)) {
        UndrawCellBorder(cw,prev);
        DrawCellBorder(cw,new);
    }
    set_box_dimensions( Pwindo,
                        cw->plot.cell[Plot_info.curr_cell].rect_cm_x,
                        cw->plot.cell[Plot_info.curr_cell].rect_cm_y,
                        cw->plot.cell_width_in_cm,
                        cw->plot.cell_height_in_cm );
    if (Plot_info.inp_file != NULL && prev != new) {
        printf("\nCell number: %i\n",new);
        init_box(Pwindo,&Settings);
    }
}

void DrawCellBorder(w, number)
Widget w;
int number;
{
    PlotWidget cw = (PlotWidget) w;
    /* update Plot_info */
    Plot_info.curr_x = Plot_info.origin_x =
                        (float) cw->plot.cell[number].rect.x;
    Plot_info.curr_y = Plot_info.origin_y =
                        (float) cw->plot.cell[number].rect.y;

    DrawBorder(cw->plot.draw_gc, DRAW, cw, number);
}

void UndrawCellBorder(w, number)
Widget w;
int number;
{
    PlotWidget cw = (PlotWidget) w;
    DrawBorder(cw->plot.undraw_gc, UNDRAW, cw, number);
}

void UpdateCell(w, number, mode)
Widget w;
int number,mode;
{
/*  XExposeEvent fake_event; */
    PlotWidget cw = (PlotWidget) w;

/*
    fake_event.x = cw->plot.cell[number].rect.x+1;
    fake_event.y = cw->plot.cell[number].rect.y+1;
    fake_event.width = cw->plot.cell[number].rect.width-2;
    fake_event.height = cw->plot.cell[number].rect.height-2;
*/
    Redisplay((Widget)cw, NULL);
}

void DeleteCell(number)
int number;
{
	XButtonEvent fake_event;
    extern PlotWidget plotArea;
    PlotWidget cw = plotArea;

    fake_event.x = cw->plot.cell[number].rect.x+1;
    fake_event.y = cw->plot.cell[number].rect.y+1;
    UndrawCell(plotArea,&fake_event,NULL,0);
}

void DeleteAllCells()
{
    extern PlotWidget plotArea;
    PlotWidget cw = plotArea;

    /* doesn't clear the title area */
    Clear_all_saved_data(cw);
    XFillRectangle(XtDisplay(cw), cw->plot.big_picture,
        cw->plot.undraw_gc, 0, Plot_info.title_offset, 
        cw->plot.pixmap_width_in_pixels + 2,
        cw->plot.pixmap_height_in_pixels + 2);
    Redisplay((Widget)cw, NULL);
}

/* BLANK cell*/
static void UndrawCell(w, event, params, num_params)
Widget w;
XEvent *event;
String *params;
Cardinal *num_params;
{
    PlotWidget cw = (PlotWidget) w;
    XExposeEvent fake_event;

    int newx = (cw->plot.cur_x + ((XButtonEvent *)event)->x) / 
             (int)cw->plot.cell_width_in_pixels;
    int newy = (cw->plot.cur_y + ((XButtonEvent *)event)->y) / 
             (int)cw->plot.cell_height_in_pixels;
    if (((XButtonEvent *)event)->y>Plot_info.title_offset) {
        if ((cw->plot.pixmap_width_in_cells*newy+newx)==Plot_info.curr_cell) {
            Clear_saved_data( cw,Plot_info.curr_cell );
            XFillRectangle(XtDisplay(cw),
                        cw->plot.big_picture, cw->plot.undraw_gc,
                        cw->plot.cell_width_in_pixels*newx + 1,
                        cw->plot.cell_height_in_pixels*newy + 1 + 
                            Plot_info.title_offset, 
                        (unsigned int)cw->plot.cell_width_in_pixels - 2, 
                        (unsigned int)cw->plot.cell_height_in_pixels - 2);
            fake_event.x = cw->plot.cell[Plot_info.curr_cell].rect.x-1;
            fake_event.y = cw->plot.cell[Plot_info.curr_cell].rect.y-1;
            if (fake_event.x<0) fake_event.x=0;
            if (fake_event.y<0) fake_event.y=0;
            fake_event.width = cw->plot.cell[Plot_info.curr_cell].rect.width+2;
            fake_event.height = cw->plot.cell[Plot_info.curr_cell].rect.height+2;
            Redisplay((Widget)cw, &fake_event);
        }
    }
    else {
        free(Plot_info.title);
        Plot_info.title = NULL;
        XFillRectangle(XtDisplay(cw), cw->plot.big_picture, cw->plot.undraw_gc,
                        1, 1, 
                        (unsigned int)cw->plot.pixmap_width_in_pixels - 1, 
                        (unsigned int)Plot_info.title_offset - 1);
        fake_event.x = 1;
        fake_event.y = 1;
        fake_event.width = cw->plot.pixmap_width_in_pixels;
        fake_event.height = Plot_info.title_offset;
        Redisplay((Widget)cw, &fake_event);
    }
}

void SetClipToPlot()
{
    extern PlotWidget plotArea;
    PlotWidget cw = plotArea;
    XRectangle clip_rect;
    int x, y, width, height, limitx, limity;

    /*
     * get limits of allowable plotting area
     * mesh limits are given by XCMAX, YCMAX values
     */
    cms_to_pixels(Pwindo[XMIN],Pwindo[YMAX],
                  &(clip_rect.x),&(clip_rect.y),
                  a00000_.psca_x,a00000_.psca_y, (int)a00000_.iyo);
    cms_to_pixels(Pwindo[XMAX],Pwindo[YMIN],&x,&y,
                      a00000_.psca_x,a00000_.psca_y, (int)a00000_.iyo);
    /*
    cms_to_pixels(Pwindo[XMIN],Pwindo[YMAX],
                  &(clip_rect.x),&(clip_rect.y),
                  a00000_.psca_x*Screen_scale_factor,
                  a00000_.psca_y*Screen_scale_factor,
                  (int)a00000_.iyo);
    cms_to_pixels(Pwindo[XMAX],Pwindo[YMIN],&x,&y,
                  a00000_.psca_x*Screen_scale_factor,
                  a00000_.psca_y*Screen_scale_factor,
                  (int)a00000_.iyo);
     */
    height = (uint)(y-clip_rect.y);
    width = (uint)((Pwindo[XCMAX]-Pwindo[XCMIN])/
                            (Pwindo[XCMAX]-Pwindo[XCMIN])*
                            a00000_.psca_x/a00000_.psca_y*height);
    limitx = cw->plot.cell[Plot_info.curr_cell].rect.x +
              cw->plot.cell[Plot_info.curr_cell].rect.width;
    limity = cw->plot.cell[Plot_info.curr_cell].rect.y +
              cw->plot.cell[Plot_info.curr_cell].rect.height;
    if (height>(limity-clip_rect.y)) height = limity-clip_rect.y;
    if (width>(limitx-clip_rect.x)) width = limitx-clip_rect.x;
    clip_rect.height = height;
    clip_rect.width = width;

    XSetClipRectangles(XtDisplay(cw), cw->plot.draw_gc, 0, 0,
                   &clip_rect, 1, YSorted);
}

void UnsetClipToPlot()
{
    extern PlotWidget plotArea;
    PlotWidget cw = plotArea;
    XSetClipMask(XtDisplay(cw), cw->plot.draw_gc, None);
}

#if XY
static void ToggleCell(w, event)
Widget w;
XEvent *event;
{
    PlotWidget cw = (PlotWidget) w;
    static int oldx = -1, oldy = -1;
    GC gc;
    int mode;
    int newx, newy;

    /* This is strictly correct, but doesn't
     * seem to be necessary */
    if (event->type == ButtonPress) {
        newx = (cw->plot.cur_x + ((XButtonEvent *)event)->x) / 
        cw->plot.cell_width_in_pixels;
        newy = (cw->plot.cur_y + ((XButtonEvent *)event)->y) / 
        cw->plot.cell_height_in_pixels;
    }
    else  {
        newx = (cw->plot.cur_x + ((XMotionEvent *)event)->x) / 
        cw->plot.cell_width_in_pixels;
        newy = (cw->plot.cur_y + ((XMotionEvent *)event)->y) / 
        cw->plot.cell_height_in_pixels;
    }


    if ((mode = cw->plot.cell[newx+newy*cw->plot.pixmap_width_in_cells].drawn)
                                                            == DRAWN) {
        gc = cw->plot.undraw_gc;
        mode = UNDRAW;
    }
    else {
        gc = cw->plot.draw_gc;
        mode = DRAW;
    }

    if (oldx != newx || oldy != newy) {
        oldx = newx;
        oldy = newy;
        DrawPixmaps(gc, mode, cw, event);
    } 
}
#endif

static void DrawBorder(gc, mode, w, number)
GC gc;
int mode;
Widget w;
int number;
{
    XGCValues valuesold, valuesnew;
    XtGCMask mask = GCLineStyle|GCForeground;
    int newx, newy;

    PlotWidget cw = (PlotWidget) w;
/*  XExposeEvent fake_event;   */

    /* otherwise, draw or undraw */
    if (mode==UNDRAW) {
        /* erase outline for old cell */
        XGetGCValues(XtDisplay(cw),cw->plot.border_gc,mask,&valuesold);
        valuesnew.line_style = LineSolid;
        valuesnew.foreground = cw->core.background_pixel;
        XChangeGC(XtDisplay(cw),cw->plot.border_gc,mask,&valuesnew);
        XDrawRectangle(XtDisplay(cw), XtWindow(cw),
            cw->plot.border_gc,
            (int)cw->plot.cell[number].rect.x,
            (int)cw->plot.cell[number].rect.y,
            (unsigned int)cw->plot.cell[number].rect.width,
            (unsigned int)cw->plot.cell[number].rect.height);
        XChangeGC(XtDisplay(cw),cw->plot.border_gc,mask,&valuesold);

/*
        fake_event.x = cw->plot.cell[number].rect.x;
        fake_event.y = cw->plot.cell[number].rect.y;
        fake_event.width = cw->plot.cell[number].rect.width+1;
        fake_event.height = cw->plot.cell[number].rect.height+1;
        Redisplay(cw, &fake_event);
*/
    }
    else {
        /* draw outline for new cell */
        XDrawRectangle(XtDisplay(cw), XtWindow(cw),
            cw->plot.border_gc,
            (int)cw->plot.cell[number].rect.x,
            (int)cw->plot.cell[number].rect.y,
            (unsigned int)cw->plot.cell[number].rect.width,
            (unsigned int)cw->plot.cell[number].rect.height);

    }


    newx = cw->plot.cell[number].col;
    newy = cw->plot.cell[number].row;
    cw->plot.cell[number].drawn = mode;
    info.mode = mode;
    info.newx = newx;
    info.newy = newy;
}

/*
 * draws label at x,y (user co-ords)
 * y_offset<0:  draw label below this y value
 * y_offset=0:  draw label centred on y value
 * y_offset>0:  draw label above this y value
 * x_offset>0:  draw label at this x value
 * x_offset=0:  draw label centred on this x value
 * x_offset<0:  draw label ending on this x value
 */
void drawlabel_(x, y, x_offset, y_offset, mode, label, len)
float *x,*y;
int *mode,*x_offset,*y_offset,*len;
char *label;
{
    XCharStruct overall;
    extern PlotWidget plotArea;

    int width=0,ascent,descent,direction,state;
    float xtmp, ytmp;
    PlotWidget cw = plotArea;

/*
    if (mode==UNDRAW)
        XFillRectangle(XtDisplay(cw), cw->plot.big_picture, gc,
                        cw->plot.cell_width_in_pixels*newx + 1,
                        cw->plot.cell_height_in_pixels*newy + 1,
                        (unsigned int)cw->plot.cell_width_in_pixels - 2,
                       (unsigned int)cw->plot.cell_height_in_pixels -
2);
    else {
    }
*/

    width = cw->plot.pixmap_width_in_pixels;
    descent = ascent = 20;
    if (Plot_info.font!=NULL) {
        XSetFont(XtDisplay(cw),cw->plot.draw_gc,Plot_info.font->fid);
        XTextExtents(Plot_info.font,
                        label,strlen(label),
                        &direction,&ascent,&descent,&overall);
        width = overall.width;
    }
    state=PEN_UP;
    xtmp= *x; ytmp= *y;
    plotu_(&xtmp,&ytmp,&state);
    xtmp=a00000_.iox; ytmp=a00000_.ioy;
    if (*y_offset==J_CENTRE) ytmp += ((float)(ascent+descent)/2.0);
    else if (*y_offset==J_TOP) ytmp += ((float)ascent+1.0);
    else if (*y_offset==J_BASE) ytmp -= ((float)descent+1);
    if (*x_offset==J_RIGHT) xtmp -= ((float)width+1);
    else if (*x_offset==J_CENTRE) {
        width /= 2;
        xtmp -= (float)width;
    }
    /*
     * keep the label on the page
     */
    if ((int)xtmp+width > (int)(cw->plot.pixmap_width_in_pixels))
        xtmp = (float)(cw->plot.pixmap_width_in_pixels-width-1);
    if ((int)ytmp > (int)(cw->plot.pixmap_height_in_pixels))
        ytmp = (float)(cw->plot.pixmap_height_in_pixels-1);

    if (Settings.text_bg)
        XDrawImageString(XtDisplay(cw),cw->plot.big_picture,
                        cw->plot.draw_gc,
                        (int)xtmp, (int)ytmp,
                        label,strlen(label));
    else
        XDrawString(XtDisplay(cw),cw->plot.big_picture,
                        cw->plot.draw_gc,
                        (int)xtmp, (int)ytmp,
                        label,strlen(label));

    a00000_.iox = xtmp+width;
    a00000_.ioy = ytmp;
}
 
/*
 * label max and min on the colour bar
 */
void labelcolourbar_(float* min, float* max, float* scale, float* barlen, 
                     float* xlevel, float* ylevel, int* vertical)
{
    char number[MAX_NUM_LEN];
    int current_pen;
    int format,len,mode,x_offset,y_offset;
    float scale_min, scale_max,x,y;
    XFontStruct *curr_font;

    current_pen=FOREGROUND;
    setpencolor_(&current_pen);
    curr_font = Plot_info.font;
    Plot_info.font = Plot_info.dflt_font;
    scale_min = *min * *scale;
    scale_max = *max * *scale;
    format = MAX(fabs((double)scale_min),fabs((double)scale_max))
                        >=10.0 ? 60 : 63;
    formatnumber_(number,&scale_min,&format);
    if (*vertical) {
        x_offset = J_CENTRE;
        y_offset = J_TOP;
    }
    else {
        x_offset = J_RIGHT;
        y_offset = J_TOP;
    }
    mode = 1; len = format/10;
    drawlabel_(xlevel,ylevel,&x_offset,&y_offset,&mode,number,&len);
    formatnumber_(number,&scale_max,&format);
    if (*vertical) {
        x_offset = J_CENTRE;
        y_offset = J_BASE;
        x = *xlevel;
        y = *ylevel + *barlen;
    }
    else {
        x_offset = J_LEFT;
        y_offset = J_TOP;
        x = *xlevel + *barlen;
        y = *ylevel;
    }
    mode = 1; len = format/10;
    drawlabel_(&x,&y,&x_offset,&y_offset,&mode,number,&len);
    Plot_info.font = curr_font;
}

void drawautolabels_(label,labellen,value,format,pwindo,current_pen)
char *label;
int *format, *labellen;
float *value,*pwindo;
int *current_pen;
{
    char number[MAX_NUM_LEN];
    XFontStruct *curr_font;
    int old_pen,x_offset,y_offset,mode,len;

    strcpy(number,"");
    old_pen = *current_pen;
    *current_pen=FOREGROUND;
    setpencolor_(current_pen);
    curr_font = Plot_info.font;
    Plot_info.font = Plot_info.dflt_font;
    x_offset = J_LEFT; y_offset = J_TOP; mode = 1;
    drawlabel_(&pwindo[XCMIN],&pwindo[YCMIN],&x_offset,&y_offset,&mode,
                                                       label,labellen);
    if (*format) {
        formatnumber_(number,value,format);
        /* draw label at centre
        x_offset = J_LEFT; y_offset = J_TOP; mode = 1; len = strlen(number);
        drawlabel_((pwindo[YCMAX]-pwindo[YCMIN])/2,pwindo[XCMIN],
                                   &x_offset,&y_offset,&mode,number,&len);
        */
        x_offset = J_RIGHT; y_offset = J_TOP; mode = 1; len = strlen(number);
        drawlabel_(&pwindo[XCMAX],&pwindo[YCMIN],
                                  &x_offset,&y_offset,&mode,number,&len);
    }
    *current_pen=old_pen;
    setpencolor_(current_pen);
    Plot_info.font = curr_font;
}

int drawlabelcm_(x,y,offset_x,offset_y,mode,label,nchar)
char *label;
int *offset_x,*offset_y,*mode,*nchar;
float *x,*y;
{
    XButtonEvent fake_event;
    extern PlotWidget plotArea;

    cms_to_pixels(*x,*y,&(fake_event.x),&(fake_event.y),
                      a00000_.psca_x*Screen_scale_factor,
                      a00000_.psca_y*Screen_scale_factor,
                      (int)a00000_.iyo);
    LocateOptionalLabel(plotArea,&fake_event);
    return(0);
}

void DrawTitle(char* title)
{
    int width=0,ascent,descent,direction;
    extern PlotWidget plotArea;
    PlotWidget cw = plotArea;
    XFontStruct *tmp;
    XExposeEvent fake_event;
    XCharStruct overall;
/*
 * replace with call to OptionalLabel
 */

    tmp = Plot_info.font;
    Plot_info.font = Plot_info.large_font;
    XSetFont(XtDisplay(cw),cw->plot.draw_gc,Plot_info.font->fid);
    XTextExtents(Plot_info.font,
                        title,strlen(title),
                        &direction,&ascent,&descent,&overall);
    width = overall.width;
    XDrawImageString(XtDisplay(cw),cw->plot.big_picture,cw->plot.draw_gc,
                        Plot_info.title_offset,
                        Plot_info.title_offset*2/3,
                        title,strlen(title));
    fake_event.x = Plot_info.title_offset;
    fake_event.width = width;
    fake_event.height = Plot_info.title_offset;
    fake_event.y = 0;
    Plot_info.font = tmp;
    Redisplay((Widget)cw, &fake_event);
}

static void OptionalLabel(w, event, params, num_params)
Widget w;
XButtonEvent *event;
String *params;
Cardinal *num_params;
{
    extern Widget LabelPrompt;

    if (!(LabelPrompt && XtIsManaged(LabelPrompt))) return;
    ReadLabel(LabelPrompt, Plot_info.curr_label);
    LocateOptionalLabel(w,event);
}

void LocateOptionalLabel(w, event)
Widget w;
XButtonEvent *event;
{
    XCharStruct overall;

    int width=0,ascent,descent,direction;
    unsigned long size;
    PlotWidget cw = (PlotWidget) w;
    int newx,newy, err;
/*
    int newx = (cw->plot.cur_x + event->x) / 
             (int)(cw->plot.cell_width_in_pixels);
    int newy = (cw->plot.cur_y + event->y) / 
             (int)(cw->plot.cell_height_in_pixels);
*/
    XExposeEvent fake_event;


    /* this fn should call other fns if allowing draw or undraw */
    int mode=DRAW;

    if (strlen(Plot_info.curr_label)==0) return;
    width = cw->plot.pixmap_width_in_pixels;
    descent = ascent = Plot_info.title_offset;
    if (Plot_info.font!=NULL) {
        XSetFont(XtDisplay(cw),cw->plot.draw_gc,Plot_info.font->fid);
        XTextExtents(Plot_info.font,
                        Plot_info.curr_label,strlen(Plot_info.curr_label),
                        &direction,&ascent,&descent,&overall);
        width = overall.width;
    }
    if (Settings.text_bg)
        XDrawImageString(XtDisplay(cw),
                        cw->plot.big_picture,cw->plot.draw_gc,
                        event->x, event->y,
                        Plot_info.curr_label,strlen(Plot_info.curr_label));
    else
        XDrawString(XtDisplay(cw),
                        cw->plot.big_picture,cw->plot.draw_gc,
                        event->x, event->y,
                        Plot_info.curr_label,strlen(Plot_info.curr_label));
    /* store the label info */
    XGetFontProperty(Plot_info.font,XA_POINT_SIZE,&size);
    store_label(w,Plot_info.curr_label,Plot_info.fontname,size,
                       event->x,event->y,width,ascent+descent,Settings.fg);

    if (!(err = FindCell(event->x,event->y,&newx,&newy,cw))) {
        cw->plot.cell[newx + newy*cw->plot.pixmap_width_in_cells].drawn
                                                                 = mode;
        info.mode = mode;
        info.newx = newx;
        info.newy = newy;
    }

    fake_event.x = event->x - cw->plot.cur_x;
    fake_event.y = event->y - cw->plot.cur_y - ascent;
    fake_event.width = width;
    fake_event.height = ascent + descent;
    Redisplay(w, &fake_event);
/*
    XtCallCallbacks(cw, XtNtoggleCallback, &info);
*/
}

/*
 * returns 0 if cell found
 * else 1
 */
int FindCell(x,y,thisx,thisy,cw)
int x,y;
int *thisx,*thisy;
PlotWidget cw;
{
    int i,j,found,foundx,foundy,indx;
    int wdth = cw->plot.cell_width_in_pixels;
    int hgt = cw->plot.cell_height_in_pixels;

    for (i=0,found=0,indx=0;i<cw->plot.pixmap_height_in_cells && !found;i++) {
        for (j=0;j<cw->plot.pixmap_width_in_cells && !found;j++,indx++) {
            if ((foundx = (x>=cw->plot.cell[indx].rect.x &&
                           x<=cw->plot.cell[indx].rect.x+wdth))
               && (foundy = (y>=cw->plot.cell[indx].rect.y &&
                           y<=cw->plot.cell[indx].rect.y+hgt))) {
                *thisx = j;
                *thisy = i;
                found = 1;
            }
        }
    }
    return((found)?0:1);
}

/*
 * calculates x,y in cms
 */
void pixels_to_cms(x,y,cm_x,cm_y,cnvrt_x,cnvrt_y,hgt)
int x,y,hgt;
float cnvrt_x,cnvrt_y,*cm_x,*cm_y;
{

/*
    *cm_x = (float)x / cnvrt_x / Screen_scale_factor;
    *cm_y = (float)(hgt-y) / cnvrt_y / Screen_scale_factor;
*/
    *cm_x = (float)x / cnvrt_x;
    *cm_y = (float)(hgt-y) / cnvrt_y;
}

/*
 * calculates x,y (cms) in pixels
 */
void cms_to_pixels(x,y,pix_x,pix_y,cnvrt_x,cnvrt_y,hgt)
int hgt,*pix_x,*pix_y;
float x,y,cnvrt_x,cnvrt_y;
{

/*
    *pix_x = (int)(x * cnvrt_x * Screen_scale_factor);
    *pix_y = hgt - (int)(y * cnvrt_y * Screen_scale_factor);
*/
    *pix_x = (int)(x * cnvrt_x);
    *pix_y = hgt - (int)(y * cnvrt_y);
}

/*
 * Save label data -
 * font and ptsize
 * x,y as fraction of page width and height
 * colour as index
 */
int store_label( w,label,fontname,size,x,y,wdth,hgt,fg_colour )
Widget w;
char *label,*fontname;
int fg_colour;
int x, y, wdth, hgt;
unsigned long size;
{
    float xsave, ysave;
    PlotWidget cw = (PlotWidget) w;
    label_node *node,**p;

    pixels_to_cms(x,y,&xsave,&ysave,
                   a00000_.psca_x*Screen_scale_factor,
                   a00000_.psca_y*Screen_scale_factor,
                   cw->plot.pixmap_height_in_pixels);
/* UPD label coords in fraction of page dim
    if (Settings.orient==PORTRAIT) {
        xsave /= DEFAULTWIDTH/10.0;
        ysave /= DEFAULTHEIGHT/10.0;
    }
    else if (Settings.orient==LANDSCAPE) {
        xsave /= DEFAULTHEIGHT/10.0;
        ysave /= DEFAULTWIDTH/10.0;
    }
*/
    if  ((node = (label_node *)malloc(sizeof(label_node)))!=NULL) {
        strcpy(node->label.label_text,label);
        strcpy(node->label.fontname,fontname);
        node->label.font_pointsize = (int)size/10;
        node->label.fg = fg_colour;
        node->label.x = xsave;
        node->label.y = ysave;
        node->label.box.x = (short)x;
        node->label.box.y = (short)(y - hgt);
        node->label.box.width = (unsigned short)(wdth);
        node->label.box.height = (unsigned short)(hgt);
        node->label.text_bg = Settings.text_bg;
        node->next = NULL;
    }
    else return(MALLOC_ERR);

    /* append to the plot's list */
    p = &(cw->plot.labels); 
    while (*p) p = &((*p)->next);
    *p = node;
    return(0);
}

void CreateBigPixmap(w)
Widget w;
{
    PlotWidget cw = (PlotWidget) w;
    extern Drawable DrawArea;

    cw->plot.big_picture = XCreatePixmap(XtDisplay(cw),
            RootWindow(XtDisplay(cw), DefaultScreen(XtDisplay(cw))),
            cw->plot.pixmap_width_in_pixels + 2, 
            cw->plot.pixmap_height_in_pixels + 2,
            XDefaultDepth(XtDisplay(cw),DefaultScreen(XtDisplay(cw))));
    DrawArea = (Drawable)( cw->plot.big_picture );
}

void InitBigPixmap(w)
Widget w;
{
    PlotWidget cw = (PlotWidget) w;

    XFillRectangle(XtDisplay(cw), cw->plot.big_picture,
        cw->plot.undraw_gc, 0, 0, 
        cw->plot.pixmap_width_in_pixels + 2,
        cw->plot.pixmap_height_in_pixels + 2);
}

#if XY
/* A Public function, not static */
char * PlotGetArray(w, width_in_cells, height_in_cells)
Widget w;
int *width_in_cells, *height_in_cells;
{
    PlotWidget cw = (PlotWidget) w;

    *width_in_cells = cw->plot.pixmap_width_in_cells;
    *height_in_cells = cw->plot.pixmap_height_in_cells;
    return (cw->plot.cell);
}
#endif

/* ARGSUSED */
static void Resize(w)
Widget w;
{
#if XY
    PlotWidget cw = (PlotWidget) w;
    /* resize does nothing unless new size is bigger than entire pixmap */
    if ((cw->core.width > cw->plot.pixmap_width_in_pixels) &&
            (cw->core.height > cw->plot.pixmap_height_in_pixels)) {
        /* 
         * Calculate the maximum cell size that will allow the
         * entire bitmap to be displayed.
         */
        Dimension w_temp_cell_width_in_pixels, h_temp_cell_height_in_pixels;
        Dimension new_cell_width_in_pixels;
        Dimension new_cell_height_in_pixels;
    
        w_temp_cell_width_in_pixels = cw->core.width / cw->plot.pixmap_width_in_cells;
        h_temp_cell_height_in_pixels = cw->core.height / cw->plot.pixmap_height_in_cells;
    
/**************fix this *********************/
        if (w_temp_cell_width_in_pixels < h_temp_cell_height_in_pixels)
            new_cell_size_in_pixels = w_temp_cell_size_in_pixels;
        else
            new_cell_size_in_pixels = h_temp_cell_size_in_pixels;
    
        /* if size change mandates a new pixmap, make one */
        if (new_cell_size_in_pixels != cw->plot.cell_size_in_pixels)
            ChangeCellSize(cw, new_cell_width_in_pixels,
                                    new_cell_width_in_pixels);
    }
#endif
    Redisplay(w,NULL);
}

static void ChangeCellSize(w, new_cell_width, new_cell_height, reorient)
Widget w;
int new_cell_width;
int new_cell_height;
int reorient;
{
#if XY
    PlotWidget cw = (PlotWidget) w;
    int x, y;

    cw->plot.cell_width_in_pixels = new_cell_width;
    cw->plot.cell_height_in_pixels = new_cell_height;

    /* recalculate variables based on cell size */
    cw->plot.pixmap_width_in_pixels = 
             cw->plot.pixmap_width_in_cells * 
             cw->plot.cell_width_in_pixels;

    cw->plot.pixmap_height_in_pixels = 
             cw->plot.pixmap_height_in_cells * 
             cw->plot.cell_height_in_pixels;
    
    /* destroy old and create new pixmap of correct size */
    XFreePixmap(XtDisplay(cw), cw->plot.big_picture);
    CreateBigPixmap(cw);
    
    /* draw lines into new pixmap */
    DrawIntoBigPixmap(cw);
    
    /* draw current cell array into pixmap */
    for (x = 0; x < cw->plot.pixmap_width_in_cells; x++) {
        for (y = 0; y < cw->plot.pixmap_height_in_cells; y++) {
            if (cw->plot.cell[x + (y * cw->plot.pixmap_width_in_cells)] == DRAWN)
                DoCell(cw, x, y, cw->plot.draw_gc);
            else
                DoCell(cw, x, y, cw->plot.undraw_gc);
        }
    }
#endif
}
/*   routine not used for now; commented out

static void DoCell(w, x, y, gc)
Widget w;
int x, y;
GC gc;
{
    PlotWidget cw = (PlotWidget) w;      */
        /* otherwise, draw or undraw */
/*  XFillRectangle(XtDisplay(cw), cw->plot.big_picture, gc,
             cw->plot.cell_width_in_pixels * x + 2,
             cw->plot.cell_height_in_pixels * y + 2,
             (unsigned int)cw->plot.cell_height_in_pixels - 3,
             (unsigned int)cw->plot.cell_height_in_pixels - 3);

}    */

void DrawPreviewLine( w,x1,y1,x2,y2 )
Widget w;
int x1,y1,x2,y2;
{
    extern PlotWidget plotArea;
    extern GC Draw_gc;

    XDrawLine(XtDisplay(w), XtWindow(plotArea), Draw_gc,
                x1,y1,x2,y2);
}

void DrawRectangle(w,xp1,yp1,xp2,yp2)
Widget w;
int xp1,yp1,xp2,yp2;
{
    extern GC Draw_gc;
    PlotWidget cw = (PlotWidget) w;

    XDrawRectangle(XtDisplay(cw),cw->plot.big_picture,Draw_gc,
                                    xp1, yp1, xp2-xp1, yp2-yp1);
}
 
static XtGeometryResult QueryGeometry(w, proposed, answer)
Widget w;
XtWidgetGeometry *proposed, *answer;
{
    PlotWidget cw = (PlotWidget) w;
    int width, height;

        /***FIX these assignments if not keeping showall=True ***/
    width = cw->plot.pixmap_width_in_pixels;
    height = cw->plot.pixmap_height_in_pixels;

    /* set fields we care about */
    answer->request_mode = CWWidth | CWHeight;

    /* initial width and height */
    if (cw->plot.show_all == True)
        answer->width = cw->plot.pixmap_width_in_pixels;
    else
        answer->width = (cw->plot.pixmap_width_in_pixels > 
                (Dimension)width) ? (Dimension)width : 
                cw->plot.pixmap_width_in_pixels;

    if (cw->plot.show_all == True)
        answer->height = cw->plot.pixmap_height_in_pixels;
    else
        answer->height = (cw->plot.pixmap_height_in_pixels > 
                (Dimension)height) ? (Dimension)height : 
                cw->plot.pixmap_height_in_pixels;

    if (  ((proposed->request_mode & (CWWidth | CWHeight))
            == (CWWidth | CWHeight)) &&
            proposed->width == answer->width &&
            proposed->height == answer->height)
        return XtGeometryYes;
    else if (answer->width == cw->core.width &&
            answer->height == cw->core.height)
        return XtGeometryNo;
    else
        return XtGeometryAlmost;
}

int Write_changed_option(fp,option)
FILE *fp;
opt_data *option;
{
    char tmpstr[MAX_OPTION_NAME*2];

    if (id_match(option_terms,option->id,tmpstr)) {
        fprintf(fp,"%s%s",tmpstr,"=");
        switch(option->type) {
        case SYB_INT  : if (option->id==O_CNTRPLT) {
                            strcpy(tmpstr,"");
                            if (option->int_val&SHADE)
                                strcpy(tmpstr,plot_types[0].name);
                            if (option->int_val&LINES) {
                                strcat(tmpstr,"+");
                                strcat(tmpstr,plot_types[1].name);
                            }
                            fprintf(fp,"%s ",tmpstr);
                        }
                        else if (option->id==O_LINESTYLE) {
                            id_match(linestyle_terms,option->int_val,tmpstr);
                            fprintf(fp,"%s ",tmpstr);
                        }
                        else if (option->id==O_COLBAR) {
                            id_match(colourbar_terms,option->int_val,tmpstr);
                            fprintf(fp,"%s ",tmpstr);
                        }
                        else if (option->id==O_TIC) {
                            id_match(ticmark_terms,option->int_val,tmpstr);
                            fprintf(fp,"%s ",tmpstr);
                        }
                        else fprintf(fp,"%3d ",option->int_val);
                        break;
        case SYB_FLOAT: fprintf(fp,"%6.2f ",option->float_val);
                        break;
        case SYB_BOOL : fprintf(fp,"%2d ",(int)(option->bool_val));
                        break;
        default:        /*return(TYPE_ERR);*/
                        break;
        }
    }
    return(0);
}

int Save_option(XtPointer val, int id, int type)
{
    opt_node *node, **p_o;
    save_node *p;
    extern PlotWidget plotArea;
    PlotWidget cw = plotArea;

    p = cw->plot.plots;
    if (p==NULL) Change_initial_options(id,val);
    else {
        if  ((node = (opt_node *)malloc(sizeof(opt_node)))!=NULL) {
            node->opt.id = id;
            node->opt.type = type;
            node->opt.int_val = 0;
            node->opt.float_val = 0;
            node->opt.bool_val = 0;
            switch(type) {
            case SYB_INT  : node->opt.int_val = *(int *)val;
                            break;
            case SYB_FLOAT: node->opt.float_val = *(float *)val;
                            break;
            case SYB_BOOL : node->opt.bool_val = *(unsigned char *)val;
                            break;
            default:        XtFree((XtPointer)node);
                            return(TYPE_ERR);
                            break;
            }
        
            node->next = NULL;
        }
        else return(MALLOC_ERR);
        /* append to changed option list in last plot */
        while (p->next) p = p->next;
        p_o = &(p->changed_opts); 
        while (*p_o) p_o = &((*p_o)->next);
        *p_o = node;
    }
    return(0);
    
}

void Change_initial_options(id,val)
int id;
XtPointer val;
{
    switch(id) {
        case O_NX3: Initial_Settings.plot_opts.nx3 = *(int *)val;
                 break;
        case O_MP: Initial_Settings.plot_opts.mp = *(int *)val;
                 break;
        case O_MPE: Initial_Settings.plot_opts.mpe = *(int *)val;
                 break;
        case O_NP: Initial_Settings.plot_opts.np = *(int *)val;
                 break;
        case O_CNTRPLT: Initial_Settings.plot_opts.contour_plot = *(int *)val;
                 break;
        case O_PRFLPTS: Initial_Settings.plot_opts.profile_pts = *(int *)val;
                 break;
        case O_FLIP: Initial_Settings.plot_opts.flip = *(int *)val;
                 break;
        case O_DBLE : Initial_Settings.plot_opts.dble = *(int *)val;
                 break;
        case O_LABEL: Initial_Settings.plot_opts.label = *(int *)val;
                 break;
        case O_STIP: Initial_Settings.plot_opts.stipple = *(int *)val;
                 break;
        case O_SOLNROT: Initial_Settings.plot_opts.solution_rot = *(int *)val;
                 break;
        case O_FOREGRND: Initial_Settings.fg = *(int *)val;
                 break;
        case O_LINESTYLE: Initial_Settings.linestyle = *(int *)val;
                 break;
        case O_LINEWDTH: Initial_Settings.linewidth = *(float *)val;
                 break;
        case O_COLBAR: Initial_Settings.plot_opts.colour_bar = *(int *)val;
                 break;
        case O_TIC: Initial_Settings.plot_opts.ticmark = *(int *)val;
                 break;
        case O_ZOOM: Initial_Settings.zoom = *(float *)val;
                 break;
        case O_CLIPTOCELL: Initial_Settings.clip_to_cell = *(unsigned char *)val;
                 break;
        case O_XCENTRE: Initial_Settings.xcentre = *(float *)val;
                 break;
        case O_YCENTRE: Initial_Settings.ycentre = *(float *)val;
                 break;
        case O_VERBOSE: Initial_Settings.verbose = *(int *)val;
                 break;
        case O_MARKCELL: Initial_Settings.mark_cell = *(unsigned char *)val;
                 break;
        case O_FILEPATH: Initial_Settings.filepath = *(unsigned char *)val;
                 break;
    }
}

int Save_plot_data( cntr_vals,prfl_vals )
float *cntr_vals, *prfl_vals;
{
    int j;
    char *name = NULL;
    save_node *node,**p;
    extern PlotWidget plotArea;
    PlotWidget cw = plotArea;
    if  ((node = (save_node *)malloc(sizeof(save_node)))!=NULL) {
        Init_save_node( node );
        if (Plot_info.inp_file) {
          /*
           * if not saving full path, find solution directory
           * FD.sols in the path or revert to full path
           */
           if (!Settings.filepath) {
              name = strstr(Plot_info.inp_file->fname,"FD.sols");
	   }
           if (name==NULL) {
               name = Plot_info.inp_file->fname;
           }
          /*strcpy(node->data.filename,Plot_info.inp_file->fname);*/
          strcpy(node->data.filename,name);
          node->data.record_num = Plot_info.inp_file->rec_curr;
          node->data.reference_num = Plot_info.inp_file->ref_curr;
        }
        if (Plot_info.elle_file)
          strcpy(node->data.ellefilename,Plot_info.elle_file->fname);
        strcpy(node->data.plot_param,Plot_info.variable);
        strcpy(node->data.dflt_label1,Plot_info.dflt_label1);
        strcpy(node->data.dflt_label2,Plot_info.dflt_label2);
        node->data.cell_num = Plot_info.curr_cell;
        node->data.plot_type = Plot_info.plot_type;
        node->data.plot_description = Plot_info.plot_description;
        node->data.option_num = Plot_info.var_num;
        node->data.stipple_type = Settings.plot_opts.stipple;
        node->data.colour_bar = Settings.plot_opts.colour_bar;
        if (cntr_vals != NULL)
            for (j=0;j<MAXCNTRVALS;j++) 
                node->data.contour_vals[j] = cntr_vals[j];
        if (prfl_vals != NULL)
            for (j=0;j<MAXPRFLVALS;j++) 
                node->data.profile_vals[j] = prfl_vals[j];
    }
    else return(MALLOC_ERR);

    /* append to list of plots */
    p = &(cw->plot.plots); 
    while (*p) p = &((*p)->next);
/* if inserting, use extra condition in while then
    node->next = *p;
*/
    *p = node;
    return(0);
}

void Init_save_node( node )
save_node *node;
{
    int j;

    strcpy(node->data.filename,"");
    strcpy(node->data.ellefilename,"");
    strcpy(node->data.plot_param,"");
    strcpy(node->data.dflt_label1,"");
    strcpy(node->data.dflt_label2,"");
    for (j=0;j<MAXCNTRVALS;j++) 
        node->data.contour_vals[j] = 0.0;
    for (j=0;j<MAXCNTRVALS;j++) 
        node->data.profile_vals[j] = 0.0;
    node->changed_opts = NULL;
    node->next = NULL;
}

int Write_plot_data( w,fp )
Widget w;
FILE *fp;
{
    unsigned char last_bg;
    char commentsfile[SYB_FILENAME_MAX+1],*dir;
    char tmpstr[MAX_OPTION_NAME];
    int i, j, k, err=0;
    save_node *s_node;
    opt_node *o_node;
    opt_data opt;
    label_node *l_node;
    FILE * fp_intro;
    PlotWidget cw = (PlotWidget)w;

      /* log the options */
    fprintf(fp,"%-10s","Options");
    if ((err=write_options(fp,&Initial_Settings))) return(err);
    fprintf(fp,"\n");
      /* log the title */
    if (Plot_info.title && strlen(Plot_info.title))
        fprintf(fp,"%-10s%s%s%s\n","Title","\"",Plot_info.title,"\"");
    s_node=cw->plot.plots;
    if (s_node) {
      /* log the plots */
        while (s_node) {
            fprintf(fp,"%-10s%d\t%d\n","Locate",
                        cw->plot.cell[s_node->data.cell_num].row,
                        cw->plot.cell[s_node->data.cell_num].col);
            if (strlen(s_node->data.filename)) {
                fprintf(fp,"%-10s%s\n","File",s_node->data.filename);
                fprintf(fp,"%-10s%d\t%d\n","Data_In",s_node->data.record_num,
                                            s_node->data.reference_num);
            }
            switch(s_node->data.plot_type) {
            case ARROWS: fprintf(fp,"%-10s","Arrow");
                         fprintf(fp,"%s\n",s_node->data.plot_param);
                         break;
            case CNTRS:  fprintf(fp,"%-10s","Contour");
                         fprintf(fp,"%s\n",s_node->data.plot_param);
                         if (s_node->data.plot_description&SHADE) {
                           fprintf(fp,"%10s%s"," ","cntr.colour");
                           for (j=0;j<3;j++) 
                             fprintf(fp," %s=%f",colour_vals[j],
                                                 s_node->data.contour_vals[j]);
                           if (id_match(colourbar_terms,s_node->data.colour_bar,
                                                 tmpstr)) 
                             fprintf(fp," %s=%s",colour_vals[j],tmpstr);
                           fprintf(fp,"\n");
                         }
                         if (s_node->data.plot_description&LINES) {
                           fprintf(fp,"%10s%s"," ","cntr.line");
                           for (j=0;j<3;j++) 
                             fprintf(fp," %s=%f",line_vals[j],
                                               s_node->data.contour_vals[j+3]);
                           fprintf(fp,"\n");
                         }
                    /* this line is not used 
                         fprintf(fp,"%10s%s"," ","cntr.limits");
                         for (j=0;j<2;j++) 
                           fprintf(fp," %s=%f",prfl_limits[j],
                                                 s_node->data.contour_vals[j]);
                         fprintf(fp,"\n");
                     */
                         if (Settings.plot_opts.stipple) {
                           fprintf(fp,"%10s%s %s=%d"," ","cntr.stipple",
                              stipple_vals[0],
                              s_node->data.stipple_type);
                           fprintf(fp,"\n");
                         }
                         break;
            case PRFLMRK:fprintf(fp,"%-10s","Profile");
                         fprintf(fp,"%s\n",profileplots[3].name);
                         fprintf(fp,"%10s%s"," ","prfl.pts");
                         for (j=0,k=PRFL_X1;j<4;j++,k++) 
                           fprintf(fp," %s=%g",prfl_pts[j],
                                           s_node->data.profile_vals[k]);
                         fprintf(fp,"\n");
                         break;
            case PRFL1D: fprintf(fp,"%-10s","Profile");
                         fprintf(fp,"%s%s%s\n",profileplots[0].name,".",
                                           s_node->data.plot_param);
                         fprintf(fp,"%10s%s"," ","prfl.line");
                         for (j=0;j<4;j++) 
                           fprintf(fp," %s=%g",line_terms[j],
                                             s_node->data.profile_vals[j]);
                         fprintf(fp,"\n");
                         fprintf(fp,"%10s%s"," ","prfl.pts");
                         for (j=0,k=PRFL_X1;j<4;j++,k++) 
                           fprintf(fp," %s=%g",prfl_pts[j],
                                           s_node->data.profile_vals[k]);
                         fprintf(fp,"\n");
                         break;
            case PRFL2D: fprintf(fp,"%-10s","Profile");
                   /*    fprintf(stdout,"(*s_node).data.plot_description = %i\n",
                                         (*s_node).data.plot_description); */
                         if (s_node->data.plot_description&DIM2_X)
                           fprintf(fp,"%s%s%s\n",profileplots[1].name,".",
                                           s_node->data.plot_param);
                         if (s_node->data.plot_description&DIM2_Y)
                           fprintf(fp,"%s%s%s\n",profileplots[2].name,".",
                                           s_node->data.plot_param);
                         fprintf(fp,"%10s%s"," ","prfl.line");
                         for (j=0;j<3;j++) 
                           fprintf(fp," %s=%g",line_terms[j],
                                             s_node->data.profile_vals[j]);
                         fprintf(fp,"\n");
                         fprintf(fp,"%10s%s"," ","prfl.pts");
                         for (j=0,k=PRFL_X1;j<4;j++,k++) 
                             fprintf(fp," %s=%g",prfl_pts[j],
                                             s_node->data.profile_vals[k]);
                         fprintf(fp,"\n");
                         fprintf(fp,"%10s%s"," ","prfl.limits");
                         fprintf(fp," %s=%g",prfl_limits[0],
                                        s_node->data.profile_vals[LWR_LIM]);
                         fprintf(fp," %s=%g",prfl_limits[1],
                                        s_node->data.profile_vals[UPR_LIM]);
                         fprintf(fp,"\n");
                         break;
            case XYPLOT: fprintf(fp,"%-10s","XYPlot");
                         fprintf(fp,"%s",s_node->data.plot_param);
                         if (name_match(s_node->data.plot_param,xyplots)
                                     == ELLE)
                             fprintf(fp," %s",s_node->data.ellefilename);
                         fprintf(fp,"\n");
                         break;
            default:     break;
            }
            if (s_node->changed_opts) {
                fprintf(fp,"%-10s","Options");
                o_node = s_node->changed_opts;
                do {
                    Write_changed_option(fp,&(o_node->opt));
                    o_node = o_node->next;
                } while (o_node); 
                fprintf(fp,"\n");
            }
            s_node = s_node->next;
        }
    }
    l_node = cw->plot.labels;
    if (l_node)
    last_bg = Initial_Settings.text_bg;
    while (l_node) {
        if (l_node->label.text_bg!=last_bg) {
            opt.id = O_TEXTBG;
            opt.type = SYB_BOOL;
            opt.bool_val = l_node->label.text_bg;
            fprintf(fp,"%-10s","Options");
            Write_changed_option(fp,&opt);
            fprintf(fp,"\n");
            last_bg = l_node->label.text_bg;
        }
        fprintf(fp,"%-10s%c%s%c","Label",LABEL_MRK,
                                      l_node->label.label_text,LABEL_MRK);
        fprintf(fp," %s", l_node->label.fontname);
        fprintf(fp," %d", l_node->label.font_pointsize);
        fprintf(fp," %d", l_node->label.fg);
        fprintf(fp," %f", l_node->label.x);
        fprintf(fp," %f", l_node->label.y);
        fprintf(fp,"\n");
        l_node = l_node->next;
    }
      /* log the explanatory comments */
    if ((dir = getenv("BASILPATH"))!=NULL) {
        strcpy(commentsfile,dir);
        strcat(commentsfile,"/docs/");
        strcat(commentsfile,CMNTS_FILE);
        if ((fp_intro=fopen(commentsfile,"r"))!=NULL) {
            while ((i=fgetc(fp_intro))!=EOF) fputc(i,fp);
            fclose(fp_intro);
            fflush( fp );
            fprintf(fp,"\n");
        }
    }
    return(0);
}

/*
 * clears plots for all cells
 */
int Clear_all_saved_data( cw )
PlotWidget cw;
{
    save_node *s_node, *s_tmp;
    opt_node *o_node, *o_tmp;
    label_node *l_node, *l_tmp;

    s_node = cw->plot.plots;
    while (s_node) {
        o_node=s_node->changed_opts;
        if (o_node) {
            while (o_node) {
                o_tmp = o_node->next;
                free(o_node);
                o_node = o_tmp;
            }
        }
        s_tmp = s_node->next;
        free(s_node);
        s_node = s_tmp;
    }
    cw->plot.plots = NULL;
    l_node = cw->plot.labels;
    while (l_node) {
        l_tmp = l_node->next;
        free(l_node);
        l_node = l_tmp;
    }
    cw->plot.labels = NULL;
    return(0);
}

/*
 * clears plots for a particular cell
 * does not clear labels
 */
int Clear_saved_data( cw,cell_num )
int cell_num;
PlotWidget cw;
{
    XtPointer val;
    save_node **s_addr, *s_tmp;
    opt_node *o_node, *o_tmp, **o_addr;

    s_addr = &(cw->plot.plots);
    while (*s_addr) {
        if ((*s_addr)->data.cell_num==cell_num) {
            o_node=(*s_addr)->changed_opts;
            if (o_node) {
                if (s_addr!=&(cw->plot.plots)) { /* append to previous */
                    while (*o_addr) o_addr = &((*o_addr)->next);
                    *o_addr = o_node;
                }
                else do { /* first plot - change initial settings */
                    switch(o_node->opt.type) {
                    case SYB_INT  : val=(XtPointer)&(o_node->opt.int_val);
                                    break;
                    case SYB_FLOAT: val=(XtPointer)&(o_node->opt.float_val);
                                    break;
                    case SYB_BOOL : val=(XtPointer)&(o_node->opt.bool_val);
                                    break;
                    }
                    Change_initial_options(o_node->opt.id,val);
                    o_tmp = o_node;
                    o_node = o_node->next;
                    free(o_tmp);
                } while (o_node!=NULL);
            }
            s_tmp = *s_addr;
            *s_addr = (*s_addr)->next;
            free(s_tmp);
        }
        else {
            o_addr = &((*s_addr)->changed_opts);
            s_addr = &((*s_addr)->next);
        }
    }
    return(0);
}

/*
 * If id matches a member of set,
 * copy the corresponding name to the string name and return 1
 * else name has length zero and return 0
 */
int id_match(set,id,name)
valid_terms set[];
int id;
char *name;
{
    int found=0;
    int i=0;

    strcpy(name,"");
    while (set[i].name && !found) {
        if (id==set[i].id) found=1;
        else i++;
    }
    if (found) strcpy(name,set[i].name);
    return(found);
}

int write_options(fp,opts)
FILE *fp;
input_options *opts;
{
    char tmpstr[MAX_OPTION_NAME*2];
    int i=0, charcnt = 0, len, err=0;

    strcpy(tmpstr,"");
    while (option_terms[i].name != NULL) {
        if (option_terms[i].id<O_RESCALE ) {
            len = strlen(option_terms[i].name);
            if ((charcnt+len+12) > 60) {
                fprintf(fp,"\n%-10s"," ");
                charcnt = 0;
            }
            fprintf(fp,"%s%s",option_terms[i].name,"=");
            charcnt += len+1;
        }
        switch(option_terms[i].id) {
        case  O_ORIENT:
                 if (opts->orient==LANDSCAPE) fprintf(fp,"%s ","LANDSCAPE");
                 else fprintf(fp,"%s ","PORTRAIT");
                 charcnt += 10;
                 break;
        case O_PAPERSZ:
                 switch(opts->paper_size) {
                 case A4_PAPER:fprintf(fp,"%s ",paper_terms[0].name);
                               charcnt += strlen(paper_terms[0].name)+1;
                               break;
                 case US_PAPER: fprintf(fp,"%s ",paper_terms[1].name);
                               charcnt += strlen(paper_terms[1].name)+1;
                               break;
                 }
                 break;
        case  O_MARKCELL: if (opts->mark_cell) fprintf(fp,"%1d ",1);
                 else fprintf(fp,"%1d ",0);
                 charcnt += 2;
                 break;
        case  O_ROWS: fprintf(fp,"%2d ",opts->rows);
                 charcnt += 3;
                 break;
        case  O_COLS: fprintf(fp,"%2d ",opts->columns);
                 charcnt += 3;
                 break;
        case  O_FONTHGT: fprintf(fp,"%2d ",opts->dfltfonthgt);
                 charcnt += 3;
                 break;
        case  O_XMARG: fprintf(fp,"%6.2f ",opts->xmargin);
                 charcnt += 3;
                 break;
        case  O_YMARG: fprintf(fp,"%6.2f ",opts->ymargin);
                 charcnt += 7;
                 break;
        case  O_PG_XMARG: fprintf(fp,"%6.2f ",opts->page_xmargin);
                 charcnt += 3;
                 break;
        case  O_PG_YMARG: fprintf(fp,"%6.2f ",opts->page_ymargin);
                 charcnt += 7;
                 break;
        case  O_FILEPATH: fprintf(fp,"%2d ",opts->filepath);
                 charcnt += 7;
                 break;
        case  O_ZOOM: fprintf(fp,"%6.2f ",opts->zoom);
                 charcnt += 7;
                 break;
        case  O_XCENTRE: fprintf(fp,"%6.2f ",opts->xcentre);
                 charcnt += 7;
                 break;
        case  O_YCENTRE: fprintf(fp,"%6.2f ",opts->ycentre);
                 charcnt += 7;
                 break;
        case  O_LINEWDTH: fprintf(fp,"%6.2f ",opts->linewidth);
                 charcnt += 7;
                 break;
        case  O_PGSCL: fprintf(fp,"%6.2f ",opts->page_scale);
                 charcnt += 7;
                 break;
        case  O_TL0: fprintf(fp,"%6.2f ",opts->physical.tl0);
                 charcnt += 7;
                 break;
        case  O_RHOC: fprintf(fp,"%6.3f ",opts->physical.rhoc);
                 charcnt += 7;
                 break;
        case O_RHOM: fprintf(fp,"%6.3f ",opts->physical.rhom);
                 charcnt += 7;
                 break;
        case O_Z0: fprintf(fp,"%6.3f ",opts->physical.z0);
                 charcnt += 7;
                 break;
        case O_ELEV0: fprintf(fp,"%6.3f ",opts->physical.elev0);
                 charcnt += 7;
                 break;
        case O_TMAX: fprintf(fp,"%6.3f ",opts->physical.tmax);
                 charcnt += 7;
                 break;
        case O_BGAM0: fprintf(fp,"%6.3f ",opts->physical.bgam0);
                 charcnt += 7;
                 break;
        case O_BGAM1: fprintf(fp,"%6.3f ",opts->physical.bgam1);
                 charcnt += 7;
                 break;
        case O_NDIVR: fprintf(fp,"%5.1f ",opts->contour.ndivr);
                 charcnt += 6;
                 break;
        case O_S0: fprintf(fp,"%5.1f ",opts->contour.s0);
                 charcnt += 6;
                 break;
        case O_DELS: fprintf(fp,"%6.3f ",opts->contour.dels);
                 charcnt += 7;
                 break;
        case O_DELOM: fprintf(fp,"%6.3f ",opts->contour.delom);
                 charcnt += 7;
                 break;
        case O_NX3: fprintf(fp,"%5d ",opts->plot_opts.nx3);
                 charcnt += 6;
                 break;
        case O_MP: fprintf(fp,"%2d ",opts->plot_opts.mp);
                 charcnt += 3;
                 break;
        case O_MPE: fprintf(fp,"%2d ",opts->plot_opts.mpe);
                 charcnt += 3;
                 break;
        case O_NP: fprintf(fp,"%2d ",opts->plot_opts.np);
                 charcnt += 3;
                 break;
        case O_CNTRPLT:
                 if (opts->plot_opts.contour_plot&SHADE)
                      strcpy(tmpstr,plot_types[0].name);
                 if (opts->plot_opts.contour_plot&LINES) {
                     strcat(tmpstr,"+");
                     strcat(tmpstr,plot_types[1].name);
                 }
                 fprintf(fp,"%s ",tmpstr);
                 charcnt += strlen(tmpstr)+1;
                 break;
        case O_PRFLPTS: fprintf(fp,"%4d ",opts->plot_opts.profile_pts);
                 charcnt += 5;
                 break;
        case O_FLIP: fprintf(fp,"%2d ",opts->plot_opts.flip);
                 charcnt += 3;
                 break;
        case O_DBLE: fprintf(fp,"%2d ",opts->plot_opts.dble);
                 charcnt += 3;
                 break;
        case O_LABEL: fprintf(fp,"%2d ",opts->plot_opts.label);
                 charcnt += 3;
                 break;
        case O_STIP: fprintf(fp,"%2d ",opts->plot_opts.stipple);
                 charcnt += 3;
                 break;
        case O_SOLNROT: fprintf(fp,"%2d ",opts->plot_opts.solution_rot);
                 charcnt += 3;
                 break;
        case O_FOREGRND: fprintf(fp,"%3d ",opts->fg);
                 charcnt += 4;
                 break;
        case O_LINESTYLE:
                 switch(opts->linestyle) {
                 default:
                 case SOLID: fprintf(fp,"%s ",linestyle_terms[0].name);
                             charcnt += strlen(linestyle_terms[0].name)+1;
                             break;
                 case DASH:  fprintf(fp,"%s ",linestyle_terms[1].name);
                             charcnt += strlen(linestyle_terms[1].name)+1;
                             break;
                 case DOT :  fprintf(fp,"%s ",linestyle_terms[2].name);
                             charcnt += strlen(linestyle_terms[2].name)+1;
                             break;
                 }
                 break;
        case O_COLMAP:
                 switch(opts->colourmap) {
                 case GREY_MAP:fprintf(fp,"%s ",colmap_terms[0].name);
                               charcnt += strlen(colmap_terms[0].name)+1;
                               break;
                 case STD_MAP: fprintf(fp,"%s ",colmap_terms[1].name);
                               charcnt += strlen(colmap_terms[1].name)+1;
                               break;
                 case ABS_MAP: fprintf(fp,"%s ",colmap_terms[2].name);
                               charcnt += strlen(colmap_terms[2].name)+1;
                               break;
                 }
                 break;
        case O_COLBAR:
                 switch(opts->plot_opts.colour_bar) {
                 case CB_VERT:fprintf(fp,"%s ",colourbar_terms[0].name);
                               charcnt += strlen(colourbar_terms[0].name)+1;
                               break;
                 case CB_HORIZ: fprintf(fp,"%s ",colourbar_terms[1].name);
                               charcnt += strlen(colourbar_terms[1].name)+1;
                               break;
                 case CB_NONE: fprintf(fp,"%s ",colourbar_terms[2].name);
                               charcnt += strlen(colourbar_terms[2].name)+1;
                               break;
                 }
                 break;
        case  O_VISCMIN: fprintf(fp,"%6.2f ",opts->viscmin);
                 charcnt += 7;
                 break;
        case  O_VISCMAX: fprintf(fp,"%6.2f ",opts->viscmax);
                 charcnt += 7;
                 break;
        case  O_SEMIN: fprintf(fp,"%6.2f ",opts->semin);
                 charcnt += 7;
                 break;
        case  O_SEMAX: fprintf(fp,"%6.2f ",opts->semax);
                 charcnt += 7;
                 break;
        case O_CLIPTOCELL: fprintf(fp,"%2d ",opts->clip_to_cell);
                 charcnt += 3;
                 break;
        case O_RESCALE: /*fprintf(fp,"%2d ",opts->rescale);
                 charcnt += 3;*/
                 break;
        case O_VERBOSE: fprintf(fp,"%2d ",opts->verbose);
                 charcnt += 3;
                 break;
        case O_TEXTBG:
                 break;
        case O_ELLE:
                 if (opts->elle) {
                     len = strlen(option_terms[i].name);
                     if ((charcnt+len+12) > 60) {
                         fprintf(fp,"\n%-10s"," ");
                         charcnt = 0;
                     }
                     fprintf(fp,"%s%s",option_terms[i].name,"=");
                     charcnt += len+1;
                     fprintf(fp,"%2d ",(int)opts->elle);
                 }
                 break;
        case O_TIC:
                 switch(opts->plot_opts.ticmark) {
                 case J_CENTRE: fprintf(fp,"%s ",ticmark_terms[0].name);
                               charcnt += strlen(ticmark_terms[0].name)+1;
                               break;
                 case J_BASE  : fprintf(fp,"%s ",ticmark_terms[1].name);
                               charcnt += strlen(ticmark_terms[1].name)+1;
                               break;
                 case J_TOP   :fprintf(fp,"%s ",ticmark_terms[2].name);
                               charcnt += strlen(ticmark_terms[2].name)+1;
                               break;
                 case TIC_NONE: fprintf(fp,"%s ",ticmark_terms[3].name);
                               charcnt += strlen(ticmark_terms[3].name)+1;
                               break;
                 }
                 break;
        case O_AXISXMIN:/* fprintf(fp,"%6.3f ",Pwindo[ULEFT]);
                 charcnt += 7;*/
                 break;
        case O_AXISXMAX: /*fprintf(fp,"%6.3f ",Pwindo[URGHT]);
                 charcnt += 7;*/
                 break;
        case O_AXISYMIN: /*fprintf(fp,"%6.3f ",Pwindo[UBASE]);
                 charcnt += 7;*/
                 break;
        case O_AXISYMAX: /*fprintf(fp,"%6.3f ",Pwindo[UTOP]);
                 charcnt += 7;*/
                 break;
        default: err = VAR_ERR;
                 break;
      }
      i++;
    }
    return(err);
}
