/*--------------------------------------------------------------------
 *    Basil / Sybil:   types.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#ifndef _newsyb_h
#define _newsyb_h
#include "defines.h"
#include "cmndefs.h"
#include "filedat.h"
#include "plotdefs.h"

struct phys_params {
    float tl0;
    float rhoc,rhom;
    float z0;
    float tmax;
    float elev0;
    float bgam0;
    float bgam1;
};
    
struct cntr_params {
    float ndivr;
    float s0;
    float dels;
    float delom;
};

struct plot_params {
    int nx3;
    int mp;
    int mpe;
    int np;
    int contour_plot;
    int colour_bar;
    int ticmark;
    int profile_pts;
    int flip;
    int dble;
    int label;
    int stipple;
    int solution_rot;
};

typedef struct {
    unsigned char orient;
    unsigned char mark_cell;
    unsigned char rescale;
    unsigned char clip_to_cell;
    unsigned char elle;
    unsigned char text_bg;
    unsigned char filepath;
    int verbose;
    int paper_size;
    int rows;
    int columns;
    int dfltfonthgt;
    int colourmap;
    int fg; /* index in Colours[] */
    int linestyle; /* -1=solid, 1=dashed */
    float xmargin,ymargin;
    float page_xmargin,page_ymargin;
    float page_scale;
    float prev_zoom;
    float zoom;
    float xcentre;
    float ycentre;
    float linewidth;
    float viscmin,viscmax;
    float semin,semax;
    struct phys_params physical;
    struct cntr_params contour;
    struct plot_params plot_opts;
} input_options;

typedef struct {
    unsigned char display_log;
    char *title;
    char fontname[MAXNAME];
    char variable[MAXNAME*2];
    char dflt_label1[MAXNAME];
    char dflt_label2[MAXNAME];
    char curr_label[MAX_LABEL_LEN];
    int curr_cell;
    int max_row, max_col;
    int title_offset;
    int max_colours;
    int plot_type; /* arrows,cntrs,profile,xy_plot */
    int plot_description; /* lines,shading */
    int var_num;
    float curr_x, curr_y;
    float origin_x, origin_y;
    float pixels_per_cm;
    fontdat font;
    fontdat dflt_font;
    fontdat large_font;
    FILE *log_file_ptr;
    file_data *elle_file;
    file_data *inp_file;
} plot_data;

/*
 * equivalent to Xlib.h XRectangle
 */
typedef struct {
    short x, y;
    unsigned short width, height;
} rectangle;

/* data for each label drawn in a cell */
typedef struct {
    char label_text[MAX_LABEL_LEN];
    char fontname[MAXNAME];
    int font_pointsize; /* GetFontProperty XA_POINT_SIZE */
    int fg;
    float x,y; /* in page cms, 0,0 at bottom left */
    rectangle box; /* label bounding box, in win pixels */
    unsigned char text_bg;
} label_data;

struct l_node{
    label_data label;
    struct l_node *next;
} ;
typedef struct l_node label_node;

typedef struct {
    int id;
    int type;
    int int_val;
    float float_val;
    unsigned char bool_val;
} opt_data;

struct o_node{
    opt_data opt;
    struct o_node *next;
} ;
typedef struct o_node opt_node;

typedef struct {
    char filename[SYB_FILENAME_MAX+1];
    char ellefilename[SYB_FILENAME_MAX+1];
    char plot_param[MAXNAME];
    char dflt_label1[MAXNAME];
    char dflt_label2[MAXNAME];
    int cell_num;
    int record_num;
    int reference_num;
    int plot_type;
    int plot_description;
    int option_num;
    int stipple_type;
    int colour_bar;
    float contour_vals[MAXCNTRVALS];
    float profile_vals[MAXPRFLVALS];
} save_data;

struct s_node{
    save_data data;
    opt_node *changed_opts;
    struct s_node *next;
};
typedef struct s_node save_node;

/* data for each cell */
typedef struct {
    char drawn;
    int row,col;
    rectangle rect;
    float rect_cm_x;
    float rect_cm_y;
    save_node *plots;
} cell_data;
 
#endif

extern colourdat Colours[];
extern plot_data Plot_info;
extern input_options Settings, Initial_Settings;
