/*--------------------------------------------------------------------
 *    Basil / Sybil:   page.c  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "types.h"
#include "globals.h"
#include "mainps.h"
#include "pref77.h"
#include "errnum.h"
#include "error.h"
#include "routines.h"

#define CM_PER_INCH    2.54
#define PTS_PER_CM    72.0/CM_PER_INCH 
/* 0.25in margin around A4 page not used by printer */
#define POINT_OFFSET    72/4
#define CM_OFFSET       CM_PER_INCH/4
#define PS_COL_MIN      500
#define PS_COL_MAX     1000
/*
#define PS_COL_MIN      8
#define PS_COL_MAX     71
*/

char LogFile[SYB_FILENAME_MAX+1],PsFile[SYB_FILENAME_MAX+1];
int File_no;
page_data Page;
void ChangeCell(),Syntax();
extern void set_box_dimensions(),clear_arrays();

void SetClipToPlot();
void UnsetClipToPlot();

int Set_Colours(red,green,blue,irgbv,steps)
int *red,*green,*blue,*irgbv,steps;
{
    int ncolors, *limits;
/*
    Colour_range[0] = PS_COL_MIN;
    Colour_range[1] = PS_COL_MAX;
*/
    if ((limits = (int *)malloc(steps*sizeof(int)))==NULL) return(1);
    ncolors = PLOT_COLOURS;
    Colour_range[0] = USEFL_SIZE;
    Colour_range[1] = MAX_COLOURS-1;
    Arrow_colours[0] = 4; /* blue */
    Arrow_colours[1] = 2; /* red */
    Plot_info.max_colours = MAX_COLOURS-1;
    createcolormap_(red,green,blue,irgbv,limits,&steps,&Colour_range[0],
                            &ncolors );
    free(limits);
    return(0);
}
int Init_App(int argc, char** argv, FILE** optsfp)
{
    char input_str[SYB_FILENAME_MAX+1],*tmpstr;
    int i,j,indx,err=0;
    int hgt,len;
    int iro,init=1,ils=1;
    float pt_wdth, pt_hgt, width, height, factor, tmp;
    float x_offset, y_offset;

    LogFile[0] = PsFile[0] = '\0';

    strncpy(LogFile,DEFAULT_LOG_FILE,SYB_FILENAME_MAX);
    if (argc>1) {
        if (ParseOptions(argc,argv,LogFile,PsFile)) {
            Syntax(argc,argv);
            return(1);
        }
    }
        /* open files */
    i = strlen(PsFile);
    if (i==0) {  /* no output filename specified */
       if ((tmpstr = strrchr(LogFile,'.'))!=NULL) {
         while (&LogFile[i]!=tmpstr) {PsFile[i] = LogFile[i]; i++; }
         PsFile[i] = '\0';
       }
       else strcpy(PsFile,LogFile);
       strcat(PsFile,DEFAULT_OUTPUT_EXT);
       i = strlen(PsFile);
    }
    fileopen_(PsFile,&i,&File_no);
    if (!File_no) error_msg(OPEN_ERR,PsFile);
    if ((*optsfp = fopen(LogFile,"r"))==NULL) error_msg(OPEN_ERR,LogFile);

        /* read the initial necessary options */
    if ((err = initial_options(*optsfp,&Settings,input_str)))
        error_msg(err,input_str);

        /* initialise the postscript file */
    iro = (int)Settings.orient;
    hplots_(&init,&iro,&File_no,&ils);
        /* initialise the font */
    hgt = DEFAULT_FONTHGT;
    len = strlen("Helvetica");
    setfont_("Helvetica",&len,&hgt);

        /* set up the page */
    if (Settings.paper_size==US_PAPER) {
        height = USLETTERHEIGHT;
        width = USLETTERWIDTH;
    }
    else {
        height = DEFAULTHEIGHT;
        width = DEFAULTWIDTH;
    }
    if (Settings.orient!=PORTRAIT) {
        tmp = height;
        height = width;
        width = tmp;
    }
    pt_wdth = width/10.0*PTS_PER_CM;
    pt_hgt = height/10.0*PTS_PER_CM;
    Page.width_in_pts = pt_wdth;
    Page.height_in_pts = pt_hgt;
    Page.width_in_cells = Settings.columns;
    Page.height_in_cells = Settings.rows;
    x_offset = Page.width_in_pts * Settings.page_xmargin;
    y_offset = Page.height_in_pts * Settings.page_ymargin;
   /*
    * allow for main title at top of page
    */
    Plot_info.title_offset += y_offset;
    Page.cell_height_in_pts = (int)( (pt_hgt-(float)Plot_info.title_offset
                                    - y_offset)/
                                    (float)Page.height_in_cells);
    Page.cell_width_in_pts = (int)((pt_wdth - x_offset*2)/
                                    (float)Page.width_in_cells);
    /*
     * cell dimensions in cm
     */
    factor = PTS_PER_CM;
    Page.cell_height_in_cm = (float)Page.cell_height_in_pts/factor;
    Page.cell_width_in_cm = (float)Page.cell_width_in_pts/factor;
 
    if (Page.cell == NULL)
        if ((Page.cell = (cell_data *)calloc(
                Page.width_in_cells *
                Page.height_in_cells, sizeof(cell_data))) == NULL)
            return(MALLOC_ERR);

    for (i=0,indx=0;i<Page.height_in_cells;i++) {
        for (j=0;j<Page.width_in_cells;j++,indx++) {
            Page.cell[indx].drawn=0;
            Page.cell[indx].row=i;
            Page.cell[indx].col=j;
            Page.cell[indx].rect.x=j*Page.cell_width_in_pts + x_offset;
            Page.cell[indx].rect.y=i*Page.cell_height_in_pts +
                                       Plot_info.title_offset;
            Page.cell[indx].rect_cm_x=(float)Page.cell[indx].rect.x/factor;
            Page.cell[indx].rect_cm_y=(float)(Page.height_in_pts -
                                          Page.cell[indx].rect.y)/factor;
            Page.cell[indx].rect.width=Page.cell_width_in_pts-1;
            Page.cell[indx].rect.height=Page.cell_height_in_pts-1;
        }
    }
    indx = Plot_info.curr_cell;
    set_box_dimensions(Pwindo,
                       Page.cell[indx].rect_cm_x,
                       Page.cell[indx].rect_cm_y,
                       Page.cell_width_in_cm,Page.cell_height_in_cm);
    return(0); 
}

/* 
 * Dummy function - not used in ps version
 */
int GetUserVals(float* user_vals, int flag, unsigned char prfl_flag, 
                int* plot, char* labels)
{
    return(0);
}

void DrawTitle(char* title)
{
    int x_offset,y_offset,mode,len,hgt;
    float x,y,factor;

    factor = PTS_PER_CM;
    hgt = DEFAULT_LGE_FONTHGT;
    len = strlen("Helvetica");
    setfont_("Helvetica",&len,&hgt);
    x = (float)(Settings.xmargin * Page.cell_width_in_pts +
               Settings.page_xmargin * Page.width_in_pts)/factor;
    y = (float)(Page.height_in_pts-Plot_info.title_offset)/factor;
/*
    x = XLMARG * Page.cell_width_in_cm;
    y = Page.height_in_cells * Page.cell_height_in_cm + 0.1;
*/
    x_offset = J_LEFT; y_offset = J_BASE; mode = 1; len = strlen(title);
    drawlabelcm_(&x,&y,&x_offset,&y_offset,&mode, title,&len);
    hgt = DEFAULT_FONTHGT;
    len = strlen("Helvetica");
    setfont_("Helvetica",&len,&hgt);
}
 
/*
 * label max and min on the colour bar
 */
void labelcolourbar_(float* min, float* max, float* scale, float* barlen,
                     float* xlevel, float* ylevel, int* vertical)
{
    char number[MAX_NUM_LEN];
    int i;
    int hgt,current_pen;
    int format,len,mode,x_offset,y_offset;
    float scale_min, scale_max,x,y;

    current_pen=FOREGROUND;
    setpencolor_(&current_pen);
    hgt = DEFAULT_FONTHGT;
    len = strlen("Helvetica");
    setfont_("Helvetica",&len,&hgt);
    scale_min = *min * *scale;
    scale_max = *max * *scale;
    format = MAX(fabs((double)scale_min),fabs((double)scale_max))
                        >10.0 ? 60 : 62;
    formatnumber_(number,&scale_min,&format);
    if (*vertical) {
        x_offset = J_CENTRE;
        y_offset = J_TOP;
    }
    else {
        x_offset = J_RIGHT;
        y_offset = J_CENTRE;
    }
    mode = 1; len = strlen(number);
    i=0;
    while(number[i]==' ' && i < MAX_NUM_LEN-1) i++;
    while(number[i]!=' ' && i < MAX_NUM_LEN-1) i++;
    len=i;
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
        y_offset = J_CENTRE;
        x = *xlevel + *barlen;
        y = *ylevel;
    }
    mode = 1; len = strlen(number);
    i=0;
    while(number[i]==' ' && i < MAX_NUM_LEN-1) i++;
    while(number[i]!=' ' && i < MAX_NUM_LEN-1) i++;
    len=i;
    drawlabel_(&x,&y,&x_offset,&y_offset,&mode,number,&len);
}

void ChangeCell(prev,new)
int prev,new;
{
    set_box_dimensions( Pwindo,
                        Page.cell[new].rect_cm_x,Page.cell[new].rect_cm_y,
                        Page.cell_width_in_cm,Page.cell_height_in_cm);
    if (Plot_info.inp_file != NULL && prev != new) {
        printf("\nCell number: %i\n",new);
        init_box(Pwindo,&Settings);
    }
}

void SetClipToPlot()
{
    float x,y,wdth,hgt;

    x = Pwindo[XCMIN];
    y = Pwindo[YCMIN];
    wdth = Pwindo[XCMAX]-Pwindo[XCMIN];
    hgt = Pwindo[YCMAX]-Pwindo[YCMIN];
    clipregionu_(&x,&y, &wdth, &hgt);
}

void UnsetClipToPlot()
{
    unsetclipregionu_();
}

void Syntax(argc,argv)
int argc;
char **argv;
{
    fprintf(stderr, "options:\n" );
    fprintf(stderr, "   -i   name of log file\n" );
    fprintf(stderr, "   -o   name for Postscript output file\n" );
    fprintf(stderr, "   -h   Print this message\n" );
    fprintf(stderr, "Example:  sybilps -i logfile -o psfile\n" );
}

void warning_msg(err_num,str)
int err_num;
char *str;
{
    char message[SYB_FILENAME_MAX];

    strcpy(message,"");
    strncpy(message,str,sizeof(message)-21);
    CreateErrorMessage(message,err_num);
    fprintf(stderr,"%s\n",message);
}

void Exit_App(err)
int err;
{
    int init=0,iro=0,ils=1;

    if (!err) hplots_(&init,&iro,&File_no,&ils);
    if (Plot_info.inp_file!=NULL)
        if (Plot_info.inp_file->fp!=NULL) fclose(Plot_info.inp_file->fp);
    clear_arrays();
    exit(err);
}

int Run_App(FILE* logfp)
{
    int err=0, dum=1;

    setpencolor_(&Settings.fg);
    setlinewidth_(&Settings.linewidth);
    dashln_(&Settings.linestyle,&dum);
    err = process_log_file(logfp,&Settings);
    Exit_App(err);
    return(err);//will not reach here if Exit_App exits
}

int Save_plot_data( cntr_vals,prfl_vals )
float *cntr_vals, *prfl_vals;
{
    return(0);
}
