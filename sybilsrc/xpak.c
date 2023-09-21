/*--------------------------------------------------------------------
 *    Basil / Sybil:   xpak.c   1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

/*      adapted from routines obtained from B.L.N. Kennett at ANU     */

#include <stdio.h>
#include <stdlib.h>
#include <Xm/Xm.h>
#include "cmndefs.h"
#include "plot.h"
#include <math.h>

#define XY      0
/* Program macros and parameter constants */

#define ABS(x) (((x) < 0) ? -(x) :(x))


void xpak_init(int *iro,int *lpl,int *ils,Display *dsply);
int CreateColormap( Display *dsply,Screen *scrn,
                    int *red,int *green,int *blue,int *irgbv,int steps );
void setpencolor_(int *ipen);
void  setlinewidth_(float *wdth);
void  origin_(float *xmin,float *px1,float *ymin,float *py1);
void  scale_(float *xmin,float *xmax,float *px1,float *px2,
             float *ymin,float *ymax,float *py1,float *py2);
void  plot_(float *x,float *y,int *i);
void  plotu_(float *x,float *y,int *i);
void symbol_(float *x,float *y,float *size,char *iword,float *angl,
             int *nchar);
void number_(float *x,float *y,float *size,float *rn,float *angl,int *nsf);
void  drawrectangleu_(float *x1,float *y1,float *x2,float *y2);
void  fillpoly_(float *xf,float *yf,int *nf,int *colour,int *outline);
int dashln_( int *dash,int *pattern );
int convertutoxpts_(float *val,float *pts);
int convertutoypts_(float *val,float *pts);
int setfont_(char *fontname,int *nchar,int *fonthgt);
int drawticmarks_(float *size,float *stop,float *val,float *start,
                  float *interval,int *justify,int *horizontal);
int drawcircle_(float *xmin,float *ymin,float *diameter,
                int *colour,int *outline);
void  Translate_u_to_win(float xu,float yu, int *xwin, int *ywin);
void formatnumber_(char *str,float *num,int *format);
extern void FindFont();

/* Global variables */

Drawable        DrawArea;
Display        *Dsply=NULL;
GC              Draw_gc, Undraw_gc, Border_gc;
Colormap        cmap;
XColor          Colours[MAX_COLOURS];

#if XY
#define STEPS        9
#define INTERVALS    (STEPS-1)
#define MIN_REQ      2*INTERVALS+USEFL_SIZE
/*
 * greyscale
 */
static int       Red0[] = { 0,255, 35, 70,105,140,175,210};
static int     Green0[] = { 0,255, 35, 70,105,140,175,210};
static int      Blue0[] = { 0,255, 35, 70,105,140,175,210};

/*
 * Useful colours
 * black,white,red,grn,blue,cyn,yellow,magenta
 */
static int       Red1[] = {0,255,255,  0,  0,255,255,  0};
static int     Green1[] = {0,255,  0,255,  0,  0,255,255};
static int      Blue1[] = {0,255,  0,  0,255,255,  0,255};

static int      Rgb0[STEPS*3]={    32767,32767,32767,
                                   36863,36863,36863,
                                   40959,40959,40959,
                                   45055,45055,45055,
                                   49151,49151,49151,
                                   53247,53247,53247,
                                   57343,57343,57343,
                                   61439,61439,61439,
                                   65535,65535,65535
                                };
static int      Rgb1[STEPS*3]={        0,    0,65535,
                                   32767,32767,65535,
                                   32767,65535,65535,
                                       0,65535,    0,
                                   65535,65535,32767,
                                   65535,65535,    0,
                                   65535,32767,32767,
                                   65535,    0,    0,
                                   65535,32767,65535
                                };
#endif

struct {int lplot,irot,il34; float a,b,c,d,asp,thet;} p00000_;
struct {int xorig,yorig;} p00001_;
struct {float psca_x, psca_y, ixo, iyo, iox, ioy;} a00000_;
/***************************************************************
 * initialize xpak globals
 ***************************************************************/
void xpak_init(iro,lpl,ils,dsply)
int  *iro,*lpl,*ils;
Display *dsply;
{
    /*Dsply = dsply;*/
    p00000_.lplot = *lpl;
    p00000_.irot = *iro;
    p00000_.il34 = *ils;
    p00000_.a = 1.0;
    p00000_.b = 0.0;
    p00000_.c = 1.0;
    p00000_.d = 0.0;
    p00000_.asp = 0.66666;

    p00001_.xorig = 0;
    p00001_.yorig = 0;
}
/*
 * This function gets the default colour map of the screen and
 * attempts to set MAX_COLOURS entries
 * The colours set are determined by the irgbv array
 * Parameters: X display opened by calling program
 *                Screen 
 * Return:       number of colours successfully allocated or
               zero if the number was <= USEFL_SIZE
 */

int CreateColormap( dsply,scrn,red,green,blue,irgbv,steps )
Display *dsply;
Screen *scrn;
int *red,*green,*blue,*irgbv;
int steps;
{
    Status result;
    int ncolors = PLOT_COLOURS;
    int *limits, intervals, max;
    float stride;
    int icolr,icolg,icolb,r_inc,g_inc,b_inc,colr,colg,colb;
    int i,k,nc0,nc1;
    int col_start=0, col_end=0;
    unsigned long pixels[MAX_COLOURS];

    cmap=DefaultColormapOfScreen(scrn);

    for ( i=0;i<USEFL_SIZE;i++ ) {
        Colours[i].red  =  red[i] * 256;
        Colours[i].green=  green[i] * 256;
        Colours[i].blue =  blue[i] * 256;
        if ((result=XAllocColor( dsply,cmap,&Colours[i] ))==0) {
            fprintf(stderr,
                "Could not match or allocate color %d\n",i);
            fprintf(stderr,"If displaying on a PC screen, the display should be set to 256 colours\n");
            return( 0 );
        }
    }

    /*
     * set up the Colour array
     */
    intervals = steps-1;

    if ((limits = (int *)malloc(steps*sizeof(int)))==NULL) {
        fprintf(stderr, "malloc failed in CreateColourmap\n");
        return(0);
    }

    limits[0] = USEFL_SIZE;
    stride = ncolors/intervals;
    for (k = 1; k < steps; k++) {
        limits[k]=(int)(k*stride)+USEFL_SIZE;
    }
    limits[steps-1] = max = ncolors+USEFL_SIZE-1;
    for (k = 0; k < intervals; k++) {
        nc0=limits[k];
        nc1=limits[k+1];
        if(nc0 < 0) nc0 = 0;
        if(nc1 > max) nc1=max;
        if (k==intervals-1 && nc1!=max) nc1=max;
        if (k==0) col_start=nc0;

        icolr=k*3;
        icolg=k*3+1;
        icolb=k*3+2;
        r_inc = (irgbv[icolr+3]-irgbv[icolr])/(nc1-nc0);
        g_inc = (irgbv[icolg+3]-irgbv[icolg])/(nc1-nc0);
        b_inc = (irgbv[icolb+3]-irgbv[icolb])/(nc1-nc0);

        colr = irgbv[icolr]; colg=irgbv[icolg]; colb = irgbv[icolb];
        for ( i=nc0;i<=nc1;i++,colr+=r_inc,
                               colg+=g_inc,colb+=b_inc,col_end++ ) {
            Colours[i].red = colr * 256;
            Colours[i].green = colg * 256;
            Colours[i].blue = colb * 256;
            Colours[i].flags=  DoRed | DoGreen | DoBlue;
        }
    }
    /*
     * try to match colours in the default map
     */
    result=1;
    for (i=col_start;i<=col_end && result!=0;i++) 
        result=XAllocColor( dsply,cmap,&Colours[i] );
    if (result==0) {
    /*
     * create new colormap
     */
        while (ncolors&&!XAllocColorCells( dsply,cmap,False,NULL,
                                                  0,pixels,ncolors))
            ncolors--;
        if (ncolors<PLOT_COLOURS) {
            cmap = XCopyColormapAndFree(dsply,cmap);
            ncolors = PLOT_COLOURS;
            while (ncolors && !XAllocColorCells( dsply,cmap,
                                     False,NULL,0,pixels,ncolors))
                ncolors--;
            if (ncolors<(2*intervals+USEFL_SIZE)) {
                fprintf(stderr,"Not enough colours available for plotting\n");
                fprintf(stderr,"If displaying on a PC screen, the display should be set to 256 colours\n");
                    return( 0 );
            }
        }
        for (i=col_start;i<=col_end;i++) 
            Colours[i].pixel=  pixels[i-USEFL_SIZE];
        XStoreColors(dsply,cmap,&Colours[USEFL_SIZE],ncolors);
    }
    
    if (limits!=0) free(limits);
    return( max ); /* max colour index */
}

/*  this version of CreateColourMap creates colours in unused cells in
the Default colour map

CreateColormap( dsply,scrn )
Display *dsply;
Screen *scrn;
{
    Status result;
    int i,c;

    cmap=DefaultColormapOfScreen(scrn);
    if (!(result=XAllocColorCells(dsply,cmap,False,NULL,0,pixels,
                            USEFL_SIZE))) {
        fprintf (stderr," Error in allocating the color cells...\n");
        return(1);
    }

    for ( i=MAP_START,c=0;i<MAP_START+USEFL_SIZE;i++,c++ ) {
        Colours[i].red  =  red[c] * 256;
        Colours[i].green=  green[c] * 256;
        Colours[i].blue =  blue[c] * 256;
        Colours[i].flags=  DoRed | DoGreen | DoBlue;
        Colours[i].pixel=  pixels[c];
    }
    XStoreColors (dsply,cmap,&Colours[MAP_START],USEFL_SIZE);
    return( 0 );
} */
/*--------------------------------------------------*/
/* This function is not called by Sybil
 * setpencolor_() and setlinewidth_() perform these functions
 */
void pen_(ipen,ivel)
int  *ipen, *ivel;

{
  int col, linestyle, linewidth;
  unsigned long valuemask;
  XGCValues values;

  valuemask = GCForeground;
  col = *ipen;
  if (col==1) col = FOREGROUND;
  values.foreground = Colours[col].pixel;

  XChangeGC(Dsply,Draw_gc,valuemask,&values);

/*
 * positive value for *ivel - brush width is *ivel pixels
 * negative value for *ivel - brush width is *ivel/10 cm
 */
  if(*ivel != 0){
    if (*ivel > 0) linewidth= *ivel;
    else linewidth= (int)((0-*ivel)*(a00000_.psca_y/10+0.5));

    linestyle=LineSolid;
    if (Draw_gc) XSetLineAttributes (Dsply,Draw_gc,
                                linewidth,linestyle, CapButt,JoinBevel);
  }
}

/***********************************************************
 * Change the pen colour
 *
 **********************************************************/
void setpencolor_(ipen)
int  *ipen;
{
    int col;
    unsigned long valuemask;
    XGCValues values;

    valuemask = GCForeground;
    col = *ipen;
    if (col==1) col = FOREGROUND;
    values.foreground = Colours[col].pixel;

    if (Draw_gc) XChangeGC(Dsply,Draw_gc,valuemask,&values);
}

/***********************************************************
 * Change the line width
 *
 **********************************************************/
void  setlinewidth_(wdth)
float *wdth;
{
    int linewidth;
    unsigned long valuemask;
    XGCValues values;

    /*
      * positive value for *wdth - brush width is *wdth pixels
      * negative value for *wdth - brush width is *wdth/10 cm
      * rounded up to 1 if less than 1
      */
    if(*wdth != 0){
        if (*wdth > 0) linewidth= (int)(*wdth);
        else linewidth= (int)((0.0 - *wdth)*(a00000_.psca_y/10) + 0.5);
        if (linewidth < 1) linewidth = 1;
    }
    else linewidth = 1;
    valuemask = GCLineWidth;
    values.line_width = linewidth;

    if (Draw_gc) XChangeGC(Dsply,Draw_gc,valuemask,&values);
/*
    linestyle=LineSolid;
    if (Draw_gc) XSetLineAttributes (Dsply,Draw_gc,
                                linewidth,linestyle, CapButt,JoinBevel);
*/
}

/******************************************************************
 *  Set new origin co_ords
 ******************************************************************/

void  origin_(xmin,px1,ymin,py1)
float *xmin,*ymin;
float *px1,*py1;

{
    p00000_.b = *px1 - p00000_.a * (*xmin);
    p00000_.d = *py1 - p00000_.c * (*ymin);
}
/******************************************************************
 *  Set up scale factors used to go between user and plotter coordinates
 ******************************************************************/

void  scale_(xmin,xmax,px1,px2,ymin,ymax,py1,py2)
float *xmin,*xmax,*ymin,*ymax;
float *px1,*px2,*py1,*py2;

{
    p00000_.a = (*px2 - *px1) / (*xmax - *xmin);
    p00000_.b = *px1 - p00000_.a * (*xmin);
    p00000_.c = (*py2 - *py1) / (*ymax - *ymin);
    p00000_.d = *py1 - p00000_.c * (*ymin);
}
/*****************************************************************
 *  Raises (i=3) or lowers (i=2) pen and moves to coordinates (x,y)
 * if i > 0 or to current position plus (x,y) if i < 0
 *****************************************************************/

void  plot_(x,y,i)
float *x,*y;
int   *i;
{
    float            xp, yp;
    float            ixv, iyv;

    xp = *x;
    yp = *y;

/*
    if(p00000_.irot != 0){
        yp = *x;
        xp = -*y;
        if(p00000_.il34 == 0){
            if(*i > 0) xp = 27.2 - *y;
        }
        else if(p00000_.il34 == 1){
            if(*i > 0) xp = 40.1 - *y;
        }
    }
*/
    ixv = xp *  a00000_.psca_x;
    iyv = yp *  a00000_.psca_y;

    if (*i < 0) {
        ixv = a00000_.iox + ixv;
        iyv = a00000_.ioy - iyv;
    }
    else {
        ixv = ixv + p00001_.xorig;
        iyv = a00000_.iyo - p00001_.yorig - iyv;
    }

    if(ABS(*i) == 2) {
        XDrawLine (Dsply,DrawArea,Draw_gc,
             (int)a00000_.iox,(int)a00000_.ioy,
                (int)ixv,(int)iyv);
    }

    /* Update the current position */

    a00000_.iox = ixv;
    a00000_.ioy = iyv;
}

/**************************************************************
 * scale user coordinates to plot coordinates - call plot_
 **************************************************************/
void  plotu_( x,y,i )
float *x, *y;
int   *i;
{
    float            xp, yp;
    float            xpro, ypro;

    xpro = *x;
    ypro = *y;
    xp = xpro * p00000_.a;
    yp = ypro * p00000_.c;

    if(*i > 0) {
        xp = xp + p00000_.b;
        yp = yp + p00000_.d;
    }

    plot_(&xp,&yp,i);

}
/*------------------------------------------------------*/
void symbol_(x,y,size,iword,angl,nchar)

/* Writes a Hollerwith string on the plot (plotter units) */

float  *x,*y,*size,*angl;
char   *iword;
int    *nchar;

{
float  xp, yp, ixv, iyv;

xp = *x;
yp = *y;
/*
if(p00000_.irot != 0){
  yp = *x;
  xp = -*y;
  if(p00000_.il34 == 0){
    if(*size > 0) xp = 27.2 - *y;}
  else if(p00000_.il34 == 1){
    if(*size > 0) xp = 40.1 - *y;}}
*/

ixv = xp *  a00000_.psca_x;
iyv = yp *  a00000_.psca_y;

if (*size < 0.0) {
  ixv = a00000_.iox + ixv;
  iyv = a00000_.ioy - iyv;}
else {
  ixv = p00001_.xorig + ixv;
  iyv = a00000_.iyo - p00001_.yorig - iyv;}

#if XY
if(ABS(*size) != xsiz){
  xsiz=ABS(*size);
  if(xsiz <=  0.25)     currfont = font1;
  else if(xsiz <= 0.35) currfont = font2;
  else if(xsiz <= 0.45) currfont = font3;
  else if(xsiz <= 0.55) currfont = font4;
  else if(xsiz <= 0.65) currfont = font5;
  else                  currfont = font8;
}

XSetFont(mydisplay,mygc,currfont->fid);
/*XDrawString(mydisplay,mywindow,mygc,
            ixv,iyv,iword,*nchar);*/
/* XDrawImageString clears bg rect before drawing string */
XDrawImageString(mydisplay,mywindow,mygc,
            (int)ixv,(int)iyv,iword,*nchar);

/* Update the current position */
len = XTextWidth(currfont,iword,*nchar);
drawlabel_(plotArea,(int)ixv,(int)iyv,&len,&ascent,&descent,1,iword);
a00000_.iox = ixv+len;
a00000_.ioy = iyv;
#endif

}
/*-------------------------------------------------------------*/
void number_(x,y,size,rn,angl,nsf)

/* Writes a number on the plot: if nsf=klm, format is fkl.m (FORTRAN),
                                if nsf = -lm, rn is fixed to an integer
                                   and format is ilm (FORTRAN) */
float *x,*y,*size,*rn,*angl;
int *nsf;

{
int  itot,idpl;
int  num=0;
char iform[8];
char iword[30];

if(*nsf > 0) {
  itot = *nsf / 10;
  if (itot>30) itot = 30;
  idpl = (*nsf)%10;
  num = snprintf(iform,7,"%%%d.%df",itot,idpl);
  num = snprintf(iword,29,iform,*rn);
}
else {
  snprintf(iword,29,"%d",(int)(*rn));
  itot = strlen( iword );
}
symbol_(x,y,size,iword,angl,&itot);
}

/*---------------------------------------------------
 * assumes x and y data are user co-ords
 */
void  drawrectangleu_(x1,y1,x2,y2)
float *x1,*y1,*x2,*y2;
{                  
    int state=3;

    plotu_(x1,y1,&state);
    state = 2;
    plotu_(x2,y1,&state);
    plotu_(x2,y2,&state);
    plotu_(x1,y2,&state);
    plotu_(x1,y1,&state);
}

/*---------------------------------------------------
 * assumes xf and yf data are user co-ords
 * *outline = 0  fill only
 * *outline = 1  fill + outline polygon
 * *outline = 2  outline only
 */
void  fillpoly_(xf,yf,nf,colour,outline)
float *xf, *yf;
int   *nf, *colour, *outline;
{                  
    XPoint *points;
    int ixv, iyv;
    int i, pen_state;
    float xp,yp;
    unsigned long valuemask;
    XGCValues values, oldvalues;

    valuemask = GCForeground;
    XGetGCValues(Dsply,Draw_gc,valuemask,&oldvalues);
    values.foreground = Colours[*colour].pixel;

    XChangeGC(Dsply,Draw_gc,valuemask,&values);

/*
    xpro = *x;
    ypro = *y;
    xp = xpro * p00000_.a;
    yp = ypro * p00000_.c;

    if(*i > 0) {
        xp = xp + p00000_.b;
        yp = yp + p00000_.d;
    ixv = xp *  a00000_.psca_x;
    iyv = yp *  a00000_.psca_y;

    else {
        ixv = ixv + p00001_.xorig;
        iyv = a00000_.iyo - p00001_.yorig - iyv;
    }
*/

    if (*outline < 2) {
        points = (XPoint *) malloc(*nf * sizeof(XPoint));
        for(i = 0; i < *nf; i++){
            xp = *(xf+i) * p00000_.a + p00000_.b;
            yp = *(yf+i) * p00000_.c + p00000_.d;
            ixv = (int)((a00000_.psca_x * xp) + p00001_.xorig +0.5);
            iyv = (int)((a00000_.iyo - yp * a00000_.psca_y) - p00001_.yorig
                                                                     +0.5);
            points[i].x = ixv;
            points[i].y = iyv;
        }

        XFillPolygon(Dsply,DrawArea,Draw_gc,points,
                                        *nf,Complex,CoordModeOrigin);
        free(points);
    }

    XChangeGC(Dsply,Draw_gc,valuemask,&oldvalues);
    if (*outline) {
        /*values.foreground = Colours[FOREGROUND].pixel;*/
        pen_state = 3;
        plotu_(&(xf[0]),&(yf[0]),&pen_state);
        pen_state = 2;
        for (i=1;i<*nf;i++) plotu_(&(xf[i]),&(yf[i]),&pen_state);
        plotu_(&(xf[0]),&(yf[0]),&pen_state);
    }
}
/*
 * Routine to change linestyle
 * sets linestyle to solid line 1 pixel wide if dash<0 or dash>6
 * sets dash pattern according to dash parameter
 * pattern parameter not used
 */
 
int dashln_( dash,pattern )
int *dash, *pattern;
{
/*
    static char      dash1[2]={3,3};
    static char      dash2[2]={6,3};
    static char      dash3[2]={9,3};
    static char      dash4[4]={9,3,3,3};
    static char      dash5[4]={9,3,6,3};
    static char      dash6[6]={9,3,6,3,6,3};
    static char      dash7[2]={1,1};
*/
    /* changed to match lppak */
    static char      dash1[2]={9,3};
    static char      dash2[2]={3,3};
    static char      dash3[2]={4,1};
    static char      dash4[4]={9,3,2,3};
    static char      dash5[4]={9,3,6,3};
    static char      dash6[6]={9,3,6,3,6,3};
    static char      dash7[2]={1,1};
    char *dashstyle;
    int style=LineOnOffDash, width=1, capstyle;
    int len;
    unsigned long valuemask;
    XGCValues values;

    valuemask = GCLineStyle|GCCapStyle;
    switch (*dash) {
    case 0:
    case 1: dashstyle = dash1;
            capstyle=CapButt;
            len=2;
            break;
    case 2: dashstyle = dash2;
            capstyle=CapNotLast;
            len=2;
            break;
    case 3: dashstyle = dash3;
            capstyle=CapProjecting;
            len=2;
            break;
    case 4: dashstyle = dash4;
            capstyle=CapProjecting;
            len=4;
            break;
    case 5: dashstyle = dash5;
            capstyle=CapButt;
            len=4;
            break;
    case 6: dashstyle = dash6;
            capstyle=CapButt;
            len=6;
            break;
    case 7: dashstyle = dash7;
            capstyle=CapButt;
            len=2;
            break;
    default:dashstyle = NULL;
            capstyle=CapButt;
            style = LineSolid;
            break;
    }
    values.line_style = style;
    values.cap_style = capstyle;
    if (dashstyle!=NULL&&Draw_gc) XSetDashes(Dsply,Draw_gc,0,dashstyle,len);
    if (Draw_gc) 
        XSetLineAttributes(Dsply,Draw_gc,width,style,capstyle,JoinBevel);
/*
    if (Draw_gc) XChangeGC(Dsply,Draw_gc,valuemask,&values);
*/
    return(0);
}

/*****************************************************************
 *  Translates user cms to pixels
 *****************************************************************/
int convertutoxpts_(val,pts)
float *val,*pts;
{
    float val_cm;
    val_cm = *val * p00000_.a;
    *pts = a00000_.psca_x * val_cm;
    return(0);
}

int convertutoypts_(val,pts)
float *val,*pts;
{
    float val_cm;
    val_cm = *val * p00000_.c;
    *pts = a00000_.psca_y * val_cm;
    return(0);
}

/*****************************************************************
 *  Translates user coordinates to window coordinates
 *****************************************************************/
void  Translate_u_to_win(float xu,float yu,int* xwin,int* ywin)
/* float xu,yu;
   int   *xwin,*ywin; */
{
    float     tmpx,tmpy;

    tmpx = xu * p00000_.a + p00000_.b;
    tmpy = yu * p00000_.c + p00000_.d;

/*
    if(p00000_.irot != 0){
        tmpy = xu;
        tmpx = -yu;
        if(p00000_.il34 == 0){
            tmpx = 27.2 - yu;
        }
        else if(p00000_.il34 == 1){
            tmpx = 40.1 - yu;
        }
    }
*/
    tmpx *=  a00000_.psca_x;
    tmpy *=  a00000_.psca_y;

    tmpx += p00001_.xorig;
    tmpy = a00000_.iyo - p00001_.yorig - tmpy;

    *xwin = (int)tmpx;
    *ywin = (int)tmpy;
    /* Update the current position */

    a00000_.iox = tmpx;
    a00000_.ioy = tmpy;
}
#if XY
/*
 * create a new colormap for the window
 */
  cmapnew=XCreateColormap(mydisplay,mywindow,myvisual,AllocNone);
  result=XAllocColorCells(mydisplay,cmapnew,False,NULL,0,
pixels,MAX_COLOURS);
  if (result==0) printf ("2 Error in allocating the color cells...\n");

  for (k = 0; k < *nstep; k++) {
    nc0=pixv[k];
    nc1=pixv[k+1];
    if(nc0 < 0) { nc0 = 0;}
    if(nc1 > ncolors-1) { nc1=ncolors-1;}
    ncr1=nc1-nc0;

    if(ncr1 > 0) {
      ip=nc0;
      icolr=k*3;
      icolg=k*3+1;
      icolb=k*3+2;
      color[ip].red = rgbv[icolr];
      color[ip].green = rgbv[icolg];
      color[ip].blue = rgbv[icolb];
      color[ip].flags=  DoRed | DoGreen | DoBlue;
      color[ip].pixel=  pixels[ip];

      if(ncr1 > 2) {
        for (i = 1;i <= ncr1; i++) {
          ip=i+nc0;
          color[ip].red = (rgbv[icolr+3] * i + rgbv[icolr] *
(ncr1-i))/ncr1;
          color[ip].green = (rgbv[icolg+3] * i + rgbv[icolg] *
(ncr1-i))/ncr1;
          color[ip].blue = (rgbv[icolb+3] * i + rgbv[icolb] *
(ncr1-i))/ncr1;
          color[ip].flags=  DoRed | DoGreen | DoBlue;
          color[ip].pixel=  pixels[ip];
        }
      }
    }
  }

  for (i = 0;i < 16; i++)
  {
  color[i].red  =  0;
  color[i].green=  0;
  color[i].blue =  0;
  color[i].flags=  DoRed | DoGreen | DoBlue;
  color[i].pixel=  pixels[i];
  }
/*
 * Useful colours
 * colmap_ reads them into indices 16-31
  for (i = 16,k=0;i < 32; i++,k++)
  {
  color[i].red  =  red[k] * 256;
  color[i].green=  green[k] * 256;
  color[i].blue =  blue[k] * 256;
  color[i].flags=  DoRed | DoGreen | DoBlue;
  color[i].pixel=  pixels[i];
 */
/*
 * Useful colours
 * colmap_ reads them into indices 8-15
 * no greyscale
 */
  for (i = 8,k=0;i < 16; i++,k++)
  {
  color[i].red  =  red[k] * 256;
  color[i].green=  green[k] * 256;
  color[i].blue =  blue[k] * 256;
  color[i].flags=  DoRed | DoGreen | DoBlue;
  color[i].pixel=  pixels[i];
  }
  XStoreColors(mydisplay,cmapnew,color,ncolors);
  XSetWindowColormap(mydisplay,mywindow,cmapnew);
}
#endif

int setfont_(fontname,nchar,fonthgt)
char *fontname;
int *nchar,*fonthgt;
{
    FindFont(fontname,*fonthgt);
    return(0);
}

int drawticmarks_(size,stop,val,start,interval,justify,horizontal)
int *justify, *horizontal;
float *size;
float *val,*start,*stop,*interval;

{
    int pen_state;
    float offset,xval,yval,p1,p2,pint;
    float tmpsize;

/* amended to allow start > stop; previously not working in that case */

    if(*start <= *stop) {p1=*start; p2=*stop;}
    else                {p1=*stop; p2=*start;}
    pint=fabs(*interval);

    tmpsize = *size; /* don't change passed value */
    if (*justify==J_CENTRE) offset = -(*size)/2;
    else offset = 0;
    if (*justify==J_TOP || *justify==J_RIGHT) tmpsize = -tmpsize;
    if (*horizontal) {
        for (xval = p1; xval < p2; xval += pint) {
            pen_state = PEN_UP;
            yval = *val + offset;
  /*        fprintf(stdout,"xval =%f yval = %f\n",xval,yval); */
            plotu_(&xval,&yval,&pen_state);
            pen_state = PEN_DN;
            yval = *val+tmpsize;
            plotu_(&xval,&yval,&pen_state);
        }
    }
    else {
        for (yval = p1; yval < p2; yval += pint) {
            pen_state = PEN_UP;
            xval = *val + offset;
   /*       fprintf(stdout,"xval =%f yval = %f\n",xval,yval);a */
            plotu_(&xval,&yval,&pen_state);
            pen_state = PEN_DN;
            xval = *val+tmpsize;
            plotu_(&xval,&yval,&pen_state);
        }
    }
    return(0);
}

#if XY
#endif
int drawcircle_(xmin,ymin,diameter,colour,outline)
float *xmin,*ymin;
float *diameter;
int *colour, *outline;
{
    int ixv, iyv, idiam;
    float xp,yp,xp1;
    unsigned long valuemask;
    XGCValues values, oldvalues;

    valuemask = GCForeground;
    XGetGCValues(Dsply,Draw_gc,valuemask,&oldvalues);
    values.foreground = Colours[*colour].pixel;

    XChangeGC(Dsply,Draw_gc,valuemask,&values);

/*
    xpro = *x;
    ypro = *y;
    xp = xpro * p00000_.a;
    yp = ypro * p00000_.c;

    if(*i > 0) {
        xp = xp + p00000_.b;
        yp = yp + p00000_.d;
    ixv = xp *  a00000_.psca_x;
    iyv = yp *  a00000_.psca_y;

    else {
        ixv = ixv + p00001_.xorig;
        iyv = a00000_.iyo - p00001_.yorig - iyv;
    }
*/

    if (*outline < 2) {
        xp = *xmin * p00000_.a + p00000_.b;
        yp = *ymin * p00000_.c + p00000_.d;
        ixv = (int)((a00000_.psca_x * xp) + p00001_.xorig +0.5);
        iyv = (int)((a00000_.iyo - yp * a00000_.psca_y) - p00001_.yorig
                                                                     +0.5);
        xp1 = (*xmin + *diameter) * p00000_.a + p00000_.b;
        idiam = (int)((a00000_.psca_x * (xp1-xp)) + p00001_.xorig +0.5);
        XFillArc(Dsply,DrawArea,Draw_gc,ixv,iyv,idiam,idiam,
                                                         0,64*360);
    }

    XChangeGC(Dsply,Draw_gc,valuemask,&oldvalues);
    if (*outline) {
        /*values.foreground = Colours[FOREGROUND].pixel;*/
        /*XChangeGC(Dsply,Draw_gc,valuemask,&values);*/
        XFillArc(Dsply,DrawArea,Draw_gc,ixv,iyv,idiam,idiam,
                                                         0,64*360);
    }
    return(0);
}

void formatnumber_(str,num,format)
char *str;
float *num;
int *format;
{
    char form[10];
    int total,dec_places,i,len;

    /* check range as form cannot have > 2 digits */
    total = (*format/10)%100;
    dec_places = *format%10;
    snprintf(form,9,"%%%d.%df",total,dec_places);
    snprintf(str,MAX_NUM_LEN,form,*num);
    /* strip leading and trailing spaces */
    while (str[0]==' ') {
        len = strlen(str)+1;
        for (i=0;i<len;i++) str[i] = str[i+1];
    }
    len = strlen(str);
    while (str[len]==' ') {
        str[len] = '\0';
        len--;
    }
}
