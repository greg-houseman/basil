/*--------------------------------------------------------------------
 *    Basil / Sybil:   sybmain.c   1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef XSYB
#include <X11/Xlib.h>  /* for type XFontStruct */
#endif

#include "types.h"
#include "data.h"
#include "errnum.h"
#include "error.h"
#include "colours.h"
#include "arrays.h"
#include "routines.h"

float Contour_vals[MAXCNTRVALS], Profile_vals[MAXPRFLVALS];
plot_data Plot_info;
input_options Settings, Initial_Settings;
int Mesh[MESH_ENTRIES];
float Pwindo[PWINDO_ENTRIES];

void init_arrays(),clear_mesh_arrays(),MeshParams();
void Init_Settings(), Init_Plot(),Init_Data();
void set_box_dimensions(),init_box();
void PrintComments();

int main(int argc, char** argv)
{
    int err=0;
    int *red,*green,*blue,*irgbv;

    Init_Data();

    if (Init_App(argc,argv,&Plot_info.log_file_ptr)) exit(1);

    set_max_cntrs(Contour_vals, Settings.contour.ndivr);

    switch(Settings.colourmap) {
    case GREY_MAP: red=Red0;
                   green=Green0;
                   blue=Blue0;
                   irgbv=Rgb0;
                   break;
    case STD_MAP:  red=Red1;
                   green=Green1;
                   blue=Blue1;
                   irgbv=Rgb1;
                   break;
    case ABS_MAP:  red=Red1;
                   green=Green1;
                   blue=Blue1;
                   irgbv=Rgb2;
                   break;
    default:       fprintf(stderr,"Settings.colourmap = %i\n",Settings.colourmap);
                   red=Red1;
                   green=Green1;
                   blue=Blue1;
                   irgbv=Rgb1;
                   break;
    }
    if ((err = Set_Colours(red,green,blue,irgbv,STEPS))) return(1);

    /* save the initial options settings */
    Initial_Settings = Settings;

    err = Run_App( Plot_info.log_file_ptr );
    return(err);
} 

void Init_Data()
{
    int i;

    Init_Plot();
    Init_Settings();

    for (i=0;i<PWINDO_ENTRIES;i++) Pwindo[i]=0.0;
    Pwindo[XCMAX] = 1.0;
    Pwindo[YCMAX] = 1.0;
    Pwindo[XCMAXREF] = 1.0;
    Pwindo[YCMAXREF] = 1.0;

    for (i=0;i<MAXCNTRVALS;i++) Contour_vals[i]=0.0;
    Contour_vals[CNTR_SCL]=1.0;
    Contour_vals[MAX_CNTRS]=Settings.contour.ndivr;
    Contour_vals[NUM_CNTRS]=10;

    for (i=0;i<MAXPRFLVALS;i++) Profile_vals[i]=0.0;
    Profile_vals[PRFL_SCL]=1.0;

    init_arrays();

    for (i=0;i<REC_1_CNT;i++) String_vars[i] = ' ';
}

void Init_Plot()
{
    Plot_info.curr_cell = 0;
    Plot_info.max_col = DEFAULT_COLUMNS-1;
    Plot_info.max_row = DEFAULT_ROWS-1;
    Plot_info.curr_x = Plot_info.curr_y = 0.0;
    Plot_info.origin_x = Plot_info.origin_y = 0.0;
    Plot_info.title_offset = DEFAULT_LGE_FONTHGT * 1.5;
    Plot_info.plot_type = -1;
    Plot_info.plot_description = 0;
    Plot_info.var_num = 0;
    Plot_info.title = NULL;
    strcpy(Plot_info.fontname,"Helvetica");
    strcpy(Plot_info.variable,"");
    strcpy(Plot_info.dflt_label1,"");
    strcpy(Plot_info.dflt_label2,"");
    strcpy(Plot_info.curr_label,"");
    Plot_info.display_log = 0;
    Plot_info.log_file_ptr = NULL;
    Plot_info.inp_file = NULL;
    Plot_info.elle_file= NULL;
}

void Init_Settings()
{
    Settings.orient = LANDSCAPE;
    Settings.paper_size = A4_PAPER;
    Settings.mark_cell = 1;
    Settings.rescale = 1;
    Settings.filepath = 1;
    Settings.verbose = 0;
    Settings.elle = 0;
    Settings.clip_to_cell = 0;
    Settings.text_bg = 1;
    Settings.rows = DEFAULT_ROWS;
    Settings.columns = DEFAULT_COLUMNS;
    Settings.dfltfonthgt = DEFAULT_FONTHGT;
    Settings.fg = FOREGROUND;
    Settings.linestyle = SOLID;
    Settings.page_scale = 1.0;
    Settings.prev_zoom = 1.0;
    Settings.zoom = 1.0;
    Settings.xcentre = 0.0;
    Settings.ycentre = 0.0;
    Settings.linewidth = 1.0;
    Settings.xmargin = XLMARG;
    Settings.ymargin = YMARG;
    Settings.page_xmargin = PG_XMARG;
    Settings.page_ymargin = PG_YMARG;
    Settings.viscmin = 0.0;
    Settings.viscmax = 0.0;
    Settings.semin = 0.0;
    Settings.semax = 0.0;
    Settings.colourmap = STD_MAP;

    Settings.physical.tl0=100;
    Settings.physical.rhoc=2.9;
    Settings.physical.rhom=3.3;
    Settings.physical.z0=160;
    Settings.physical.tmax=2.5;
    Settings.physical.elev0=0.25;
    Settings.physical.bgam0=0.0725;
    Settings.physical.bgam1=0.139;

    Settings.contour.ndivr=20.0;
    Settings.contour.s0=0.0;
    Settings.contour.dels=0.025;
    Settings.contour.delom=10.0;

    Settings.plot_opts.nx3 = NNX3_SIZE;
/*  Settings.plot_opts.mp=4;
    Settings.plot_opts.np=4;  */
    Settings.plot_opts.mp=(NNX3_SIZE-3)/20;
    Settings.plot_opts.np=(NNX3_SIZE-3)/20;
    Settings.plot_opts.mpe=6;
    Settings.plot_opts.contour_plot = LINES | SHADE;
    Settings.plot_opts.ticmark = J_TOP;
    Settings.plot_opts.colour_bar=CB_VERT;
    Settings.plot_opts.profile_pts=51;
    Settings.plot_opts.flip=0;
    Settings.plot_opts.dble=0;
    Settings.plot_opts.label=1;
    Settings.plot_opts.stipple=0;
    Settings.plot_opts.solution_rot=0;
}

void clear_mesh_arrays()
{
    int i;

    for (i=PLOT_INT_ARRAYS;i<MAX_INT_ARRAYS;i++)
        if (Data_arrays_int[i]!=NULL) {
            free(Data_arrays_int[i]);
            Data_arrays_int[i]=NULL;
        }
    for (i=PLOT_FL_ARRAYS;i<MAX_FL_ARRAYS;i++)
        if (Data_arrays_fl[i]!=NULL) {
            free(Data_arrays_fl[i]);
            Data_arrays_fl[i]=NULL;
        }
}

/*
 * if a file has been opened, clear the mesh arrays and reset
 * the interpolation mesh dimensions
 */
void NewMesh()
{
    if (Plot_info.inp_file!=NULL) {
        Mesh[NNY3] = 0;
        MeshParams( Pwindo[YCMAX]-Pwindo[YCMIN],Pwindo[XCMAX]-Pwindo[XCMIN],
                                                                 Mesh );
    }
}

void MeshParams( hgt,wdth,mesh )
float hgt, wdth;
int *mesh;
{
    int prev_nny3,prev_nnx3;
 /* 
  * put the specified resolution on the longest dimension 
  */
    prev_nnx3=mesh[NNX3];
    prev_nny3=mesh[NNY3];
    if (hgt < wdth) {
        mesh[NNX3] = Settings.plot_opts.nx3;
        mesh[NNY3] = (int)((mesh[NNX3]-3)*hgt/wdth)+3;
    }
    else {
        mesh[NNY3] = Settings.plot_opts.nx3;
        mesh[NNX3] = (int)((mesh[NNY3]-3)*wdth/hgt)+3;
    }
 /*
  * set array size for work arrays using value in cmndefs.h
  */
    if((mesh[NNX3]!=prev_nnx3)||(mesh[NNY3]!=prev_nny3)) {
        clear_mesh_arrays();
        mesh[NP3] = mesh[NNX3]*mesh[NNY3];
        mesh[ICON] = mesh[NP3]/7;
        mesh[NX3] = mesh[NNX3];
        mesh[NY3] = mesh[NNY3];
        mesh[NX2]=mesh[NX3]-1;
        mesh[NX1]=mesh[NX2]-1;
        mesh[NX0]=mesh[NX1]-1;
        mesh[NY2]=mesh[NY3]-1;
        mesh[NY1]=mesh[NY2]-1;
        mesh[NY0]=mesh[NY1]-1;
        mesh[NXY]=mesh[NX3]*mesh[NY3];
    }
}

void set_box_dimensions( pwin,x,y,width,height )
float *pwin,x,y,width,height;
{
    /* Pwindo values for plotting routines */
    pwin[XMIN] = x + width*Settings.xmargin;
    pwin[YMIN] = y - height*(1.0-Settings.ymargin);
    pwin[XMAX] = pwin[XMIN] + width*(1.0-2.0*Settings.xmargin);
    pwin[YMAX] = pwin[YMIN] + height*(1.0-2.0*Settings.ymargin);
//  printf("x = %f y = %f width = %f height = %f Settings.xmargin = %f Settings.ymargin = %f\n",
//          x,y,width,height,Settings.xmargin,Settings.ymargin);
//  printf("pwin[XMIN,XMAX,YMIN,YMAX = %f %f %f %f\n",pwin[XMIN],pwin[XMAX],pwin[YMIN],pwin[YMAX]);
}
 
void init_box(float* pwin, input_options* opts)
{
    int pen_state,ncomp;
    int err=0,len=1;
    float plotwdth,plothght,wdth,hght,urangex,urangey;
    float pcentrex,pcentrey,ucentrex,ucentrey,phwdthx,phhghty;
    float zoomnow,scalex,scaley,poffsetx,poffsety,toler;
    float xref,yref;
    float xledge,xredge,ybedge,ytedge;
    ncomp = Data_vars_int[NCOMP];
    wdth = pwin[XMAX] - pwin[XMIN];
    hght = pwin[YMAX] - pwin[YMIN];
    toler = hght;
    if(wdth > toler)toler=wdth;
    toler=0.005*toler;
/*
 *  on entry, XCMINREF etc are user-unit limits on solution domain
 *    whereas XCMIN etc are user-unit limits on active plotting
 *  If clip_to_cell is activated plotting will be restricted to
 *  domain constrained by XCMIN etc redefined here.
 */
    urangex = pwin[XCMAXREF]-pwin[XCMINREF];
    urangey = pwin[YCMAXREF]-pwin[YCMINREF];
    plotwdth = wdth;
    plothght = plotwdth*urangey/urangex;
    if (plothght > hght) {
        plothght = hght;
        plotwdth = plothght*urangex/urangey;
    }
/*
 *  set scale and centre the plot in the available space
 */
    scalex=scaley=urangex/plotwdth;
    pwin[XMARGE] = 0.5*(wdth-plotwdth);
    pwin[YMARGE] = 0.5*(hght-plothght);
/*
 * scalars that start with p are in plotter units, with u in user units
 */
    ucentrex = 0.5*(pwin[XCMINREF]+pwin[XCMAXREF]);
    ucentrey = 0.5*(pwin[YCMINREF]+pwin[YCMAXREF]);
//  printf("wdth = %f hght = %f urangex = %f urangey = %f scalex = %f ucentrex = %f ucentrey = %f\n",
//          wdth,hght,urangex,urangey,scalex,ucentrex,ucentrey);
    phwdthx = 0.5*wdth;
    phhghty = 0.5*hght;
    pcentrex = 0.5*(pwin[XMAX]+pwin[XMIN]);
    pcentrey = 0.5*(pwin[YMAX]+pwin[YMIN]);
    zoomnow=1.0;
/*
 *   by using XCMINREF, YCMINREF we retain the current scale unless
 *   changed to current solution domain limits by a rescale command.  
 *   Origin may be shifted to another cell, in the same relative location.
 */
    pwin[ULEFT] = pwin[XCMINREF] - pwin[XMARGE]*scalex;
    pwin[URGHT] = pwin[XCMAXREF] + pwin[XMARGE]*scalex;
    pwin[UBASE] = pwin[YCMINREF] - pwin[YMARGE]*scaley;
    pwin[UTOP]  = pwin[YCMAXREF] + pwin[YMARGE]*scaley;
/*
 *    limits of plot in user coordinates, are set by CheckCurrentRec
 *    and modified below if plot is clipped or zoomed */
/*
 *    the values in xcentre, ycentre are used to center the plot on the page
 *    only if zoom != 1. Otherwise plot centre is basec on XCMINREF, etc
 */
    zoomnow = opts->zoom;
    if(fabs(zoomnow-1.0)>1.e-4) {
        ucentrex = opts->xcentre;
        ucentrey = opts->ycentre;
        scalex = scalex/zoomnow;
        scaley = scaley/zoomnow;
/*  project zoom centre from spherical to (x-y) if required */
        if (ncomp==-1) {               /* sinusoidal projection */
           xref=Data_vars_fl[XREFM];
           yref=Data_vars_fl[YREFM];
           if (Data_vars_int[IDEFTYP]>110) {  /* corotate applied */
              xref=0.;
              yref=0.;
           }
           projectxy_(&ucentrex,&ucentrey,
                       &xref,&yref,&len,&ncomp,&err);
        }
     }
/* set the model unit limits relative to current centre using a margin   */
     pwin[SCLX]=scalex;
     pwin[SCLY]=scaley;
     pwin[ULEFT] = ucentrex-phwdthx*scalex;
     pwin[URGHT] = ucentrex+phwdthx*scalex;
     pwin[UBASE] = ucentrey-phhghty*scaley;
     pwin[UTOP]  = ucentrey+phhghty*scaley;
/*
 *    set the scale for the graph that follows
 */
     scale_(&pwin[ULEFT],&pwin[URGHT],&pwin[XMIN],&pwin[XMAX],
           &pwin[UBASE],&pwin[UTOP],&pwin[YMIN],&pwin[YMAX]);
     pen_state=PEN_UP;
     plotu_(&pwin[ULEFT],&pwin[UBASE],&pen_state);
/*
 *    if clipping activated, plot is restricted to the domain
 *               XCMIN, XCMAX, YCMIN, YCMAX reset here
 */
     poffsetx=pcentrex - ucentrex/scalex;
     poffsety=pcentrey - ucentrey/scaley;
     if (opts->clip_to_cell) {
        xledge=pwin[XCMINREF]/scalex + poffsetx;
        if((pwin[XMIN] - xledge) > -toler) 
            pwin[XCMIN]=(pwin[XMIN]-poffsetx)*scalex;
	xredge=pwin[XCMAXREF]/scalex + poffsetx;
        if((pwin[XMAX] - xredge) < toler) 
            pwin[XCMAX]=(pwin[XMAX]-poffsetx)*scalex;
        ybedge=pwin[YCMINREF]/scaley + poffsety;
        if((pwin[YMIN] - ybedge) > -toler) 
            pwin[YCMIN]=(pwin[YMIN]-poffsety)*scaley;
        ytedge=pwin[YCMAXREF]/scaley + poffsety;
        if((pwin[YMAX] - ytedge) < toler)
            pwin[YCMAX]=(pwin[YMAX]-poffsety)*scaley;
        fprintf(stdout,"Plot clip limits are: X: %f to %f, Y: %f to %f\n",
                          pwin[XCMIN],pwin[XCMAX],pwin[YCMIN],pwin[YCMAX]);
    }
/*
 *   check the plot is not moved off the page entirely
 *   as could happen with a poor choice of (xcentre,ycentre)
 */ 
    if(((pwin[XCMAX]/scalex+poffsetx) < pwin[XMIN]) ||
       ((pwin[XCMIN]/scalex+poffsetx) > pwin[XMAX]) ||
       ((pwin[YCMAX]/scaley+poffsety) < pwin[YMIN]) ||
       ((pwin[YCMIN]/scaley+poffsety) > pwin[YMAX]) ) {
       printf("Choice of zoom centre coordinates +/- clipping has ");
       printf("moved the plot entirely off the page\n");
    }
//  fprintf(stdout,"init_box: XCMIN,XCMAX,YCMIN,YCMAX:  %f %f %f %f\n",
//          pwin[XCMIN],pwin[XCMAX],pwin[YCMIN],pwin[YCMAX]);
    opts->rescale=0;
    opts->prev_zoom = opts->zoom;
}

int CheckCurrentRec(char* msg, int* err_num)
{
    int err=0;
    float ucxmin, ucxmax, ucymin, ucymax;
    float eps=1E-4;

    if (Plot_info.inp_file==NULL || Plot_info.inp_file->fp==NULL) {
        strcpy(msg,"No solution file open");
        return(1);
    }
/*
 *  read the relevant data record if not already in memory
 */
    if (Plot_info.inp_file->rec_curr!=Plot_info.inp_file->rec_req) {
        err=read_data(Plot_info.inp_file->fp,
                        String_vars,
                        Data_vars_int,Data_vars_fl,
                        Data_arrays_int,Data_arrays_fl,
			&ucxmin,&ucxmax,&ucymin,&ucymax,
                        &Plot_info.inp_file->rec_curr,
                        Plot_info.inp_file->rec_req,
                        &Plot_info.inp_file->rec_max,
                        Settings.plot_opts.solution_rot
                        );
        Pwindo[XCMIN]=ucxmin;
        Pwindo[XCMAX]=ucxmax;
        Pwindo[YCMIN]=ucymin;
        Pwindo[YCMAX]=ucymax;
        if (err) {
            strcpy(msg,"");
            *err_num = err;
        }
/*
 *  records created before 1.6.5 will have the default VC value of zero
 *    (Data_vars_fl[63] was unused in earler versions)
 *  Set to 1.0 for use in STRAIX, STRAIL if IVIS==0 (no vhb array)
 */
        if (fabs(Data_vars_fl[IVC])<eps)
            Data_vars_fl[IVC]=1.0;
    }
    if (!err) {
      if ((Pwindo[SCLX]==0 && Pwindo[SCLY]==0) || Settings.rescale!=0) {
        Pwindo[XCMINREF]=Pwindo[XCMIN];
    	Pwindo[XCMAXREF]=Pwindo[XCMAX];
    	Pwindo[YCMINREF]=Pwindo[YCMIN];
    	Pwindo[YCMAXREF]=Pwindo[YCMAX];
      }
      Settings.rescale = 0;
    }
    init_box( Pwindo,&Settings );
   /* YCMIN, YCMAX may be reset by init_box in case zoom/clip is activated */	
    MeshParams( Pwindo[YCMAX]-Pwindo[YCMIN],
                Pwindo[XCMAX]-Pwindo[XCMIN],Mesh );
    return(err);
}
 
int CheckReferenceRec(char* msg, int* err_num)
{
    int err=0;
    if (Plot_info.inp_file->ref_curr!=Plot_info.inp_file->ref_req) {
        err=read_reference(Plot_info.inp_file->fp,
                    Data_vars_int,Data_arrays_fl,
                    Plot_info.inp_file->rec_curr,
                    &Plot_info.inp_file->ref_curr,
                    Plot_info.inp_file->ref_req,
                    Settings.plot_opts.solution_rot
                    );
        if (err) {
            strcpy(msg,"");
            *err_num = err;
        }
    }
    return(err);
}

void CreateErrorMessage( message,err_num )
char *message;
int err_num;
{
    char buf[20];

    switch (err_num) {
    case MALLOC_ERR:strcat(message," - Memory error");
                break;
    case READ_ERR:strcat(message," - Error reading file");
                break;
    case OPEN_ERR:strcat(message," - Error opening file");
                break;
    case EOF_ERR:strcat(message," - End of file");
                if (Plot_info.inp_file!=NULL)
                sprintf(buf," %2d records",Plot_info.inp_file->rec_max);
                strcat(message,buf);
                break;
    case VHB_ERR:strcat(message," - No VHB array");
                break;
    case SSQ_ERR:strcat(message," - No SSQ array");
                break;
    case FROT_ERR:strcat(message," - No FROT array");
                break;
    case DENS_ERR:strcat(message," - No DENS array");
                break;
    case ELFX_ERR:strcat(message," - No IELFIX array");
                break;
    case REF_ERR:strcat(message," - Reference error");
                break;
    case MSHTYP_ERR:strcat(message," - Calculation not valid for mesh \ngenerated from poly file");
                break;
    case ISEG_ERR:strcat(message," - No segment array");
                break;
    case LMSH_ERR:strcat(message," - No Lagrangian mesh arrays");
                break;
    case LMRK_ERR:strcat(message," - No marker arrays");
                break;
    case STR_ERR:strcat(message," - Error reading string");
                break;
    case NUM_ERR:strcat(message," - Error reading number");
                break;
    case VAR_ERR:strcat(message," - Unknown variable");
                break;
    case KEY_ERR:strcat(message," - Unknown keyword");
                break;
    case INIT_ERR:strcat(message," - Initial variable must be Options");
                break;
    case RANGE_ERR:strcat(message," - Value out of range");
                break;
    case STGS_ERR:strcat(message," - row, col or orientation mismatch");
                break;
    case PATH_ERR:strcat(message," - full path saved to log file");
                break;
    default:    break;
    }
}

void error_msg(err_num,str)
int err_num;
char *str;
{
    char message[SYB_FILENAME_MAX];

    strncpy(message,str,sizeof(message)-1);
    CreateErrorMessage(message,err_num);
    fprintf(stderr,"%s\n",message);
    Exit_App(err_num);
}

void PrintComments()
{
    char str[REC_1_CNT+1];
    int i,j;

    for (i=0,j=DATE_START;i<8;i++,j++)
        str[i]=String_vars[j];
    str[i]='\0';
    fprintf(stdout,"\nBasil run commenced %s",str);
    for (i=0;i<8;i++,j++)
        str[i]=String_vars[j];
    str[i]='\0';
    fprintf(stdout," at %s\n",str);
    for (i=0,j=COMMENTS_START;i<80;i++,j++)
        str[i]=String_vars[j];
    str[i]='\0';
    fprintf(stdout,"%s\n",str);
    fflush(stdout);
}
