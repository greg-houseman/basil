/*--------------------------------------------------------------------
 *    Basil / Sybil:   pref77.c  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#ifdef XSYB
#include <X11/Xlib.h>  /* for type XFontStruct in Plot_info - tidy */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "types.h"
#include "errnum.h"
#include "error.h"
#include "data.h"
#include "strain.h"
#include "deform.h"
#include "pref77.h"
#include "plotdefs.h"
#include "arrays.h"
#include "routines.h"

#define XY 0

int Colour_range[2], Arrow_colours[2];
/* float Calc_Step(), Calc_Level(); */

int plot_strain_markers(int* int_vars, float* fl_vars, float** arraysf, float* pwindo)
{
    int verbose;
    float xmin,xmax,ymin,ymax,tmp;

    if (arraysf[STELPX]==NULL) return(LMRK_ERR);
    if (arraysf[EXLG]!=NULL && arraysf[EYLG]!=NULL)
        find_max_min(arraysf[EXLG],arraysf[EYLG],
                     int_vars[NUL],
                     &xmin,&xmax,&ymin,&ymax);
    else find_max_min(arraysf[STELPX],arraysf[STELPY],
                     int_vars[NSM] * int_vars[NPM],
                     &xmin,&xmax,&ymin,&ymax);
        tmp=pwindo[XCMIN]; pwindo[XCMIN]=xmin; xmin=tmp;
        tmp=pwindo[XCMAX]; pwindo[XCMAX]=xmax; xmax=tmp;
        tmp=pwindo[YCMIN]; pwindo[YCMIN]=ymin; ymin=tmp;
        tmp=pwindo[YCMAX]; pwindo[YCMAX]=ymax; ymax=tmp;
    verbose=0;
    stmark_(arraysf[STELPX],arraysf[STELPY],
                  &int_vars[NSM],&int_vars[NPM],&int_vars[NRM],
                  &fl_vars[STELPR],
                  &verbose,
                  &pwindo[XCMIN],&pwindo[XCMAX],
                  &pwindo[YCMIN],&pwindo[YCMAX]
                  );
    verbose = Settings.verbose;
    Save_plot_data( NULL,NULL );
    pwindo[XCMIN]=xmin;
    pwindo[XCMAX]=xmax;
    pwindo[YCMIN]=ymin;
    pwindo[YCMAX]=ymax;
    return(0);
}

int plot_mesh(int option, int* int_vars, float* fl_vars, int** arraysi, 
              float** arraysf, float* pwindo, int* mesh, int* current_pen)
{
    int verbose;
    int opt, i, *iptr;
    int startcol,endcol,rangecol,mpe;
    int orient_vert=1;
    int reset=0;
    int ireg,isvp;
    float xmin,xmax,ymin,ymax,tmp,vmin,vmax,densreg,xref,yref;
    float xlevel,ylevel,barwdth,barlen,scale=1.0;

    if ((option==VISCMSH||option==VISCDBL||option==SEMSH||option==SEDBL)
					&&arraysf[VHB]==NULL)
        return(VHB_ERR);
    if ((option==INTLBNDRY||option==INTLBNDRYDBL)&&arraysi[ISEG]==NULL)
        return(ISEG_ERR);
    if ((option==LGBNDRY||option==LGELMNT)) {
        if (arraysf[EXLG]==NULL) return(LMSH_ERR);
        reset=1;
        find_max_min(arraysf[EXLG],arraysf[EYLG],
                         int_vars[NUL],
                         &xmin,&xmax,&ymin,&ymax);
        tmp=pwindo[XCMIN]; pwindo[XCMIN]=xmin; xmin=tmp;
        tmp=pwindo[XCMAX]; pwindo[XCMAX]=xmax; xmax=tmp;
        tmp=pwindo[YCMIN]; pwindo[YCMIN]=ymin; ymin=tmp;
        tmp=pwindo[YCMAX]; pwindo[YCMAX]=ymax; ymax=tmp;
    }
/*    XREFM and YREFM usage depends on context */
    xref=fl_vars[XREFM];
    yref=fl_vars[YREFM];
    if (int_vars[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }
    startcol = Colour_range[0];
    endcol = Colour_range[1];
    rangecol = Colour_range[1] - startcol;
    mpe = Settings.plot_opts.mpe;
    vmin = Settings.viscmin; vmax = Settings.viscmax;
    if (option==SEMSH || option==SEDBL) {
      vmin = Settings.semin; vmax = Settings.semax;
    }
    verbose = Settings.verbose;
    if ((option==LGBNDRY||option==LGELMNT)) {
      /*
       * Lagrangian plot chosen
       */
      if (arraysi[LGEM]==NULL) return(LMSH_ERR);
      if (arraysi[LNOR]==NULL) {
        if ((arraysi[LNOR]=(int *)malloc(int_vars[NUL]*sizeof(int)))==NULL)
            return(MALLOC_ERR);
        for (i=1,iptr=arraysi[LNOR]; i<=int_vars[NUL]; i++,iptr++) *iptr = i;
      }
      if (option==LGBNDRY) opt = BNDRY;
      else if (option==LGELMNT) opt = ELMNT;
      verbose=0;
      plmesh_(&opt,&int_vars[NEL],&int_vars[NUL],&int_vars[NUL],
                        &int_vars[NBL],arraysf[EXLG],arraysf[EYLG],
                        arraysi[LGEM],arraysi[LNOR],arraysi[LGIBC],
                        arraysi[LGIBCF],arraysf[VHB],
                        &fl_vars[TBXOFF],&fl_vars[TBYOFF],&int_vars[IFLT],
                        &mpe,&startcol,&rangecol,
                        &pwindo[XCMIN],&pwindo[XCMAX],
                        &pwindo[YCMIN],&pwindo[YCMAX],
                        &vmin,&vmax,&Settings.zoom,
			&int_vars[NCOMP],&xref,&yref,&verbose );
       verbose = Settings.verbose;
    }
    else if (option==INTLBNDRY) {
        verbose=0;
/*	fprintf(stdout,"calling plseg_, NSEG = %i\n",int_vars[NSEG]); */
        plseg_(&option,&int_vars[NUP], &int_vars[NBP], &int_vars[NSEG],
                        arraysf[EX],arraysf[EY],arraysi[NOR],arraysi[IBC],
                        arraysi[IBCTYP],arraysi[ISEG],&rangecol,
			&pwindo[XCMIN],&pwindo[XCMAX],
			&pwindo[YCMIN],&pwindo[YCMAX],
			&int_vars[NCOMP],&xref,&yref,&verbose );
        verbose = Settings.verbose;
    }
    else if ((option!=INTBND && option!=INTBNDDBL)) {
    
        verbose=0;
        plmesh_(&option,&int_vars[NE],&int_vars[NN],&int_vars[NUP],
                        &int_vars[NBP],arraysf[EX],arraysf[EY],
                        arraysi[LEM],arraysi[NOR],arraysi[IBC],
                        arraysi[IBCTYP],arraysf[VHB],
                        &fl_vars[TBXOFF],&fl_vars[TBYOFF],&int_vars[IFLT],
                        &mpe,&startcol,&rangecol,
                        &pwindo[XCMIN],&pwindo[XCMAX],
                        &pwindo[YCMIN],&pwindo[YCMAX],
                        &vmin,&vmax,&Settings.zoom,
			&int_vars[NCOMP],&xref,&yref,&verbose );
        verbose = Settings.verbose;
        if (option==VISCMSH||option==VISCDBL||option==SEMSH||option==SEDBL){
            /*
             * draw colourbar for vmin->vmax
             */
            if (Settings.plot_opts.colour_bar!=CB_NONE) {
                if (Settings.plot_opts.colour_bar==CB_HORIZ) orient_vert=0;
                drawcolourbar_(&startcol,&endcol,
                        &barwdth,&barlen,&xlevel,&ylevel,
                        &pwindo[XCMIN],&pwindo[XCMAX],
                        &pwindo[YCMIN],&pwindo[YCMAX],
                        &orient_vert);
            /*
             * label max and min on the colour bar
             */
                if (Settings.plot_opts.label)
                    labelcolourbar_(&vmin,&vmax,&scale,&barlen,
                                    &xlevel,&ylevel,&orient_vert);
                setlinewidth_(&Settings.linewidth);
            }
        }
    }
    else {
    /*  amendments added for plotting of internal boundaries */
        if (arraysi[IMP]==NULL)
          if ((arraysi[IMP]=(int *)malloc(2*int_vars[NMP]*sizeof(int)))==NULL)
		                  return(MALLOC_ERR);
        if (arraysf[BNDS]==NULL)
          if ((arraysf[BNDS]=(float *)malloc(4*sizeof(float)))==NULL)
		                  return(MALLOC_ERR);
/*
 *   hardwired to plot only region 2.  A dialog box to input variable ireg
 *   to allow any region number to be plotted is needed.
 */
        ireg=2;  isvp=0;
        intbnd_(arraysi[IMP],arraysi[LEM],arraysi[NOR],
		arraysi[IMAT],arraysi[IBC],arraysi[IBNGH],&int_vars[NMP],
	        &int_vars[NE],&int_vars[NN],&int_vars[NUP],&int_vars[NBP],
	        arraysf[EX],arraysf[EY],arraysf[BNDS],arraysf[DENS],
	        &densreg,&isvp,&ireg,&int_vars[NCOMP],
                &xref,&yref,&verbose);
    }

    /* in case viscosity changed pen colour */
    if (option==VISCMSH||option==VISCDBL) setpencolor_(current_pen);
    if (reset) {
        pwindo[XCMIN]=xmin;
        pwindo[XCMAX]=xmax;
        pwindo[YCMIN]=ymin;
        pwindo[YCMAX]=ymax;
    }
    Save_plot_data( NULL,NULL );

    if (arraysi[LNOR]!=0) { free(arraysi[LNOR]); arraysi[LNOR]=0; } 
    if (arraysi[IMP]!=0) { free(arraysi[IMP]); arraysi[IMP]=0; } 
    if (arraysf[BNDS]!=0) { free(arraysf[BNDS]); arraysf[BNDS]=0; } 

    return(0);
}

int plot_velocity(int* int_vars, float* fl_vars, int** arraysi, float** arraysf,
                  char* name, int option, int plot_type, int plot_descrip, 
                  float* cntr_vals, float* profile_vals, float* pwindo,
                  int* mesh, int* current_pen, int* line)
{
    int y_direction;
    char plot_label[5];
    int i, format, len, meshnum=1;
    int namelen=0;
    float ss0,uarrowf,*u,*v, *data=0, *fptr;
    float xmid, ymid, xx, yy, xref, yref;
    float lower,upper;
    float vmin,vmax;
    float voffset = 0;
    double uarrow,sp,theta;
    int m1, n1, mp, np, labels, npts, *iptr;
    int startcol,rangecol;
    int err=0, verbose, mpe=1;

    if (Settings.plot_opts.dble) {
        meshnum=2;
        pwindo[XCMAX] += fl_vars[TBXOFF];
        pwindo[YCMAX] += fl_vars[TBYOFF];
        Settings.plot_opts.nx3 *= 2;
        NewMesh();
    }
    if (arraysf[SNTRP]==NULL)
        if ((arraysf[SNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))
                             ==NULL) return(MALLOC_ERR);
    if (arraysf[BNTRP]==NULL)
        if ((arraysf[BNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))
                             ==NULL) return(MALLOC_ERR);
    if (arraysi[IHELP]==NULL)
        if ((arraysi[IHELP]=(int *)malloc(mesh[NP3]*sizeof(int)))
                             ==NULL) return(MALLOC_ERR);
/*    XREFM and YREFM usage depends on context */
    xref=fl_vars[XREFM];
    yref=fl_vars[YREFM];
    if (int_vars[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }

    ss0=0.0;
    verbose = Settings.verbose;
    if (plot_type==ARROWS) {
        if(verbose){
            vprint_(arraysf[EX],arraysf[EY],arraysf[UVP],arraysi[NOR],
                    &int_vars[NCOMP],&int_vars[NUP],&int_vars[NN],
                    &xref);
        }
        data=u=arraysf[UVP];
        ntrplt_(data,arraysf[SNTRP],&ss0,
            &pwindo[XCMIN],&pwindo[YCMIN], &pwindo[XCMAX],&pwindo[YCMAX],
            &fl_vars[TBXOFF],&fl_vars[TBYOFF],&fl_vars[VELXO],
            arraysf[EX],arraysf[EY],
            arraysi[LEM],arraysi[NOR],arraysi[IHELP],
            &int_vars[NE],&int_vars[NUP], &mesh[NP3],
            &mesh[NX3], &mesh[NY3], &meshnum, &verbose
            );
        data=v=arraysf[UVP]+int_vars[NUP];
        ntrplt_(data,arraysf[BNTRP],&ss0,
            &pwindo[XCMIN],&pwindo[YCMIN], &pwindo[XCMAX],&pwindo[YCMAX],
            &fl_vars[TBXOFF],&fl_vars[TBYOFF],&fl_vars[VELYO],
            arraysf[EX],arraysf[EY],
            arraysi[LEM],arraysi[NOR],arraysi[IHELP],
            &int_vars[NE],&int_vars[NUP], &mesh[NP3],
            &mesh[NX3], &mesh[NY3], &meshnum, &verbose
            );

        for (i=0,uarrow=0.0;i<(int)int_vars[NUP];i++) {
            sp = sqrt((double)(u[i]*u[i] + v[i]*v[i]));
            uarrow = MAX(uarrow,sp);
        }
        uarrowf = (float)uarrow;
        m1=n1=2;
        mp = Settings.plot_opts.mp;
        np = Settings.plot_opts.np;
        strcpy(plot_label,"velo");
        verbose = Settings.verbose;
        arrow_(arraysf[SNTRP],arraysf[BNTRP],&uarrowf,&mesh[NXY],&m1,
            &mesh[NX2],&mp,&n1,&mesh[NY2],&np,&mesh[NX3],&mesh[NY2],
            &pwindo[XCMIN],&pwindo[XCMAX],&pwindo[YCMIN],&pwindo[YCMAX],
            &xref,&yref,&verbose,&int_vars[NCOMP]);
/* if verbose print arraysf[SNTRP], arraysf[BNTRP]  used to be in ARROW */
        if (verbose!=0) {
            err=DataPrintMeshVelo(plot_label,arraysf[SNTRP],arraysf[BNTRP],
                           m1,mesh[NX2],mp,n1,mesh[NY2],np,
                           mesh[NX3],mesh[NY3],pwindo,
                           xref,yref,int_vars[NCOMP],verbose);
        }
        labels= (Settings.plot_opts.dble  ? BNDRYDBL : BNDRY);
        startcol = Colour_range[0];
        rangecol = Colour_range[1] - startcol;
        /*
         * if zoom factor not set, draw the mesh boundary
         */
        if ((fabs(Settings.zoom-1.0)>1.e-4)) {
            verbose = 0;
            plmesh_(&labels,&int_vars[NE],&int_vars[NN],&int_vars[NUP],
                        &int_vars[NBP],arraysf[EX],arraysf[EY],
                        arraysi[LEM],arraysi[NOR],arraysi[IBC],
                        arraysi[IBCTYP],arraysf[VHB],
                        &fl_vars[TBXOFF],&fl_vars[TBYOFF],&int_vars[IFLT],
                        &mpe,&startcol,&rangecol,
                        &pwindo[XCMIN],&pwindo[XCMAX],&pwindo[YCMIN],&pwindo[YCMAX],
                        &vmin,&vmax,&Settings.zoom,&int_vars[NCOMP],&xref,&yref,
			&verbose );
            verbose = Settings.verbose;
        }
        if (Settings.plot_opts.label) {
            format = 63; len = strlen(plot_label);
            drawautolabels_(plot_label,&len,&uarrowf,&format,pwindo,current_pen);
        }
    }
    if (plot_type==CNTRS||plot_type==PRFL1D||plot_type==PRFL2D) {
        xmid = (pwindo[XCMAX] - pwindo[XCMIN])/2;
        ymid = (pwindo[YCMAX] - pwindo[YCMIN])/2;
        if ((option==UR||option==UTH||option==UM) && arraysf[AMESH]==NULL)
            if ((arraysf[AMESH]=(float *)
                malloc(int_vars[NUP]*sizeof(float)))==NULL) return(MALLOC_ERR);
        switch(option) {
        case UX  : data=arraysf[UVP];
                   voffset = fl_vars[VELXO];
                    break;
        case UY  : data=arraysf[UVP] + int_vars[NUP];
                   voffset = fl_vars[VELYO];
                    break;
        case UM  :  for (i=0,fptr=arraysf[AMESH],
                        u=arraysf[UVP],v=arraysf[UVP]+int_vars[NUP];
                               i<int_vars[NUP];i++,fptr++) {
                        *fptr = (float)sqrt((double)(u[i]*u[i]+v[i]*v[i]));
                    }
                    data = arraysf[AMESH];
                    voffset = fl_vars[VELYO];
                    break;
        case UR  : for (i=0,iptr=arraysi[NOR],fptr=arraysf[AMESH],
                        u=arraysf[UVP],v=arraysf[UVP]+int_vars[NUP];
                               i<int_vars[NUP];i++,iptr++,fptr++) {
                       xx = arraysf[EX][*iptr] - xmid;
                       yy = arraysf[EY][*iptr] - ymid;
                       if (xx==0.0 && yy==0.0) theta=0.0;
                       else theta = atan2((double)yy,(double)xx);
                       *fptr = u[i] * (float)cos(theta) + 
                                               v[i]*(float)sin(theta);
                    }
                    data = arraysf[AMESH];
                    break;
        case UTH : for (i=0,iptr=arraysi[NOR],fptr=arraysf[AMESH],
                        u=arraysf[UVP],v=arraysf[UVP]+int_vars[NUP];
                               i<int_vars[NUP];i++,iptr++,fptr++) {
                       xx = arraysf[EX][*iptr] - xmid;
                       yy = arraysf[EY][*iptr] - ymid;
                       if (xx==0.0 && yy==0.0) theta=0.0;
                       else theta = atan2((double)yy,(double)xx);
                       *fptr = v[i] * (float)cos(theta) - 
                                               u[i]*(float)sin(theta);
                    }
                    data = arraysf[AMESH];
                    break;
        default  :  break;
        }
        if (plot_type==CNTRS) {
            verbose=0;
            ntrplt_(data,arraysf[SNTRP],&ss0,
                        &pwindo[XCMIN],&pwindo[YCMIN],
                        &pwindo[XCMAX],&pwindo[YCMAX],
                        &fl_vars[TBXOFF],&fl_vars[TBYOFF],&voffset,
                        arraysf[EX],arraysf[EY],
                        arraysi[LEM],arraysi[NOR],arraysi[IHELP],
                        &int_vars[NE],&int_vars[NUP], &mesh[NP3],
                        &mesh[NX3], &mesh[NY3], &meshnum, &verbose
                        );
            verbose = Settings.verbose;
/* if verbose dataprint arraysf[SNTRP] */
            err = contour(int_vars,fl_vars,arraysi,arraysf,pwindo,mesh,
                    plot_descrip,cntr_vals,profile_vals,name,current_pen,
                    &Settings.verbose);
        }
        if (plot_type==PRFL1D) {
            namelen = strlen(name);
            npts=profile_vals[NUM_PTS] = Settings.plot_opts.profile_pts;
            verbose=0;
            ntrpln_(data,arraysf[SNTRP],&npts,&ss0,
                    &profile_vals[PRFL_X1],&profile_vals[PRFL_Y1],
                    &profile_vals[PRFL_X2],&profile_vals[PRFL_Y2],
                    &fl_vars[TBXOFF],&fl_vars[TBYOFF],&voffset,
                    arraysf[EX],arraysf[EY],&xref,&yref,
                    arraysi[LEM],arraysi[NOR],arraysi[IHELP],
                    &int_vars[NE],&int_vars[NUP], &mesh[NP3],
                    &int_vars[NCOMP],&verbose,name, &namelen
                    );
            verbose = Settings.verbose;
            err=Profile(arraysf[SNTRP],arraysi[IHELP],npts,profile_vals,
                int_vars[NCOMP],name,pwindo,current_pen,line);
        }
        if (plot_type==PRFL2D) {
            npts=profile_vals[NUM_PTS] = Settings.plot_opts.profile_pts;
            if (plot_descrip==DIM2_Y) {
                profile_vals[PRFL_X1] = 0.0;
                profile_vals[PRFL_Y1] = lower = profile_vals[LWR_LIM];
                profile_vals[PRFL_X2] = 0.0;
                profile_vals[PRFL_Y2] = upper = profile_vals[UPR_LIM];
                Plot_info.plot_description=DIM2_Y;
                y_direction = 1;
            }
            else /*plot_descrip==DIM2_X*/ {
                profile_vals[PRFL_X1] = lower = profile_vals[LWR_LIM];
                profile_vals[PRFL_Y1] = 0.0;
                profile_vals[PRFL_X2] = upper = profile_vals[UPR_LIM];
                profile_vals[PRFL_Y2] = 0.0;
                Plot_info.plot_description=DIM2_X;
                y_direction = 0;
            }
      /*    fprintf(stdout," calling profile_2D from plot_velocity \n");   */
            profile_2D(name,y_direction,arraysf,arraysi,data,arraysf[SNTRP],
                              mesh,npts,int_vars,fl_vars,&lower,&upper);
            err=Profile(arraysf[SNTRP],arraysi[IHELP],npts,profile_vals,
               int_vars[NCOMP],name,pwindo,current_pen,line);
        }
    }
    if (!err) Save_plot_data( cntr_vals,profile_vals );
    if (Settings.plot_opts.dble) {
        pwindo[XCMAX] -= fl_vars[TBXOFF];
        pwindo[YCMAX] -= fl_vars[TBYOFF];
        Settings.plot_opts.nx3 /= 2;
    }
    return((err==USER_CANCEL)? 0 : err);
}

int plot_gravity(int* int_vars, float* fl_vars, int** arraysi, float** arraysf,
                 float* pwindo, int* mesh, int plot_type, int plot_descrip, 
                 float* cntr_vals, float* profile_vals, char* label, int option,
                 int* current_pen, int* line)
{
    int i, meshnum=1, y_direction=0;
    float rmaxl,zmeas,proflen,zsurf,xmin,xmax,sfac;
    float densreg,zmin,zmax,gload,delr,xref,yref;
    int npts;
    int err=0, verbose, opt, len=0;
    int ireg, isvp, ir, ibout, ibin;
    float *profiltmp=0, *profiltmp2=0, *gravprof=0, *edgepts=0, *pts=0;
    int *help2=0;
/*
 *    the following parameters define the length scale, and density scale
 *    used to scale up the dimensionless stress, and the surface density
 *    which are required to convert dimensionless stress to surface deflection
 *    (hardwired here for the time being - maybe put in the options list later)
 */
    float dscale=1.0e+5,rhscale=30.0,sdens=2.6e+3;

    verbose = Settings.verbose;
    npts=profile_vals[NUM_PTS] = Settings.plot_opts.profile_pts;
    proflen = profile_vals[PRFL_X2]-profile_vals[PRFL_X1];
    rmaxl = pwindo[XCMAX];
    zmeas = profile_vals[PRFL_Y1];
/*    XREFM and YREFM usage depends on context */
    xref=fl_vars[XREFM];
    yref=fl_vars[YREFM];
    if (int_vars[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }

    if (arraysf[SNTRP]==NULL)
        if ((arraysf[SNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))
                             ==NULL) return(MALLOC_ERR);
    if (arraysi[IHELP]==NULL)
        if ((arraysi[IHELP]=(int *)malloc(mesh[NP3]*sizeof(int)))
                             ==NULL) return(MALLOC_ERR);
    if ((gravprof=(float *)malloc(npts*sizeof(float)))==NULL)
            return(MALLOC_ERR);
	for (i=0;i<npts;i++) gravprof[i]=0.0;
/*
 *  for gravity anomaly, compute the surface integrals around the material regions
 *  the gravity calculation requires material region numbers have been set in basil
 *  check that these region numbers are available. 
 */
    if(option==GRBOUG||option==GRFREE)
    {

        if (arraysi[IMP]==NULL)
          if ((arraysi[IMP]=(int *)malloc(2*int_vars[NMP]*sizeof(int)))==NULL)
		                  return(MALLOC_ERR);
        isvp=int_vars[NMP];
        if (arraysf[BNDS]==NULL)
          if ((arraysf[BNDS]=(float *)malloc(4*isvp*sizeof(float)))==NULL)
		                  return(MALLOC_ERR);
/*   
 *    assume for the time being that there are only two distinct
 *    internal regions
 *    may be modified by changing the number in the next line,
 *    to include more regions
 */
        ireg=8;
        gload=0.0;
        for (ir=1; ir< ireg+1; ir++) {
            /*
             *    first locate the internal boundary, defined by its segments,
             *    for each region
             */
            isvp=int_vars[NMP];
            intbnd_(arraysi[IMP],arraysi[LEM],arraysi[NOR],
		    arraysi[IMAT],arraysi[IBC],arraysi[IBNGH],&int_vars[NMP],
                    &int_vars[NE],&int_vars[NN],&int_vars[NUP],&int_vars[NBP],
                    arraysf[EX],arraysf[EY],arraysf[BNDS],arraysf[DENS],
                    &densreg,&isvp,&ir,&int_vars[NCOMP],&xref,&yref,&verbose);
            /*
             *    if there is at least one segment on the boundary,
             *    add in its contribution to the gravity
             *    also compute the total buoyant load (not now used to get
             *    a zero level on topography below)
             */
            if(isvp>=1 || ir==1){
                gravcp_(&int_vars[NCOMP],&ir,&densreg,&zmin,&zmax,
                         arraysf[BNDS],&isvp,gravprof,
                         arraysf[SNTRP],arraysi[IHELP],&npts,&rmaxl,
                        &proflen,&zmeas,&dscale,&rhscale,&err);
                gload=gload+(zmax-zmin)*densreg;
            }
            else{
                fprintf(stdout,"Warning: gravity calculation, region %i has %i segments\n",ir,isvp);
            }
        }
        if(err)return(err);
    }
    if(option==TOPOG||option==GRFREE){
    /*
     *    now calculate the stress on the upper surface (y=1 by default)
	 *    using stmesh
     */
    /*    if (arraysf[AMESH]==NULL)
            if ((arraysf[AMESH]=(float *)
                malloc(int_vars[NUP]*sizeof(float)))==NULL) return(MALLOC_ERR);
        if (arraysf[BMESH]==NULL)
            if ((arraysf[BMESH]=(float *)
                malloc(int_vars[NUP]*sizeof(float)))==NULL) return(MALLOC_ERR);
        if (arraysf[CMESH]==NULL)
            if ((arraysf[CMESH]=(float *)
                malloc(int_vars[NUP]*sizeof(float)))==NULL) return(MALLOC_ERR);
        if (arraysf[DMESH]==NULL)
            if ((arraysf[DMESH]=(float *)
                malloc(int_vars[NUP]*sizeof(float)))==NULL) return(MALLOC_ERR);
        if (arraysf[EMESH]==NULL)
            if ((arraysf[EMESH]=(float *)
                malloc(int_vars[NUP]*sizeof(float)))==NULL) return(MALLOC_ERR); */

        if (arraysf[BNTRP]==NULL)
            if ((arraysf[BNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))==NULL) 
                return(MALLOC_ERR);

        if ((help2=(int *)malloc(npts*sizeof(int)))==NULL)
            return(MALLOC_ERR);

	for (i=0;i<mesh[NP3];i++) arraysf[BNTRP][i]=0.0;
        /*
         *  ambiguity of profile level parameters exists, since input level 
         *  is used to define gravity measurement level,
         *  and can't simultaneously be used for defining topo level
         *  without modification of parameter input window.
         */
         /* zsurf=profile_vals[PRFL_Y1]; */
        zsurf=1.0;
        xmin=profile_vals[PRFL_X1];
        xmax=profile_vals[PRFL_X2];
        opt = SIYY;
        len = strlen(label);
        stmesh_(&opt,&plot_type,arraysf[SNTRP],arraysf[BNTRP],
                    arraysf[CNTRP],arraysf[DNTRP],
                    /* arraysf[AMESH],arraysf[BMESH],
                    arraysf[CMESH],arraysf[DMESH],arraysf[EMESH],*/
                    profiltmp,profiltmp2,edgepts,pts,
                    arraysf[EX],arraysf[EY],arraysf[VHB],arraysf[UVP],
                    arraysi[LEM],arraysi[NOR],arraysi[IHELP],arraysi[IBNGH],
                    arraysi[IBC],help2,arraysf[SSQ],arraysi[IELFIX],
                    int_vars,fl_vars,
                    &int_vars[NE],&int_vars[NUP],&int_vars[NUVP],
                    &int_vars[NFP],&int_vars[NN],&int_vars[NBP],
                    &mesh[NX3],&mesh[NY3],&mesh[NP3],&npts,
                    &y_direction,&meshnum,&xmin,&zsurf,&xmax,&zsurf,
                    &pwindo[XCMIN],&pwindo[YCMIN],&pwindo[XCMAX],&pwindo[YCMAX],
                    &verbose,label,&len
                    );
/*
 *   scale the dimensionless surface stress into surface deflection in m 
 *   a reference surface load is also subtracted so that surface level at 
 *   the extreme right of the solution is zero.  This is a bit arbitrary,
 *   but so is almost any other reference level.
 */
        gload=arraysf[SNTRP][npts-1];
        if(option==TOPOG) {
            sfac=-dscale*rhscale/sdens;
            for (i=0; i<npts; i++)
                   gravprof[i]=sfac*(arraysf[SNTRP][i] - gload);
        }
        else {
            sfac=-rhscale/sdens;
            for (i=0; i<npts; i++) 
                   arraysf[SNTRP][i]=sfac*(arraysf[SNTRP][i] - gload);
        }
    }
/*
 *   now calclate the free-air anomaly if required, by constructing a new region
 *   boundary that contains the topographic load
 */
    if(option==GRFREE){
        delr=rmaxl/(float)(npts-1);
        isvp=2*(npts-1);
        for (i=1,ibout=0,ibin=4; i<npts; i++,ibout=ibout+8,ibin=ibin+8) {
           arraysf[BNDS][ibout]=arraysf[BNDS][ibin+2]=delr*(i-1);
           arraysf[BNDS][ibout+1]=1.0;
           arraysf[BNDS][ibout+2]=arraysf[BNDS][ibin]=delr*i;
           arraysf[BNDS][ibout+3]=1.0;
           arraysf[BNDS][ibin+1]=zsurf+arraysf[SNTRP][i];
           arraysf[BNDS][ibin+3]=zsurf+arraysf[SNTRP][i-1];
       }
 /*     for (i=0,ibout=0; i<isvp; i++,ibout=ibout+4) {
           fprintf(stdout," i = %i x1,y1 = %f %f, x2,y2 = %f %f\n",i,
             arraysf[BNDS][ibout],arraysf[BNDS][ibout+1],arraysf[BNDS][ibout+2],arraysf[BNDS][ibout+3]); } */
/*
 *   the gravity contribution from the surface layer is accummulated
 *   into gravprof on top of the Bouguer anomaly calculated above,
 *   then copied into SNTRP for plotting
 */
        densreg=sdens/rhscale;
        ir=ireg+1;
        gravcp_(&int_vars[NCOMP],&ir,&densreg,&zmin,&zmax,
                  arraysf[BNDS],&isvp,gravprof,
                  arraysf[SNTRP],arraysi[IHELP],&npts,&rmaxl,
                  &proflen,&zmeas,&dscale,&rhscale,&err);
    }
	/*
	 *   for all options, now plot the profile stored in gravprof
	 */
    err=Profile(gravprof,arraysi[IHELP],npts,profile_vals,
                 int_vars[NCOMP],label,pwindo,current_pen,line);
    if (!err) Save_plot_data( NULL,profile_vals );
    if (gravprof!=0) { free(gravprof); gravprof=0; }
    if (help2!=0) { free(help2); help2=0; }
    return((err==USER_CANCEL)? 0 : err);
}

int plot_rotation(int* int_vars, float* fl_vars, int** arraysi, float** arraysf,
                  char* name, int option, int plot_type, int plot_descrip, 
                  float* cntr_vals, float* profile_vals, float* pwindo, int* mesh,
                  int* current_pen, int* line)
{
    char plot_label[5];
    int format, len, meshnum=1;
    float ss0,frmin,*data,rmax,baseval;
    float vmin,vmax,xref,yref;
    float value;
    int m1, n1, mp, np, labels;
    int startcol,rangecol;
    int err=0, verbose, mpe=1;

    if (arraysf[FROT]==NULL) return(FROT_ERR);
    if (Settings.plot_opts.dble) {
        meshnum=2;
        pwindo[XCMAX] += fl_vars[TBXOFF];
        pwindo[YCMAX] += fl_vars[TBYOFF];
        Settings.plot_opts.nx3 *= 2;
        NewMesh();
    }
    if (arraysf[SNTRP]==NULL)
        if ((arraysf[SNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))
                             ==NULL) return(MALLOC_ERR);
    if (arraysi[IHELP]==NULL)
        if ((arraysi[IHELP]=(int *)malloc(mesh[NP3]*sizeof(int)))
                             ==NULL) return(MALLOC_ERR);
/*    XREFM and YREFM usage depends on context */
    xref=fl_vars[XREFM];
    yref=fl_vars[YREFM];
    if (int_vars[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }

    ss0=0.0;
    frmin=5.0;
    rmax=1.0;
    verbose = Settings.verbose;
    if (plot_type==ARROWS) {
        data=arraysf[FROT];
        verbose=0;
        ntrplt_(data,arraysf[SNTRP],&ss0,
            &pwindo[XCMIN],&pwindo[YCMIN], &pwindo[XCMAX],&pwindo[YCMAX],
            &fl_vars[TBXOFF],&fl_vars[TBYOFF],&fl_vars[VELXO],
            arraysf[EX],arraysf[EY],
            arraysi[LEM],arraysi[NOR],arraysi[IHELP],
            &int_vars[NE],&int_vars[NUP], &mesh[NP3],
            &mesh[NX3], &mesh[NY3], &meshnum, &verbose
            );
        verbose = Settings.verbose;

        m1=n1=2;
        mp = Settings.plot_opts.mp;
        np = Settings.plot_opts.np;
        baseval = -90.0 + Settings.plot_opts.solution_rot*90.0;
        strcpy(plot_label,"FROT");
        rotor_(arraysf[SNTRP],&frmin,&baseval,&rmax,&mesh[NXY],
            &m1,&mesh[NX2],&mp,&n1,&mesh[NY2],&np,
            &mesh[NX3],
            &pwindo[XCMIN],&pwindo[XCMAX],
            &pwindo[YCMIN],&pwindo[YCMAX]);
        labels= (Settings.plot_opts.dble  ? BNDRYDBL : BNDRY);
        startcol = Colour_range[0];
        rangecol = Colour_range[1] - startcol;
        /*
         * if zoom factor not set, draw the mesh boundary
         */
        if ((fabs(Settings.zoom-1.0)>1.e-4)) {
            verbose = 0;
            plmesh_(&labels,&int_vars[NE],&int_vars[NN],&int_vars[NUP],
                        &int_vars[NBP],arraysf[EX],arraysf[EY],
                        arraysi[LEM],arraysi[NOR],arraysi[IBC],
                        arraysi[IBCTYP],arraysf[VHB],
                        &fl_vars[TBXOFF],&fl_vars[TBYOFF],&int_vars[IFLT],
                        &mpe,&startcol,&rangecol,
                        &pwindo[XCMIN],&pwindo[XCMAX],
                        &pwindo[YCMIN],&pwindo[YCMAX],
                        &vmin,&vmax,&Settings.zoom,
			&int_vars[NCOMP],&xref,&yref,&verbose );
            verbose = Settings.verbose;
        }
        value=0.0;
        if (Settings.plot_opts.label) {
            format = 63; len = strlen(plot_label);
            drawautolabels_(plot_label,&len,&value,&format,pwindo,current_pen);
        }
    }
    if (!err) Save_plot_data( cntr_vals,profile_vals );
    if (Settings.plot_opts.dble) {
        pwindo[XCMAX] -= fl_vars[TBXOFF];
        pwindo[YCMAX] -= fl_vars[TBYOFF];
        Settings.plot_opts.nx3 /= 2;
    }
    return((err==USER_CANCEL)? 0 : err);
}

int plot_strain(int* int_vars, float* fl_vars, int** arraysi, float** arraysf,
                float* pwindo, int* mesh, int option, float* cntr_vals,
                float* profile_vals, char* label, int plot_type, int plot_descrip,
                int* current_pen, int* line, int* colours)
{
    int opt,npts;
    int y_direction=0;
    int err=0, meshnum=1, verbose;
    int len=0;
    int *help2=0;
    float xmin, xmax, ymin, ymax;
    float *profiltmp=0, *profiltmp2=0, *edgepts=0, *pts=0;

    opt=(int)option;
    if (Settings.plot_opts.dble) {
        meshnum=2;
        pwindo[XCMAX] += fl_vars[TBXOFF];
        pwindo[YCMAX] += fl_vars[TBYOFF];
        Settings.plot_opts.nx3 *= 2;
        NewMesh();
    }
    if (arraysf[BNTRP]==NULL)
        if ((arraysf[BNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    if (arraysf[SNTRP]==NULL)
        if ((arraysf[SNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    if (arraysf[CNTRP]==NULL)
        if ((arraysf[CNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    if (arraysf[DNTRP]==NULL)
        if ((arraysf[DNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    if (arraysi[IHELP]==NULL)
        if ((arraysi[IHELP]=(int *)malloc(mesh[NP3]*sizeof(int)))==NULL)
            return(MALLOC_ERR);
    if (plot_type==PRFL1D) {
        npts=profile_vals[NUM_PTS] = Settings.plot_opts.profile_pts;
        if ((profiltmp=(float *)malloc(npts*sizeof(float)))==NULL)
            return(MALLOC_ERR);
        xmin = profile_vals[PRFL_X1];
        ymin = profile_vals[PRFL_Y1];
        xmax = profile_vals[PRFL_X2];
        ymax = profile_vals[PRFL_Y2];
    }
    else if (plot_type==PRFL2D) {
        npts=profile_vals[NUM_PTS] = Settings.plot_opts.profile_pts;

        if ((profiltmp=(float *)malloc(npts*sizeof(float)))==NULL)
            return(MALLOC_ERR);
        if ((profiltmp2=(float *)malloc(npts*sizeof(float)))==NULL)
            return(MALLOC_ERR);
        if ((edgepts=(float *)malloc(4*npts*sizeof(float)))==NULL)
            return(MALLOC_ERR);
        if ((pts=(float *)malloc(int_vars[NBP]*sizeof(float)))==NULL)
            return(MALLOC_ERR);

        if (plot_descrip==DIM2_Y) {
            profile_vals[PRFL_X1] = xmin = 0.0;
            profile_vals[PRFL_Y1] = ymin = profile_vals[LWR_LIM];
            profile_vals[PRFL_X2] = xmax = 0.0;
            profile_vals[PRFL_Y2] = ymax = profile_vals[UPR_LIM];
            Plot_info.plot_description=DIM2_Y;
            y_direction = 1;
        }
        else /*plot_descrip==DIM2_X*/ {
            profile_vals[PRFL_X1] = xmin = profile_vals[LWR_LIM];
            profile_vals[PRFL_Y1] = ymin = 0.0;
            profile_vals[PRFL_X2] = xmax = profile_vals[UPR_LIM];
            profile_vals[PRFL_Y2] = ymax = 0.0;
            Plot_info.plot_description=DIM2_X;
            y_direction = 0;
        }
    }
    else npts = mesh[NX3];
    if ((help2=(int *)malloc(npts*sizeof(int)))==NULL)
        return(MALLOC_ERR);
    verbose = Settings.verbose;
    len = strlen(label);
    verbose = Settings.verbose;
    stmesh_(&opt,&plot_type,arraysf[SNTRP],arraysf[BNTRP],
            arraysf[CNTRP],arraysf[DNTRP],
            profiltmp,profiltmp2,edgepts,pts,
            arraysf[EX],arraysf[EY],arraysf[VHB],arraysf[UVP],
            arraysi[LEM],arraysi[NOR],arraysi[IHELP],arraysi[IBNGH],
            arraysi[IBC],help2,arraysf[SSQ],arraysi[IELFIX],
            int_vars,fl_vars,
            &int_vars[NE],&int_vars[NUP],&int_vars[NUVP],
            &int_vars[NFP],&int_vars[NN],&int_vars[NBP],
            &mesh[NX3],&mesh[NY3],&mesh[NP3],&npts,
            &y_direction,&meshnum,&xmin,&ymin,&xmax,&ymax,
            &pwindo[XCMIN],&pwindo[YCMIN], &pwindo[XCMAX],&pwindo[YCMAX],
            &verbose,label,&len
            );
    verbose = Settings.verbose;
    free(help2);
    if (plot_type==PRFL2D) {
        free(pts);
        free(edgepts);
        free(profiltmp2);
    }
    if (plot_type==CNTRS) {
        err = contour(int_vars,fl_vars,arraysi,arraysf,pwindo,mesh,
                    plot_descrip,cntr_vals,profile_vals,label,current_pen,
                    &Settings.verbose);
    }
    if (plot_type==ARROWS) {
        plot_arrows(opt,label,arraysf,arraysi,fl_vars,int_vars,pwindo,mesh,
                                                current_pen,colours);
    }
    else if (plot_type==PRFL1D||plot_type==PRFL2D) {
        err = Profile(profiltmp,arraysi[IHELP],npts,profile_vals,
              int_vars[NCOMP],label,pwindo, current_pen,line);
        free(profiltmp);
    }
    if (!err) Save_plot_data( cntr_vals,profile_vals );
    if (Settings.plot_opts.dble) {
        pwindo[XCMAX] -= fl_vars[TBXOFF];
        pwindo[YCMAX] -= fl_vars[TBYOFF];
        Settings.plot_opts.nx3 /= 2;
    }
    return((err==USER_CANCEL)? 0 : err);
}

int plot_density(int* int_vars, float* fl_vars, int** arraysi, float** arraysf, 
                 float* pwindo, int* mesh, char* name, int plot_type, 
                 int plot_descrip, float* cntr_vals, float* profile_vals,
                 int* current_pen, int* line)
{
    int y_direction;
    int npts,npts1,meshnum=1;
//  int npts2;
    float lower,upper;
    float bgdens;
    float *data;
    float xmin, xmax, ymin, ymax;
    int err = 0, verbose;

    verbose = Settings.verbose;
    if (arraysf[DENS]==NULL) return(DENS_ERR);
    if (arraysf[SNTRP]==NULL)
        if ((arraysf[SNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    if (arraysi[IHELP]==NULL)
        if ((arraysi[IHELP]=(int *)malloc(mesh[NP3]*sizeof(int)))==NULL)
            return(MALLOC_ERR);

    bgdens = 0.0;
    if (Settings.plot_opts.dble) {
        meshnum=2;
        pwindo[XCMAX] += fl_vars[TBXOFF];
        pwindo[YCMAX] += fl_vars[TBYOFF];
        Settings.plot_opts.nx3 *= 2;
        NewMesh();
    }
    if (plot_type==CNTRS ) {
       
        verbose = 0;
        ntrpld_(arraysf[DENS],arraysf[SNTRP],&bgdens,
            &pwindo[XCMIN],&pwindo[YCMIN],&pwindo[XCMAX],&pwindo[YCMAX],
            &fl_vars[TBXOFF],&fl_vars[TBYOFF],&fl_vars[VELXO],
            arraysf[EX],arraysf[EY],
            arraysi[LEM],arraysi[NOR],arraysi[IHELP],
            &int_vars[NE],&int_vars[NUP], &mesh[NP3],
            &mesh[NX3], &mesh[NY3], &meshnum, &verbose);
            npts1=mesh[NX3]; 
//	    npts2=mesh[NY3];
        verbose = Settings.verbose;
        err = contour(int_vars,fl_vars,arraysi,arraysf,pwindo,mesh,
                          plot_descrip,cntr_vals,profile_vals,name,current_pen,
                          &Settings.verbose);
    }
    if (plot_type==PRFL1D) {
        npts1=profile_vals[NUM_PTS] = Settings.plot_opts.profile_pts;
        xmin = profile_vals[PRFL_X1];
        ymin = profile_vals[PRFL_Y1];
        xmax = profile_vals[PRFL_X2];
        ymax = profile_vals[PRFL_Y2];
  /*    fprintf(stdout,"plot_density, xmin=%f ymin=%f, xmax=%f, ymax=%f\n",
         xmin,ymin,xmax,ymax);  */
        verbose = 0;
        ntrplnd_(arraysf[DENS],arraysf[SNTRP],&npts1,&bgdens,
            &xmin,&ymin,&xmax,&ymax,
            &fl_vars[TBXOFF],&fl_vars[TBYOFF],&fl_vars[VELXO],
            arraysf[EX],arraysf[EY],
            arraysi[LEM],arraysi[NOR],arraysi[IHELP],
            &int_vars[NE],&int_vars[NUP], &mesh[NP3], &verbose);
        verbose = Settings.verbose;
        err = Profile(arraysf[SNTRP],arraysi[IHELP],npts1,profile_vals,
              int_vars[NCOMP],name,pwindo, current_pen,line);
    }
    if (plot_type==PRFL2D) {
        npts=profile_vals[NUM_PTS] = Settings.plot_opts.profile_pts;
        if (plot_descrip==DIM2_Y) {
            profile_vals[PRFL_X1] = 0.0;
            profile_vals[PRFL_Y1] = lower = profile_vals[LWR_LIM];
            profile_vals[PRFL_X2] = 0.0;
            profile_vals[PRFL_Y2] = upper = profile_vals[UPR_LIM];
            Plot_info.plot_description=DIM2_Y;
            y_direction = 1;
        }
        else /*plot_descrip==DIM2_X*/ {
            profile_vals[PRFL_X1] = lower = profile_vals[LWR_LIM];
            profile_vals[PRFL_Y1] = 0.0;
            profile_vals[PRFL_X2] = upper = profile_vals[UPR_LIM];
            profile_vals[PRFL_Y2] = 0.0;
            Plot_info.plot_description=DIM2_X;
            y_direction = 0;
        }
  /*    fprintf(stdout," calling profile_2D from plot_density \n");   */
  /*    fprintf(stdout,"plot_density, name = %s lower = %f, upper = %f\n",
                                      name, lower, upper);  */
        data=arraysf[DENS];
        profile_2D(name,y_direction,arraysf,arraysi,data,arraysf[SNTRP],
                              mesh,npts,int_vars,fl_vars,&lower,&upper);
        err=Profile(arraysf[SNTRP],arraysi[IHELP],npts,profile_vals,
               int_vars[NCOMP],name,pwindo, current_pen,line);
    }
    if (!err) Save_plot_data( cntr_vals,profile_vals );
    if (Settings.plot_opts.dble) {
        pwindo[XCMAX] -= fl_vars[TBXOFF];
        pwindo[YCMAX] -= fl_vars[TBYOFF];
        Settings.plot_opts.nx3 /= 2;
    }
    return((err==USER_CANCEL)? 0 : err);
}

int plot_layer(int* int_vars, float* fl_vars, int** arraysi, float** arraysf, 
               float* pwindo, int* mesh, char* name, int plot_type, 
               int plot_descrip, float* cntr_vals, float* profile_vals, 
               int topo, int* current_pen, int* line)
{
    int y_direction;
    char namex[5];
    int npts1,npts2,meshnum=1;
    float logss0,scale_factor,zero_shift;
    float upper,lower,xref,yref;
//  float thickm;
    int i, err = 0, verbose;
    int len=0;

    strcpy(namex,"");
    if (arraysf[SSQ]==NULL) return(SSQ_ERR);
    if (topo && arraysi[IELFIX]==NULL) return(ELFX_ERR);
    if (arraysf[SNTRP]==NULL)
        if ((arraysf[SNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    if (arraysi[IHELP]==NULL)
        if ((arraysi[IHELP]=(int *)malloc(mesh[NP3]*sizeof(int)))==NULL)
            return(MALLOC_ERR);
    if (arraysf[AMESH]==NULL)
        if ((arraysf[AMESH]=(float *)
            malloc(int_vars[NUP]*sizeof(float)))==NULL) return(MALLOC_ERR);
    if (arraysf[BMESH]==NULL)
        if ((arraysf[BMESH]=(float *)
            malloc(int_vars[NUP]*sizeof(float)))==NULL) return(MALLOC_ERR);
/*    XREFM and YREFM usage depends on context */
    xref=fl_vars[XREFM];
    yref=fl_vars[YREFM];
    if (int_vars[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }

    logss0 = -(float)log((double)fl_vars[HLENSC]);
    if(fl_vars[BDEPSC] == 0.0)  fl_vars[BDEPSC] = 0.35;
    if (cntr_vals[CNTR_LVL]==0.0)  cntr_vals[CNTR_LVL]=fl_vars[BDEPSC];
    if (Settings.plot_opts.dble) {
        meshnum=2;
        pwindo[XCMAX] += fl_vars[TBXOFF];
        pwindo[YCMAX] += fl_vars[TBYOFF];
        Settings.plot_opts.nx3 *= 2;
        NewMesh();
    }
    /*
     * crustal thickness (topo = THICKNESS);
     * Gravitational potential energy (topo=GRAVPE);
     * elevation (topo=TOPOG)
     */

//  fprintf(stdout,"in plot_layer, BDEPSC = %f, HLENSC = %f, RISOST = %f, REFLEV = %f\n",
//                       fl_vars[BDEPSC],fl_vars[HLENSC],fl_vars[RISOST],fl_vars[REFLEV]);

/*   the following segment used to take account of convective thinning - needs rethink and testing

          tl0 = 100.0;
          zero_shift = -tl0*Settings.physical.bgam0/fl_vars[HLENSC];
          lnbgam0 = log(tl0*Settings.physical.bgam0);
          lnbgam1 = log(tl0*Settings.physical.bgam1);
          for (i=0; i < int_vars[NUP]; i++) {
              arraysf[BMESH][i] = lnbgam0 + arraysf[SSQ][i];
              if(arraysi[IELFIX][i]) arraysf[BMESH][i] = arraysf[BMESH][i]+
                                                       lnbgam1 - lnbgam0;
        }   */

/*  print SSQ by node if required, rescaling to TOPO or GRAVPE should be done externally 
    verbose=0;
      fprintf(stdout, "ln(SSQ) values for vertex nodes follow\n");
      lprint_(arraysf[EX],arraysf[EY],arraysf[SSQ],arraysi[NOR],
                &int_vars[NCOMP],&int_vars[NUP],&int_vars[NN],
                &xref);
  */
    if (plot_type==CNTRS) {
        verbose = 0;
        ntrplt_(arraysf[SSQ],arraysf[SNTRP],&logss0,
            &pwindo[XCMIN],&pwindo[YCMIN], &pwindo[XCMAX],&pwindo[YCMAX],
            &fl_vars[TBXOFF],&fl_vars[TBYOFF],&fl_vars[VELXO],
            arraysf[EX],arraysf[EY],
            arraysi[LEM],arraysi[NOR],arraysi[IHELP],
            &int_vars[NE],&int_vars[NUP], &mesh[NP3],
            &mesh[NX3], &mesh[NY3], &meshnum, &verbose
            );
        verbose = Settings.verbose;
        npts1=mesh[NX3]; npts2 = mesh[NY3];
    }
    if (plot_type==PRFL1D) {
        npts1=profile_vals[NUM_PTS] = Settings.plot_opts.profile_pts;
        npts2 = 1;
        len = strlen(name);
        verbose = 0;
        ntrpln_(arraysf[SSQ],arraysf[SNTRP],&npts1,&logss0,
            &profile_vals[PRFL_X1],&profile_vals[PRFL_Y1],
            &profile_vals[PRFL_X2],&profile_vals[PRFL_Y2],
            &fl_vars[TBXOFF],&fl_vars[TBYOFF],&fl_vars[VELXO],
            arraysf[EX],arraysf[EY],&xref,&yref,
            arraysi[LEM],arraysi[NOR],arraysi[IHELP],
            &int_vars[NE],&int_vars[NUP], &mesh[NP3],
            &int_vars[NCOMP],&verbose,name,&len
            );
        verbose = Settings.verbose;
    }
    else if (plot_type==PRFL2D) {
        npts1=profile_vals[NUM_PTS] = Settings.plot_opts.profile_pts;
        npts2 = 1;
        if (plot_descrip==DIM2_Y) {
            profile_vals[PRFL_X1] = 0.0;
            profile_vals[PRFL_Y1] = lower = profile_vals[LWR_LIM];
            profile_vals[PRFL_X2] = 0.0;
            profile_vals[PRFL_Y2] = upper = profile_vals[UPR_LIM];
            Plot_info.plot_description=DIM2_Y;
            y_direction = 1;
        }
        else /*plot_descrip==DIM2_X*/ {
            profile_vals[PRFL_X1] = lower = profile_vals[LWR_LIM];
            profile_vals[PRFL_Y1] = 0.0;
            profile_vals[PRFL_X2] = upper = profile_vals[UPR_LIM];
            profile_vals[PRFL_Y2] = 0.0;
            Plot_info.plot_description=DIM2_X;
            y_direction = 0;
        }
  /*    fprintf(stdout," calling profile_2D from plot_layer \n");   */
        profile_2D(name,y_direction,arraysf,arraysi,arraysf[SSQ],arraysf[SNTRP],
                               mesh,npts1,int_vars,fl_vars,&lower,&upper);
    }

/*  at this stage the array to be plotted is still in the form of log(thickness) */

//  thickm=100.0;         /* nominal lithosphere thickness in km */
    scale_factor = fl_vars[HLENSC]*fl_vars[BDEPSC];

    if(topo==THICKNESS) {
/*    for units of layer thickness, for km uncomment the following line */
/*    scale_factor=scale_factor*thickm;   */
      zero_shift = 0.0;
    }
    else if (topo==TOPOG) {
/*    for units in km, assumes lithosphere thickness = 100 km (included in RISOST factor) */
/*    in current state, RISOST = (rhom/(rhmo-rhoc))/thickm;  and REFLEV in m;  These should 
      be later changed for self-consistency    */
      scale_factor=scale_factor/fl_vars[RISOST];
      zero_shift=0.001*fl_vars[REFLEV] - fl_vars[BDEPSC]/fl_vars[RISOST];
    }
/*    GRAVPE is proportional to crustal thickness squared  */
    else if (topo==GRAVPE) {
        for (i=0; i < npts1*npts2; i++) {
            arraysf[SNTRP][i] = 2.0*arraysf[SNTRP][i];
        }
        scale_factor = 0.5*fl_vars[ARGANP]*fl_vars[BDEPSC]*fl_vars[BDEPSC]*
                       (pow(fl_vars[HLENSC],(-fl_vars[SE])));
        zero_shift = -scale_factor/(fl_vars[HLENSC]*fl_vars[HLENSC]);
/*   If convective thinning is required, need to replace ARGANP with BRGANP
     according to IELFIX: not implemented yet    */
    }

/*  now convert interpolated array to linear scale with appropriate units  */
    snorm_(arraysf[SNTRP],&npts1,&npts2,&scale_factor,&zero_shift);

    if (plot_type==CNTRS) {
        err = contour(int_vars,fl_vars,arraysi,arraysf,pwindo,mesh,
                plot_descrip,cntr_vals,profile_vals,name,current_pen,
                &Settings.verbose);
    }
    else if (plot_type==PRFL1D) {
        err = Profile(arraysf[SNTRP],arraysi[IHELP],npts1,profile_vals,
                      int_vars[NCOMP],name,pwindo,
                      current_pen,line);
    }
    if (!err) Save_plot_data( cntr_vals,profile_vals );
    if (Settings.plot_opts.dble) {
        pwindo[XCMAX] -= fl_vars[TBXOFF];
        pwindo[YCMAX] -= fl_vars[TBYOFF];
        Settings.plot_opts.nx3 /= 2;
    }
    return((err==USER_CANCEL)? 0 : err);
}

int plot_deformation(int* int_vars, float* fl_vars, int** arraysi, float** arraysf, 
                     float* pwindo, int* mesh, int option, int plot_descrip,
                     float* cntr_vals, float* profile_vals, char* name, int* current_pen)
{
    int format,len;
    int opt,mpe,indx_offset;
    int err=0, verbose;
//  int meshnum=1;
    float value;

    if (int_vars[IMSH]==3) return(MSHTYP_ERR);
    if (arraysf[AMESH]==NULL)
      if ((arraysf[AMESH]=(float *)malloc(int_vars[NUP]*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    if (Settings.plot_opts.dble) {
//      meshnum=2;
        pwindo[XCMAX] += fl_vars[TBXOFF];
        pwindo[YCMAX] += fl_vars[TBYOFF];
        Settings.plot_opts.nx3 *= 2;
        NewMesh();
    }
    opt = (int)option;
    verbose = Settings.verbose;
    if (option!=ROTA) {
        /*mpe = (option>2?1:6);*/
        mpe = Settings.plot_opts.mpe;
        verbose=0;
        meshdf_(&opt,&mpe,name,arraysf[EX],arraysf[EY],arraysf[EXREF],
                arraysf[EYREF],arraysi[LEM],arraysi[NOR],arraysf[AMESH],
                &int_vars[NE],&int_vars[NUP],&int_vars[NY],
                &pwindo[YCMIN], &pwindo[YCMAX], &verbose
                );
        verbose = Settings.verbose;
        if (Settings.plot_opts.label) {
            value = 0.0; format = 0; len = strlen(name);
            drawautolabels_(name,&len,&value,&format,pwindo,current_pen);
        }
    }
    else {
        if (arraysf[SNTRP]==NULL)
            if ((arraysf[SNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))
                ==NULL) return(MALLOC_ERR);
        indx_offset = 2;
        cntr_vals[CNTR_LVL]=0;
        fdput_(&indx_offset,arraysf[SNTRP],arraysf[AMESH],
            arraysf[EX],arraysf[EY],arraysi[NOR],arraysi[LEM],
            &int_vars[NE],&int_vars[NUP],&mesh[NP3],
            &mesh[NX0],&mesh[NX2],&mesh[NX3],
            &mesh[NY2],&mesh[NY3]
            );
        err = contour(int_vars,fl_vars,arraysi,arraysf,pwindo,mesh,
                plot_descrip,cntr_vals,profile_vals,name,current_pen,
                &Settings.verbose);
    }
    if (!err) Save_plot_data( cntr_vals,NULL );
    if (Settings.plot_opts.dble) {
        pwindo[XCMAX] -= fl_vars[TBXOFF];
        pwindo[YCMAX] -= fl_vars[TBYOFF];
        Settings.plot_opts.nx3 /= 2;
    }
    return((err==USER_CANCEL)? 0 : err);
}

int plot_arrows(opt,label,arraysf,arraysi,fl_vars,int_vars,pwindo,mesh,
                            current_pen,colours)
float *fl_vars,**arraysf;
int *int_vars, **arraysi;
int opt;
char *label;
float *pwindo;
int *mesh;
int *current_pen,*colours;
{
    int i_smin,j_smin,i_smax,j_smax,i_bmin,j_bmin,i_bmax,j_bmax;
    int m1,n1,mp,np,cols[2], format, len;
    int verbose;
    int err=0;
    float smin,smax,bmin,bmax,rmin,rmax,xref,yref,sumint;
    double maxs,mins;

/*    XREFM and YREFM usage depends on context */
    xref=fl_vars[XREFM];
    yref=fl_vars[YREFM];
    if (int_vars[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }
    verbose = Settings.verbose;
    m1=n1=2;
    mp = Settings.plot_opts.mp;
    np = Settings.plot_opts.np;
    rangxy_(arraysf[SNTRP],arraysi[IHELP],&mesh[NXY],&m1,&mesh[NX2],&mp,
            &n1,&mesh[NY2],&np,&smax,&i_smax,&j_smax,&smin,
            &i_smin,&j_smin,&mesh[NX3],&sumint);
    rangxy_(arraysf[BNTRP],arraysi[IHELP],&mesh[NXY],&m1,&mesh[NX2],&mp,
            &n1,&mesh[NY2],&np,&bmax,&i_bmax,&j_bmax,&bmin,
            &i_bmin,&j_bmin,&mesh[NX3],&sumint);

    cols[0] = colours[0];
    cols[1] = colours[1];
    maxs = MAX((double)bmax,(double)smax);
    mins = MIN((double)bmin,(double)smin);
    rmax = (float)(MAX(maxs,mins));
    rmin = (float)(MIN(maxs,mins));
    verbose = 0;
    arwtwo_(arraysf[CNTRP],arraysf[DNTRP],arraysf[SNTRP],arraysf[BNTRP],
            &rmax,&rmin,&mesh[NXY],&m1,&mesh[NX2],
            &mp,&n1,&mesh[NY2],&np,&mesh[NP3],
            &mesh[NY2],&mesh[NX3],cols,
            &pwindo[XCMIN],&pwindo[YCMIN],&pwindo[XCMAX],&pwindo[YCMAX],
            &xref,&yref,&int_vars[NCOMP],label,&verbose
            );
    verbose = Settings.verbose;
    if (verbose!=0) {
        err = DataPrintArrow(label,arraysf[CNTRP],arraysf[DNTRP],
                             arraysf[SNTRP], arraysf[BNTRP],
                             m1,mesh[NX2],mp,n1,mesh[NY2],np,mesh[NX3],
                             pwindo,xref,yref,int_vars[NCOMP],verbose);
    }
    if (Settings.plot_opts.label) {
        format = 63; len = strlen(label);
        drawautolabels_(label,&len,&rmax,&format,pwindo,current_pen);
    }
//  return(0);
    return(err);
}

int profile_2D( char* name, int ydir, float** arraysf, int** arraysi,
                float* data, float* res, int* mesh, int npts, int* int_vars,
                float* fl_vars, float* lower, float* upper)
{
    int i, ix1, ix2, iy1, iy2, len=0;
    int *help2;
    int verbose=0;
//  int meshnum=1;
    float *profiltmp,*edgepts,*pts,*xvals,*yvals;
    float x, y, logss0, area, hlen, xref, yref;
    if ((profiltmp=(float *)malloc(npts*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    if ((edgepts=(float *)malloc(4*npts*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    if ((pts=(float *)malloc(int_vars[NBP]*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    if ((help2=(int *)malloc(npts*sizeof(int)))==NULL)
            return(MALLOC_ERR);
/*    XREFM and YREFM usage depends on context */
    xref=fl_vars[XREFM];
    yref=fl_vars[YREFM];
    if (int_vars[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }
    if (ydir) {
        xvals = arraysf[EX];
        yvals = arraysf[EY];
        ix1=0;
        ix2=2;
        iy1=1;
        iy2=3;
    }
    else {
        xvals = arraysf[EY];
        yvals = arraysf[EX];
        ix1=1;
        ix2=3;
        iy1=0;
        iy2=2;
    }
    fprintf(stdout,"profile_2D: lower =%f upper = %f ydir = %i\n",*lower,*upper,ydir);
    edges_(edgepts,lower,upper,xvals,yvals,pts,
           arraysi[IBNGH],arraysi[IBC], arraysi[NOR],
           &int_vars[NUP],&int_vars[NE],&int_vars[NBP],&npts);
    logss0 = 0.0;
//  if (Settings.plot_opts.dble) meshnum=2;
    for (i=0;i<npts;i++) {
  /*    if (arraysi[IHELP][i]) {  
          this array may be modified by call to ntrpln or ntrplnd below 
          so if condition is switched out - blank lines now flagged
          by two end points being identical            
  */
  /*    fprintf(stdout,"profile_2D: i,edgepts = %i %f %f %f %f\n",i,&edgepts[ix1+4*i],
		&edgepts[iy1+4*i],&edgepts[ix2+4*i],&edgepts[iy2+4*i]);  */
            /* If name is not Density  */
            if (strcmp(name, "Density")) {
               len=strlen(name);
               verbose = 0;
               ntrpln_(data,profiltmp,&npts,&logss0,
                  &edgepts[ix1+4*i],&edgepts[iy1+4*i],
                  &edgepts[ix2+4*i],&edgepts[iy2+4*i],
                  &fl_vars[TBXOFF],&fl_vars[TBYOFF],&fl_vars[VELXO],
                  arraysf[EX],arraysf[EY],&xref,&yref,
                  arraysi[LEM],arraysi[NOR],arraysi[IHELP],
                  &int_vars[NE],&int_vars[NUP], &mesh[NP3],
                  &int_vars[NCOMP],&verbose,name,&len
               );  
               verbose = Settings.verbose;
            }
            /* If name is Density  */
            else {
               verbose = 0;
               ntrplnd_(data,profiltmp,&npts,&logss0,
                  &edgepts[ix1+4*i],&edgepts[iy1+4*i],
                  &edgepts[ix2+4*i],&edgepts[iy2+4*i],
                  &fl_vars[TBXOFF],&fl_vars[TBYOFF],&fl_vars[VELXO],
                  arraysf[EX],arraysf[EY],
                  arraysi[LEM],arraysi[NOR],arraysi[IHELP],
                  &int_vars[NE],&int_vars[NUP], &mesh[NP3], &verbose);
               verbose = Settings.verbose;
            }
            x = edgepts[ix2+4*i] - edgepts[ix1+4*i];
            y = edgepts[iy2+4*i] - edgepts[iy1+4*i];
            hlen = (float)sqrt((double)(x*x + y*y));
            integrattrpz_(profiltmp,arraysi[IHELP],&npts,&hlen,
                                                     &area,NULL,&len);
            res[i] = area;
   /*   } */
    }
    free(help2);
    free(pts);
    free(edgepts);
    free(profiltmp);
    return(0);
}

int contour(int_vars,fl_vars,arraysi,arraysf,pwindo,mesh,
            plot_descrip,cntr_vals,prfl_vals,label,current_pen,
            verbose)
int *int_vars, **arraysi;
float *fl_vars, **arraysf;
float *pwindo;
int *mesh, plot_descrip;
float *cntr_vals, *prfl_vals;
char *label;
int *current_pen, *verbose;
{
 /* char number[MAX_NUM_LEN];  */
    char dfltlabels[2][MAXNAME];
    int old_pen,format,len;
    float hmesh,vmesh,xref,yref,sumint;
    float *wkarea,xlevel,ylevel,barwdth,barlen,con[3];
    float vmin,vmax,zcmin,zcmax,interval;
    int m1,n1,nnx,nny,startcol,endcol;
    int num_cntrs,labels;
    int rangecol, orient_vert=1, mpe=1;
    int err=0,verbound=0;     // verbound = verbose flag for boundary

    if (arraysf[CNTRP]==NULL)
        if ((arraysf[CNTRP]=(float *)malloc(mesh[NP3]*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    if (arraysi[IHELP]==NULL)
        if ((arraysi[IHELP]=(int *)malloc(mesh[NP3]*sizeof(int)))==NULL)
            return(MALLOC_ERR);
    if (arraysi[IWORK]==NULL)
        if ((arraysi[IWORK]=(int *)malloc(mesh[NP3]*sizeof(int)))==NULL)
            return(MALLOC_ERR);
    m1=n1=2;
//  if (Settings.plot_opts.dble) meshnum=2;
    old_pen = *current_pen;
/*    XREFM and YREFM usage depends on context */
    xref=fl_vars[XREFM];
    yref=fl_vars[YREFM];
    if (int_vars[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }

    cntour_(arraysf[SNTRP],arraysf[CNTRP],arraysi[IWORK],arraysi[IHELP],
            &mesh[NXY],&m1,&mesh[NX2],&n1,&mesh[NY2],&mesh[NX3],
            &nnx,&nny,&hmesh,&vmesh,&zcmin,&zcmax,
            &pwindo[XCMIN],&pwindo[XCMAX],&pwindo[YCMIN],&pwindo[YCMAX],
            &sumint);
    fprintf(stdout,"%s range: Min = %g  Max = %g  Int = %g\n",
                                            label,zcmin,zcmax,sumint);
    startcol=Colour_range[0]; endcol=Colour_range[1];
    rangecol = endcol - startcol;

/*  propose contour min, max, step and level based on present zcmin, zcmax */

    ticstep(zcmin,zcmax,&interval,&cntr_vals[CNTR_MIN],&cntr_vals[CNTR_MAX]);
    cntr_vals[CNTR_STEP]=interval;
    if (cntr_vals[CNTR_MIN]*cntr_vals[CNTR_MAX] <= 0)cntr_vals[CNTR_LVL] = 0.0;
    else cntr_vals[CNTR_LVL] = ((int) (cntr_vals[CNTR_MIN]/interval))* interval;

  /*    cntr_vals[CNTR_STEP]=Calc_Step(cntr_vals[CNTR_MIN],
                                cntr_vals[CNTR_MAX], cntr_vals[NUM_CNTRS]);  
        cntr_vals[CNTR_LVL]=Calc_Level(cntr_vals[CNTR_MIN],
                                cntr_vals[CNTR_MAX], cntr_vals[CNTR_STEP]); */

    if (cntr_vals[CNTR_MIN_U]==0.0 && cntr_vals[CNTR_MAX_U]==0.0) {
        cntr_vals[CNTR_MIN_U] = cntr_vals[CNTR_MIN];
        cntr_vals[CNTR_MAX_U] = cntr_vals[CNTR_MAX];
        cntr_vals[CNTR_LVL_U] = cntr_vals[CNTR_LVL];
    }
    if (cntr_vals[CNTR_STEP_U]==0.0) cntr_vals[CNTR_STEP_U] = cntr_vals[CNTR_STEP];
/*
    strncpy(&dfltlabels[0][0],label,MAXNAME);
    formatnumber_(number,&cntr_vals[CNTR_STEP],&format);
    strncpy(&dfltlabels[1][0],number,MAXNAME);
*/
/*     allow user to adjust values or cancel 
       The variables that end with _U (e.g CNTR_MIN_U) are confirmed
       by GetUserVals, then used below. The values without _U are just
       the recommendation for the plot, chosen by ticset */

    if (GetUserVals( cntr_vals,CNTRVALS,0,
                    NULL, (char*) dfltlabels )==USER_CANCEL) return(USER_CANCEL);

/*      make the colour plot       */

    if (plot_descrip&SHADE) {
        c3code_(arraysf[CNTRP],arraysi[IWORK],&mesh[NY2],&mesh[NX2],
            &nny,&nnx,&startcol,&endcol,&hmesh,&vmesh,
            &cntr_vals[CNTR_MIN_U],&cntr_vals[CNTR_MAX_U],
            &pwindo[XCMIN],&pwindo[XCMAX],
            &pwindo[YCMIN],&pwindo[YCMAX]
            );
        /*
         * draw the the colour bar
         */
        if (Settings.plot_opts.colour_bar!=CB_NONE) {
            if (Settings.plot_opts.colour_bar==CB_HORIZ) orient_vert=0;
            drawcolourbar_(&startcol,&endcol,
                &barwdth,&barlen,&xlevel,&ylevel,
                &pwindo[XCMIN],&pwindo[XCMAX],
                &pwindo[YCMIN],&pwindo[YCMAX],
                &orient_vert);
        /*
         * label max and min on the colour bar
         */
            if (Settings.plot_opts.label)
                labelcolourbar_(&cntr_vals[CNTR_MIN_U],&cntr_vals[CNTR_MAX_U],
                       &cntr_vals[CNTR_SCL],&barlen,&xlevel,&ylevel,
                       &orient_vert);
        }
        setlinewidth_(&Settings.linewidth);
    }

    if (plot_descrip&LINES) {
        if (plot_descrip&SHADE) {
               /* black contour lines over screen colour shading */
            *current_pen=BLACK_PEN;
            setpencolor_(current_pen);
        }
        if ((wkarea=(float *)malloc(mesh[ICON]*sizeof(float)))==NULL)
            return(MALLOC_ERR);
        num_cntrs = (int)cntr_vals[MAX_CNTRS];
        con[0] = cntr_vals[CNTR_LVL_U];
        con[1] = cntr_vals[CNTR_STEP_U];
  /*    con[2] = cntr_vals[MAX_CNTRS];
  fprintf(stdout," LVL %f STEP %f MAX %f\n",con[0],con[1],con[2]); */
        linecntour_(arraysf[SNTRP],arraysf[CNTRP],arraysi[IWORK],
            &mesh[ICON], &mesh[NX2],&mesh[NY2],
            &nnx,&nny, &hmesh,&vmesh,
            &cntr_vals[CNTR_MIN_U],&cntr_vals[CNTR_MAX_U],
            &pwindo[XCMIN], &pwindo[YCMIN], &Settings.plot_opts.stipple,
            &num_cntrs,con,wkarea
            );
        free(wkarea);
        cntr_vals[CNTR_LVL_U]=con[0];
        cntr_vals[CNTR_STEP_U]=con[1];
        /*
         * if zoom factor not set, draw the mesh boundary
         */
        labels= (Settings.plot_opts.dble  ? BNDRYDBL : BNDRY);
        if ((fabs(Settings.zoom-1.0)<1.e-4))
            plmesh_(&labels,&int_vars[NE],&int_vars[NN],&int_vars[NUP],
                        &int_vars[NBP],arraysf[EX],arraysf[EY],
                        arraysi[LEM],arraysi[NOR],arraysi[IBC],
                        arraysi[IBCTYP],arraysf[VHB],
                        &fl_vars[TBXOFF],&fl_vars[TBYOFF],&int_vars[IFLT],
                        &mpe,&startcol,&rangecol,
                        &pwindo[XCMIN],&pwindo[XCMAX],
                        &pwindo[YCMIN],&pwindo[YCMAX],
                        &vmin,&vmax,&Settings.zoom,
			&int_vars[NCOMP],&xref,&yref,&verbound);
        if ((plot_descrip&SHADE) && Settings.plot_opts.colour_bar!=CB_NONE) {
            markcolourbar_( &cntr_vals[CNTR_MIN_U],&cntr_vals[CNTR_MAX_U],
                        &num_cntrs,&cntr_vals[CNTR_STEP_U],
                        &cntr_vals[CNTR_LVL_U],
                        &orient_vert,
                        &xlevel,&ylevel,&barwdth,&barlen);
        }
    }
    if (Settings.plot_opts.label) {
        if (plot_descrip&LINES) format = 63;
        else format = 0;
        len = strlen(label);
        drawautolabels_(label,&len,&cntr_vals[CNTR_STEP_U],&format,
                                                  pwindo,current_pen);
    }

    *current_pen = old_pen;
    setpencolor_(current_pen);

    sybflush_();

    if(*verbose!=0){
        err = DataPrintMesh(label,arraysf[SNTRP],arraysi[IHELP],
                            pwindo,mesh[NX3],mesh[NY3],
                            xref,yref,int_vars[NCOMP],*verbose);
    }
    return(err);
}

     /*  name of profile changed to Profile to facilitate finding */

int Profile(data,flags,npts,profile_vals,
            ncomp,label,pwindo,current_pen,line)
/* char label[5]; */
char *label;
float *data,*profile_vals,*pwindo;
int *current_pen,*line;
int ncomp;
int *flags,npts;
{
    char number[MAX_NUM_LEN];
    int old_pen, old_line, format, len, idir;
    int mode, x_offset, y_offset;
    int m1,n1,mps,nps,i_min,i_max,j_min,j_max;
    int verbose,tic,retval,err=0,leng=1;
    float scale_min, scale_max;
    float line_len,line_min,line_max,ave_lat,ave_lon;
    float xrng, yrng, xref, yref;
    float *data2;
    float xmin, xmax, ymin, ymax, inc, area, sumint;
    float iv_ticint, iv_min, iv_max, dv_ticint, dv_min, dv_max;
    float scfactor, ticint, ticstart, size ;
/*  float rtodeg=57.2957795;  */
    float xtmp,ytmp;
//  float dtmp,xltmp,xinc,yinc,dinc,line_inc;

    old_line = *line; old_pen = *current_pen;

/*    XREFM and YREFM usage depends on context */
    xref=Data_vars_fl[XREFM];
    yref=Data_vars_fl[YREFM];
    if (Data_vars_int[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }
    m1=n1=mps=nps=1;
    /* find the data max and min */
    rangxy_(data,flags,&npts,&m1,&npts,&mps,&m1,&n1,&nps,
                &profile_vals[PRFL_MAX],&i_max,&j_max,
                &profile_vals[PRFL_MIN],&i_min,&j_min,&npts,&sumint);
    fprintf(stdout,"Profile values for %s:  Min = %f  Max = %f\n",label,
                      profile_vals[PRFL_MIN],profile_vals[PRFL_MAX]);

    /* propose plotting limits and tick interval based on rounded numbers
        this for dependent variable: dv_min -> dv_max         */

    if (ticstep(profile_vals[PRFL_MIN],profile_vals[PRFL_MAX],
                     &dv_ticint,&dv_min,&dv_max)) {
            warning_msg(RANGE_ERR,"tic max = min - not drawing tics"); 
    }
    profile_vals[PRFL_MIN] = dv_min;
    profile_vals[PRFL_MAX] = dv_max;
    if (profile_vals[PRFL_MIN_U]==0.0 && profile_vals[PRFL_MAX_U]==0.0) {
        profile_vals[PRFL_MIN_U] = profile_vals[PRFL_MIN];
        profile_vals[PRFL_MAX_U] = profile_vals[PRFL_MAX];
    }

    /* allow user to adjust values or cancel, then re-scale and 
       reset ticint etc.  profile_vals[PRFL_MIN_U] and [PRFL_MAX_U]
       will be the limiting values on the graph axis (user units) */

    retval=GetUserVals( profile_vals,PRFLVALS,0,NULL,NULL );
    if (retval==USER_CANCEL)    return(USER_CANCEL);
    if ((data2=(float *)malloc(npts*sizeof(float)))==NULL)
            return(MALLOC_ERR);
    for (m1=0; m1<npts; m1++) {
              data[m1]*=profile_vals[PRFL_SCL];
              data[m1]+=profile_vals[PRFL_SHFT];
              }
    if (ticstep(profile_vals[PRFL_MIN_U],profile_vals[PRFL_MAX_U],
                     &dv_ticint,&dv_min,&dv_max)) {
            warning_msg(RANGE_ERR,"tic max = min - not drawing tics");
    }
    dv_min=profile_vals[PRFL_MIN_U];
    dv_max=profile_vals[PRFL_MAX_U];
 /* fprintf(stdout,"dv_min = %f dv_max =%f\n",dv_min,dv_max); */

/* (PRFL_X1, PRFL_Y1) and (PRFL_X2, PRFL_Y2) define the coordinates of profile
    end points in model space, PRFL_MIN_U and PRFL_MAX_U define lower and
    upper limits to be used for plotting the dependent variable */

    line_len = (float)sqrt((double)(
                               (profile_vals[PRFL_X2]-profile_vals[PRFL_X1])
                              *(profile_vals[PRFL_X2]-profile_vals[PRFL_X1])
                             + (profile_vals[PRFL_Y2]-profile_vals[PRFL_Y1])
                              *(profile_vals[PRFL_Y2]-profile_vals[PRFL_Y1])));
    len = strlen(label);
    integrattrpz_(data,flags,&npts,&line_len,&area,label,&len);
    fprintf(stdout,"integrattrpz area %f\n",area);
    /*
    integrat_(data,flags,&npts,&line_len,&area,label,&len);
    fprintf(stdout,"integrat area %f\n",area);
    */

    verbose = Settings.verbose;

/* choose plot orientation depending whether X-range > Y-range in model space 
   probably choice is ok for vertical or horizontal profile:
     idir=0: independent variable is plotted horizontal
     idir=1: independent variable is plotted vertical    
   same test works for 1D and 2D plots                   */

    idir=0;
    xrng=profile_vals[PRFL_X2]-profile_vals[PRFL_X1];
    yrng=profile_vals[PRFL_Y2]-profile_vals[PRFL_Y1];
    if (fabs(xrng) < fabs(yrng)) idir = 1;

/*  to scale x-axis  as for contour plots use this block in next section
        xmin = pwindo[XCMINREF];
        xmax = xmin + line_len;
        ymin = pwindo[UBASE];
        ymax = pwindo[UTOP];
        inc = (line_len)/(npts-1);   */

/*  this block for plotting independent variable horizontal  */

    xmin = pwindo[ULEFT];
    xmax = pwindo[URGHT];
    scale_min = profile_vals[PRFL_MIN_U];
    scale_max = profile_vals[PRFL_MAX_U];
/*    max < min to reverse the axis of the dependent variable */
    if (scale_min>scale_max) dv_ticint = -dv_ticint;

    if(idir==0) {
        ymin = pwindo[UBASE];
        ymax = pwindo[UTOP];
        inc = (xmax-xmin)/(npts-1);
        for (m1=0,xtmp=xmin; m1<npts; m1++,xtmp+=inc) data2[m1]=xtmp;
        calcprofil_(data2,data,flags,&npts,
            &xmin,&xmax,
            &scale_min,
            &scale_max,
            &xmin,&ymin,&xmax,&ymax,&verbose
            );
    }
/* this block for plotting independent variable vertical
    data and data2 are swapped in the call to calcprofil */
    else {
        ymin = pwindo[YCMIN];
        ymax = pwindo[YCMAX];
        inc = (ymax-ymin)/(npts-1);
        for (m1=0,ytmp=ymin; m1<npts; m1++,ytmp+=inc) data2[m1]=ytmp;
        calcprofil_(data,data2,flags,&npts,
            &scale_min,&scale_max,&ymin,&ymax,
            &xmin,&ymin,&xmax,&ymax,&verbose
            );
    }

/* at this stage (xmin, ymin) and (xmax, ymax) are corners of plot
     region in model coordinates          */

/*    independent variable for graph is the coordinate
      that is changing most across the profile.  For spherical
      projection, convert distance from radians to degrees */

    if(idir==0){
       line_min=profile_vals[PRFL_X1];
       line_max=profile_vals[PRFL_X2];
       if(ncomp==-1){
         ave_lat=0.5*(profile_vals[PRFL_Y1]+profile_vals[PRFL_Y2]);
         projectdeg_(&line_min,&ave_lat,&xref,&yref,&leng,&ncomp,&err);
         ave_lat=0.5*(profile_vals[PRFL_Y1]+profile_vals[PRFL_Y2]);
         projectdeg_(&line_max,&ave_lat,&xref,&yref,&leng,&ncomp,&err);
       }
    }
    else{
       line_min=profile_vals[PRFL_Y1];
       line_max=profile_vals[PRFL_Y2];
       if(ncomp==-1){
         ave_lon=0.5*(profile_vals[PRFL_X1]+profile_vals[PRFL_X2]);
         projectdeg_(&ave_lon,&line_min,&xref,&yref,&leng,&ncomp,&err);
         ave_lon=0.5*(profile_vals[PRFL_X1]+profile_vals[PRFL_X2]);
         projectdeg_(&ave_lon,&line_max,&xref,&yref,&leng,&ncomp,&err);
       }
    }
    line_len=line_max-line_min;

/*  if verbose flag is set, output profile values in (x,y,value,xlabel,distance)
    sets.  If spherical thin sheet, x, y, distance are converted to degrees           */

    if(verbose!=0){
        err = DataPrintLine(label,data,profile_vals,line_len,line_min,
                        xref,yref,npts,ncomp,verbose);
    }
/*#if XY
    if(verbose!=0){
  //   fprintf(stdout,"(x,y,value,xlabel,distance) for %s profile follow\n",label);
       xtmp0=profile_vals[PRFL_X1];
       ytmp0=profile_vals[PRFL_Y1];
       xinc=(profile_vals[PRFL_X2]-xtmp0)/(float)(npts-1);
       yinc=(profile_vals[PRFL_Y2]-ytmp0)/(float)(npts-1);
       line_inc=line_len/(float)(npts-1);
       dinc=rtodeg*sqrt((double)xinc*(double)xinc+(double)yinc*(double)yinc);   
       for (m1=0; m1<npts; m1++){
          xtmp=xtmp0+xinc*(float)m1;
          ytmp=ytmp0+yinc*(float)m1;
          if(ncomp==-1){
            projectdeg_(&xtmp,&ytmp,&xref,&yref,&leng,&ncomp,&err);
          };
          dtmp=dinc*(float)m1;
          xltmp=line_min+line_inc*(float)m1;   //  already in degrees
   //     fprintf(stdout," %f %f %f %f %f\n",xtmp,ytmp,data[m1],xltmp,dtmp);
       }
    }
#endif */

/*     draw external boundary in model units  */

      *current_pen=FOREGROUND;
      setpencolor_(current_pen);
      *line = -1; mode = 1;
      dashln_(line,&mode);
      drawrectangleu_(&xmin,&ymin,&xmax,&ymax);

/*     Write title on graph */

      if (Settings.plot_opts.label) {
        xtmp = xmin+(xmax-xmin)/2.0;
        ytmp = ymax+(ymax-ymin)/50.;
        len=strlen(label);
        x_offset = J_CENTRE; y_offset = J_BASE; mode = 1;
        drawlabel_(&xtmp,&ytmp,&x_offset,&y_offset,&mode, label,&len);
      }

/* tick marks, then label on vertical axis */

      tic = Settings.plot_opts.ticmark;
      size = (float)sqrt((double)((ymax-ymin)*(xmax-xmin)))/30.0;
      if (idir == 0) {    /*   dependent variable in user units  */
        scale_min = profile_vals[PRFL_MIN_U];
        scale_max = profile_vals[PRFL_MAX_U];
        scfactor = (ymax-ymin)/(scale_max-scale_min);
        ticint = dv_ticint*scfactor;  /* convert to plotter units */
        ticstart = ymin + (dv_min - scale_min)*scfactor;
      }
      else {           /*  independent variable in model units  */
        scale_min = line_min;
        scale_max = line_min+line_len;
        scfactor = (ymax-ymin)/(scale_max-scale_min);
        ticstep(scale_min,scale_max,&iv_ticint,&iv_min,&iv_max);
        ticint = iv_ticint*scfactor;  /* convert to plotter units */
        ticstart = ymin + (iv_min - scale_min)*scfactor;
      }
      mode = 0;
      drawticmarks_(&size,&ymax,&xmin,&ticstart,&ticint,&tic,&mode);

            /* leave gap next to axis */
    if (Settings.plot_opts.label) {
      xtmp = xmin-(xmax-xmin)/30.0; 
      format     = MAX(fabs((double)scale_min),fabs((double)scale_max))
                        >10.0     ? 60 : 63;
      formatnumber_(number,&scale_min,&format);
      x_offset = J_RIGHT; y_offset = J_BASE; mode = 1;
      len=strlen(number); if (len > (MAX_NUM_LEN-1)) len = MAX_NUM_LEN-1;
      drawlabel_(&xtmp,&ymin,&x_offset,&y_offset,&mode, number,&len);
      formatnumber_(number,&scale_max,&format);
      x_offset = J_RIGHT; y_offset = J_TOP; mode = 1;
      len=strlen(number); if (len > (MAX_NUM_LEN-1)) len = MAX_NUM_LEN-1;
      drawlabel_(&xtmp,&ymax,&x_offset,&y_offset,&mode, number,&len);
      if (idir != 0) {
         ytmp = ymin+(ymax-ymin)/2.0;
         strcpy(number,"Y");
         len=strlen(number);
         x_offset = J_RIGHT; y_offset = J_CENTRE; mode = 1;
         drawlabel_(&xtmp,&ytmp,&x_offset,&y_offset,&mode,number,&len);
      }
    }

        /* tick marks then label on horizontal axis */

    if (idir!=0) {         /*   dependent variable in user units */
       scale_min = profile_vals[PRFL_MIN_U];
       scale_max = profile_vals[PRFL_MAX_U];
       scfactor = (xmax-xmin)/(scale_max-scale_min);
       ticint = dv_ticint*scfactor;  /* convert to plotter units */
       ticstart = xmin + (dv_min - scale_min)*scfactor;
    }
    else {           /*   independent variable in model units  */
       scale_min = line_min;
       scale_max = line_min+line_len;
       ticstep(scale_min,scale_max,&iv_ticint,&iv_min,&iv_max);
       scfactor = (xmax-xmin)/(scale_max-scale_min);
       ticint = iv_ticint*scfactor;  /* convert to plotter units */
       ticstart = xmin + (iv_min - scale_min)*scfactor;
        }
    mode = 1;
    drawticmarks_(&size,&xmax,&ymin,&ticstart,&ticint,&tic,&mode);

            /* leave gap next to axis */
    if (Settings.plot_opts.label) {
        ytmp = ymin-(ymax-ymin)/30.0;
        format     = MAX(fabs((double)scale_min),fabs((double)scale_max))
                          >10.0     ? 60 : 63;
        formatnumber_(number,&scale_min,&format);
 /*     x_offset = J_LEFT; y_offset = J_TOP; mode = 1; */
        x_offset = J_CENTRE; y_offset = J_TOP; mode = 1;
        len=strlen(number); if (len > (MAX_NUM_LEN-1)) len = MAX_NUM_LEN-1;
        drawlabel_(&xmin,&ytmp,&x_offset,&y_offset,&mode,number,&len);
        formatnumber_(number,&scale_max,&format);
 /*     x_offset = J_RIGHT; y_offset = J_TOP; mode = 1; */
        len=strlen(number); if (len > (MAX_NUM_LEN-1)) len = MAX_NUM_LEN-1;
        drawlabel_(&xmax,&ytmp,&x_offset,&y_offset,&mode,number,&len);
        if (idir == 0) {
          xtmp = xmin+(xmax-xmin)/2.0;
          strcpy(number,"X");
          len=strlen(number);
          x_offset = J_CENTRE; y_offset = J_TOP; mode = 1;
          drawlabel_(&xtmp,&ytmp,&x_offset,&y_offset,&mode,number,&len);
        }
    }

/*   if in spherical coordinates, re-project to lat-long before saving */

    if(ncomp==-1){
        projectdeg_(&profile_vals[PRFL_X1],&profile_vals[PRFL_Y1],
                    &xref,&yref,&leng,&ncomp,&err); 
        projectdeg_(&profile_vals[PRFL_X2],&profile_vals[PRFL_Y2],
                    &xref,&yref,&leng,&ncomp,&err);

/*  following store operation now done in GetProfileVals */
/*      for (i=0;i<MAXPRFLVALS;i++)
                     current_vals[i]=profile_vals[i]; */
    }

    *current_pen=old_pen;
    setpencolor_(current_pen);
    *line=old_line; mode = 1;
    dashln_(line,&mode);
    free(data2);
    sybflush_();
    return( 0 );
}

int DrawProfileLine(float* prfl_vals, float* pwindo, int ncomp)
{
    int pen_state;
    int leng=1,err=0;
    float xref,yref;
/*    XREFM and YREFM usage depends on context */
    xref=Data_vars_fl[XREFM];
    yref=Data_vars_fl[YREFM];
    if (Data_vars_int[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }

    /*
     * if the user requests a single point profile,
     * draw an isosceles triangle with the smallest angle
     * apex at the requested point
     */
    if ((prfl_vals[PRFL_X1]==prfl_vals[PRFL_X2]) &&
        (prfl_vals[PRFL_Y1]==prfl_vals[PRFL_Y2])) {
        int num=3; 
        float len=(pwindo[XCMAX]-pwindo[XCMIN])*0.05;
        float xp[3], yp[3]; 
        int outline=2; 
        int col=Settings.fg;
        xp[0]=prfl_vals[PRFL_X1]; yp[0]=prfl_vals[PRFL_Y1];
        xp[1]=xp[0]-len*0.3; yp[1]=yp[0]-len;
        xp[2]=xp[0]+len*0.3; yp[2]=yp[1];
        fillpoly_(xp,yp,&num, &col,&outline);
    }

/*  or else draw the line from (X1,Y1) to (X2, Y2)  */

    else {
    setlinewidth_(&Settings.linewidth);
    pen_state = PEN_UP;
    plotu_(&prfl_vals[PRFL_X1],&prfl_vals[PRFL_Y1],&pen_state);
    pen_state = PEN_DN;
    plotu_(&prfl_vals[PRFL_X2],&prfl_vals[PRFL_Y2],&pen_state);
    }

/*   if in spherical coordinates, re-project to lat-long before saving */
 /*   fprintf(stdout," DrawProfileLine: ncomp = %i\n",ncomp);  */
    if(ncomp==-1){
        projectdeg_(&prfl_vals[PRFL_X1],&prfl_vals[PRFL_Y1],
                    &xref,&yref,&leng,&Data_vars_int[NCOMP],&err);
        projectdeg_(&prfl_vals[PRFL_X2],&prfl_vals[PRFL_Y2],
                    &xref,&yref,&leng,&Data_vars_int[NCOMP],&err);
    }

    Save_plot_data( NULL,prfl_vals );
    return( 0 );
}

int DrawBndBox(float* pwin, float* prfl_vals)
{
    drawrectangleu_(&pwin[XCMIN],&pwin[YCMIN],&pwin[XCMAX],&pwin[YCMAX]);
    Save_plot_data( NULL,prfl_vals );
    return(0);
}

int set_max_cntrs(float* cntr_vals, float value)
{
    cntr_vals[MAX_CNTRS] = value;
    return(0);
}

/*     ticstep chooses sensible limits (g_min, g_max) and tick intervals 
       (interval) for variable values that are in the range (min, max) */

int ticstep(float min,float max,float *interval,float* g_min,float* g_max)
//int ticstep(min,max,interval,g_min,g_max)
//float min,max,*interval,*g_min,*g_max;
{
    float p1,stp1,eps=-0.00001;
    int bnd1,bnd2,icnt=0;

    stp1 = fabs (max-min);
    if(stp1 == 0.){
      fprintf(stdout," max = min = %f not permissible for profiles\n",max);
      return(1);
    }
    p1=stp1/8.0;

/* use factors of 10 to put the step size in range 1-10 */

    while ( (abs(icnt)< 20) && (p1 <= 1.0) )
                      { p1=p1*10. ; icnt=icnt-1; }
    while ( (abs(icnt)< 20) && (p1 >= 10.0) )
                      { p1=p1*0.1 ; icnt=icnt+1; }

/* if still outside sensible range, write warning and set arbitrary value */

    if((p1 < 1.0)||(p1 > 10)){
       stp1=0.01;
       fprintf(stdout,"Could not choose sensible scale after %i steps\n",icnt);
    }

/* truncate and re-scale by icnt factors of 10  */

    if(p1 <= 1.5) p1=1.0;
    if((p1 <= 3.0)&&(p1 > 1.5)) p1 = 2.0;
    if((p1 <= 7.0)&&(p1 > 3.0)) p1 = 5.0;
    if(p1 > 7.0) p1=10.0;
    p1=p1 * pow (10.,icnt);
    
    *interval = p1;

/*  set g_min and g_max as multiples of p1 - user units 
    casting to integer truncates, different for (+, -)  
    eps allows round-off, to get ticks on end-points */

    if(min >= 0){ bnd1 =  (int) ((min/p1)*(1.0+eps));}
    else       { bnd1 = ((int) ((min/p1)*(1.0-eps))) - 1;}
    *g_min = bnd1*p1;
    if(max <= 0){ bnd2 =  (int) ((max/p1)*(1.0-eps)); }
    else       { bnd2 = ((int) ((max/p1)*(1.0+eps))) + 1;}
    *g_max = bnd2*p1;
/*  fprintf(stdout,"ticstep g_min= %f g_max=%f\n",*g_min,*g_max); */

    return(0);
}

int DataPrintLine(char *label,float *data,float *profile_vals,float line_len,
              float line_min, float xref, float yref, int npts, int ncomp,
              int verbose)
/* Pass xmin,ymin,xmax,ymax? */
{
    char *fname;
    char datafile[SYB_FILENAME_MAX+1];
    char header1[SYB_FILENAME_MAX+12], header2[MAX_LABEL_LEN+1];
    int leng=1, m1;
    int err=0;
//  int len=0;
    float line_inc;
    float xtmp0,xtmp,ytmp0,ytmp,xltmp,xinc,yinc;
//  float rtodeg=57.2957795;
//  float dtmp,dinc;
    FILE *fp=0;

    /* using global Plot_info */
    if ((fname=strrchr(Plot_info.inp_file->fname,'/'))==NULL)
        fname=Plot_info.inp_file->fname;
    else fname++;
//  len = strlen(fname);
    snprintf(header1,SYB_FILENAME_MAX+12,"%s record %02d\n",
                     fname,Plot_info.inp_file->rec_curr);
    if (ncomp==-1)
        snprintf(header2,MAX_LABEL_LEN,"%s\t%s\t%s\t%s\n",
                                          "X","Y",label,"dist");    
    else
        snprintf(header2,MAX_LABEL_LEN,"%s\t%s\t%s\n","X","Y",label);    
    snprintf(datafile,SYB_FILENAME_MAX,"%s.%s","sybil",label);
    fprintf(stdout,"%s %s to %s; columns: %s\n",
                   "Writing interpolated values of",label,datafile,header2);
    if (verbose==2) fp = stdout;
    else if ((fp=fopen(datafile,"w"))==NULL) err=OPEN_ERR;
    if (!err) {
        fputs(header1,fp);
        fputs(header2,fp);
        xtmp0=profile_vals[PRFL_X1];
        ytmp0=profile_vals[PRFL_Y1];
        xinc=(profile_vals[PRFL_X2]-xtmp0)/(float)(npts-1);
        yinc=(profile_vals[PRFL_Y2]-ytmp0)/(float)(npts-1);
        line_inc=line_len/(float)(npts-1);
//      dinc=rtodeg*sqrt((double)xinc*(double)xinc+(double)yinc*(double)yinc);
        for (m1=0; m1<npts; m1++){
            xtmp=xtmp0+xinc*(float)m1;
            ytmp=ytmp0+yinc*(float)m1;
            if(ncomp==-1){
              projectdeg_(&xtmp,&ytmp,&xref,&yref,&leng,&ncomp,&err);
            };
//          dtmp=dinc*(float)m1;
            xltmp=line_min+line_inc*(float)m1;   /*  already in degrees  */
            if(ncomp==-1)
              fprintf(fp," %f %f %e %f\n",xtmp,ytmp,data[m1],xltmp);
            else
              fprintf(fp," %f %f %e\n",xtmp,ytmp,data[m1]);
        }
        if (verbose!=2) fclose(fp);
    }
    return(err);
}

int DataPrintMesh(char *label,float *data,int *flags,float *pwindo,
                  int nx3,int ny3,float xref,float yref,int ncomp,
                  int verbose)
{
    char *fname;
    char datafile[SYB_FILENAME_MAX+1];
    char header1[SYB_FILENAME_MAX+12], header2[MAX_LABEL_LEN+1];
    int err=0;
    int i, j, ij;
    int leng=1;
//  int len;
    float ybp;
    float xtmp0,xtmp,ytmp0,ytmp,xinc,yinc;
    FILE *fp=0;

    /* using global Plot_info */
    if ((fname=strrchr(Plot_info.inp_file->fname,'/'))==NULL)
        fname=Plot_info.inp_file->fname;
    else fname++;
//  len = strlen(fname);
    snprintf(header1,SYB_FILENAME_MAX+12,"%s record %02d\n",
                     fname,Plot_info.inp_file->rec_curr);
    snprintf(header2,MAX_LABEL_LEN,"     %s         %s     %s\n",
                                          "X","Y",label);    
    snprintf(datafile,SYB_FILENAME_MAX,"%s.%s","sybil",label);
    printf("%s %s to %s; columns: %s\n",
                   "Writing interpolated values of",label,datafile,header2);
    if (verbose==2) fp = stdout;
    else if ((fp=fopen(datafile,"w"))==NULL) err=OPEN_ERR;
    if (!err) {
        fputs(header1,fp);
        fputs(header2,fp);
	xtmp0=pwindo[XCMIN];
        ytmp0=pwindo[YCMIN];
        xinc=(pwindo[XCMAX]-xtmp0)/(float)(nx3-3);
        yinc=(pwindo[YCMAX]-ytmp0)/(float)(ny3-3);
        for (j=1; j<(ny3-1); j++) {
            ybp=ytmp0+yinc*(j-1);
            ytmp=ybp;
            for (i=1; i<(nx3-1); i++) {
                xtmp=xtmp0+xinc*(i-1);
                ij= j*nx3+i;
                if (flags[ij]!=0) {
                    if(ncomp==-1){
                        ytmp=ybp;
                        projectdeg_(&xtmp,&ytmp,&xref,&yref,&leng,&ncomp,&err);
                    }
                    fprintf(fp," %f %f %e\n",xtmp,ytmp,data[ij]);
                }
            }
        }
        if (verbose!=2) fclose(fp);
    }
    return(err);
}

int DataPrintArrow(char *label,float *c2t,float *s2t, float *wx,float *wy,
                   int m1,int m2,int mp,int n1,int n2,int np,int nx3,
                   float *pwindo, float xref, float yref, int ncomp,
                   int verbose)
{
    char *fname;
    char datafile[SYB_FILENAME_MAX+1];
    char header1[SYB_FILENAME_MAX+12], header2[MAX_LABEL_LEN+1];
    int err=0;
    int i, j, ij;
//  int len;
    int leng=1;
    float xtmp0,xtmp,ytmp0,ytmp,xinc,yinc;
    float s1,s2,cc,theta,thetadeg;
    FILE *fp=0;

    /* using global Plot_info */
    if ((fname=strrchr(Plot_info.inp_file->fname,'/'))==NULL)
        fname=Plot_info.inp_file->fname;
    else fname++;
//  len = strlen(fname);
    snprintf(header1,SYB_FILENAME_MAX+12,"%s record %02d\n",
                     fname,Plot_info.inp_file->rec_curr);
    snprintf(header2,MAX_LABEL_LEN,"%8s%9s%13s%13s%13s  %s\n",
                     "centreX","centreY","axis1","axis2","axisC","orientation");
    snprintf(datafile,SYB_FILENAME_MAX,"%s.%s","sybil",label);
    if (verbose==2) fp = stdout;
    else if ((fp=fopen(datafile,"w"))==NULL) err=OPEN_ERR;
    if (!err) {
        fputs(header1,fp);
        fputs(header2,fp);
        xtmp0=pwindo[XCMIN];
        ytmp0=pwindo[YCMIN];
        xinc=(pwindo[XCMAX]-xtmp0)/(float)(m2-m1);
        yinc=(pwindo[YCMAX]-ytmp0)/(float)(n2-n1);
        xtmp = xtmp0;
        ytmp = ytmp0;
        for (j=n1; j<=n2; j+=np) {
            for (i=m1; i<=m2; i+=mp) {
                xtmp=xtmp0+xinc*(i-m1);
                ytmp=ytmp0+yinc*(j-n1);
                ij= (j-1)*nx3+i-1;
                s1 = wx[ij];
                s2 = wy[ij];
                if ((s1!=0.0)||(s2!=0.0)) {
                    if ((s2t[ij]==0.0) && (c2t[ij]==0.0))
                        theta = 0.0;
                    else
                        theta = atan2(s2t[ij],c2t[ij]);
                    cc = -s1-s2;
                    thetadeg = theta*180.0/3.141592654;
                    if(ncomp==-1){
                        projectdeg_(&xtmp,&ytmp,&xref,&yref,&leng,&ncomp,&err);
                    }
                    fprintf(fp," %f %f %e %e %e %e\n",
                                    xtmp,ytmp,s1,s2,cc,thetadeg);
                 }
             }
        }
        if (verbose!=2) fclose(fp);
    }
    return(err);
}

int DataPrintMeshVelo(char *label, float *data1, float *data2,
                      int m1, int m2, int mp, int n1, int n2, int np,
                      int nx3, int ny3,
                      float *pwindo, float xref, float yref, int ncomp,
                      int verbose)
{
    char *fname;
    char datafile[SYB_FILENAME_MAX+1];
    char header1[SYB_FILENAME_MAX+12], header2[SYB_FILENAME_MAX+12];
    int err=0;
    int i, j, ij;
/*  int len;  */
    int leng=1;
    float xtmp0,xtmp,ytmp0,ytmp,xinc,yinc,s1,s2;
    FILE *fp=0;

    header1[0] = header2[0] = '\0';
    /* using global Plot_info */
    if ((fname=strrchr(Plot_info.inp_file->fname,'/'))==NULL)
        fname=Plot_info.inp_file->fname;
    else fname++;
 /* len = strlen(fname);  */
    snprintf(header1,SYB_FILENAME_MAX+12,"%s record %02d\n",
                     fname,Plot_info.inp_file->rec_curr);
    snprintf(header2,MAX_LABEL_LEN,"%6s%6s%6s%6s\n",
                                          "X","Y","Ux","Uy");    
    snprintf(datafile,SYB_FILENAME_MAX,"%s.%s","sybil",label);
    fprintf(stdout,"%s %s to %s; columns: %s\n",
                   "Writing interpolated values of",label,datafile,header2);
    if (verbose==2) fp = stdout;
    else if ((fp=fopen(datafile,"w"))==NULL) err=OPEN_ERR;
    if (!err) {
        fputs(header1,fp);
        fputs(header2,fp);
        xtmp0=pwindo[XCMIN];
        ytmp0=pwindo[YCMIN];
        xinc=(pwindo[XCMAX]-xtmp0)/(float)(m2-m1);
        yinc=(pwindo[YCMAX]-ytmp0)/(float)(n2-n1);
        for (j=n1; j<=n2; j+=np) {
            for (i=m1; i<=m2; i+=mp) {
                xtmp=xtmp0+xinc*(i-m1);
                ytmp=ytmp0+yinc*(j-n1);
                ij= (j-1)*nx3+i-1;
		s1=data1[ij];
		s2=data2[ij];
		if ((s1!=0.0)||(s2!=0.0)) {
                   if(ncomp==-1){
                      projectdeg_(&xtmp,&ytmp,&xref,&yref,&leng,&ncomp,&err);
                   }
                   fprintf(fp," %f %f %e %e\n",xtmp,ytmp,s1,s2);
		}
            }
        }
        if (verbose!=2) fclose(fp);
    }
    return(err);
}
