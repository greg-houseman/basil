/*--------------------------------------------------------------------
 *    Basil / Sybil:   log.c    1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#ifdef XSYB
#include <X11/Xlib.h>  /* for type XFontStruct in Plot_info - tidy */
#include <X11/Intrinsic.h>  /* for type XtPointer - tidy */
#endif

#ifdef SUN
#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "types.h"
#include "mainps.h"
#include "globals.h"
#include "strain.h"
#include "deform.h"
#include "string_utils.h"
#include "log.h"
#include "errnum.h"
#include "error.h"
#include "routines.h"

/* int Save_option(XtPointer, int, int); */
int Save_option();
extern int Arrow_colours[];
extern void ChangeCell(),clear_arrays();

int num_items(set)
char *set[];
{
    int i=0;

    while (set[i]) i++;
    return(i);
}

/*
 * searches the given set of valid_terms and returns
 * the associated id if a name match is found
 * else returns -1
 */
int name_match(str,set)
char *str;
valid_terms set[];
{
    unsigned char found;
    int i=0;

    found=0;
    while (set[i].name && !found) {
        if (!strcmp(str,set[i].name)) found=1;
        else i++;
    }
    return( found ? set[i].id : -1 );
}

/*
 * searches the given set of names and returns
 * the index if a name match is found
 * else returns -1
 */
int word_match(str,set)
char *str, *set[];
{
    unsigned char found;
    int i=0;

    found=0;
    while (set[i] && !found) {
        if (!strcmp(str,set[i])) found=1;
        else i++;
    }
    return( found ? i : -1 );
}

/*
int valid_key( key )
int key;
{
    if (key==OPTIONS || key==FILE_OPT || key==TITLE ||
        key==DATA_IN || key==XY_PLOT || key==ARROWS ||
        key==PRFL || key==CNTRS || key==LABEL || key==LOCATE)
        return(1);
    return(0);
}
*/

int validate( str,key,set )
char *str;
int *key;
valid_terms set[];
{
    int err=0;

    if (str[0]=='#') *key = COMMENT;
    else *key = name_match(str,set);
    return( err );
}

int initial_options(FILE* logfp, input_options* opts, char* input_str)
{
    unsigned char finished=0, initial=1;
    int err = 0, key;
    int num=0;
    long pos;

    input_str[0] = '\0';
    if (fscanf(logfp,"%s",input_str)!=1) err = READ_ERR;
    else {
      while (!feof(logfp) && !finished) {
        validate(input_str,&key,keywords);
        switch( key ) {
        case OPTIONS:  err=read_options(logfp,input_str,opts,initial);
                       if (!err && !feof(logfp)) {
                         validate(input_str,&key,keywords);
                         if (key==-1) err = VAR_ERR;
                         else {
                           pos = ftell(logfp); pos -= strlen(input_str);
                           fseek(logfp,pos,SEEK_SET);
                         }
                       }
                       finished = 1;
                       break;
        case COMMENT:  dump_comments(logfp);
                       num=fscanf(logfp,"%[^A-Za-z#]",input_str);
                       num=fscanf(logfp,"%s",input_str);
                       break;
        default:       err = INIT_ERR;
                       finished = 1;
                       break;
        }
      }
    }
    return( err );
}

int process_log_file(FILE* logfp, input_options* opts)
{
    unsigned char finished=0, initial=0;
    int err = 0, key, err_num=0;
    int num=0;
    char input_str[SYB_FILENAME_MAX], current_soln_file[SYB_FILENAME_MAX];

    input_str[0] = current_soln_file[0] = '\0';
    if (fscanf(logfp,"%s",input_str)!=1) {
        if (feof(logfp)) return(0);
        else error_msg(READ_ERR,input_str);
    }
    while (!feof(logfp) && !finished) {
      validate(input_str,&key,keywords);
      switch( key ) {
      case FILE_OPT: if ((err=read_filename(logfp,input_str)))
                        error_msg(err,input_str);
                     if (strcmp(current_soln_file,input_str)) {
                         if (Plot_info.inp_file!=NULL) {
                             if (Plot_info.inp_file->fp!=NULL) {
                                 fclose(Plot_info.inp_file->fp);
                                 clear_arrays();
                             }
                         }
                         if ((err=file_open( &(Plot_info.inp_file),
                                          input_str,SYB_FILENAME_MAX,"rb" )))
                           error_msg(err,input_str);
                           fprintf(stdout,"Opened file: %s\n",input_str);

                         strncpy(current_soln_file,input_str,SYB_FILENAME_MAX);
                         /*
                          * the DATA_IN statement will check the record
                          * (else zoom done twice etc)
                          */
                         /*if (err=CheckCurrentRec(Err_str,&err_num))*/
                             /*error_msg(err,input_str);*/
                     }
                     num=fscanf(logfp,"%[^A-Za-z#]",input_str);
                     num=fscanf(logfp,"%s",input_str);
                     break;
      case DATA_IN:  if ((err=set_records(logfp,Plot_info.inp_file)))
                        error_msg(err,input_str);
                 /*  fprintf(stdout,"CheckCurrentRec from DATA_IN\n"); */
                     if ((err=CheckCurrentRec(Err_str,&err_num)))
                         error_msg(err,input_str);
                     num=fscanf(logfp,"%[^A-Za-z#]",input_str);
                     num=fscanf(logfp,"%s",input_str);
                     break;
      case XY_PLOT:  if ((err=create_xy_plot(logfp,input_str)))
                        error_msg(err,input_str);
                     num=fscanf(logfp,"%[^A-Za-z#]",input_str);
                     num=fscanf(logfp,"%s",input_str);
                     break;
      case ARROWS:   if ((err=create_arrow_plot(logfp,input_str)))
                        error_msg(err,input_str);
                     num=fscanf(logfp,"%[^A-Za-z#]",input_str);
                     num=fscanf(logfp,"%s",input_str);
                     break;
      case PRFL:     if ((err=create_profile_plot(logfp,input_str)))
                        error_msg(err,input_str);
                     break;
      case CNTRS:    if ((err=create_contour_plot(logfp,input_str)))
                        error_msg(err,input_str);
                     break;
      case LABEL:    if ((err=create_optional_label(input_str,logfp)))
                        error_msg(err,input_str);
                     num=fscanf(logfp,"%[^A-Za-z#]",input_str);
                     num=fscanf(logfp,"%s",input_str);
                     break;
      case LOCATE:   if ((err=set_cell(logfp)))
                        error_msg(err,input_str);
                     num=fscanf(logfp,"%[^A-Za-z#]",input_str);
                     num=fscanf(logfp,"%s",input_str);
                     break;
      case TITLE:    if ((err=read_title(logfp,input_str,&Plot_info.title)))
                        error_msg(err,input_str);
                     if (Plot_info.title) DrawTitle(Plot_info.title);
                     num=fscanf(logfp,"%[^A-Za-z#]",input_str);
                     num=fscanf(logfp,"%s",input_str);
                     break;
      case OPTIONS:  if ((err=read_options(logfp,input_str,opts,initial)))
                        error_msg(err,input_str);
                     /* if rescale set then process in DataIn */
                     /*if (opts->rescale!=0)
                         err=CheckCurrentRec(Err_str,&err_num);*/
                     break;
      case COMMENT:  dump_comments(logfp);
                     num=fscanf(logfp,"%[^A-Za-z#]",input_str);
                     num=fscanf(logfp,"%s",input_str);
                     break;
      default:       err = KEY_ERR;
                     error_msg(err,input_str);
                     break;
      }
      if (err) finished = 1;
    }
    return( 0 );
}

int create_xy_plot(FILE* logfp, char* input_str)
{
    fontdat curr_font;
    int err=0,plot_parameter= -1,err_num=0;
    int num=0;
    int len, hgt=8;
	unsigned char tmp;

/*    fprintf(stdout,"log.c create_xy_plot input_str = %s\n",input_str); */
    if ((err=read_xyplot_data(logfp,input_str,&plot_parameter)))
        error_msg(err,input_str);
    if (plot_parameter!=ELLE)
        if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
           strcpy(Err_str,"No solution file open");
           error_msg(0,Err_str);
        }
        /*if (err=CheckCurrentRec(Err_str,&err_num)) return(err);*/
    switch( plot_parameter ) {
    case LGMESH:    /* fall through */
    case MESH: if (Settings.plot_opts.dble) {
                   switch(Plot_info.var_num) {
                   case ELMNT     : Plot_info.var_num = ELMNTDBL;
                                    break;
                   case BNDRY     : Plot_info.var_num = BNDRYDBL;
                                    break;
		   case INTBND    : Plot_info.var_num = INTBNDDBL;
                                    break;
		   case INTLBNDRY : Plot_info.var_num = INTLBNDRYDBL;
                                    break;
                   case VISCMSH:    Plot_info.var_num = VISCDBL;
                                    break;
                   case SEMSH:      Plot_info.var_num = SEDBL;
                                    break;
                   default:         break;
                   }
      fprintf(stdout,"log.c case MESH Plot_info.var_num = %i\n",
        Plot_info.var_num);
               }
               if (Plot_info.var_num==ELMNTNUM) {
                        /* set font to small Helvetica for element labels */
                   curr_font = Plot_info.font;
                   len = strlen("Helvetica");
                   setfont_("Helvetica",&len,&hgt);
               }
               err=plot_mesh(Plot_info.var_num,Data_vars_int,Data_vars_fl,
                        Data_arrays_int,Data_arrays_fl,
                        Pwindo,Mesh,&Settings.fg
                        );
               if (Plot_info.var_num==ELMNTNUM) Plot_info.font = curr_font;
               if (err) error_msg(err,input_str);
               break;
    case BBOX: DrawBndBox(Pwindo,Profile_vals);
               break;
    case STRN_MRKR:
               plot_strain_markers(Data_vars_int,Data_vars_fl,
                                               Data_arrays_fl,Pwindo);
               break;
    case DEFORM:if ((err=CheckReferenceRec(Err_str,&err_num)))
                   error_msg(err,input_str);
               else if ((err=plot_deformation(Data_vars_int,Data_vars_fl,
                   Data_arrays_int,Data_arrays_fl,
                   Pwindo,Mesh,Plot_info.var_num,Plot_info.plot_description,
                   Contour_vals,Profile_vals,Plot_info.dflt_label1,
                   &Settings.fg)))
                       error_msg(err,"deform");
               break;
    case ELLE: num=fscanf(logfp,"%s",input_str);
               if ((err=file_open(&(Plot_info.elle_file),input_str,
                            SYB_FILENAME_MAX,"r"))) error_msg(err,input_str);
               strncpy(Plot_info.elle_file->fname,input_str,SYB_FILENAME_MAX);
               if (Plot_info.inp_file==NULL) {
				   tmp = Settings.rescale;
				   Settings.rescale = 1;
                   init_box(Pwindo,&Settings);
				   Settings.rescale = tmp;
               }
               if ((err=plot_elle_regions(Pwindo)))
                   error_msg(err,"elle");
               break;
    default:   error_msg(VAR_ERR,input_str);
               break;
    }
    return(0);
}

int create_arrow_plot(FILE* logfp, char* input_str)
{
    int err=0,plot_parameter= -1;
/*  int err_num=0; */

    if ((err=read_arrow_data(logfp,input_str,&plot_parameter)))
        error_msg(err,input_str);
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
       strcpy(Err_str,"No solution file open");
       error_msg(0,Err_str);
    }
    /*if (err=CheckCurrentRec(Err_str,&err_num)) return(err);*/
    switch( plot_parameter ) {
    case SYB_VELOCITY: plot_velocity(Data_vars_int,Data_vars_fl,
                        Data_arrays_int,Data_arrays_fl,
                        Plot_info.dflt_label1,Plot_info.var_num,
                        Plot_info.plot_type,Plot_info.plot_description,
                        Contour_vals,Profile_vals,
                        Pwindo,Mesh,&Settings.fg,&Settings.linestyle
                        );
               break;
    case SYB_STRAIN:
    case SYB_STRESS: if ((err=plot_strain(Data_vars_int,Data_vars_fl,
                        Data_arrays_int,Data_arrays_fl,
                        Pwindo,Mesh,Plot_info.var_num,
                        Contour_vals,Profile_vals,Plot_info.dflt_label1,
                        Plot_info.plot_type,Plot_info.plot_description,
                        &Settings.fg,
                        &Settings.linestyle,Arrow_colours
                        )))
                     error_msg(err,input_str);
               break;
    case ROTATION:    plot_rotation(Data_vars_int,Data_vars_fl,
                        Data_arrays_int,Data_arrays_fl,
                        Plot_info.dflt_label1,Plot_info.var_num,
                        Plot_info.plot_type,Plot_info.plot_description,
                        Contour_vals,Profile_vals,
                        Pwindo,Mesh,&Settings.fg,&Settings.linestyle
                        );
               break;
    default:   error_msg(VAR_ERR,input_str);
               break;
    }
    return(0);
}

int create_contour_plot(FILE* logfp, char* input_str)
{
    int err=0,plot_parameter= -1;
/*  int err_num=0; */


    if ((err=read_contour_data(logfp,input_str,&plot_parameter)))
        error_msg(err,input_str);
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
       strcpy(Err_str,"No solution file open");
       error_msg(0,Err_str);
    }
    /*if (err=CheckCurrentRec(Err_str,&err_num)) return(err);*/
    switch( plot_parameter ) {
    case SYB_VELOCITY: plot_velocity(Data_vars_int,Data_vars_fl,
                        Data_arrays_int,Data_arrays_fl,
                        Plot_info.dflt_label1,Plot_info.var_num,
                        Plot_info.plot_type,Plot_info.plot_description,
                        Contour_vals,Profile_vals,
                        Pwindo,Mesh,&Settings.fg,&Settings.linestyle
                        );
               break;
    case SYB_STRAIN:
    case SYB_STRESS: if ((err=plot_strain(Data_vars_int,Data_vars_fl,
                        Data_arrays_int,Data_arrays_fl,
                        Pwindo,Mesh,Plot_info.var_num,
                        Contour_vals,Profile_vals,Plot_info.dflt_label1,
                        Plot_info.plot_type,Plot_info.plot_description,
                        &Settings.fg,
                        &Settings.linestyle,Arrow_colours
                        )))
                     error_msg(err,input_str);
               break;
    case LAYER:    if ((err = plot_layer(Data_vars_int,Data_vars_fl,
                    Data_arrays_int,Data_arrays_fl,Pwindo,Mesh,
                    Plot_info.variable,
                    Plot_info.plot_type,Plot_info.plot_description,
                    Contour_vals,Profile_vals,Plot_info.var_num,
                    &Settings.fg,&Settings.linestyle
                    )))
                     error_msg(err,input_str);
               break;
    case DENSITY:  if ((err = plot_density(Data_vars_int,Data_vars_fl,
                    Data_arrays_int,Data_arrays_fl,Pwindo,Mesh,
                    Plot_info.variable,
                    Plot_info.plot_type,Plot_info.plot_description,
                    Contour_vals,Profile_vals,
                    &Settings.fg,&Settings.linestyle
                    )))
                     error_msg(err,input_str);
               break;
    default:   error_msg(VAR_ERR,input_str);
               break;
    }
    return(0);
}

int create_profile_plot(FILE* logfp, char* input_str)
{
    int err=0,plot_parameter= -1;

    if ((err=read_profile_data(logfp,input_str,&plot_parameter)))
        error_msg(err,input_str);
    if ((Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) &&
        (Plot_info.elle_file==NULL||Plot_info.elle_file->fp==NULL)) {
       strcpy(Err_str,"No solution file open");
       error_msg(0,Err_str);
    }
    /*if (err=CheckCurrentRec(Err_str,&err_num)) return(err);*/
    switch( plot_parameter ) {
    case SYB_VELOCITY: plot_velocity(Data_vars_int,Data_vars_fl,
                        Data_arrays_int,Data_arrays_fl,
                        Plot_info.dflt_label1,Plot_info.var_num,
                        Plot_info.plot_type,Plot_info.plot_description,
                        Contour_vals,Profile_vals,
                        Pwindo,Mesh,&Settings.fg,&Settings.linestyle
                        );
               break;
    case SYB_STRAIN:
    case SYB_STRESS: if ((err=plot_strain(Data_vars_int,Data_vars_fl,
                        Data_arrays_int,Data_arrays_fl,
                        Pwindo,Mesh,Plot_info.var_num,
                        Contour_vals,Profile_vals,Plot_info.dflt_label1,
                        Plot_info.plot_type,Plot_info.plot_description,
                        &Settings.fg,
                        &Settings.linestyle,Arrow_colours
                        )))
                     error_msg(err,input_str);
               break;
    case LAYER:   if ((err = plot_layer(Data_vars_int,Data_vars_fl,
                    Data_arrays_int,Data_arrays_fl,Pwindo,Mesh,
                    Plot_info.variable,
                    Plot_info.plot_type,Plot_info.plot_description,
                    Contour_vals,Profile_vals,Plot_info.var_num,
                    &Settings.fg,&Settings.linestyle
                    )))
                     error_msg(err,input_str);
               break;
    case DENSITY:  if ((err = plot_density(Data_vars_int,Data_vars_fl,
                    Data_arrays_int,Data_arrays_fl,Pwindo,Mesh,
                    Plot_info.variable,
                    Plot_info.plot_type,Plot_info.plot_description,
                    Contour_vals,Profile_vals,
                    &Settings.fg,&Settings.linestyle
                    )))
                     error_msg(err,input_str);
               break;
    case GRAVITY:  if ((err = plot_gravity(Data_vars_int,Data_vars_fl,
                    Data_arrays_int,Data_arrays_fl,Pwindo,Mesh,
                    Plot_info.plot_type,Plot_info.plot_description,
                    Contour_vals,Profile_vals,
					Plot_info.dflt_label1,Plot_info.var_num,
                    &Settings.fg,&Settings.linestyle
                    )))
                     error_msg(err,input_str);
               break;
    case PRFLMRK:  DrawProfileLine(Profile_vals,Pwindo,Data_vars_int[NCOMP]);
               break;
    default:   error_msg(VAR_ERR,input_str);
               break;
    }
    return(0);
}

int read_xyplot_data(FILE* fp, char* str, int* param)
{
    int indx,len;
    int num=0;

    num=fscanf(fp,"%[^A-Za-z]",str);
/*    fprintf(stdout,"log.c read_xyplot_data, param = %i\n",*param); */
    if (fscanf(fp,"%[^ .\r\n\t]",str)!=1) return(STR_ERR);
    if ((indx = name_match(str,xyplots))==-1)
        error_msg(VAR_ERR,str);
    *param = indx;
    len = strlen(str);
    switch(indx) {
    case LGMESH:if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
               if (str[len]!='.') return (SYNTAX_ERR);
               if (fscanf(fp,"%s",&str[len+1])!=1) return(STR_ERR);
               if ((indx = name_match(&str[len+1],lgmeshplots))==-1)
                   error_msg(VAR_ERR,str);
               break;
    case MESH :if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
               if (str[len]!='.') return (SYNTAX_ERR);
               if (fscanf(fp,"%s",&str[len+1])!=1) return(STR_ERR);
               if ((indx = name_match(&str[len+1],meshplots))==-1)
                   error_msg(VAR_ERR,str);
               break;
    case DEFORM:if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
               if (str[len]!='.') return (SYNTAX_ERR);
               if (fscanf(fp,"%s",&str[len+1])!=1) return(STR_ERR);
               if ((indx = name_match(&str[len+1],deformplots))==-1)
                   error_msg(VAR_ERR,str);
               break;
    case STRN_MRKR:
               break;
    case BBOX:
               break;
    case ELLE: 
               break;
    default:   error_msg(VAR_ERR,str);
               break;
    }
    strncpy(Plot_info.variable,str,MAXNAME*2-1);
    strncpy(Plot_info.dflt_label1,&str[len+1], MAXNAME-1);
    Plot_info.var_num= indx;
    Plot_info.plot_type= XYPLOT;
    return(0);
}

int read_arrow_data(FILE* fp, char* str, int* param)
{
    int indx,len;
    int num=0;

    num=fscanf(fp,"%[^A-Za-z]",str);
    if (fscanf(fp,"%[^. \r\n\t]",str)!=1) return(STR_ERR);
    if ((indx = name_match(str,arrowplots))==-1)
        error_msg(VAR_ERR,str);
    *param = indx;
    if (*param != SYB_VELOCITY && *param != ROTATION) {
        len = strlen(str);
        if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
        if (str[len]!='.') return (SYNTAX_ERR);
        if (fscanf(fp,"%s",&str[len+1])!=1) return(STR_ERR);
        if ((indx = name_match(&str[len+1],strain_arrow_plots))==-1)
            error_msg(VAR_ERR,str);
        strncpy(Plot_info.dflt_label1,&str[len+1], MAXNAME-1);
    }
    strncpy(Plot_info.variable,str,MAXNAME*2-1);
    Plot_info.var_num= indx;
    Plot_info.plot_type= ARROWS;
    return(0);
}

int read_contour_data(FILE* fp, char* str, int* param)
{
    char tmpstr[MAXNAME];
    int indx,indx2,len,err=0;
    int num=0;
    valid_terms *set;

    set = NULL;
      /* read the first descriptor */
    num=fscanf(fp,"%[^A-Za-z]",str);
    if (fscanf(fp,"%[^. \r\n\t]",str)!=1) return(STR_ERR);
    if ((indx = name_match(str,contourplots))==-1)
          error_msg(VAR_ERR,str);
    *param = indx2 = indx;
    switch(indx) {
    case SYB_VELOCITY:
            set = velocityplots;
            break;
    case SYB_STRAIN:   /* fall though */
    case SYB_STRESS:
            set = strain_cntr_plots;
            break;
    case LAYER:
            set = layerplots;
            break;
    default:break;
    }
    if (set != NULL) {
        len = strlen(str);
            /* read the fullstop */
        if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
        if (str[len]!='.') return (SYNTAX_ERR);
              /* find the parameter name for plot label */
        if (fscanf(fp,"%s",&str[len+1])!=1) return(STR_ERR);
        len++;
        if ((indx = name_match(&str[len],set))==-1)
              error_msg(VAR_ERR,str);
        strncpy(Plot_info.dflt_label1,&str[len], MAXNAME-1);
    }

    /* 
     * don't break old log files which use Thickness instead of Layer
     */
    if (indx2==THICKNESS) {
        strncpy(tmpstr,"thickness", MAXNAME-1);
        strncpy(Plot_info.dflt_label1,tmpstr, MAXNAME-1);
        sprintf(Plot_info.variable,"%s%s%s","Layer",".",tmpstr);
        indx = name_match("thickness",layerplots);
        *param = name_match("Layer",contourplots);
    }
    else strncpy(Plot_info.variable,str,MAXNAME*2-1);

    Plot_info.var_num= indx;
    Plot_info.plot_type = CNTRS;
    Plot_info.plot_description = 0;
         /* read contour vals */
    err = read_contour_values(fp,str);
    return( err );;
}

int read_contour_values(FILE* fp, char* str)
{
    unsigned char finished;
    int len,indx,i,j,tmp_int;
    int num=0;
    float tmp;
    long pos;

    finished = 0;
    while (!finished) {
      num=fscanf(fp,"%[^A-Za-z#]",str);
      while ((str[0] = getc(fp))=='#') {
          dump_comments(fp);
          num=fscanf(fp,"%[^A-Za-z#]",str);
      }
      if (str[0] !='c') {
        ungetc(str[0],fp);
        num=fscanf(fp,"%s",str);
        finished = 1;
      }
      else {
        if (fscanf(fp,"%[^. \r\n\t]",&str[1])!=1) return(STR_ERR);
        if (!strcmp(str,"cntr")) {
          len = strlen(str);
            /* read the fullstop */
          if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
            /* find the type of contour values on this line */
          if (fscanf(fp,"%s",&str[len+1])!=1) return(STR_ERR);
          len++;
          if ((indx = word_match(&str[len],cntr_terms))==-1)
              error_msg(VAR_ERR,str);
          switch(indx) {
            /* colour */
          case 0:  for (j=0;j<num_items(colour_vals);j++) {
                     num=fscanf(fp,"%[^A-Za-z]",str);
                     if (fscanf(fp,"%[^= \r\n\t]",str)!=1) return(STR_ERR);
                     if ((i = word_match(str,colour_vals))==-1) {
                         pos = ftell(fp); pos -= strlen(str);
                         fseek(fp,pos,SEEK_SET);
                     }
                     else {
                       len = strlen(str);
                       if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
                       if (str[len] != '=') return(STR_ERR);
                       switch( i ) {
                       case 0: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                               Contour_vals[CNTR_MIN_U] = tmp;
                               break;
                       case 1: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                               Contour_vals[CNTR_MAX_U] = tmp;
                               break;
                       case 2: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                               Contour_vals[CNTR_SCL] = tmp;
                               break;
                       case 3: if (fscanf(fp,"%s",str)!=1) return(STR_ERR);
                               if ((tmp_int=name_match(str,colourbar_terms))
                                    ==-1) return(RANGE_ERR);
                               Settings.plot_opts.colour_bar = tmp_int;
                               break;
                       }
                     }
                   }
                   Plot_info.plot_description =
                       Plot_info.plot_description|SHADE;
                   break;
            /* line */
          case 1:  for (j=0;j<num_items(line_vals);j++) {
                     num=fscanf(fp,"%[^A-Za-z]",str);
                     if (fscanf(fp,"%[^= \r\n\t]",str)!=1) return(STR_ERR);
                     if ((i = word_match(str,line_vals))==-1)
                         error_msg(VAR_ERR,str);
                     len = strlen(str);
                     if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
                     if (str[len] != '=') return(STR_ERR);
                     if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                     switch( i ) {
                     case 0: Contour_vals[CNTR_LVL_U] = tmp;
                             break;
                     case 1: Contour_vals[CNTR_STEP_U] = tmp;
                             break;
                     case 2: Contour_vals[MAX_CNTRS] = tmp;
                             break;
                     }
                   }
                   Plot_info.plot_description =
                       Plot_info.plot_description|LINES;
                   break;
            /* limits */
          case 2:  for (j=0;j<num_items(prfl_limits);j++) {
                     num=fscanf(fp,"%[^A-Za-z]",str);
                     if (fscanf(fp,"%[^= \r\n\t]",str)!=1) return(STR_ERR);
                     if ((i = word_match(str,prfl_limits))==-1)
                       error_msg(VAR_ERR,str);
                     len = strlen(str);
                     if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
                     if (str[len] != '=') return(STR_ERR);
                     if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
            /*
             * this line is no longer in use so just read it but
             * take no action
                     switch( i ) {
                     case 0: Contour_vals[CNTR_MIN_U] = tmp;
                             break;
                     case 1: Contour_vals[CNTR_MAX_U] = tmp;
                             break;
                     }
             */
                   }
                   break;
            /* stipple */
          case 3:  for (j=0;j<num_items(stipple_vals);j++) { 
                     num=fscanf(fp,"%[^A-Za-z]",str);
                     if (fscanf(fp,"%[^= \r\n\t]",str)!=1) return(STR_ERR);
                     if ((i = word_match(str,stipple_vals))==-1)
                         error_msg(VAR_ERR,str);
                     len = strlen(str);
                     if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
                     if (str[len] != '=') return(STR_ERR);
                     if (fscanf(fp,"%d",&tmp_int)!=1) return(NUM_ERR);
                     switch( i ) {
                     case 0: Settings.plot_opts.stipple = tmp_int;
                             break;
                     }
                   }
                   break;
          default: break;
          }
        }
        else return(STR_ERR);
      }
    }
    return(0);
}

int read_profile_data(FILE* fp,char* str,int* param)
{
    char tmpstr[MAXNAME];
    int indx,indx2,len,err=0;
    int num=0;
    valid_terms *set;

    set = NULL;
    num=fscanf(fp,"%[^A-Za-z0-9]",str);
       /* read the profile type */
    if (fscanf(fp,"%[^. \r\n\t]",str)!=1) return(STR_ERR);
    if ((indx = name_match(str,profileplots))==-1)
        error_msg(VAR_ERR,str);
    Plot_info.plot_description = indx;
    if (indx==PRFLMRK) Plot_info.plot_type = *param = indx;
    else {
        if (indx==DIM1) Plot_info.plot_type = PRFL1D;
        else Plot_info.plot_type = PRFL2D;
        
        len = strlen(str);
        if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
           /* read the parameter */
        if (fscanf(fp,"%[^. \r\n\t]",str)!=1) return(STR_ERR);
        if ((indx = name_match(str,contourplots))==-1)
            error_msg(VAR_ERR,str);
        *param = indx2 = indx;
        switch(indx) {
        case SYB_VELOCITY:
                set = velocityplots;
                break;
        case SYB_STRAIN:
        case SYB_STRESS:
                set = strain_cntr_plots;
                break;
        case LAYER:
                set = layerplots;
                break;
        case GRAVITY:
                set = gravityplots;
                break;
        default:break;
        }
        if (set != NULL) {
            len = strlen(str);
            if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
            if (str[len]!='.') return (SYNTAX_ERR);
                /* find the parameter name for plot label */
            if (fscanf(fp,"%s",&str[len+1])!=1) return(STR_ERR);
            len++;
            if ((indx = name_match(&str[len],set))==-1)
                error_msg(VAR_ERR,str);
            strncpy(Plot_info.dflt_label1,&str[len], MAXNAME-1);
        }
        /* 
         * don't break old log files which use Thickness instead of Layer
         */
        if (indx2==THICKNESS) {
            strncpy(tmpstr,"thickness", MAXNAME-1);
            strncpy(Plot_info.dflt_label1,tmpstr, MAXNAME-1);
            sprintf(Plot_info.variable,"%s%s%s","Layer",".",tmpstr);
            indx = *param= name_match("Layer",contourplots);
        }
        else strncpy(Plot_info.variable,str,MAXNAME*2-1);

        Plot_info.var_num= indx;
        if (Plot_info.plot_description==DIM1) Plot_info.plot_type= PRFL1D;
        else Plot_info.plot_type= PRFL2D;
    }
        /* read profile vals */
    err = read_profile_values(fp,str,Data_vars_int[NCOMP]);
    return( err );;
}

int read_profile_values(FILE* fp, char* str, int ncomp)
{
    unsigned char finished;
    int len,indx,i,j;
    int num=0;
    float tmp,xref,yref;
    long pos;
    int leng=1,err=0;

/*    XREFM and YREFM usage depends on context */
    xref=Data_vars_fl[XREFM];
    yref=Data_vars_fl[YREFM];
    if (Data_vars_int[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }

    finished = 0;
    while (!finished) {
      num=fscanf(fp,"%[^A-Za-z#]",str);
      if ((str[0] = getc(fp))!='p') {
        ungetc(str[0],fp);
        num=fscanf(fp,"%s",str);
        finished = 1;
      }
      else {
        if (fscanf(fp,"%[^. \r\n\t]",&str[1])!=1) return(STR_ERR);
        if (!strcmp(str,"prfl")) {
          len = strlen(str);
            /* read the fullstop */
          if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
            /* find the type of contour values on this line */
          if (fscanf(fp,"%s",&str[len+1])!=1) return(STR_ERR);
          len++;
          if ((indx = word_match(&str[len],prfl_terms))==-1)
              error_msg(VAR_ERR,str);
          switch(indx) {
            /* line - min,max,scale */
          case 0:  for (j=0;j<num_items(line_terms);j++) {
                     num=fscanf(fp,"%[^A-Za-z]",str);
                     if (fscanf(fp,"%[^= \r\n\t]",str)!=1) return(STR_ERR);
                     if ((i = word_match(str,line_terms))==-1) {
                         pos = ftell(fp); pos -= strlen(str);
                         fseek(fp,pos,SEEK_SET);
                     }
                     else {
                       len = strlen(str);
                       if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
                       if (str[len] != '=') return(STR_ERR);
                       if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                       switch( i ) {
                       case 0: Profile_vals[PRFL_MIN_U] = tmp;
                               break;
                       case 1: Profile_vals[PRFL_MAX_U] = tmp;
                               break;
                       case 2: Profile_vals[PRFL_SCL] = tmp;
                               break;
                       case 3: Profile_vals[PRFL_SHFT] = tmp;
                               break;
                       }
                     }
                   }
                   break;
            /* points */
          case 1:  for (j=0;j<num_items(prfl_pts);j++) { 
                     num=fscanf(fp,"%[^A-Za-z]",str);
                     if (fscanf(fp,"%[^= \r\n\t]",str)!=1) return(STR_ERR);
                     if ((i = word_match(str,prfl_pts))==-1)
                         error_msg(VAR_ERR,str);
                     len = strlen(str);
                     if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
                     if (str[len] != '=') return(STR_ERR);
                     if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                     switch( i ) {
                     case 0: Profile_vals[PRFL_X1] = tmp;
                             break;
                     case 1: Profile_vals[PRFL_Y1] = tmp;
                             break;
                     case 2: Profile_vals[PRFL_X2] = tmp;
                             break;
                     case 3: Profile_vals[PRFL_Y2] = tmp;
                             break;
                     }
                   }

/*   if in spherical coordinates, project to x-y for plotting  */
                   if(ncomp==-1){
                      projectxy_(&Profile_vals[PRFL_X1],&Profile_vals[PRFL_Y1],
                                 &xref,&yref,&leng,&ncomp,&err);
                      projectxy_(&Profile_vals[PRFL_X2],&Profile_vals[PRFL_Y2],
                                 &xref,&yref,&leng,&ncomp,&err);
                   }

                   break;
            /* limits */
          case 2:  for (j=0;j<num_items(prfl_limits);j++) { 
                     num=fscanf(fp,"%[^A-Za-z]",str);
                     if (fscanf(fp,"%[^= \r\n\t]",str)!=1) return(STR_ERR);
                     if ((i = word_match(str,prfl_limits))==-1)
                         error_msg(VAR_ERR,str);
                     len = strlen(str);
                     if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
                     if (str[len] != '=') return(STR_ERR);
                     if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                     switch( i ) {
                     case 0: Profile_vals[LWR_LIM] = tmp;
                             break;
                     case 1: Profile_vals[UPR_LIM] = tmp;
                             break;
                     }
                   }

/*   conversion needed for spherical coordinates, but conversion
     depends on whether limits provided in longitude or latitude */

                   break;
          default: break;
          }
        }
        else return(STR_ERR);
      }
    }
   /* this print statement placed to catch an occasional read failure 
      caused by who knows what.  With the print statement here, it 
      seems not to fail

      fprintf(stdout,"log.c profile X1=%f Y1=%f, X2=%f, Y2=%f\n",
         Profile_vals[PRFL_X1],Profile_vals[PRFL_Y1],
         Profile_vals[PRFL_X2],Profile_vals[PRFL_Y2]); */
    return(0);
}

int create_optional_label(char* str, FILE* fp)
{
    char fontname[80];
    int fonthgt,len,x_offset,y_offset,mode,colour,err;
    float x,y;
#ifndef XSYB
    int bg;
    float width=0,hgt,xtmp,ytmp,margin;
#endif

    if (!(err = read_label(str,fontname,&fonthgt,&colour,&x,&y,fp))) {
        strcpy(Plot_info.fontname,fontname);
/* UPD label coords in fraction of page dim
        if (Settings.orient==PORTRAIT) {
            x *= DEFAULTWIDTH/10.0;
            y *= DEFAULTHEIGHT/10.0;
        }
        else if (Settings.orient==LANDSCAPE) {
            x *= DEFAULTHEIGHT/10.0;
            y *= DEFAULTWIDTH/10.0;
        }
*/
        len = strlen(fontname);
        setfont_(fontname,&len,&fonthgt);
        x_offset = J_LEFT; y_offset = J_BASE; mode = 1;
        len = strlen(Plot_info.curr_label);
#ifndef XSYB
        if (Settings.text_bg) {
            /*
             * clear text background before drawing
             */
            bg = BACKGROUND;
            setpencolor_(&bg);
            margin = (float)fonthgt*0.5*2.54/72.0;
            width = (float)(fonthgt*len)*2.54/72.0/2.0 + margin;
            hgt = (float)fonthgt*2.54/72.0 /* + margin*/;
            xtmp = x - margin/2;
            ytmp = y - margin/2;
            shadrt_(&xtmp,&ytmp,&width,&hgt);
/*
            labelbg_(&x,&y,&Plot_info.curr_label,&len);
*/
        }
#endif
        Settings.fg = colour;
        setpencolor_(&colour);
        drawlabelcm_(&x,&y,&x_offset,&y_offset,&mode,
                                   Plot_info.curr_label,&len);
    }
    return(0);
}

int read_label(char* str,char* fntname,int* fnthgt,int* fg,float* x,float* y,FILE* fp)
{
    unsigned char valid=1;
    char name[80];
    int num=0;

    num=fscanf(fp,"%[^\"#\r\n]",str);
    if (fscanf(fp,"%1s",str)!=1) return(STR_ERR);
    if (str[0]== '"') {
        num=fscanf(fp,"%[^\"#\r\n]",name);
        if (fscanf(fp,"%1s",str)!=1) return(STR_ERR);
        if (str[0]== '"') {
            strncpy(Plot_info.curr_label,name,MAX_LABEL_LEN-1);
            if (fscanf(fp,"%s",str)!=1) return(STR_ERR);
            strcpy(fntname,str);
            if (fscanf(fp,"%d",fnthgt)!=1) return(NUM_ERR);
            if (fscanf(fp,"%d",fg)!=1) return(NUM_ERR);
            if (fscanf(fp,"%f %f",x,y)!=2) return(NUM_ERR);
        }
        else valid=0;
    }
    else valid=0;
    if (!valid) ungetc(str[0],fp);
/*
    if (fscanf(fp,"%s",str)!=1) return(STR_ERR);
    strncpy(Plot_info.curr_label,str,MAX_LABEL_LEN-1);
    if (fscanf(fp,"%f %f",x,y)!=2) return(NUM_ERR);
    if (fscanf(fp,"%s",str)!=1) return(STR_ERR);
    strcpy(fntname,str);
    if (fscanf(fp,"%d",fnthgt)!=1) return(NUM_ERR);
*/
    return( 0 );
}

int read_title(FILE* fp, char* str, char** title)
{
    unsigned char valid=1;
    char name[80];
    int num=0;

    num=fscanf(fp,"%[^\"#\r\n]",str);
    if (fscanf(fp,"%1s",str)!=1) return(STR_ERR);
    if (str[0]== '"') {
        num=fscanf(fp,"%[^\"#\r\n]",name);
        if (fscanf(fp,"%1s",str)!=1) return(STR_ERR);
        if (str[0]== '"') {
            if ((*title=(char *)malloc((strlen(name)+1)*sizeof(char)))
                                                           ==NULL)
                return( MALLOC_ERR );
            strcpy(*title,name);
        }
        else valid=0;
    }
    else valid=0;
    if (!valid) ungetc(str[0],fp);
    return( 0 );
}

int set_records(FILE* fp, file_data* data)
{
    if (fscanf(fp,"%d",&(data->rec_req))!=1) return(NUM_ERR);
    if (fscanf(fp,"%d",&(data->ref_req))!=1) /* may not be present */
        printf("reference rec not found, using default value of 1\n");
    return(0);
}

int read_filename(FILE* fp, char* name)
{
    if (fscanf(fp,"%s",name)!=1) return(STR_ERR);
    return( 0 );
}

int set_cell(FILE* fp)
{
    int row,col,prev;

    if (fscanf(fp,"%d",&row)!=1) return(NUM_ERR);
    if (fscanf(fp,"%d",&col)!=1) return(NUM_ERR);
    if (row<0 || row>Plot_info.max_row || col<0 || col>Plot_info.max_col)
        return(RANGE_ERR);
    prev = Plot_info.curr_cell;
    Plot_info.curr_cell = row*(Plot_info.max_col+1) + col;
    ChangeCell(prev,Plot_info.curr_cell);
    return( 0 );
}

int read_plot_type(fp,iptr)
FILE *fp;
int *iptr;
{
    char tmpstr[80];
    int tmp;
    int num=0;

    *iptr = 0;
    do {
        num=fscanf(fp,"%[^A-Za-z]",tmpstr);
        num=fscanf(fp,"%[^+ \r\n\t]",tmpstr);
        if ((tmp = name_match(tmpstr,plot_types))==-1) return( 1 );
        *iptr = *iptr | tmp;
    } while ((tmpstr[0] = getc(fp))=='+');
    return( 0 );        
}

int read_options(FILE* fp, char* str, input_options* opts, unsigned char initial)
{
    unsigned char finished=0, bool_val=0;
    int indx,len,i,dum,err=0,opt_type=0;
    int num=0;
    float tmp;

    while (!finished) {
      num=fscanf(fp,"%[^A-Za-z#]",str);
      if ((str[0] = getc(fp))=='#') dump_comments(fp);
      else {
        num=fscanf(fp,"%[^= \r\n\t]",&str[1]);
        /*
         * if an option is not recognized, read over it and continue
         */
        /*if ((indx = name_match(str,option_terms))==-1) finished = 1;*/
        indx = name_match(str,option_terms);
        if (indx==-1 && ((indx = name_match(str,keywords))!=-1||feof(fp)))
            finished = 1;
        else {
          len = strlen(str);
          if (fscanf(fp,"%1s",&str[len])!=1) return(STR_ERR);
          if (str[len] != '=') return(STR_ERR);
          opt_type = 0;
          switch( indx ) {
            case O_ORIENT: num=fscanf(fp,"%s",str);
                    if (!strcmp(str,"PORTRAIT")) opts->orient = PORTRAIT;
                    else if (!strcmp(str,"LANDSCAPE"))
                                       opts->orient = LANDSCAPE;
                    else return(RANGE_ERR);
                    break;
            case O_PAPERSZ: num=fscanf(fp,"%s",str);
                    if ((i = name_match(str,paper_terms))==-1)
                         return(RANGE_ERR);
                    opts->paper_size = i;
                    break;
            case O_MARKCELL: num=fscanf(fp,"%d",&i);
                    if (i==1) opts->mark_cell = 1;
                    else if (i==0) opts->mark_cell = 0;
                    else return(RANGE_ERR);
                    break;
            case O_FILEPATH: num=fscanf(fp,"%d",&i);
                    if (i==1) opts->filepath = 1;
                    else if (i==0) opts->filepath = 0;
                    else return(RANGE_ERR);
                    break;
            case O_ROWS: if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    if (i<MINHEIGHT || i>MAXHEIGHT ) return(RANGE_ERR);
                    opts->rows = i;
                    Plot_info.max_row = i-1;
                    break;
            case O_COLS: if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    if (i<MINWIDTH || i>MAXWIDTH ) return(RANGE_ERR);
                    opts->columns = i;
                    Plot_info.max_col = i-1;
                    break;
            case O_FONTHGT: if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    if (i<0) return(RANGE_ERR);
                    opts->dfltfonthgt = i;
                    opt_type = SYB_INT;
#ifdef XSYB
                    if (!initial) {
                      FindFont("Helvetica",Settings.dfltfonthgt);
                      Plot_info.dflt_font = Plot_info.font;
                    }
#endif
                    break;
            case O_XMARG: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->xmargin = tmp;
                    break;
            case O_YMARG: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->ymargin = tmp;
                    break;
            case O_PG_XMARG: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->page_xmargin = tmp;
                    break;
            case O_PG_YMARG: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->page_ymargin = tmp;
                    break;
            case O_ZOOM: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->prev_zoom = opts->zoom;
                    opts->zoom = tmp;
                    opt_type = SYB_FLOAT;
                    break;
            case O_XCENTRE: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->xcentre = tmp;
                    opt_type = SYB_FLOAT;
                    break;
            case O_YCENTRE: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->ycentre = tmp;
                    opt_type = SYB_FLOAT;
                    break;
            case O_LINEWDTH: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->linewidth = tmp;
                    if (!initial) setlinewidth_(&tmp);
                    opt_type = SYB_FLOAT;
                    break;
            case O_PGSCL: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->page_scale = tmp;
                    break;
            case O_TL0: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->physical.tl0 = tmp;
                    break;
            case O_RHOC: if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->physical.rhoc = tmp;
                    break;
            case O_RHOM:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->physical.rhom = tmp;
                    break;
            case O_Z0:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->physical.z0 = tmp;
                    break;
            case O_ELEV0:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->physical.elev0 = tmp;
                    break;
            case O_TMAX:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->physical.tmax = tmp;
                    break;
            case O_BGAM0:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->physical.bgam0 = tmp;
                    break;
            case O_BGAM1:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->physical.bgam1 = tmp;
                    break;
            case O_NDIVR:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->contour.ndivr = tmp;
                    break;
            case O_S0:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->contour.s0 = tmp;
                    break;
            case O_DELS:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->contour.dels = tmp;
                    break;
            case O_DELOM:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->contour.delom = tmp;
                    break;
            case O_NX3:if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    if (i<0) return(RANGE_ERR);
                    opts->plot_opts.nx3 = i;
                    opt_type = SYB_INT;
                    break;
            case O_MP:if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->plot_opts.mp = i;
                    opt_type = SYB_INT;
                    break;
            case O_MPE:if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->plot_opts.mpe = i;
                    opt_type = SYB_INT;
                    break;
            case O_NP:if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->plot_opts.np = i;
                    opt_type = SYB_INT;
                    break;
            case O_CNTRPLT:if (read_plot_type(fp,&i)) return(VAR_ERR);
                    opts->plot_opts.contour_plot = i;
                    opt_type = SYB_INT;
                    break;
            case O_PRFLPTS:if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->plot_opts.profile_pts = i;
                    opt_type = SYB_INT;
                    break;
            case O_FLIP:if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->plot_opts.flip = i;
                    opt_type = SYB_INT;
                    bool_val = (i==1);
                    break;
            case O_DBLE:if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->plot_opts.dble = i;
                    opt_type = SYB_INT;
                    bool_val = (i==1);
                    break;
            case O_LABEL:if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->plot_opts.label = i;
                    opt_type = SYB_INT;
                    bool_val = (i==1);
                    break;
            case O_STIP:if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->plot_opts.stipple = i;
                    opt_type = SYB_INT;
                    break;
            case O_SOLNROT:if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->plot_opts.solution_rot = i;
                    opt_type = SYB_INT;
                    break;
            case O_COLBAR:if (fscanf(fp,"%s",str)!=1) return(NUM_ERR);
                    if ((i = name_match(str,colourbar_terms))==-1)
                         return(RANGE_ERR);
                    opts->plot_opts.colour_bar = i;
                    opt_type = SYB_INT;
                    break;
            case O_FOREGRND:if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->fg = i;
                    if (!initial) setpencolor_(&i);
                    opt_type = SYB_INT;
                    break;
            case O_LINESTYLE:if (fscanf(fp,"%s",str)!=1) return(NUM_ERR);
                    if ((i = name_match(str,linestyle_terms))==-1)
                         return(RANGE_ERR);
                    opts->linestyle = i;
                    /* X window may not exist yet */
                    dum=1;
                    if (!initial) dashln_(&opts->linestyle,&dum);
                    opt_type = SYB_INT;
                    break;
            case O_COLMAP:if (fscanf(fp,"%s",str)!=1) return(NUM_ERR);
                    if ((i = name_match(str,colmap_terms))==-1)
                         return(RANGE_ERR);
                    opts->colourmap = i;
                    opt_type = SYB_INT;
                    break;
            case O_VISCMIN:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->viscmin = tmp;
                    opt_type = SYB_FLOAT;
                    break;
            case O_VISCMAX:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->viscmax = tmp;
                    opt_type = SYB_FLOAT;
                    break;
            case O_SEMIN:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->semin = tmp;
                    opt_type = SYB_FLOAT;
                    break;
            case O_SEMAX:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    opts->semax = tmp;
                    opt_type = SYB_FLOAT;
                    break;
            case O_CLIPTOCELL: if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->clip_to_cell = i;
                    opt_type = SYB_INT;
                    break;
            case O_RESCALE: if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->rescale = i;
                    opt_type = SYB_BOOL;
                    bool_val = (i==1);
                    break;
            case O_VERBOSE: if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->verbose = i;
                    opt_type = SYB_INT;
                    break;
            case O_ELLE: if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->elle = i;
                    opt_type = SYB_BOOL;
                    bool_val = (i==1);
                    break;
            case O_TEXTBG: if (fscanf(fp,"%d",&i)!=1) return(NUM_ERR);
                    opts->text_bg = i;
                    opt_type = SYB_BOOL;
                    bool_val = (i==1);
                    break;
            case O_TIC:if (fscanf(fp,"%s",str)!=1) return(NUM_ERR);
                    if ((i = name_match(str,ticmark_terms))==-1)
                         return(RANGE_ERR);
                    opts->plot_opts.ticmark = i;
                    opt_type = SYB_INT;
                    break;
            case O_AXISXMIN:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    Pwindo[ULEFT] = tmp;
                    if (opts->rescale) opts->rescale= 0;
                    break;
            case O_AXISXMAX:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    Pwindo[URGHT] = tmp;
                    if (opts->rescale) opts->rescale= 0;
                    break;
            case O_AXISYMIN:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    Pwindo[UBASE] = tmp;
                    if (opts->rescale) opts->rescale= 0;
                    break;
            case O_AXISYMAX:if (fscanf(fp,"%f",&tmp)!=1) return(NUM_ERR);
                    Pwindo[UTOP] = tmp;
                    if (opts->rescale) opts->rescale= 0;
                    break;
            /*default:return(VAR_ERR);*/
            default:warning_msg(VAR_ERR,str);
                    num=fscanf(fp,"%[^A-Za-z0-9]",str);
                    num=fscanf(fp,"%[^ \r\n\t]",str);
                    break;
          }
#ifdef XSYB
          if (!initial && indx<=O_RESCALE && indx>O_ORIENT) 
            switch(opt_type) {
            case SYB_INT  : Save_option((XtPointer)&i,indx,SYB_INT);
                            break;
            case SYB_FLOAT: Save_option((XtPointer)&tmp,indx,SYB_FLOAT);
                            break;
            case SYB_BOOL : Save_option((XtPointer)&bool_val,indx,SYB_BOOL);
                            break;
            default:        break;
            }
#endif
        }
      }
    }
    return( err );
}

int dump_comments(FILE *fp)
{
    unsigned char done=0;
    int c;

    while (!done) {
        c = getc(fp);
        done = (c=='\n'||c==EOF);
    }
    return(0);
}
