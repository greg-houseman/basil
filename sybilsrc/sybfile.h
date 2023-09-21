/*--------------------------------------------------------------------
 *    Basil / Sybil:   file.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
#ifndef _sybfile_h
#define _sybfile_h

#ifdef SUN
int file_open();
int read_data();
#else
#include "filedat.h"
#ifdef __cplusplus
extern "C" {
#endif
int file_open( file_data **data, char *name, int max, char *mode );
int read_data( FILE *fp,char *string_vars,int *int_vars,
               float *fl_vars, int **arraysi,float **arraysf,
               float *xmin,float *xmax, float *ymin,float *ymax,
               int *current, int req,int *maxrec,int rotate );
#ifdef __cplusplus
}
#endif
#endif

#endif

