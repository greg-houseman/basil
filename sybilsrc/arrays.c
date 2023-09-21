/*--------------------------------------------------------------------
 *    Basil / Sybil:  arrays.c  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "cmndefs.h"
#include "data.h"
#include "arrays.h"

char Err_str[80], String_vars[REC_1_CNT];
int Data_vars_int[REC_2_CNT], *Data_arrays_int[MAX_INT_ARRAYS];
float Data_vars_fl[REC_3_CNT], *Data_arrays_fl[MAX_FL_ARRAYS];

void find_max_min(xvec,yvec,dim,xmin,xmax,ymin,ymax)
int dim;
float *xvec,*yvec,*xmin,*xmax,*ymin,*ymax;
{
    /*
     * find the min and max for Pwindo
     */
    int i;
    float min,max;

    min = max = xvec[0];
    for (i=1;i<dim;i++) {
        min=MIN(min,xvec[i]);
        max=MAX(max,xvec[i]);
    }
    *xmin = min; *xmax = max;
    min = max = yvec[0];
    for (i=1;i<dim;i++) {
        min=MIN(min,yvec[i]);
        max=MAX(max,yvec[i]);
    }
    *ymin = min; *ymax = max;
}

void init_arrays()
{
    int i;

    for (i=0;i<MAX_INT_ARRAYS;i++) Data_arrays_int[i]=NULL;
    for (i=0;i<MAX_FL_ARRAYS;i++) Data_arrays_fl[i]=NULL;
    for (i=0;i<REC_2_CNT;i++) Data_vars_int[i]=0;
    for (i=0;i<REC_3_CNT;i++) Data_vars_fl[i]=0;
}

void clear_arrays()
{
    int i;

    for (i=0;i<MAX_INT_ARRAYS;i++)
        if (Data_arrays_int[i]!=NULL) {
            free(Data_arrays_int[i]);
            Data_arrays_int[i]=NULL;
        }
    for (i=0;i<MAX_FL_ARRAYS;i++)
        if (Data_arrays_fl[i]!=NULL) {
            free(Data_arrays_fl[i]);
            Data_arrays_fl[i]=NULL;
        }
    for (i=0;i<REC_2_CNT;i++) Data_vars_int[i]=0;
    for (i=0;i<REC_3_CNT;i++) Data_vars_fl[i]=0;
}

void clear_alloc_arrays()
{
    int i;

    for (i=0;i<MAX_INT_ARRAYS;i++)
        if (Data_arrays_int[i]!=NULL) {
            free(Data_arrays_int[i]);
            Data_arrays_int[i]=NULL;
        }
    for (i=0;i<MAX_FL_ARRAYS;i++)
        if (Data_arrays_fl[i]!=NULL) {
            free(Data_arrays_fl[i]);
            Data_arrays_fl[i]=NULL;
        }
}
