/*--------------------------------------------------------------------
 *    Basil / Sybil:   arrays.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
#ifndef _S_arrays_h
#define _S_arrays_h

#define FNAME_W_START   0
#define FNAME_W_SIZE    16
#define DATE_START      16
#define DATE_SIZE       16
#define FNAME_R_START   32
#define FNAME_R_SIZE    16
#define COMMENTS_START  48
#define COMMENTS_SIZE   32

void init_arrays();
void clear_arrays();
void clear_alloc_arrays();
void find_max_min(float *, float *, int, float *, float *, float *, float *);

extern char Err_str[80], String_vars[REC_1_CNT];
extern int Data_vars_int[REC_2_CNT], *Data_arrays_int[MAX_INT_ARRAYS];
extern float Data_vars_fl[REC_3_CNT], *Data_arrays_fl[MAX_FL_ARRAYS];
#endif
