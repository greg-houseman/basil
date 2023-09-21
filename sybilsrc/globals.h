/*--------------------------------------------------------------------
 *    Basil / Sybil:   globals.h   1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#include "data.h"
#include "cmndefs.h"

extern int Mesh[MESH_ENTRIES];
extern float Pwindo[PWINDO_ENTRIES];
extern int Data_vars_int[REC_2_CNT], *Data_arrays_int[MAX_INT_ARRAYS];
extern float Data_vars_fl[REC_3_CNT], *Data_arrays_fl[MAX_FL_ARRAYS];
extern float Contour_vals[MAXCNTRVALS], Profile_vals[MAXPRFLVALS];
extern char Err_str[256];
