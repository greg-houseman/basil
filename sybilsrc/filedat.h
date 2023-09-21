/*--------------------------------------------------------------------
 *    Basil / Sybil:   filedat.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
#ifndef _filedat_h
#define _filedat_h
#include "cmndefs.h"

typedef struct {
	char fname[SYB_FILENAME_MAX+1];
	FILE *fp;
	int rec_req;
	int rec_curr;
	int ref_req;
	int ref_curr;
	int rec_max;
} file_data;

#endif
