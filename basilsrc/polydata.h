#ifndef _polydata_h
#define _polydata_h
#include "triangle.h"

typedef struct{
    int numpts;
    double *ppts;
} polydata;

#ifdef __cplusplus
extern "C" {
#endif
int getpolydata(struct triangulateio *io,polydata *polys,int nr);
#ifdef __cplusplus
}
#endif

#endif
