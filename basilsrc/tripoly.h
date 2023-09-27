#ifndef _tripoly_h
#define _tripoly_h

/* #define SINGLE */

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include "triangle.h"

struct flagvals {
    unsigned char midptnodes;
    unsigned char bndpts;
    int quality;
    float area;
};

#ifndef SUN

#ifdef __cplusplus
extern "C" {
#endif
int tripoly(struct triangulateio *out, struct triangulateio *inp, char *infile,
            struct flagvals *flags, char *optionaltxt,
            int *npolyp,int **ipoly);
int tripolypts(struct triangulateio *out,
            struct flagvals *flags,
            int num, float *xvals, float *yvals );
int tripolydata( struct triangulateio *out,struct triangulateio *inp,
               struct flagvals *flags,
                 int num, int numseg, int num_pt_attr, int *segs,
                 float *xvals, float *yvals,
                 float *pt_attributes,
                 float *rgnattribs, int numrgns );
int readpolyinput( struct triangulateio *inp, char *infile,
             struct flagvals *flags, char *optionaltxt,
             int *maxindx, int **ielle,int *nbm,int *internal_bnd);
int readelleinput( char *infile, char *optionaltxt,
                   int *maxindx, int **ielle);
void initio(struct triangulateio *io);
void cleanio(struct triangulateio *io);
int dump_comments(FILE *fp);
int WritePoly(struct triangulateio *out,char *name);
#ifdef __cplusplus
}
#endif
#else /* SUN */
int tripoly();
int tripolypts();
int tripolydata();
int readpolyinput();
int readelleinput();
void initio();
void cleanio();
#endif /* SUN */

#endif
