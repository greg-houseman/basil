/*--------------------------------------------------------------------
 *    Basil / Sybil:   trimesh.c  1.1  13 May 1997
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

/*                                                                           */
/*  (trimesh.c)                                                              */
/*     program to read input file of nodes and polygon                       */
/*        and create a delauney mesh                                         */
/*                                                                           */
/*****************************************************************************/

/* If SINGLE is defined when triangle.o is compiled, it should also be       */
/*   defined here.  If not, it should not be defined here.                   */

/* #define SINGLE */

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "triangle.h"
#include "trimesh.h"
#include "tripoly.h"
#include "poly.h"
#include "polydata.h"

#ifndef _STDLIB_H_
extern void *malloc();
extern void free();
extern void exit(int);
extern int WritePoly(struct triangulateio *out,char *name);
#endif /* _STDLIB_H_ */

struct node {
    int id;
    struct node *next;
};

struct angle_node {
    int id;
    double ang;
    struct angle_node *next;
};

struct gen_fault {
    int id;
    int end1;
    int end2;
    struct node *ordered_nodes;
    struct gen_fault *next;
};

int attribute_set(char *infile,char *key,int max,
                  float *rgnattribs,float *dflt);
int unode_attribute(char *infile,char *key,
                        int maxunodes,float **unodeattribs,float *dflt);
int unode_data(char *infile,char *key,
                        int *maxunodes,float **unodedata);
int vcfrompoly_( float *vhb,float *ex,float *ey,int *lem,int *nor,
             int *ivv,int *nup,int *mesh_ne,
             char *polyfile,char *comments,int *ierr);
void  print_triangulate_data(struct triangulateio *out);
struct gen_fault *newfltnode(int id);
unsigned char in_flt_list(struct gen_fault *list, int id);
void insert_in_flt_list(struct gen_fault **list,struct gen_fault *node);
void free_flt_list(struct gen_fault *list);
int findfltpts(struct gen_fault **list,struct triangulateio *out,
               int *nfp,int *newvert,int *newseg,int *newibc,int *ibcmoffv);
void countpbpts(int *nn,int *countpb,int pbxy,float *ex,float *ey,float pdist);
unsigned char is_bnd_pt(struct triangulateio *out,int index);
unsigned char is_inex_bnd_pt(int id,int *iseg,int *label,int nseg);
unsigned char is_flt_intersect_pt(int id,int fltid,struct gen_fault *list);
unsigned char in_list(struct node *list, int id);
struct node *newnode(int id);

void bndextlimits_(int *,int *);
void bndfltlimits_(int *,int *);
int  CartesianToPolar(REAL, REAL, REAL, double *, double *);
int createafault(struct triangulateio *out,
                 int *mpstart, int *vstart, int *fstart,
                 int *lem, float *ex, float *ey,
                 int *ne, int *ibc, int *ibcngh, int *ibctyp,
                 int *ifbc1,int *ifbc2,int *iseg,
                 int nbpfinal,int *nbpext,int *nseg,int *vxtotal,int *mptotal,
                 struct gen_fault *flt);
int  crossingstest_(double *,int *, double *, int *);
void findNeighbourPts(int node, int prev, int next,
                      struct triangulateio *out,
                        int *cnt, int *nbnodes, int maxnbs);
void findTriangles(int start,int *nbs,int cnt,int *lem,int ne,
                                   int *tri_nbs, int *tricnt);
int findTriangV(int count, int start, int end, int *lem,
                  int ne, int *tri_1, int *tri_2);
int findTriangM(int count, int midp, int nabr, int *lem,
                  int ne, int *tri_3);
void findbbox_(float *,float *, float *, int *);
int  find_maxnbs(float *ex,float *ey,int *lem,int ne,int nn);
void insert_node(struct node *nnode,struct node **list,REAL *coords);
void print_lists(struct gen_fault *list);
int  reallocateArrays();
/* these are defined and used in basmain.c - make new header?  */
extern int AllocIntArray(int cnt, int **addr, int dflt);
extern int AllocFlArray(int cnt, float **addr, float dflt);
extern int AllocDblArray(int cnt, double **addr, double dflt);
#ifndef MALLOC_ERR
#define MALLOC_ERR = -1
#endif

void showbasilvar();

/*****************************************************************************/
/*                                                                           */
/*  trimesh()   Create and refine a mesh.                                       */
/*                                                                           */
/*****************************************************************************/

int trimesh( int *iflt, int *nup,int *nn,int *nmp,int *nbp, int *nnnof, int *nmpnof,
             int  *ne,int *nfp,int *nfpf3,int *nseg,int *ivis,int *maxnbs,
             int *nellep,int **ielle_adr,int **ipolyn_adr,
             int **lem_adr,int **lemneighbor_adr,
             int **ibc_adr,int **ibcnb_adr,int **ibctyp_adr,int **ifbc1_adr,int **ifbc2_adr,
             int **iseg_adr,
             float **ex_adr,float **ey_adr,float **vhb_adr,
             float *xlen, float *ylen,
             float area,float quality,float *visdflt,float *se,int pbxy,int *nbpext,
             float pdist,char *infile, char *comments)
{
  struct triangulateio out, inp;
  struct flagvals triflags;
  struct gen_fault *flts=NULL;
  int err=0;
  int i,j,k,k1;
//int *ncnt=NULL;
  int *lem,*lemneighbor,*ibc,*ibcngh,*ibctyp,*ifbc1,*ifbc2,*polyn,*iseg;
  float *ex,*ey;
  float *vhb=NULL;
  extern int pntonsegment_();
  float eps;
  int *numnbs=NULL;
//int b,c,d,mpstart;
  int a;
  int newvert, newseg, newibc, ibcmoffv, newbpts, fvert, nbpfinal;
  int bvert,bmidpts;
  int aj,ajac,mac,macnv,lemnela,mpk,mpk0;
  float dflt_visc,dflt_sexp,dflt_uvisc,curr;
  float *unodedata=NULL,*unodeattribs=NULL,*rgnattribs=NULL;
  float min_dist, dist, c_x, c_y, val;
  int maxrgns=0,maxunodes, dist_set=0, vxtotal=0, mptotal=0, vnofault, mnofault, ibcadded;
  int clrptr=0,countpb=0;

//int visc_var=0, sexp_var=0;

  /* default area is  0.025*0.025 (area of tri with minNodeSep sides) */
  triflags.area = area; /* area of tri with minNodeSep sides */
  /*triflags.area = 0;*/
  triflags.quality = quality;
  /*triflags.quality = 5;*/
  triflags.midptnodes = 1;

//  triflags.bndpts=0 suppresses all changes to boundary nodes in remesh operation
  if(pbxy) {
     triflags.bndpts = 0;
     if(pbxy==1)
     printf("Periodic boundaries (X) are enabled. pdist = %f\n",pdist);
     if(pbxy==2)
     printf("Periodic boundaries (Y) are enabled. pdist = %f\n",pdist);
  }
  else triflags.bndpts = 1;

  initio( &inp );
  initio( &out );

  if (tripoly(&out,&inp,infile,&triflags,comments,nellep,ielle_adr)!=0)
      return(1);
//WritePoly(&inp, "in.poly");

  /**********************************/
  /******* Setup Basil Arrays *******/
  /**********************************/

  *nbp = 0;
  *nseg = 0;
  *nfp = 0;
  *nup=out.numberofpoints;
  *nn=out.numberofpoints-out.numberofedges;
  *nmp=out.numberofedges;
  newvert=newseg=newibc=newbpts=ibcmoffv=0;
  bvert=bmidpts=fvert=0;

//  add up space needed for fault points
  if (*iflt!=0 && out.numberofsegments>0) {
    findfltpts(&flts,&out,nfp,&newvert,&newseg,&newibc,&ibcmoffv); 
    if(!(*nfp)){
        printf("warning:fault type is non-zero but no internal faults found, iflt reset to 0\n");
        *iflt=0;
    }
  }
    *nfpf3=*nfp;
    newbpts=newseg+newvert;
 // printf("trimesh: ex,ey array dimensions nup = %i newbpts = %i nupfinal = %i\n",*nup,newbpts,*nup+newbpts);

 // establish the EX and EY arrays
  if ((err=AllocFlArray((*nup+newbpts),ex_adr,0.0))!=0){
      printf("Problem allocating EX array\n");
      return(err);
  }
  ex = *ex_adr;
  if ((err=AllocFlArray((*nup+newbpts),ey_adr,0.0))!=0){
      printf("Problem allocating EY array\n");
      return(err);
  }
  ey = *ey_adr;

  eps = 1.0e-6;
  for (i = 0; i < *nup; i++) {
    j=i;
    if (j>=*nn) j+=newvert;
    ex[j] = (float) out.pointlist[i*2];
    ey[j] = (float) out.pointlist[i*2+1];
    if( is_bnd_pt(&out,i) && out.pointmarkerlist[i] != 99) {
      (*nbp)++;
      if (i<*nn) bvert++;
      else bmidpts++;
    }
  }
  //  bvert and bmidpts = no of vertex and midpoint nodes on ext. boundary
  *nbpext=bvert+bmidpts;  

  // if periodic boundaries in use, count the nodes prior to allocation of ifbc1 etc
  if(pbxy>0) countpbpts(nn,&countpb,pbxy,ex,ey,pdist);
  *nfp = *nfp + countpb;    // allow space in fault arrays for pb points

//nbpfinal includes all nodes on external boundary and internal faults
  nbpfinal = *nbp+newibc;
//printf("trimesh: ibc array dimensions nbp = %i newibc = %i nbpfinal = %i\n",*nbp,newibc,nbpfinal);

  if ((err=AllocIntArray(nbpfinal,ibc_adr,-1))!=0){
      printf("Problem allocating IBC array\n");
      return(err);
  }
  ibc = *ibc_adr;
  if ((err=AllocIntArray((nbpfinal*2),ibcnb_adr,0))!=0){
      printf("Problem allocating IBCNGH array\n");
      return(err);
  } ibcngh = *ibcnb_adr;

   /* prefer default for zero stress is -1 */
  if ((err=AllocIntArray((nbpfinal*2),ibctyp_adr,-1))!=0){
      printf("Problem allocating IBCTYP array\n");
      return(err);
  } ibctyp = *ibctyp_adr;

  ifbc1=ifbc2=0;
  if (*iflt>0) {
    if ((err=AllocIntArray(*nfp,ifbc1_adr,0))!=0){
        printf("Problem allocating IFBC array\n");
        return(err);
    } ifbc1 = *ifbc1_adr;

    if ((err=AllocIntArray(*nfp,ifbc2_adr,0))!=0){
        printf("Problem allocating IFBC2 array\n");
        return(err);
    } ifbc2 = *ifbc2_adr;
  }

  /***** count the number of nodes connected to each vertex ******/
  numnbs = (int *) malloc((*nn+newvert) * sizeof(int));
  for (i = 0; i < *nn+newvert; i++) numnbs[i]=0;

  *ne=out.numberoftriangles;
  if ((err=AllocIntArray((*ne*6),lem_adr,0))!=0){
      printf("Problem allocating LEM array\n");
      return(err);
  }
  lem = *lem_adr;
  if ((err=AllocIntArray((*ne*3),lemneighbor_adr,0))!=0){
      printf("Problem allocating LEM array\n");
      return(err);
  }
  lemneighbor = *lemneighbor_adr;
  /*
   * ivv    visc (vhb[i,1-7])    se (vhb[i,8])
   *  0       constant           constant(1)    No vhb array
   *  1       variable           constant(1)
   *  2       constant           constant(>1) 
   *  3       variable           constant(>1)
   *  4       variable           variable
   */
  if ( *ivis>0 && *vhb_adr==NULL ) {
        /*
         * allocate and initialise the vhb array
         * the input file contains VIS settings
         */
      if ((err=AllocFlArray((*ne * 8),vhb_adr,*visdflt))!=0){
          printf("Problem allocating VHB array\n");
          return(err);
      }
      vhb = *vhb_adr;
      for (i = 0; i < *ne; i++) vhb[i*8+7]=*se;
  }
  /*
   * the triangle attribute is the region identifier
   * save it and set triangle attributes, if found
   */
  if (out.triangleattributelist!=NULL) {
      if ((err=AllocIntArray(*ne,ipolyn_adr,0))!=0){
          printf("Problem allocating POLY array\n");
          return(err);
      }
      polyn = *ipolyn_adr;
      maxrgns = (int)out.triangleattributelist[0];
      for (i = 0; i < *ne; i++) {
          k = (int)out.triangleattributelist[i];
          polyn[i] = k;
          if (k > maxrgns) maxrgns = k;
      }
      if ((rgnattribs= (float *)malloc((maxrgns+1) *
                                    sizeof(float)))==NULL) return(1);
      dflt_visc = *visdflt;
      for (i=0;i<maxrgns;i++) rgnattribs[i] = dflt_visc;
      if (attribute_set(infile,VISCOSITY_KEY,
                        maxrgns+1,rgnattribs,&dflt_visc)) {
          if ( *vhb_adr==NULL ) {
            /*
             * allocate and initialise the vhb array
             */
              if ((err=AllocFlArray((*ne * 8),vhb_adr,dflt_visc))!=0){
                  printf("Problem allocating VHB array\n");
                  return(err);
              }
              vhb = *vhb_adr;
              for (i = 0; i < *ne; i++) vhb[i*8+7]=1.0;
              *ivis = 1;
          }
          curr = dflt_visc;
          for (i = 0; i < *ne; i++) {
              k = (int)out.triangleattributelist[i];
              polyn[i] = k;
              for (j = 0; j < 7; j++) {
                  vhb[i*8+j] = rgnattribs[k];
                  if (fabs((double)rgnattribs[k]-curr)>eps) {
                    curr =  rgnattribs[k];
          //        visc_var = 1;
                  }
              }
          } /* viscosity */
          *visdflt = dflt_visc;
      }
      dflt_sexp = *se;
      for (i=0;i<maxrgns;i++) rgnattribs[i] = dflt_sexp;
      if (attribute_set(infile,SE_KEY,
                        maxrgns+1,rgnattribs,&dflt_sexp)) {
          if ( *vhb_adr==NULL ) {
            /*
             * allocate and initialise the vhb array
             */
              if ((err=AllocFlArray((*ne * 8),vhb_adr,1.0))!=0){
                  printf("Problem allocating VHB array\n");
                  return(err);
              }
              vhb = *vhb_adr;
              for (i = 0; i < *ne; i++) vhb[i*8+7]=dflt_sexp;
              *ivis = 1;
          }
          curr = dflt_sexp;
          for (i = 0; i < *ne; i++) {
              k = (int)out.triangleattributelist[i];
              polyn[i] = k;
              vhb[i*8+7] = rgnattribs[k];
              if (fabs((double)rgnattribs[k]-curr)>eps) {
                curr =  rgnattribs[k];
     //         sexp_var = 1;
              }
          }
          *se = dflt_sexp;
      } /* stress/strain exponent */

      if (unode_data(infile,UNODES_KEY,
                        &maxunodes,&unodedata)) {
     //   visc_var=0;
          dflt_uvisc=1.0;
          if (unode_attribute(infile,U_VISCOSITY_KEY,
                        maxunodes,&unodeattribs,&dflt_uvisc)) {
              if ( *vhb_adr==NULL ) {
                /*
                 * allocate and initialise the vhb array
                 */
                  if ((err=AllocFlArray((*ne * 8),vhb_adr,1.0))!=0){
                      printf("Problem allocating VHB array\n");
                      return(err);
                  }
                  vhb = *vhb_adr;
                  for (i = 0; i < *ne; i++) vhb[i*8+7]=1.0;
                  *ivis = 1;
              }
              curr= dflt_uvisc;
              for (i = 0; i < *ne; i++) {
                  k = (int)out.triangleattributelist[i];
                  polyn[i] = k;
                  for (j = 0,c_x=0.0,c_y=0.0; j < 3; j++) {
                    a = out.trianglelist[i*6+j];
                    c_x += out.pointlist[a*2];
                    c_y += out.pointlist[a*2+1];
                  }
                  c_x /= 3.0;
                  c_y /= 3.0;
                  min_dist=0.0;
	//	  b=-1;
                  for (j = 0, dist_set=0,val=dflt_uvisc; j < maxunodes; j++) {
                    if (unodedata[j+maxunodes*2]==k) {
                        dist = (c_x-unodedata[j])*(c_x-unodedata[j]) +
                        (c_y-unodedata[j+maxunodes])*(c_y-unodedata[j+maxunodes]);
                        if (!dist_set || (dist<min_dist)) {
                            min_dist = dist;
                            val = unodeattribs[j];
                        }
                    }
                  }
                  for (j = 0; j<7; j++) vhb[i*8+j]=val;
                  if (fabs(val-curr)>eps) {
                    curr =  rgnattribs[k];
       //           visc_var = 1;
                  }
              }
              *visdflt = dflt_uvisc;
          }                            // unode visc
    }                                  // unodes
  }                                    // if triangle attributes
/*for (i = 0; i < *ne; i++) {
    for (j = 0; j < 3; j++) {
      a=out.trianglelist[i*6+j];
      numnbs[a]++;
      lem[i*6+j]=out.trianglelist[i*6+j] + 1;
      lemneighbor[i*3+j]=out.neighborlist[i*3+j] + 1;
      if (is_bnd_pt(&out,a)) {
        k=0;
        while (k< *nbp && ibc[k]!=-1 && ibc[k]!=a+1) k++;
        if (ibc[k] == -1) ibc[k] = a+1;
        switch (j) {
 
 //       ibcngh[k] is vertex in clockwise direction,
 //       ibcngh[k+nbp] is vertex in anticlockwise direction,
  
          case 0:
          if (is_bnd_pt(&out,out.trianglelist[i*6+1]))
            ibcngh[k+nbpfinal]=out.trianglelist[i*6+1] + 1;
          if (is_bnd_pt(&out,out.trianglelist[i*6+2]))
            ibcngh[k]=out.trianglelist[i*6+2] + 1;
          break;
          case 1:
          if (is_bnd_pt(&out,out.trianglelist[i*6+2]))
            ibcngh[k+nbpfinal]=out.trianglelist[i*6+2] + 1;
          if (is_bnd_pt(&out,out.trianglelist[i*6]))
            ibcngh[k]=out.trianglelist[i*6] + 1;
          break;
          case 2:
          if (is_bnd_pt(&out,out.trianglelist[i*6]))
            ibcngh[k+nbpfinal]=out.trianglelist[i*6] + 1;
          if (is_bnd_pt(&out,out.trianglelist[i*6+1]))
            ibcngh[k]=out.trianglelist[i*6+1] + 1;
          break;
        }
      }
    }
  }

    //  vertices have been loaded in ibc
//  printf("ibc  array follows, nbp = %i\n",*nbp);
//  for (k=0; k < *nbp; k++){printf("k = %i ibc=%i ",k,ibc[k]);}
    mpstart = bvert+fvert+newvert;
      for (i = 0; i < *ne; i++) {
    b=out.trianglelist[i*6+4];
    lem[i*6+3]=b + 1 + newvert;
    c=out.trianglelist[i*6+5];
    lem[i*6+4]=c + 1 + newvert;
    d=out.trianglelist[i*6+3];
    lem[i*6+5]=d + 1 + newvert;
    ncnt[d]=1;
    ncnt[b]=1;
    ncnt[c]=1;
    if (is_bnd_pt(&out,b)) {
        k=mpstart;
        while (k<  nbpfinal && ibc[k]!=-1 && ibc[k]!=b+1+newvert) k++;
        if (ibc[k] == -1) ibc[k] = b+1+newvert;
            if (is_bnd_pt(&out,out.trianglelist[i*6+2]))
              ibcngh[k]=out.trianglelist[i*6+2] + 1;
            if (is_bnd_pt(&out,out.trianglelist[i*6]))
              ibcngh[k+nbpfinal]=out.trianglelist[i*6] + 1;
    }
    if (is_bnd_pt(&out,c)) {
        k=mpstart;
        while (k< nbpfinal && ibc[k]!=-1 && ibc[k]!=c+1+newvert) k++;
            if (ibc[k] == -1) ibc[k] = c+1+newvert;
            if (is_bnd_pt(&out,out.trianglelist[i*6]))
              ibcngh[k]=out.trianglelist[i*6] + 1;
            if (is_bnd_pt(&out,out.trianglelist[i*6+1]))
              ibcngh[k+nbpfinal]=out.trianglelist[i*6+1] + 1;
    }
    if (is_bnd_pt(&out,d)) {
        k=mpstart;
        while (k< nbpfinal && ibc[k]!=-1 && ibc[k]!=d+1+newvert) k++;
        if (ibc[k] == -1) ibc[k] = d+1+newvert;
            if (is_bnd_pt(&out,out.trianglelist[i*6+1]))
              ibcngh[k]=out.trianglelist[i*6+1] + 1;
            if (is_bnd_pt(&out,out.trianglelist[i*6+2]))
              ibcngh[k+nbpfinal]=out.trianglelist[i*6+2] + 1;
    }
  }
*/
   //   midpoints now loaded in ibc
/*  printf("ibc  array follows, nbp = %i\n",nbpfinal);
    for (k=0; k < nbpfinal; k++){printf("k = %i ibc=%i ",k,ibc[k]);} */
/*
for (i = 0; i < *ne; i++) {
    printf("%d %d %d %d %d %d\n",out.trianglelist[i*6+0],
           out.trianglelist[i*6+1],out.trianglelist[i*6+2],
           out.trianglelist[i*6+3],out.trianglelist[i*6+4],
           out.trianglelist[i*6+5]);
}
*/

// setup new lem, ibc and ibcngh array - faults will be added by createafault
/* beware basil and triangle order the midpoint nodes differently
                      1                1
                     /|               /|
             basil 5/ |4            3/ |5 triangle
                   /  |             /  |
                 2/___|0          2/___|0
                    3                4
*/
  mpk0=mpk=bvert + ibcmoffv;                    // in ibc, midpoints stored after all vertices
  k=0;
  for (i = 0; i < *ne; i++) {
    for (j = 0; j < 3; j++) {
      aj=out.trianglelist[i*6+j];                      // node no. of lem(j,i)
      ajac=out.trianglelist[i*6+(j+1)%3];              // anticlockwise vertex
      mac=out.trianglelist[i*6+3+(j+2)%3];             // anticlockwise midpoint
      macnv=mac+newvert;                               // offset for fault vertex nodes added
      lem[i*6+j]=aj + 1;                               // node j enters lem
      lem[i*6+3+(j+1)%3]=macnv + 1;                    // anticlockwise midpoint enters lem
      lemneighbor[i*3+j]=out.neighborlist[i*3+j] + 1;  // where is this thing used ?

 // elements with a side on the external boundary have neighbourlist = -1 in triangle
      lemnela=out.neighborlist[i*3+(j+2)%3] + 1;       // is there a neighbour on ac side ?

 //   if antclockwise segment (relative to j) is on the external boundary
      if(lemnela==0){                                  // if external boundary segment
	 ibc[k] = aj+1;
 	 ibcngh[k+nbpfinal]=ajac + 1;                  // record anticlockwise neighbour
//	 printf("trimesh: i, j, k = %i %i %i aj, ajac = %i %i\n",i,j,k,aj,ajac);
	 k++;

 //    ibcngh[k] is clockwise vertex; ibcngh[k+nbpfinal] is anticlockwise
         ibc[mpk] = macnv + 1;                         // anticlockwise midpoint
	 ibcngh[mpk+nbpfinal]=ajac + 1;
	 ibcngh[mpk]=aj + 1;
// 	 printf("ac: aj = %i ajac = %i macnv = %i k = %i mpk = %i\n",aj,ajac,macnv,k,mpk);
	 mpk++;
      }
    }
  }
//printf("%i vertex nodes and %i midpoint nodes for ext boundary stored in ibc, ibcmoffv = %i\n",k,mpk-mpk0,ibcmoffv);
  if((k!=bvert)||((mpk-mpk0)!=bvert)){   //check valid condition if any node duplicated
     printf("trimesh problem in ibc: %i vertex nodes, %i midpoints, bvert = %i\n",
            k,(mpk-mpk0),bvert);
  }

//   to identify clockwise vertex neighbours, examine all of the anticlockwise segments
  for(k1=0; k1<bvert; k1++){
     aj=ibc[k1];
     ajac=ibcngh[k1+nbpfinal];
     for(k=0; k<bvert; k++){
        if(ajac==ibc[k])ibcngh[k]=aj;
     }
  }

  /***** store segment data - boundaries defined in .poly file ******/
  iseg=0;
  if (out.numberofsegments>0 && out.segmentlist) {
    *nseg = out.numberofsegments;
//  printf("trimesh: *nseg = %i *nbp = %i newseg = %i\n",*nseg,*nbp,newseg);
    if ((err=AllocIntArray((*nseg+newseg) * 3, iseg_adr, 0))!=0){
      printf("Problem allocating ISEG array\n");
      return(err);
    }
    iseg = *iseg_adr;
    for (i = 0; i < *nseg; i++) {
        iseg[i*3]=out.segmentlist[i*2]+1;
        iseg[i*3+1]=out.segmentlist[i*2+1]+1;
        iseg[i*3+2]=out.segmentmarkerlist[i];
    //  printf("i = %i iseg[i*3+ 0,1,2] = %i %i %i\n",i,iseg[i*3],iseg[i*3+1],iseg[i*3+2]);
    }
  }

  if (*iflt){
    int mpstart = out.numberofpoints+newvert-1;
    int vstart = out.numberofpoints-out.numberofedges-1;
    int fstart = -1;
    vxtotal=0;
    while (ibc[vxtotal]!=-1) vxtotal++;
    vxtotal--;
    mptotal=nbpfinal-1;
    while (ibc[mptotal]==-1) mptotal--;
    vnofault=vxtotal;
    mnofault=mptotal;
//  printf("trimesh: vnofault = %i mnofault = %i bvert = %i\n",vnofault,mnofault,bvert);
    while (flts) {
        if (createafault(&out,&mpstart,&vstart,&fstart,
                lem,ex,ey,ne,ibc,ibcngh,ibctyp,ifbc1,ifbc2,iseg,
                nbpfinal,nbpext,nseg,&vxtotal,&mptotal,flts)) {
                printf("error creating fault\n");
                return(1);
        }
        flts=flts->next;
    }
/*  printf("ibc  array follows, nbp = %i\n",nbpfinal);
    for (k=0; k < nbpfinal; k++){printf("k = %i ibc=%i ",k,ibc[k]);} */
    /***** update the size parameters to include new pts and     ******/
    /***** recount the max number of nodes connected to a vertex ******/
    *nbp=nbpfinal;
    *nnnof = *nn;
    *nn += newvert;
    *nmpnof = *nmp;
    *nmp += newseg;
    *nup += newbpts;
    ibcadded=vxtotal+mptotal-vnofault-mnofault;
//  printf("trimesh: vxtotal = %i mptotal = %i vnofault = %i mnofault = %i\n",vxtotal,mptotal,vnofault,mnofault);
    printf("Check: space used in IBC for faults %i should be newibc = %i\n",ibcadded,newibc);
    printf("Check: space used in EX and EY for fault vertices %i reaches NN = %i\n",vstart+1,*nn);
    if((ibcadded != newibc) || ((vstart+1) != *nn)){
       printf("********  Warning, check out above mismatch  *********\n");
    }
  }   // non-zero iflt

// for debugging, check ibc, ibcngh ok
//for (k=0; k< nbpfinal; k++){
//   printf("trimesh, k = %i ibc[k] = %i nac = %i ncw = %i ex = %f ey = %f \n",k,ibc[k],
//	    ibcngh[k],ibcngh[k+nbpfinal],ex[ibc[k]-1],ey[ibc[k]-1]);
// }
  
  *maxnbs=find_maxnbs(ex,ey,lem,*ne,*nn);
  printf("After mesh constructed: maxnbs=%d used to define matrix bandwidth\n",*maxnbs);

//fflush(stdout);

  /************************************/
  /******* Output for debugging *******/
  /************************************/

  /*printf("Basil Variables:\n\n");*/
  /*showbasilvar(*nup,*ne,*nbp,ex,ey,ibc,lem,lemneighbor);*/
//  printf("in trimesh, checking for zeros in ibcngh: nbpfinal = %i\n",nbpfinal);
    for (j = 0; j<nbpfinal; j++){
      if ((ibcngh[j]==0)||(ibcngh[j+nbpfinal]==0)){
         printf("trimesh WARNING: j = %i ibc[j] = %i ibcngh[j] = %i ibcngh[j+nbpfinal] = %i\n",
                     j,ibc[j],ibcngh[j],ibcngh[j+nbpfinal]);
      }
    }


  /***********************/
  /******* Cleanup *******/
  /***********************/
  /* Free all allocated arrays, including those allocated by Triangle. */
//if (ncnt) free(ncnt);
  if (numnbs) free(numnbs);
  if (rgnattribs) free(rgnattribs);
  free_flt_list(flts);
  /* triangle pts out.regionlist at inp.regionlist */
  if (inp.regionlist!=NULL && out.regionlist==inp.regionlist) clrptr=1;
  cleanio( &inp );
  if (clrptr) out.regionlist=0;
  cleanio( &out );

  return 0;

}    // end routine trimesh

//*******************************************************************

int tri_remesh( int *iflt,int *nup,int *nn,int *nmp,int *nbp,int *nnnof,
             int *nmpnof,int *ne,int *nfp,int *nfpf3,int *nseg,int *ivv,int *maxnbs,
             int *nellep,int **ielle_adr,int **ipolyn_adr,
             int **lem_adr,int **lemneighbor_adr,
             int **ibc_adr,int **ibcnb_adr,int **ibctyp_adr,int **ifbc1_adr,int **ifbc2_adr,
             int **ifeqv_adr,int **jfbc1_adr,int **jfbc2_adr,int **iseg_adr,
             int *isegnew, int newnseg, int *ibp, int newpts, int num_pt_attr,
             int **imat_adr, int mreg,
             float **ex_adr,float **ey_adr,float **vhb_adr,float **qbnd_adr,
             float *exbnd, float *eybnd, float *pt_attributes,
             float *xlen, float *ylen,
             float area, float quality, float *visdflt,float *se,int pbxy,int *nbpext,
             int npbpt)
{
  struct triangulateio out, inp;
  struct flagvals triflags;
  struct gen_fault *flts=NULL;
  int i,j,k,k1;
  int *lem,*lemneighbor,*ibc,*ibcngh,*ibctyp,*imat;
  int *ifbc1,*ifbc2,*ifeqv,*jfbc1,*jfbc2;
  float *ex,*ey,*qbnd;
  extern int pntonsegment_();
  int aj,ajac,mac,macnv,lemnela,mpk;
  int newvert, newseg, newibc, ibcmoffv, newbpts, nbpfinal;
  int bvert,bmidp,nvptb,nvptf,nvpti;
  float *rgnattribs=NULL;
/*float *unodedata=NULL,*unodeattribs=NULL; */
  float c_x, c_y;
/*int maxrgns=0,maxunodes, dist_set=0; */
  int clrptr=0,icount=0,imatmax=0,vxtotal=0,mptotal=0;
  int npat,numrgns=0,vnofault,mnofault,ibcadded;
  int itp1,itp2;
  unsigned char ptflag;

  /* default area is  0.025*0.025 (area of tri with minNodeSep sides) */
  triflags.area = area; /* area of tri with minNodeSep sides */
  /*triflags.area = 0;*/
  triflags.quality = quality;
  /*triflags.quality = 5;*/
  triflags.midptnodes = 1;
//fprintf(stdout,"in tri_remesh, area = %f quality = %f\n", area, quality);
//
//  triflags.bndpts=0 suppresses all changes to boundary nodes in remesh operation
  if(pbxy) {
     triflags.bndpts = 0;
     if(pbxy==1)
     printf("Periodic boundaries (X) are enabled. npbpt = %i\n",npbpt);
     if(pbxy==2)
     printf("Periodic boundaries (Y) are enabled. npbpt = %i\n",npbpt);
  }
  else triflags.bndpts = 1;

  initio( &inp );
  initio( &out );

/*    new section to identify material numbers   */

if  (mreg>0) {
   imat = *imat_adr;
   for (i=0; i<*ne; i++) {
       if (imat[i]>numrgns) numrgns=imat[i];
   }
   if ((rgnattribs= (float *)malloc((numrgns) * 4 *
                                    sizeof(float)))==NULL) return(1);
   for (i=0; i<numrgns; i++) {
       unsigned char found=0;
       for (j=0; j<*ne && !found; j++) {
         if (imat[j]==(i+1)) {
              for (k = 0,c_x=0.0,c_y=0.0; k < 3; k++) {
                            c_x += (*ex_adr)[(*lem_adr)[j*6+k]-1];
                            c_y += (*ey_adr)[(*lem_adr)[j*6+k]-1];
              }
              c_x /= 3.0;
              c_y /= 3.0;
              found = 1;
         }
       }
    if (found) {
      rgnattribs[i*4] = c_x;
      rgnattribs[i*4 + 1] = c_y;
      rgnattribs[i*4 + 2] = i+1;
      rgnattribs[i*4 + 3] = area;
    }
    /* report error */
   }
}
//printf("before tripolydata, newpts = %i newnseg = %i\n",newpts,newnseg);
  if (tripolydata(&out,&inp,&triflags,newpts,newnseg,num_pt_attr,
                  isegnew,exbnd,eybnd,pt_attributes,rgnattribs,numrgns)!=0)
      return(1);
//printf("after tripolydata: newpts = %i numrgns = %i\n",newpts,numrgns);
//  WritePoly(&inp, "regrid_7_5in.poly");

    if (numrgns>0) {free(rgnattribs); rgnattribs=NULL;}

    if (*lem_adr!=NULL) {free(*lem_adr); *lem_adr=NULL;}
    if (*lemneighbor_adr!=NULL) {free(*lemneighbor_adr); *lemneighbor_adr=NULL;}
    if (*ibc_adr!=NULL) {free(*ibc_adr);*ibc_adr=NULL;}
    if (*ibcnb_adr!=NULL) {free(*ibcnb_adr); *ibcnb_adr=NULL;}
    if (*ibctyp_adr!=NULL) {free(*ibctyp_adr); *ibctyp_adr=NULL;}
//  if (*ibpole_adr!=NULL) {free(*ibpole_adr); *ibpole_adr=NULL;}
    if (*iseg_adr!=NULL) {free(*iseg_adr); *iseg_adr=NULL;}
    if (*imat_adr!=NULL) {free(*imat_adr); *imat_adr=NULL;}
    if (*qbnd_adr!=NULL) {free(*qbnd_adr); *qbnd_adr=NULL;}
//  if (*taperf_adr!=NULL) {free(*taperf_adr); *taperf_adr=NULL;}
  /**********************************/
  /******* Setup Basil Arrays *******/
  /**********************************/

//  identify and count fault points in seg array, set value of nfpf3, nfp
  newvert=newseg=newbpts=newibc=ibcmoffv=0;
  *nfpf3=*nfp=0;
  if(*iflt && (out.numberofsegments>0)){
    findfltpts(&flts,&out,nfpf3,&newvert,&newseg,&newibc,&ibcmoffv);
    if(!(*nfpf3)){
        printf("Warning: fault type is non-zero but no internal faults found; iflt reset to 0\n");
	*iflt=0;
    }
  }
  *nfp=*nfpf3+npbpt;
//printf("after findfltpts: newvert, newseg, newibc, *nfp, *nfpf3 = %i %i %i %i %i\n", 
//         newvert,newseg,newibc,*nfp,*nfpf3);

// new mesh size parameters after tripolydata
  *nup=out.numberofpoints;
  *nmpnof=*nmp=out.numberofedges;
  *nnnof=*nn=*nup-(*nmp);
  *ne=out.numberoftriangles;
//printf("tri_remesh: *nup = %i *nn = %i *nmp = %i *ne = %i\n",*nup,*nn,*nmp,*ne);

//  re-establish the seg array
  *nseg = out.numberofsegments;
  nvpti=nvptf=nvptb=0;
  if(*nseg){
    *iseg_adr = isegnew = (int *) malloc((*nseg+newseg) * 3 * sizeof(int));
    if(isegnew==NULL) return(1);
    for (k = 0; k < *nseg; k++) {
        isegnew[k*3]=out.segmentlist[k*2]+1;
        isegnew[k*3+1]=out.segmentlist[k*2+1]+1;
        isegnew[k*3+2]=out.segmentmarkerlist[k];
//	printf("k = %i isegnew = %i %i %i \n",k,isegnew[k*3],isegnew[k*3+1],isegnew[k*3+2]);

//  determine number of vertex points on external boundary and on internal faults
    }
    for(i=0; i<*nup; i++){
      ptflag=is_inex_bnd_pt(i,out.segmentlist,out.segmentmarkerlist,*nseg);
      if(ptflag==1){nvptb++;}
      else if(ptflag==2){nvptf++;}
      else if(ptflag==3){nvpti++;}
    }
  }
  *nbpext=*nbp=nvptb;   //  the no. of external boundary vertex points = no. of midpoints
//printf("tri_remesh: nvptb = %i nvptf = %i nvpti = %i *nseg = %i\n",nvptb,nvptf,nvpti,*nseg);
//
//  Warning: check what happens if fault intersects boundary or other fault, or if
//  multiple internal faults. The following check provides a warning that has occurred.
//
  if((nvptb+nvptf+nvpti)!= *nseg){
      printf("Remeshed segment boundaries appear to intersect at %i points\n",
		      (*nseg-(nvptb+nvptf+nvpti)));
  }

// allocate space for ex, ey coordinates, allowing for duplicate fault nodes to be added
// newseg is number of new midpoints; newvert is number of new vertices
  newbpts=newseg+newvert;
  *ex_adr = ex = (float *) malloc((*nup+newbpts) * sizeof(float));
  *ey_adr = ey = (float *) malloc((*nup+newbpts) * sizeof(float));
  for (i = 0; i < *nup+newbpts; i++) {
    ex[i] = 0.0;
    ey[i] = 0.0;
  }
  for (i = 0; i < *nup; i++) {
    j=i;
    if (j>=*nn) j+=newvert;   	// allows new vertex nodes to be inserted ahead of midpoints
    ex[j] = (float) out.pointlist[i*2];
    ey[j] = (float) out.pointlist[i*2+1];
  }
//  printf("markerlist %d\n",*nbp);
//
//  space for ibc = boundary vertices + boundary midpoints + 2*(fault vertices + fault midpts)
//  following should work for closed faults that don't intersect.  Other cases may need further thought
//
  bvert=nvptb + 2*newvert;
  bmidp=nvptb + 2*newseg;
  nbpfinal = bvert+bmidp;
//printf("tri_remesh: bvert = %i bmidp = %i nbpfinal = %i\n",bvert,bmidp,nbpfinal);
  *ibc_adr = ibc = (int *) malloc(nbpfinal * sizeof(int));
  *ibcnb_adr = ibcngh = (int *) malloc(nbpfinal * 2 * sizeof(int));
  *ibctyp_adr = ibctyp = (int *) malloc(nbpfinal * 2 * sizeof(int));
  *qbnd_adr = qbnd = (float *) malloc(nbpfinal * 2 * sizeof(float));
  *ifbc1_adr = ifbc1 = (int *) malloc(*nfp * sizeof(int));
  *ifbc2_adr = ifbc2 = (int *) malloc(*nfp * sizeof(int));
  *ifeqv_adr = ifeqv = (int *) malloc(*nfp * sizeof(int));
  *jfbc1_adr = jfbc1 = (int *) malloc(*nfp * sizeof(int));
  *jfbc2_adr = jfbc2 = (int *) malloc(*nfp * sizeof(int));
  for (i=0; i < nbpfinal; i++) ibc[i]=-1;
  for (i=0; i < nbpfinal * 2; i++) ibcngh[i]=0;
  for (i=0; i < nbpfinal * 2; i++) ibctyp[i]=-1;
  for (i=0; i < nbpfinal * 2; i++) qbnd[i]=0.0;
  for (i=0; i < (*nfp); i++) ifbc1[i]=0;
  for (i=0; i < (*nfp); i++) ifbc2[i]=0;
  for (i=0; i < (*nfp); i++) jfbc1[i]=0;
  for (i=0; i < (*nfp); i++) jfbc2[i]=0;

//  re-establish material numbers from the region identifier
  if (mreg>0) {
     *imat_adr = imat= (int *) malloc(*ne * sizeof(int));
     icount=0;
     imatmax=0;
     if (imat!=NULL && out.triangleattributelist!=NULL) {
        for (i = 0; i < *ne; i++) {
           imat[i] = (int)out.triangleattributelist[i];
           if(imat[i]<=0) {imat[i]=1; icount++;}
           if(imat[i]>imatmax) imatmax=imat[i];
        }
      printf ("Prior to interpolation, %i unset material numbers set to 1, maximum material number = %i\n",
                      icount,imatmax);
     }
  }

  /*
   * ivv    visc (vhb[i,1-7])    se (vhb[i,8])
   *  0       constant(1)        constant(1)    No vhb array
   *  1       variable           constant(1)
   *  2       constant           constant(>1)     No vhb array if visc=1
   *  3       variable           constant(>1)
   *  4       variable           variable
   */
/*  if ( *ivv>0 && *vhb_adr==NULL ) {
 
        // allocate and initialise the vhb array
        // the input file contains VIS settings
  
      if ((*vhb_adr = vhb=
          (float *) malloc(*ne * 8 * sizeof(float)))==NULL)
            return(1);
      for (i = 0; i < *ne * 8; i++) vhb[i]=1.0;
  }
  printf("allocated and initialized vhb, *ne = %i\n");

  // the triangle attribute is the region identifier
  // save it and set triangle attributes, if found
 
  if (out.triangleattributelist!=NULL) {
      if ((*ipolyn_adr = polyn=
                  (int *) malloc(*ne * sizeof(int)))==NULL)
          return(1);
      maxrgns = (int)out.triangleattributelist[0];
      for (i = 0; i < *ne; i++) {
          k = (int)out.triangleattributelist[i];
          polyn[i] = k;
          if (k > maxrgns) maxrgns = k;
      }
      if ((rgnattribs= (float *)malloc((maxrgns+1) *
                                    sizeof(float)))==NULL) return(1);
      dflt = 1.0;
      for (i=0;i<maxrgns;i++) rgnattribs[i] = dflt;
      if (attribute_set(infile,VISCOSITY_KEY,
                        maxrgns+1,rgnattribs,&dflt)) {
          if ( *vhb_adr==NULL ) {
            
            // allocate and initialise the vhb array
         
              if ((*vhb_adr = vhb=
                  (float *) malloc(*ne * 8 * sizeof(float)))==NULL)
                    return(1);
              for (i = 0; i < *ne * 8; i++) vhb[i]=dflt;
              for (i = 0; i < *ne; i++) vhb[i*8+7]=1.0;
            }
          if (*ivv<1)  *ivv=1;
          if (dflt>1.0 && *ivv<2)  *ivv=2;
          for (i = 0; i < *ne; i++) {
              k = (int)out.triangleattributelist[i];
              polyn[i] = k;
              for (j = 0; j < 7; j++) 
                  vhb[i*8+j] = rgnattribs[k];
          }
      }
      dflt = 1.0;
      for (i=0;i<maxrgns;i++) rgnattribs[i] = dflt;
      if (attribute_set(infile,SE_KEY,
                        maxrgns+1,rgnattribs,&dflt)) {
          if ( *vhb_adr==NULL ) {
          
            // allocate and initialise the vhb array
         
              if ((*vhb_adr = vhb=
                  (float *) malloc(*ne * 8 * sizeof(float)))==NULL)
                    return(1);
              for (i = 0; i < *ne * 8; i++) vhb[i]=1.0;
              for (i = 0; i < *ne; i++) vhb[i*8+7]=dflt;
            }
          if (fabs((double)dflt-1.0)>eps) {
            if ( *ivv<1)  *ivv=2;
            else if ( *ivv==1)  *ivv=3;
          }
          if ((fabs((double)dflt-1.0)>eps) && *ivv<2)  *ivv=2;
          curr = dflt;
          for (i = 0; i < *ne; i++) {
              k = (int)out.triangleattributelist[i];
              polyn[i] = k;
              vhb[i*8+7] = rgnattribs[k];
              if (fabs((double)rgnattribs[k]-curr)>eps) {
                curr =  rgnattribs[k];
                if (*ivv<4)  *ivv=4;
              }
          }
      }
      if (unode_data(infile,UNODES_KEY,
                        &maxunodes,&unodedata)) {
        dflt=1.0;
          if (unode_attribute(infile,U_VISCOSITY_KEY,
                        maxunodes,&unodeattribs,&dflt)) {
          if ( *vhb_adr==NULL ) {
    
            // allocate and initialise the vhb array
   
              if ((*vhb_adr = vhb=
                  (float *) malloc(*ne * 8 * sizeof(float)))==NULL)
                    return(1);
              for (i = 0; i < *ne * 8; i++) vhb[i]=dflt;
              for (i = 0; i < *ne; i++) vhb[i*8+7]=1.0;
            }
          if (*ivv<1)  *ivv=1;
          if (dflt>1.0 && *ivv<2)  *ivv=2;
          for (i = 0; i < *ne; i++) {
              k = (int)out.triangleattributelist[i];
              polyn[i] = k;
              for (j = 0,c_x=0.0,c_y=0.0; j < 3; j++) {
                a = out.trianglelist[i*6+j];
                c_x += out.pointlist[a*2];
                c_y += out.pointlist[a*2+1];
              }
              c_x /= 3.0;
              c_y /= 3.0;
              min_dist=0.0; b=-1;
              for (j = 0, dist_set=0,val=dflt; j < maxunodes; j++) {
                if (unodedata[j+maxunodes*2]==k) {
                    dist = (c_x-unodedata[j])*(c_x-unodedata[j]) +
                    (c_y-unodedata[j+maxunodes])*(c_y-unodedata[j+maxunodes]);
                    if (!dist_set || dist<min_dist) {
                        min_dist = dist;
                        val = unodeattribs[j];
                    }
                }
              }
              for (j = 0; j<7; j++) vhb[i*8+j]=val;
          }
          }
      }
  }

 */

//  re-establish lem and lemneighbor arrays
  *lem_adr = lem= (int *) malloc(*ne * 6 * sizeof(int));
  *lemneighbor_adr = lemneighbor= (int *) malloc(*ne * 3 * sizeof(int));

// setup new lem, ibc and ibcngh array - faults will be added by createafault
/* beware basil and triangle order the midpoint nodes differently
                      1                1
                     /|               /|
             basil 5/ |4            3/ |5 triangle
                   /  |             /  |
                 2/___|0          2/___|0
                    3                4
*/
  npat=out.numberofpointattributes;
  mpk=*nbpext+2*newvert;                    // in ibc, midpoints stored after all vertices
  k=0;
  for (i = 0; i < *ne; i++) {
    for (j = 0; j < 3; j++) {
      aj=out.trianglelist[i*6+j];                      // node no. of lem(j,i)
      ajac=out.trianglelist[i*6+(j+1)%3];              // anticlockwise vertex
      mac=out.trianglelist[i*6+3+(j+2)%3];             // anticlockwise midpoint
      macnv=mac+newvert;                               // offset for fault vertex nodes added
      lem[i*6+j]=aj + 1;                               // node j enters lem
      lem[i*6+3+(j+1)%3]=macnv + 1;                    // anticlockwise midpoint enters lem
      lemneighbor[i*3+j]=out.neighborlist[i*3+j] + 1;  // where is this thing used ?

 // elements with a side on the external boundary have neighbourlist = -1 in triangle
      lemnela=out.neighborlist[i*3+(j+2)%3] + 1;       // is there a neighbour on ac side ?

 //   if antclockwise segment (relative to j) is on the external boundary
      if(lemnela==0){                                  // if external boundary segment
	 ibc[k] = aj+1;
         qbnd[k]=out.pointattributelist[aj*npat];      // boundary entries for node j
         qbnd[k+nbpfinal]=out.pointattributelist[aj*npat+1];
         ibctyp[k]=(int)ceilf((float)(out.pointattributelist[aj*npat+2]));
         ibctyp[k+nbpfinal]=(int)ceilf((float)(out.pointattributelist[aj*npat+3]));
 	 ibcngh[k+nbpfinal]=ajac + 1;                  // record anticlockwise neighbour
	 k++;

 //    ibcngh[k] is clockwise vertex; ibcngh[k+nbpfinal] is anticlockwise
         ibc[mpk] = macnv + 1;                         // anticlockwise midpoint
         qbnd[mpk]=out.pointattributelist[mac*npat];
         qbnd[mpk+nbpfinal]=out.pointattributelist[mac*npat+1];
         ibctyp[mpk]=(int)ceilf((float)(out.pointattributelist[mac*npat+2]));
         ibctyp[mpk+nbpfinal]=(int)ceilf((float)(out.pointattributelist[mac*npat+3]));
	 ibcngh[mpk+nbpfinal]=ajac + 1;
	 ibcngh[mpk]=aj + 1;
// 	 printf("ac: aj = %i ajac = %i macnv = %i k = %i mpk = %i\n",aj,ajac,macnv,k,mpk);
//       printf("ex[macnv], ey[macnv] = %f %f ex[aj] ey[aj] = %f %f ex[ajac] ey[ajac] = %f %f\n",
//	         ex[macnv], ey[macnv], ex[aj], ey[aj], ex[ajac], ey[ajac]);
	 mpk++;
      }
    }
  }
  if((k!=nvptb)||((mpk-*nbpext-2*newvert)!=nvptb)){
     printf("tri_remesh problem in ibc: %i vertex nodes added, %i midpoints added. nvptb = %i\n",
            k,(mpk-*nbpext-2*newvert),nvptb);
  }

//   to identify clockwise vertex neighbours, examine all of the anticlockwise segments
  for(k1=0; k1<nvptb; k1++){
     aj=ibc[k1];
     ajac=ibcngh[k1+nbpfinal];
     for(k=0; k<nvptb; k++){
        if(ajac==ibc[k])ibcngh[k]=aj;
     }
  }
//  for (k=0; k< nbpfinal; k++){
//   printf("tri_remesh, k = %i ibc[k] = %i ex = %f ey = %f typx = %i qx = %f typy = %i qy = %f\n",k,ibc[k],
//	    ex[ibc[k]-1],ey[ibc[k]-1],ibctyp[k],qbnd[k],ibctyp[k+nbpfinal],qbnd[k+nbpfinal]);
//   }

//  re-establish the faults
  if (*iflt && *nseg) {
    int mpstart = out.numberofpoints+newvert-1;
    int vstart = out.numberofpoints-out.numberofedges-1;
    int fstart = -1;
    vxtotal=0;
    while (ibc[vxtotal]!=-1) vxtotal++;
    vxtotal--;
    mptotal=nbpfinal-1;
    while (ibc[mptotal]==-1) mptotal--;
    vnofault=vxtotal;
    mnofault=mptotal;
    while (flts) {
       if (createafault(&out,&mpstart,&vstart,&fstart,
            lem,ex,ey,ne,ibc,ibcngh,ibctyp,ifbc1,ifbc2,
            isegnew,nbpfinal,nbpext,nseg,&vxtotal,&mptotal,flts)) {
                printf("error creating fault\n");
                return(1);
        }
        flts=flts->next;
    }

// reset fault boundary conditions from previous stored attributes (defined real)
// external boundary conditions have been reset earlier before faults added
// every second point was re-initiated by createafault; no point attribute exists

    if (npat>=5) {
 // printf("vnofault = %i, newvert = %i nbpfinal = %i\n",vnofault,newvert,nbpfinal);
      for (i=vnofault+1;i<vnofault+1+2*newvert;i=i+2) {       // for fault vertices
        j=ibc[i]-1;
        qbnd[i+1]=qbnd[i]=out.pointattributelist[j*npat];
        qbnd[i+1+nbpfinal]=qbnd[i+nbpfinal]=out.pointattributelist[j*npat+1];
        itp1=(int)ceilf((float)(out.pointattributelist[j*npat+2]));
        itp2=(int)ceilf((float)(out.pointattributelist[j*npat+3]));
        ibctyp[i+1]=ibctyp[i]=itp1;
        ibctyp[i+1+nbpfinal]=ibctyp[i+nbpfinal]=itp2;
 /*     printf("tri_remesh: i = %i ibc[i] = %i ibctyp[i], [i+nbpfinal] = %i %i qbnd = %f %f\n",
        	        i,ibc[i],ibctyp[i],ibctyp[i+nbpfinal],qbnd[i],qbnd[i+nbpfinal]);
        printf("tri_remesh: i = %i ibc[i] = %i ibctyp[i], [i+nbpfinal] = %i %i qbnd = %f %f\n",
        	        i+1,ibc[i+1],ibctyp[i+1],ibctyp[i+1+nbpfinal],qbnd[i+1],qbnd[i+1+nbpfinal]); */
      }
 // printf("mnofault = %i, nbpfinal = %i newvert = %i\n",mnofault,nbpfinal,newvert);
      for (i=mnofault+1;i<nbpfinal;i=i+2) {        // for fault midpoints
        j=ibc[i]-newvert-1;
        qbnd[i+1]=qbnd[i]=out.pointattributelist[j*npat];
        qbnd[i+1+nbpfinal]=qbnd[i+nbpfinal]=out.pointattributelist[j*npat+1];
        itp1=(int)ceilf((float)(out.pointattributelist[j*npat+2]));
        itp2=(int)ceilf((float)(out.pointattributelist[j*npat+3]));
        ibctyp[i+1]=ibctyp[i]=itp1;
        ibctyp[i+1+nbpfinal]=ibctyp[i+nbpfinal]=itp2;
 /*     printf("tri_remesh: i = %i ibc[i] = %i ibctyp[i], [i+nbpfinal] = %i %i qbnd = %f %f\n",
        	        i,ibc[i],ibctyp[i],ibctyp[i+nbpfinal],qbnd[i],qbnd[i+nbpfinal]);
        printf("tri_remesh: i = %i ibc[i] = %i ibctyp[i], [i+nbpfinal] = %i %i qbnd = %f %f\n",
        	        i+1,ibc[i+1],ibctyp[i+1],ibctyp[i+1+nbpfinal],qbnd[i+1],qbnd[i+1+nbpfinal]); */
      }
    }
 //  update the size parameters to include new pts 
    *nn += newvert;
    *nmp += newseg;
    *nup += newbpts;
    ibcadded=vxtotal+mptotal-vnofault-mnofault;
    printf("Check: space used in IBC for faults %i should be newibc = %i\n",ibcadded,newibc);
    printf("Check: space used in EX and EY for fault vertices %i reaches NN = %i\n",(vstart+1),*nn);
    if((ibcadded != newibc) || ((vstart+1) != *nn)){
       printf("********  Warning, check out above mismatch  *********\n");
    }
  }   // non-zero iflt
  *nbp=nbpfinal;

//  printf("tri_remesh: k = %i ibc[k] = %i ibcngh[k], ibcngh[k+nbpfinal] = %i %i\n",
//      	        k,ibc[k],ibcngh[k],ibcngh[k+nbpfinal]); }
//int nbp2=nbpfinal*2;
//printf("after faults added with nbp2 =%i; ibctyp follows\n",nbp2);
//imatpp_(ibctyp,&nbp2);
//printf("after faults added with nbp2 =%i; qgnd follows\n",nbp2);
//smatpp_(qbnd,&nbp2);

  /***** recount the max number of nodes connected to a vertex ******/
  *maxnbs=find_maxnbs(ex,ey,lem,*ne,*nn);
  printf("After mesh constructed: maxnbs=%d used to define matrix bandwidth\n",*maxnbs);
  if(*maxnbs > 20)return 1;

  //printf("Basil Variables:\n\n");               // for debugging
  //showbasilvar(*nup,*ne,*nbp,ex,ey,ibc,lem,lemneighbor);
  //fflush(stdout);

  // Free all allocated arrays, including those allocated by Triangle.
  if (rgnattribs) free(rgnattribs);
  free_flt_list(flts);
  /* triangle pts out.regionlist at inp.regionlist */
  if (inp.regionlist!=NULL && out.regionlist==inp.regionlist) clrptr=1;
  cleanio( &inp );
  if (clrptr) out.regionlist=0;
  cleanio( &out );

  return 0;

}   //end routine tri_remesh

//*******************************************************************************

int attribute_set(char *infile,char *key,int max,
                  float *rgnattribs, float *dflt)
{
    FILE *fp=NULL;
    char buf[81], str[32];
    int found = 0, finished=0, i;
    int num=0;
    float val;

    if((fp=fopen(infile,"r"))==NULL) {
        num=fprintf(stderr,"cannot open file %s\n",infile);
        return(1);
    }
    /*
     * find keyword
     */
    while (!found && (fgets(buf,80,fp)!=NULL)) {
        if (isalpha(buf[0])) {
            num=sscanf(buf,"%s",str);
            if (!strcmp(str,key)) found = 1;
        }
    }
    if (found) {
    /*
     * set default value
     */
        if (fgets(buf,80,fp)!=NULL) {
            if (isalpha(buf[0])) {
                num=sscanf(buf,"%s",str);
                if (!strcmp(str,DEFAULT)) {
                    if ((num=sscanf(&buf[strlen(DEFAULT)],"%f",dflt))==1) {
                        for (i=0;i<max;i++) rgnattribs[i] = *dflt;
                    }
                    finished = (fgets(buf,80,fp)==NULL);
                }
                else finished = 1;
            }
        }
        while (!finished && ((isdigit(buf[0])||(buf[0]=='#')))) {
            if (isdigit(buf[0])) {
                num = sscanf(buf,"%d %f",&i,&val);
                if (num != 2) fprintf(stderr,"Error scanning %s\n",buf);
                else {
                    if (i<max) rgnattribs[i] = val;
                    else fprintf(stderr,"Invalid region id %d\n",i);
                }
                finished = (fgets(buf,80,fp)==NULL);
            }
            else if (buf[0]=='#') {
                dump_comments(fp);
                finished = (fgets(buf,80,fp)==NULL);
            }
            else finished = 1;
        }
    }
    if (fp) fclose(fp);
    return(found);
}

int unode_data(char *infile,char *key,
                        int *maxunodes,float **unodedata)
{
    FILE *fp=NULL;
    char buf[81], str[32];
    int found = 0, finished=0, i;
    int num=0;
    float *ptr_x=NULL, *ptr_y=NULL, *ptr_rgn=NULL;

    if((fp=fopen(infile,"r"))==NULL) {
        num=fprintf(stderr,"cannot open file %s\n",infile);
        return(0);
    }
    /*
     * find keyword
     */
    while (!found && (fgets(buf,80,fp)!=NULL)) {
        if (isalpha(buf[0])) {
            num=sscanf(buf,"%s",str);
            if (!strcmp(str,key)) found = 1;
        }
    }
    if (found && (fgets(buf,80,fp)!=NULL)) {
    /*
     * allocate arrays 
     */
        num=sscanf(buf,"%d",maxunodes);
        if ((num==1) && (*maxunodes>0) && (*unodedata=
                (float *) malloc(*maxunodes * 3 * sizeof(float)))!=NULL) {
            for (i = 0; i < *maxunodes*3; i++) (*unodedata)[i]=0.0;
            ptr_x=*unodedata;
            ptr_y=ptr_x + *maxunodes;
            ptr_rgn=ptr_y + *maxunodes;
        }
        else finished = 1;
    
        finished = (fgets(buf,80,fp)==NULL);
        while (!finished && (isdigit(buf[0]))) {
            num=sscanf(buf,"%d %f %f %f",&i,ptr_x,ptr_y,ptr_rgn);
            ptr_x++; ptr_y++; ptr_rgn++;
        }
    }
    if (fp) fclose(fp);
    return(found);
}

int unode_attribute(char *infile,char * key,
                        int maxunodes,float **unodeattribs,float *dflt)
{
    FILE *fp=NULL;
    char buf[81], str[32];
    int found = 0, finished=0, i;
    int num=0;
    float val;

    if((fp=fopen(infile,"r"))==NULL) {
        num=fprintf(stderr,"cannot open file %s\n",infile);
        return(0);
    }
    /*
     * find keyword
     */
    while (!found && !feof(fp)) {
        if ((fgets(buf,80,fp)!=NULL) && isalpha(buf[0])) {
            num=sscanf(buf,"%s",str);
            if ((num==1) && (!strcmp(str,key))) found = 1;
        }
    }
    if (found) {
    /*
     * allocate arrays and
     * set default value
     */
        if (maxunodes>0 &&
            (*unodeattribs=
                        (float *) malloc(maxunodes * sizeof(float)))==NULL) {
            for (i = 0; i < maxunodes ; i++) (*unodeattribs)[i]=*dflt;
        }
        else finished = 1;
        finished = (fgets(buf,80,fp)==NULL);
        if (!finished && (isalpha(buf[0]))) {
            num=sscanf(buf,"%s",str);
            if (num==1 && !strcmp(str,DEFAULT)) {
                if ((num=sscanf(&buf[strlen(DEFAULT)],"%f",dflt))==1) {
                    for (i=0;i<maxunodes;i++) (*unodeattribs)[i] = *dflt;
                }
                finished = (fgets(buf,80,fp)==NULL);
            }
            else finished = 1;
        }
    
        while (!finished && (isdigit(buf[0])) && !feof(fp)) {
            num=sscanf(buf,"%d %f",&i,&val);
            if ((num==2) && (i<maxunodes)) (*unodeattribs)[i] = val;
            else fprintf(stderr,"Invalid unode id %d\n",i);
            finished = (fgets(buf,80,fp)==NULL);
        }
    }
    if (fp) fclose(fp);
    return(found);
}

int vcfrompoly_( float *vhb,float *ex,float *ey,int *lem,int * nor,
             int *ivv,int *nup,int *mesh_ne,
             char *polyfile,char *comments,
             int *ierr)
{
  struct triangulateio out,inp;
  struct flagvals triflags;
  char optionaltxt[80];
  int i,k,n,err=0;
  int maxrgns=0,nbm=0,internal_bnd=0,found;
  int dum=0,*dum_adr=0;
  float val,dflt;
  float *rgnattribs=NULL;
  double pt[2], cpt[2];
  polydata *polys;

  initio( &inp );

  if((err = readpolyinput(&inp,polyfile,&triflags,optionaltxt,
                      &dum,&dum_adr,&nbm,&internal_bnd))!=0)
      return(err);
  if (dum_adr) free(dum_adr);
  if ((polys = (polydata *)
                malloc(sizeof(polydata)*(inp.numberofregions+1)))
      ==0) return(1);
//for (i = 0; i < inp.numberofregions; i++)
  for (i=0;i<inp.numberofregions+1;i++) { polys[i].numpts=0; polys[i].ppts=0; }
  if((err = getpolydata(&inp,polys,inp.numberofregions))!=0)
      return(err);
  /*
   * ivv    visc (vhb[i,1-7])    se (vhb[i,8])
   *  0       constant(1)        constant(1)    No vhb array
   *  1       variable           constant(1)
   *  2       constant           constant(>1)     No vhb array if visc=1
   *  3       variable           constant(>1)
   *  4       variable           variable
   */
      if ((rgnattribs= (float *)malloc((maxrgns+1) *
                                    sizeof(float)))==NULL) return(1);
      dflt = 1.0;
      for (i=0;i<maxrgns;i++) rgnattribs[i] = dflt;
      if (attribute_set(polyfile,VISCOSITY_KEY,
                        maxrgns+1,rgnattribs,&dflt)) {
          for (i = 0; i < *mesh_ne * 8; i++) vhb[i]=dflt;
          for (i = 0; i < *mesh_ne; i++) vhb[i*8+7]=1.0;
          if (*ivv<1)  *ivv=1;
          if (dflt>1.0 && *ivv<2)  *ivv=2;
          for (k = 0; k < *mesh_ne; k++) {
              cpt[0] = cpt[1] = 0.0;
              for (n = 0; n < 7; n++) {
                if (n==6) {
                    pt[0] = cpt[0]/3.0;
                    pt[1] = cpt[1]/3.0;
                }
                else {
                    pt[0] = ex[nor[lem[k*6+n]-1]-1];
                    pt[1] = ey[nor[lem[k*6+n]-1]-1];
                    if (n<3) {
                        cpt[0] += pt[0];
                        cpt[1] += pt[1];
                    }
                }
                /*
                 * for each mesh vertex, find the vc value
                 * from the poly data
                 */
                found = 0;
                val = dflt;
                  for (i = 0; i < inp.numberofregions && !found; i++) {
                    if (polys[i].numpts>0 &&
                         (found=crossingstest_(polys[i].ppts,
                                             &polys[i].numpts,pt,&dum)))
                      val = rgnattribs[(int)inp.regionlist[i*4+2]];
                  }
                if (found) vhb[k*8+n] = val;
              }
          }
      }

  /***********************/
  /******* Cleanup *******/
  /***********************/
  /* Free all allocated arrays, including those allocated by Triangle. */
  if (rgnattribs) free(rgnattribs);
  if (inp.regionlist==out.regionlist)  out.regionlist=NULL;
  for (k=0;k<inp.numberofregions;k++)
      if (polys[k].ppts) free(polys[k].ppts);
  cleanio( &inp );

  return 0;
}

#if XY
#ifdef SUN
int vcfrompoly_( vhb,ex,ey,lem,nor,ivv,nup,mesh_ne,polyfile,comments,
                  ierr)
int *nup,*mesh_ne,*ivv,*lem,*nor,*nellep,**ielle,*ierr;
float *ex, *ey, *vhb;
char *polyfile,*comments;
#else
int vcfrompoly_( float *vhb,float *ex,float *ey,int *lem,int * nor,
             int *ivv,int *nup,int *mesh_ne,
             char *polyfile,char *comments,int *ierr)
#endif
{
  struct triangulateio out,inp;
  struct flagvals triflags;
  int i,j,k,n,a;
  float val,dflt,curr;
  float eps;
  float *rgnattribs=NULL;
  double ppts[6], pt[2];
  int maxrgns,ne=0,found,num=3,dum=0;
  int *dum_adr;
/******************
this is a quick hack only using regions with constant vc.
Note area flag is <0 to override value in poly file.
Just want rough triangulation to use previous attribute code
for interpolating regions values onto basil mesh
********************/

  triflags.area = 0.025*0.025/**0.5*/; /* area of tri with minNodeSep sides */
  triflags.quality = -1;
  triflags.midptnodes = 1;
  triflags.bndpts = 1;

  initio( &inp );
  initio( &out );

  if ((*ierr=tripoly(&out,&inp,polyfile,&triflags,comments,&dum,&dum_adr))!=0)
      return(1);

  ne=out.numberoftriangles;
  eps = 1.0e-6;

  /*
   * ivv    visc (vhb[i,1-7])    se (vhb[i,8])
   *  0       constant(1)        constant(1)    No vhb array
   *  1       variable           constant(1)
   *  2       constant           constant(>1)     No vhb array if visc=1
   *  3       variable           constant(>1)
   *  4       variable           variable
   */
  if (out.triangleattributelist!=NULL) {
      maxrgns = (int)out.triangleattributelist[0];
      for (i = 0; i < ne; i++) {
          k = (int)out.triangleattributelist[i];
          if (k > maxrgns) maxrgns = k;
      }
      if ((rgnattribs= (float *)malloc((maxrgns+1) *
                                    sizeof(float)))==NULL) return(1);
      dflt = 1.0;
      for (i=0;i<maxrgns;i++) rgnattribs[i] = dflt;
      if (attribute_set(polyfile,VISCOSITY_KEY,
                        maxrgns+1,rgnattribs,&dflt)) {
          for (i = 0; i < *mesh_ne * 8; i++) vhb[i]=dflt;
          for (i = 0; i < *mesh_ne; i++) vhb[i*8+7]=1.0;
          if (*ivv<1)  *ivv=1;
          if (dflt>1.0 && *ivv<2)  *ivv=2;
          for (k = 0; k < *mesh_ne; k++) {
              for (n = 0; n < 3; n++) {
                pt[0] = ex[nor[lem[k*6+n]-1]-1];
                pt[1] = ey[nor[lem[k*6+n]-1]-1];
                /*
                 * for each mesh vertex, find the vc value
                 * from the poly data
                 */
                found = 0;
                val = dflt;
                  for (i = 0; i < ne && !found; i++) {
                      for (j = 0; j < 3; j++) {
                        a = out.trianglelist[i*6+j];
                        ppts[j*2] = out.pointlist[a*2];
                        ppts[j*2+1] = out.pointlist[a*2+1];
                      }
                    if (found=crossingstest_(ppts,&num,pt,&dum))
                      val = rgnattribs[(int)out.triangleattributelist[i]];
                  }
                if (found) vhb[k*8+n] = val;
              }
              vhb[k*8+3] = (vhb[k*8+2]+vhb[k*8])/2.0;
              vhb[k*8+4] = (vhb[k*8+1]+vhb[k*8])/2.0;
              vhb[k*8+5] = (vhb[k*8+1]+vhb[k*8+2])/2.0;
              /*
               * check this!!!
               */
              vhb[k*8+6] = (vhb[k*8+1]+vhb[k*8+2]+vhb[k*8])/3.0;
          }
      }
  }

  /***********************/
  /******* Cleanup *******/
  /***********************/
  /* Free all allocated arrays, including those allocated by Triangle. */
  if (rgnattribs) free(rgnattribs);
  if (dum_adr) free(dum_adr);
  if (inp.regionlist==out.regionlist)  out.regionlist=NULL;
  cleanio( &inp );
  cleanio( &out );

  return 0;
}
#endif

/*****************************************************************************/
/*                                                                           */
/*  showbasilvar()   Print the basil variables.                              */
/*                 for debugging                                             */
/*                                                                           */
/*****************************************************************************/

void showbasilvar(nup,ne,nbp,ex,ey,ibc,lem,lemneighbor)

int nup,ne,nbp,*ibc,*lem,*lemneighbor;
float *ex,*ey;
{
  int i;

  printf("Basil sizes: nup=%d ne=%d nbp=%d\n",nup,ne,nbp);

  printf("Basil points\n");
  for (i = 0; i < nup; i++) {
    printf(" pt=%d ex=%f ey=%f \n",i,ex[i],ey[i]);
  }

  printf("Basil boundary points\n");
  for (i = 0; i < nbp; i++) {
    printf(" i=%d  node=%d \n",i,ibc[i]);
  }

  printf("Basil elements\n");
  for (i = 0; i < ne; i++) {
      printf(" el=%d nodes=%d %d %d %d %d %d neighbors=%d %d %d \n",i,
               lem[i*6+0],
               lem[i*6+1],
               lem[i*6+2],
               lem[i*6+3],
               lem[i*6+4],
               lem[i*6+5],
               lemneighbor[i*3+0],
               lemneighbor[i*3+1],
               lemneighbor[i*3+2]);
  }
}

void  print_triangulate_data(struct triangulateio *out)
{
  int i,j,k;
  int fltMIN, fltMAX;

  /* limits for boundary id numbers for internal fault bnds */
  bndfltlimits_(&fltMIN,&fltMAX);

  printf("segments\n");
  if (out->segmentmarkerlist) {
    for (i=0;i<out->numberofsegments;i++) {
          if (out->segmentmarkerlist[i]>=fltMIN &&
              out->segmentmarkerlist[i]<=fltMAX) {
                  j=out->segmentlist[i*2];
                  k=out->segmentlist[i*2+1];
                  printf(" %d %d %f %f %d %f %f\n",i,
                         j, out->pointlist[j*2], out->pointlist[j*2+1],
                         k, out->pointlist[k*2], out->pointlist[k*2+1]);
          }
    }
  }
  printf("edges\n");
  if (out->edgemarkerlist) {
    for (i=0;i<out->numberofedges;i++) {
                  j=out->edgelist[i*2];
                  k=out->edgelist[i*2+1];
                  printf("%d %d %d %f %f %d %f %f\n",out->edgemarkerlist[i],i,
                         j, out->pointlist[j*2], out->pointlist[j*2+1],
                         k, out->pointlist[k*2], out->pointlist[k*2+1]);
    }
  }
}

  /***********************************************/
  /******* functions for generating faults *******/
  /***********************************************/


int createafault(struct triangulateio *out,
                 int *mpstart, int *vstart, int *fstart,
                 int *lem, float *ex, float *ey,int *ne,
                 int *ibc,int *ibcngh,int *ibctyp,int *ifbc1,int *ifbc2,
                 int *iseg,
                 int nbpfinal,int *nbpext,int *nseg,int *vxtotal,int *mptotal,
                 struct gen_fault *flt)
// this routine uses list of vertex fault nodes that are specific to a
// fault discontinuity to modify the existing finite element mesh by
// duplication of all nodes on the fault. A fault may be open ended or 
// close on itself. The end nodes of an open fault are not duplicated 
// unless that fault terminates on an external boundary.
// This routine replaces createfault and createclosedfault. 
// ifbc1 and ifbc2 are now established in this routine.
// ibctyp for all fault nodes is initialised to 11
{
    int fltid,  curr, next, away, first, last, curr1, curr2, next1, next2;
    int *tri_1=NULL, *tri_2=NULL, *tri_3=NULL, *midpts=NULL, *newnext=NULL;
    int nums, i, j, k, nelf, inxt, nnxt, itmp;
    int nel1, nel2, j1curr, j1next, j2curr, j2next, j2midp, midpt1;
    int tmplem, tmpnxt, tmpmid, faultclosed=0;
    int vindex,mpindex;
//  int vfirst, mfirst;
    int ift1,ift2,igt1=0,igt2=0,jrp1,jrp2,jacw,jcwn,kn,ierr=0;
    int counting=0, mounting=0;
    struct node *fltnode;

//  check sizes of relevant arrays already defined as this routine may
//  be repeatedly called to add additional fault boundaries
//  we start adding vertex nodes at vindex and midpoint nodes at mpindex
//  followingg operation done once before sequence of calls to createafault
    vindex=0;
    while (ibc[vindex]!=-1) vindex++;
    vindex--;
//  vfirst=vindex+1;
    mpindex=nbpfinal-1;
    while (ibc[mpindex]==-1) mpindex--;
//  printf("entering createafault: vindex = %i mpindex = %i nbpfinal = %i\n",vindex,mpindex,nbpfinal);
    fltid=flt->id;
    fltnode=flt->ordered_nodes;
    nums=0; 
    struct node *pos = fltnode;
    first = pos->id;
    while (pos->next) { pos=pos->next; nums++; }
//  printf(" %i",first); while (pos->next) { pos=pos->next; printf(" %i",pos->id); nums++; }
//  printf("\n");
    last = pos->id;
//  increment first and last for use with lem, ibc
    first++; last++;

// check if fault is closed and allocate space for temporary arrays needed
    faultclosed=0;
    if (first == last) faultclosed=1;
    if((tri_1 = (int *) malloc((2*nums) * sizeof(int)))==NULL) return(1);
    if((tri_2 = (int *) malloc((2*nums) * sizeof(int)))==NULL) return(1);
    if((tri_3 = (int *) malloc((2*nums) * sizeof(int)))==NULL) return(1);
    if((midpts = (int *) malloc((nums) * sizeof(int)))==NULL) return(1);
    if((newnext = (int *) malloc((2*nums) * sizeof(int)))==NULL) return(1);
    
// find the triangles associated with each segment, store element number
// and local node number of leading vertex node in tri_1 and tri_2
// tri_1 and tri_2 refer to elements on +ve and -ve side of fault respectively
// node numbers in tri_1 are retained, in tri_2 are replaced by new nodes
// for a closed fault the tri_1 elements are inside the fault boundary
// based on sequence of nodes in elements being in inti-clockwise order
//
    curr=fltnode->id+1;
    for(k=0;k<nums;k++){
        next=fltnode->next->id+1;
  //    printf("finding tris: k = %i curr = %i next = %i\n",k,curr,next);
        newnext[2*k]=next;
        nelf=findTriangV(k,curr,next,lem,*ne,tri_1,tri_2);
        if (nelf!=2)  {
            printf("Can't find element for nodes k = %i curr = %i next = %i on fault %i\n",k,curr,next,fltid);
            ierr=1;
            if(*fstart > 0){
  //    check if curr/next has been replaced by one of the newly added fault nodes; can be a problem
  //    if current fault is being joined to an earlier fault on the wrong side.
              for(kn = 0; kn < *fstart; kn++){
                 if(ifbc1[kn]==curr){
                    printf("Substituting fault pair node %i for node curr = %i\n",-ifbc1[kn+1],curr);
                    curr=-ifbc1[kn+1];
  //                printf("kn = %i curr = %i ifbc1[kn+1] = %i\n",kn,curr,ifbc1[kn+1]);
                 }
                 else if(ifbc1[kn]==next){
                    printf("Substituting fault pair node %i for node next = %i\n",-ifbc1[kn+1],next);
                    next=-ifbc1[kn+1];
  //                printf("kn = %i next = %i ifbc1[kn+1] = %i\n",kn,next,ifbc1[kn+1]);
                 }
              }
              nelf=findTriangV(k,curr,next,lem,*ne,tri_1,tri_2);
              if(nelf==2)ierr=0;
            }
            if(ierr)return ierr;
        }
        j2midp=3 + (tri_2[2*k+1] + 1) % 3;
        nel2=tri_2[2*k];
        nel1=tri_1[2*k];
  //      printf("ccf lem1 = %i %i %i %i %i %i\n",lem[6*nel1],
  //        lem[6*nel1+1],lem[6*nel1+2],lem[6*nel1+3],lem[6*nel1+4],lem[6*nel1+5]); 
  //      printf("ccf lem2 = %i %i %i %i %i %i\n",lem[6*nel2],
  //        lem[6*nel2+1],lem[6*nel2+2],lem[6*nel2+3],lem[6*nel2+4],lem[6*nel2+5]); 
        midpts[k]=lem[6*tri_2[2*k] + j2midp];
        fltnode=fltnode->next;
        curr=fltnode->id+1;
    }     // for k

    if(faultclosed){
       printf("For (closed) fault = %i first and last nodes are %i and %i with %i segments\n",
              fltid,first,last,nums);
    }
    else{
       printf("For  (open)  fault = %i first and last nodes are %i and %i with %i segments\n",
              fltid,first,last,nums);
    }

  // set ift1 and/or ift2 if fault terminates on external boundary, first and last may be
  //  reversed relative to input sequence defined in polyfile
    ift1 = 0;
    ift2 = 0;
    if (!faultclosed){
       for(k=0;k<vindex;k++) {
          if(ibc[k]==first)ift1=k;
          if(ibc[k]==last)ift2=k;
       }
  // add unattached internal fault end-points to IBC, IBCTYP arrays
       if(ift1){
          printf("     fault %i is connected to node %i on a prior boundary segment, ift1 = %i\n",
                       fltid,ibc[ift1],ift1);
       }
       else{
          vindex++;
          ibc[vindex]=first;
	  counting++;
          ibctyp[vindex]=15;
          ibctyp[vindex+nbpfinal]=15;
//	  printf("fltend, first = %i added in position vindex = %i \n",first,vindex);
       }
       if(ift2){
          printf("     fault %i is connected to node %i on a prior boundary segment, ift2 = %i\n",
                       fltid,ibc[ift2],ift2);
       }
       else{
          vindex++;
          ibc[vindex]=last;
	  counting++;
          ibctyp[vindex]=15;
          ibctyp[vindex+nbpfinal]=15;
//	  printf("fltend, last = %i added in position vindex = %i \n",last,vindex);
       }
    }

 // printf("before new vertices, *mpstart = %i *vstart = %i, vindex = %i\n",*mpstart,*vstart,vindex);
   //  work through the list of segments identified in tri_1 and tri_2
   //   in each pass of this loop nodes are identified as:
   //        curr2---------next2    new nodes added
   //        curr1---------next1    existing nodes retained
   //   'next' vertex fault node (== 'curr' fault node of the next pass)
   //   -----curr-----next----    ----  at step i
   //   -----    -----curr----next----  at step i+1
   //   for open fault, end points not duplicated unless on prior boundary (ift1 / ift2)

    for(i=0;i<nums;i++) {
	nel1=tri_1[2*i];
	j1curr=tri_1[2*i+1];
	j1next=(j1curr+1) % 3;
	curr=lem[6*nel1+j1curr];
	next=lem[6*nel1+j1next];
 //     printf("createafault fltid = %i curr = %i next = %i\n",fltid,curr,next); 
        if((i==0) && ift1){     // duplicate first point if node ift1 already in ibc
             (*fstart)++;
	     vindex++;
             ifbc1[*fstart]=ibc[ift1];
             ibctyp[ift1]=11;
             ibctyp[ift1+nbpfinal]=11;
             ifbc2[*fstart]=vindex+1;
             (*vstart)++;                 // new vertex node number of duplicate node
	     ex[*vstart]=ex[curr-1];
    	     ey[*vstart]=ey[curr-1];
	     counting++;
	     ibc[vindex]=igt1=*vstart+1;       // duplicate node added
             ibctyp[vindex]=11;
             ibctyp[vindex+nbpfinal]=11;
             (*fstart)++;
             ifbc1[*fstart]=-(*vstart+1);
             ifbc2[*fstart]=ibc[ift1];
             nel2=tri_2[0];
	     j2next=tri_2[1];
	     j2curr=(j2next+1) % 3;
	     lem[6*nel2+j2curr]=*vstart+1;  // reset lem value for 'first' duplicate node
        }
   //   printf("i = %i nel1 = %i, curr = %i, next = %i last = %i\n",i,nel1,curr,next,last);
    
   //   for the existing vertex node on the + side of the fault
	if((next != last) || faultclosed){
	     vindex++;                    // ibc address
	     counting++;
	     ibc[vindex]=next;            // add existing node to ibc; ift2 node already present
             ibctyp[vindex]=11;
             ibctyp[vindex+nbpfinal]=11;
        }
	if((next != last) || faultclosed || ift2){
             (*fstart)++;                 // add ifbc1 entries, also if boundary reached
             ifbc1[*fstart]=next;
             ifbc2[*fstart]=vindex+1;
   //   now define and add the entries for the duplicate vertex node in tri_2
   //   vertex lem entries will be reset after loop i finished in multiple elements
	     (*vstart)++;                 // new vertex node number of duplicate node
	     nnxt=newnext[2*i+1]=*vstart+1;
	     ex[*vstart]=ex[next-1];
	     ey[*vstart]=ey[next-1];
	     vindex++;                    // ibc address
	     counting++;
	     ibc[vindex]=nnxt;
	     if((next==last)&&(ift2))igt2=nnxt;
             ibctyp[vindex]=11;
             ibctyp[vindex+nbpfinal]=11;
             (*fstart)++;
             ifbc1[*fstart]=-(*vstart+1);
             ifbc2[*fstart]=vindex+1;

	   //   resetting vertex lem entries for 'next' node, find other relevant elements.
	   //   step from element to element using midpoints until fault or boundary found again
	     nel2=tri_2[2*i];
	     j2next=tri_2[2*i+1];
	     j2curr=(j2next+1) % 3;
	     j2midp=3 + ((j2next+3) % 3);
	     tmplem=nel2;
	     tmpnxt=j2next;
	     inxt=(i+1) % nums;
	     int done=0; itmp=0;
	     while (!done){
		      itmp++;
		      if(itmp>=20){
			 printf("Problem finding element that contains next = %i away = %i\n",next,away);
			 return 1;
		      }
		      tmpmid=lem[6*tmplem + j2midp];
		      lem[6*tmplem+tmpnxt]=nnxt;                    // reset lem value for 'next' node
		      if(tmpmid == midpts[inxt]){                   // now back on the fault
			 done++; 
		      }
		      else{                                         // check the adjacent triangle
			 nelf=findTriangM(i,tmpmid,tmplem,lem,*ne,tri_3);
                         if (nelf==0) {done++;}                     // must be external boundary
			 else if (nelf>1)  {
			    printf("Problem %i finding element for nodes itmp = %i next = %i tmpmid = %i\n",
			       nelf,itmp,next,tmpmid);
			    return 1;
			 }
                         else{                                      // next triangle
			    tmplem=tri_3[2*i];
			    j2midp=3 + (tri_3[2*i+1] + 2) % 3;
			    tmpnxt=j2midp-3;
                         }
		      }
	     }                                              // end of steps through adjacent elements
	}                                                   //  if((next != last) || faultclosed)
    }                                                       // end of main:loop on segment: i
    *vxtotal=vindex;
 // printf("vindex = %i vstart = %i after main vertex loop\n",vindex,*vstart);

	 
   //   treatment of midpoints separated from vertices, as midpoints used to fix vertices
 // mfirst=mpindex;       // save for use in next operation
    for(i=0;i<nums;i++) {

   //   for the existing midpoint node on the + side of the fault
        mpindex++;
        midpt1=midpts[i];
	mounting++;
        ibc[mpindex]=midpt1;
        ibctyp[mpindex]=11;
        ibctyp[mpindex+nbpfinal]=11;
        (*fstart)++;
        ifbc1[*fstart]=midpt1;
        ifbc2[*fstart]=mpindex+1;
	nel1=tri_1[2*i];
	j1curr=tri_1[2*i+1];
	j1next=(j1curr+1) % 3;
	curr=lem[6*nel1+j1curr];
	next=lem[6*nel1+j1next];
        ibcngh[mpindex]=curr;
        ibcngh[mpindex+nbpfinal]=next;
   //   printf("midpoint1: mpindex = %i midpt1 = %i ibcngh1 =%i ibcngh2 = %i\n",mpindex,midpt1,
   //                           ibcngh[mpindex],ibcngh[mpindex+nbpfinal]);
   //   now define and add the entries for the duplicate midpoint node in tri_2
        (*mpstart)++;
        ex[*mpstart]=ex[midpt1-1];
        ey[*mpstart]=ey[midpt1-1];
        mpindex++;
	mounting++;
        ibc[mpindex]=*mpstart+1;
        ibctyp[mpindex]=11;
        ibctyp[mpindex+nbpfinal]=11;
        (*fstart)++;
        ifbc1[*fstart]=-(*mpstart+1);
        ifbc2[*fstart]=mpindex+1;
        nel2=tri_2[2*i];
        j2next=tri_2[2*i+1];
        j2curr=(j2next+1) % 3;
	curr=lem[6*nel2+j2curr];
	next=lem[6*nel2+j2next];
        ibcngh[mpindex]=next;
        ibcngh[mpindex+nbpfinal]=curr;
   //   printf("midpoint2: mpindex = %i midpt1 = %i ibcngh1 =%i ibcngh2 = %i\n",mpindex,midpt1,
   //                           ibcngh[mpindex],ibcngh[mpindex+nbpfinal]);

   //   reset the lem entry for the new midpoint node  (only one, unlike vertices}
        nel2=tri_2[2*i];
        j2next=tri_2[2*i+1];
        j2midp=3 + ((j2next+4) % 3);
        lem[6*nel2+j2midp]=*mpstart+1;
   //	printf("createafault: nel2 = %i j2midp = %i *mpstart+1 = %i\n",nel2,j2midp,*mpstart+1);
   //   add in the new iseg entries - original ones still defined
        iseg[(*nseg)*3]=curr;
        iseg[(*nseg)*3+1]=next;                            // end of current segment
	iseg[(*nseg)*3+2]=fltid;                           // segment label
	(*nseg)++;
   //  if ift1 or ift2 are non-zero, existing segment end points need to be changed; see below
    }     // end of midpoint loop on segment: i
    *mptotal=mpindex;

   //   now set / reset the ibcngh arrays for the vertex nodes and
   //   note that special treatment for ends of unattached faults
   //	printf("createafault: end0 = %i end2 = %i\n",ibc[vfirst],ibc[*vxtotal]);
    for(i=0; i<nums; i++){
	nel1=tri_1[2*i];
	j1curr=tri_1[2*i+1];
	j1next=(j1curr+1) % 3;
        nel2=tri_2[2*i];
	j2next=tri_2[2*i+1];
        j2curr=(j2next+1) % 3;
	curr1=lem[6*nel1+j1curr];
	next1=lem[6*nel1+j1next];
	curr2=lem[6*nel2+j2curr];
	next2=lem[6*nel2+j2next];
   //	printf("createafault: i = %i curr1 = %i curr2 = %i next1 = %i next2 = %i\n",
   //			i,curr1,curr2,next1,next2);
   //    attached nodes should be addressed by the following sets of conditions
        for(j=0; j<*vxtotal+1; j++){                  // end nodes may be before vfirst
           if (ibc[j]==curr1){
              ibcngh[j+nbpfinal]=next1;             // anti-clockwise connection
	   }
	   else if (ibc[j]==next1){                 // clockwise connection
              ibcngh[j]=curr1;
	   }
	   else if (ibc[j]==next2){ 
              ibcngh[j]=curr2;                      // clockwise connection
	   }
	   else if (ibc[j]==curr2){
              ibcngh[j+nbpfinal]=next2;             // anti-clockwise connection
	   }
    //   connections on unattached node: curr1 and curr2 should be the same node
	   if((i==0) && (curr1==curr2) && (ibc[j]==curr1)){
              ibcngh[j+nbpfinal]=next1;             // anti-clockwise connection
              ibcngh[j]=next2;                      // clockwise connection
           }
    //   connections on unattached node: next1 and next2 should be the same node
	   if((i==(nums-1)) && (next1==next2) && (ibc[j]==next1)){
              ibcngh[j+nbpfinal]=curr2;             // anti-clockwise connection
              ibcngh[j]=curr1;                      // clockwise connection
           }
        }
    }

   // final task, locate and replace entries in ibcngh and iseg where igt1(2) should replace
   // ibc[ift1(2)]. Search all entries in ibcngh; only the ones that need to be changed will be found
   // we then look to replace the same nodes in the seg array for external boundary
   // these operations apply when fault terminates on a pre-defined structure
    if(ift1){
       jrp1=0;
       jcwn=0;
       for (j=0; j<*vxtotal+1; j++){
          if(ibc[j]==igt1) jrp1=j;                 // position in ibc of node igt1
          if(ibcngh[j]==ibc[ift1]) jcwn=j;         // position of clockwise node to ift1
       }   	  
       if((igt1*jrp1*jcwn)==0){
          printf("problem resetting node connection in createafault: ift1, igt1, jrp1,jcwn = %ii %i %i %i\n",
			  ift1,igt1,jrp1,jcwn);
       }
  //   printf("IFT1 = %i igt1 = %i jrp1 = %i jcwn = %i ibc[jrp1] = %i ibcngh[jcwn] = %i\n",
  //           ift1,igt1,jrp1,jcwn,ibc[jrp1],ibcngh[jcwn]);
       ibcngh[jcwn+nbpfinal]=igt1;                 // new anti-clockwise link to jcwn
       ibcngh[jrp1]=ibc[jcwn];                     // new clockwise link to igt1

       for (k=0; k<*nseg; k++){
          if(iseg[k*3+2]!=fltid){                  // if a prior structure
             if((iseg[k*3]==ibc[ift1])&&(iseg[k*3+1]==ibc[jcwn]))
                          iseg[k*3]=igt1;
             if((iseg[k*3+1]==ibc[ift1])&&(iseg[k*3]==ibc[jcwn]))
                          iseg[k*3+1]=igt1;
          }        // fltid ok
       }           // for k
    }              // if ift1
    if(ift2){
       jrp2=0;
       jacw=0;
       for (j=0; j<*vxtotal+1; j++){
          if(ibc[j]==igt2) jrp2=j;                    // position in ibc of node igt2
          if(ibcngh[j+nbpfinal]==ibc[ift2]) jacw=j;  // position of anti-clockwise node to ift2
       }
       if((igt2*jrp2*jacw)==0){
          printf("problem resetting node connection in createafault: ift2, igt2, jrp2, jacw = %i %i %i %i\n",
			  ift2,igt2,jrp2,jacw);
       }
   //  printf("IFT2 = %i igt2 = %i jrp2 = %i jacw = %i ibc[jrp2] = %i ibcngh[jacw] = %i\n",
   //          ift2,igt2,jrp2,jacw,ibc[jrp2],ibcngh[jacw]);
       ibcngh[jacw]=igt2;                          // new clockwise link to jacw
       ibcngh[jrp2+nbpfinal]=ibc[jacw];            // new anti-clockwise link to igt2

       for (k=0; k<*nseg; k++){
          if(iseg[k*3+2]!=fltid){                  // if a prior structure
             if((iseg[k*3]==ibc[ift2])&&(iseg[k*3+1]==ibc[jacw]))
                          iseg[k*3]=igt2;
             if((iseg[k*3+1]==ibc[ift2])&&(iseg[k*3]==ibc[jacw]))
                          iseg[k*3+1]=igt2;
          }        // fltid ok
       }           // for k
    }              // if ift2

  //printf("createafault: counting = %i mounting = %i \n",counting,mounting);
  // next block only for checking iseg if problems arise
 /* for(i=0; i<(*nseg); i++){
       printf("i = %i iseg[3*i] = %i iseg[3*i+1] = %i iseg[3*i+2] = %i\n",
               i,iseg[3*i],iseg[3*i+1],iseg[3*i+2]);
    }

  // next block only for checking ibcngh if problems arise
    printf("new vertexnodes, vfirst = %i *vxtotal = %i \n",vfirst,*vxtotal);
    for(i=vfirst; i<*vxtotal+1; i++){
       printf("i = %i ibc[i] = %i ibcngh[i] = %i ibcngh[i+nbpfinal] = %i\n",
               i,ibc[i],ibcngh[i],ibcngh[i+nbpfinal]);
    } */

 /* int newmnodes=nums*2;
    printf("new midpoint nodes : %i\n",newmnodes);
    for(i=0; i<newmnodes; i++){
       itmp=mfirst+1+i;
       printf("itmp = %i ibc[itmp] = %i ibcngh[itmp] = %i ibcngh[itmp+nbpfinal] = %i\n",
               itmp,ibc[itmp],ibcngh[itmp],ibcngh[itmp+nbpfinal]);
    } */

    if (tri_1) free(tri_1);
    if (tri_2) free(tri_2);
    if (tri_3) free(tri_3);
    if (midpts) free(midpts);
    if (newnext) free(newnext);
    return 0;
}       // end of createafault

unsigned char is_bnd_pt(struct triangulateio *out,int index)
{
    unsigned char bnd=0;
    int extMIN, extMAX, fltMIN, fltMAX;

    /* limits for boundary id numbers for external bnds */
    bndextlimits_(&extMIN,&extMAX);
    /* limits for boundary id numbers for internal fault bnds */
    bndfltlimits_(&fltMIN,&fltMAX);
    if (((out->pointmarkerlist[index]>=extMIN)
             && (out->pointmarkerlist[index]<=extMAX))
       /*|| ((out->pointmarkerlist[index]>=fltMIN)
             && (out->pointmarkerlist[index]<=fltMAX))*/)
        bnd=1;
    return(bnd);
}

#ifndef  ANSI_DECLARATORS
struct gen_fault *newfltnode(id)
int id;
#else
struct gen_fault *newfltnode(int id)
#endif
{
    struct gen_fault *node = (struct gen_fault *)
                               malloc(sizeof(struct gen_fault));
    if (node) {
        node->id = id;
        node->end1=node->end2=-1;
        node->ordered_nodes=0;
        node->next=0;
    }
    return(node);
};

unsigned char in_flt_list(struct gen_fault *list,int id)
{
    unsigned char found = 0;
    struct gen_fault *node = list;
    while (node && node->id != id) node=node->next;
    if (node) found = 1;
    return(found);
}

void print_lists(struct gen_fault *list)
{
    struct node *nnode;
    while (list) {
        printf("flt=%d %d",list->id,list->end1);
        nnode=list->ordered_nodes;
        while (nnode) {
            printf(" %d",nnode->id);
            nnode=nnode->next;
        }
        printf(" %d\n",list->end2);
        list=list->next;
    }
}

void insert_in_flt_list(struct gen_fault **list,struct gen_fault *fltnode)
{
    struct gen_fault *node = 0;
    if (list) node=*list;
    if (node && node->id<fltnode->id) {
        while (node->next && (node->next->id<fltnode->id))
                node=node->next;
        fltnode->next = node->next;
        node->next = fltnode;
    }
    else {
        *list = fltnode;
    }
}

void free_flt_list(struct gen_fault *list)
{
    struct node *tmp,*nodes;
    if (list) {
        free_flt_list(list->next);
        nodes=list->ordered_nodes;
        while (nodes){
            tmp=nodes;
            nodes=tmp->next;
            free (tmp);
        }
        free(list);
    }
}

/*
 * insert flt id in list of faults
 * insert pt in appropriate list of flt nodes, ordered on x,y coords
 */
int findfltpts(struct gen_fault **fltlist,struct triangulateio *out,
               int *nfp,int *newvert,int *newseg,int *newibc,int *ibcmoffv)
// entries produced in fltlist are all one less than node numbers used
//  in fortran entities
{
    int i,j,m,next,id,fltid,nbseg=0,err=0;
    int finished, segno;
    int fltMIN, fltMAX, extMIN, extMAX;
    int *fltsegs=NULL;
    int ift1,ift2,ipta,iptb,ka,kb,met,news,setn;
    struct gen_fault *list;
    struct gen_fault *nfnode=NULL;
    struct node *nnode=NULL, *pos;
    int found=0,fltcnt=0;

    *nfp = *newvert = *newseg = *newibc = 0;
    bndfltlimits_(&fltMIN,&fltMAX);
    bndextlimits_(&extMIN,&extMAX);
    for (i = 0; i < out->numberofsegments; i++) {
	// count the number of segments in the external boundary
        if (out->segmentmarkerlist[i] >= extMIN    &&
                out->segmentmarkerlist[i] <= extMAX &&
                    !in_flt_list(*fltlist,out->segmentmarkerlist[i])) {
	    nbseg++;
        }
        // initialise a list for each fault with distinct label
        if (out->segmentmarkerlist[i] >= fltMIN    &&
                out->segmentmarkerlist[i] <= fltMAX &&
                    !in_flt_list(*fltlist,out->segmentmarkerlist[i])) {
            if ((nfnode=newfltnode(out->segmentmarkerlist[i]))==NULL)
                    return(1);
            insert_in_flt_list(fltlist,nfnode);
        }
    }
    printf("External boundary has %i segments = no. of midpoints = no. of vertices\n",nbseg);
    list=*fltlist;
    while (list && !err) {
      
  //    for each labelled structure, count the number of segments
        fltid=list->id;
        fltcnt++;
        news=0;
        for (i=0;i<out->numberofsegments && !err;i++)
            if (out->segmentmarkerlist[i]==fltid) news++;
        if (news<3) {
            printf("Need a minimum of 3 vertex points to create fault %d\n",
                            fltid);
            err=1;
        }
  
  //    store the segment indices to help find the nodes
        if ((fltsegs = (int *) malloc(news*2 * sizeof(int)))==NULL)
          return(1);
        for (j=0;j<news*2;j++) fltsegs[j] = 0;
        for (i=0,j=0;i<out->numberofsegments && !err;i++)
            if (out->segmentmarkerlist[i]==fltid) fltsegs[j++] = i;
       
  //    find two end nodes i.e. they occur only once in the list
        ift1=0;
        ift2=0;
        for (j=0; j<news; j++) {
           for (ka=0; ka<2; ka++) {
	      ipta=out->segmentlist[fltsegs[j]*2+ka];
	      met=0;
              for (m=0; m<news; m++) {
                 if(m!=j){
	            for (kb=0; kb<2; kb++) {
		       iptb=out->segmentlist[fltsegs[m]*2+kb];
		       if(ipta==iptb)met++;
                    }
		 }
              }
	      if((met==0)&&(ift1!=0)&&(ift2==0))ift2=ipta;
	      if((met==0)&&(ift1==0))ift1=ipta;
           }
	}
  //   ift1 and ift2 here are the start and end nodes for the fault
        printf("Fault structure %i has %i segments with end nodes %i and = %i \n",fltid,news,ift1+1,ift2+1);

  //   if the fault is closed, we add 2*news new nodes, we add 4*news new
  //   entries in ibc and we add 4*news entries to the fault arrays (ifbc etc)
  //   the offset between vertex and midpoint nodes (ibcmoffv) is 2*news
  //   first point in list is arbitrary, as fault reconnects to itself
  //  
        if((ift1==0) && (ift2==0)){
           id=out->segmentlist[fltsegs[1]*2];
           *nfp = *nfp + 4*news;              // space for fault arrays
           *newibc = *newibc + 4*news;        // space for IBC 
           *newvert=*newvert+news;            // space for EX,EY
        }

  //  if fault is open and unconnected to prior structures, news+1 existing v nodes
  //  and news-1 new v nodes are added to IBC, news-1 pairs of v nodes + news pairs
  //  of midpoint nodes are added to ifbc, news-1 nodes added to EX and EY
  //
        else {
	   id = ift1;
           *nfp = *nfp + 4*news - 2;          // space for fault arrays
           *newibc = *newibc + 4*news;        // space for IBC 
           *newvert=*newvert+news-1;          // space for EX,EY
        }
	*newseg=*newseg+news;                 // space for ISEG
        *ibcmoffv = *ibcmoffv + 2*news;       // IBC offset V-M

   //  check if fault end points are also on an earlier defined boundary
   //  (i.e. lesser fltid), in which case additional two entries needed in
   //  fault arrays and for EX,EY, extra node in ibc already taken care of
        if(ift1){
           found=0;
           for (i = 0; i < out->numberofsegments && !found; i++) {
              segno=out->segmentmarkerlist[i];
              if(((segno >= extMIN) && (segno <= extMAX)) || 
                 ((segno >= fltMIN) && (segno < fltid  ))) {
                   if((out->segmentlist[i*2]  ==ift1) || 
                      (out->segmentlist[i*2+1]==ift1)) {
                      found = out->segmentmarkerlist[i];
                      printf("Fault structure %i connects to prior structure at node number %i\n",
                             fltid,ift1+1);
                   }
              }
           }
           if(found){
	      *nfp = *nfp + 2;                // extra entries for fault
	      *newvert=*newvert+1;            // extra node on join
	   }
        }

        if(ift2){
           found=0;
           for (i = 0; i < out->numberofsegments && !found; i++) {
              segno=out->segmentmarkerlist[i];
              if(((segno >= extMIN) && (segno <= extMAX)) || 
                 ((segno >= fltMIN) && (segno < fltid  ))) {
                   if((out->segmentlist[i*2]  ==ift2) || 
                      (out->segmentlist[i*2+1]==ift2)) {
                      found = out->segmentmarkerlist[i];
                      printf("Fault structure %i connects to prior structure at node number %i\n",
                             fltid,ift2+1);
                   }
              }
           }
           if(found){
	      *nfp = *nfp + 2;                // extra entries for fault
	      *newvert=*newvert+1;            // extra node on join
	   }
        }
   //	printf("findfltpts: fltid = %i *newseg = %i *newvert = %i *nfp = %i *newibc = %i\n",
   //			fltid,*newseg,*newvert,*nfp,*newibc);

   //  start the list for this fault with node ift1
        if ((nnode=newnode(id))==0) return(1);
        list->ordered_nodes = nnode;
        pos = nnode;

   //  and add points according to connections
        finished=0;
        while (finished<news && !err) {
           found=0;
           for (j=0; j<news ;j++) {
              setn=0;
              if (out->segmentlist[fltsegs[j]*2]==id){
                 next=out->segmentlist[fltsegs[j]*2+1];
                 setn=1;
              }
              if(out->segmentlist[fltsegs[j]*2+1]==id) {
                 next=out->segmentlist[fltsegs[j]*2];
                 setn=1;
              }
              if(setn==1){
    //           printf("j = %i setn = %i id = %i next = %i\n",j,setn,id,next);
                 if ((nnode=newnode(next))==NULL) return(1);
                 pos->next = nnode;
                 pos = nnode;
                 id = next;
                 finished++;
                 found=1;
              }
           }
	}
   //   add final point and close list
        if ((nnode=newnode(id))==NULL) return(1);
        pos = nnode;
        pos->next=NULL;

   // for debugging
   /*   pos = list->ordered_nodes;
        while (pos) {printf("%d ",pos->id); pos = pos->next;}
        printf("\n"); */

 // free up fltsegs for the next fault in the list
        if (fltsegs) { free(fltsegs); fltsegs = NULL; }
        list=list->next;
    }  // end of list of faults
    printf("Faults identified: %i Fault array size: = %i Places added to IBC = %i\n",fltcnt,*nfp,*newibc);

    return(0);
}   // end of findfltpts

void countpbpts(int *nn,int *countpb,int pbxy,float *ex,float *ey,float pdist)
//  this routine added to compute the contribution of periodic boundaries
//  to IFBC1, etc prior to allocation of arrays by trimesh.c
//  pdist must be correct for established mesh, or method will not work
{
    int i,j,count=0;
    float eps=1.e-5,pa1,pb1,pa2,pb2;
    if(pbxy==1){
      for (i=0;i<*nn;i++){
        pa1=ex[i];
        pb1=ey[i];
        for (j=0;j<*nn;j++){
          pa2=ex[j];
          pb2=ey[j];
          if((fabsf(pa2-pa1-pdist)<eps)&&(fabsf(pb1-pb2)<eps)) count++;
        }
      }
    }
    else if(pbxy==2){
      for (i=0;i<*nn;i++){
        pb1=ex[i];
        pa1=ey[i];
        for (j=0;j<*nn;j++){
          pb2=ex[j];
          pa2=ey[j];
          if((fabsf(pa2-pa1-pdist)<eps)&&(fabsf(pb1-pb2)<eps)) count++;
        }
      }
    }
    (*countpb) = (count+count-1)*2;   // midpoints one less than vertices
    printf("Allowing space for %i nodes in periodic boundaries\n",*countpb);
}

unsigned char is_inex_bnd_pt(int id, int *iseg, int *label, int nseg)
/*
 * this routine designed to return 1 if id is on the external boundary,
 * 2 if on a faulted boundary, 3 if on a non-fault internal boundary or 0 
 * if not on any boundary. For points that may be on more than one type of
 * boundary 1 is precedent over 2 is precedent over 3.  The
 * segment arrays returned from tripoly are used for this purpose 
 */
{
    unsigned char found;
    int i, extMIN, extMAX, fltMIN, fltMAX;

    bndextlimits_(&extMIN,&extMAX);             // get the seg label bounds from Fortran
    bndfltlimits_(&fltMIN,&fltMAX);
    found=0;
      for(i=0; i<nseg; i++) {
        if((iseg[i*2]==id) || (iseg[i*2+1]==id)) {
	  if((label[i]>extMAX)&&(label[i]<fltMIN)) found=3;  // id is an internal boundary node
	}
      }
      for(i=0; i<nseg; i++) {
        if((iseg[i*2]==id) || (iseg[i*2+1]==id)) {
	  if((label[i]>=fltMIN)) found=2;      // id is a fault node
	}
      }
      for(i=0; i<nseg; i++) {
        if((iseg[i*2]==id) || (iseg[i*2+1]==id)) {
	  if((label[i]<=extMAX)) found=1;      // id is an external boundary node
	}
      }
    return(found);
}

unsigned char is_flt_intersect_pt(int id,int fltid,struct gen_fault *list)
{
    unsigned char found=0;
    int fltMIN, fltMAX, fltid2;

    /* limits for boundary id numbers for internal fault bnds */
    bndfltlimits_(&fltMIN,&fltMAX);
    while (list && !found) {
        fltid2=list->id;
        if (fltid2!=fltid) 
            if (in_list(list->ordered_nodes,id)) found=1;
        list=list->next;
    }
    return(found);
}

struct node *newnode(int id)
{
    struct node *nnode = NULL;
    nnode = (struct node *) malloc(sizeof(struct node));
    if (nnode) {
        nnode->id = id;
        nnode->next=NULL;
    }
    return(nnode);
}

unsigned char in_list(struct node *list,int id)
{
    unsigned char found = 0;
    struct node *nnode = list;
    while (nnode && nnode->id != id) nnode=nnode->next;
    if (nnode) found = 1;
    return(found);
}

void insert_node(struct node *nnode,struct node **list,
                REAL *coords)
{
    struct node *pos = 0;
    if (list) pos=*list;
        /* insert in list odered by y then x */
#if XY
    if (pos && (coords[pos->id*2+1]<coords[nnode->id*2+1] ||
                   ( coords[pos->id*2+1]==coords[nnode->id*2+1]    &&
                    coords[nnode->id*2]>=coords[pos->id*2]))) {
        while (pos->next &&
                    coords[pos->next->id*2+1]<coords[nnode->id*2+1])
                pos=pos->next;
        while (pos->next &&
                    coords[pos->next->id*2]<=coords[nnode->id*2] &&
                        coords[pos->next->id*2+1]==coords[nnode->id*2+1])
                pos=pos->next;
        if (pos->id != nnode->id) {
            nnode->next = pos->next;
            pos->next = nnode;
        }
    }
#endif
        /* insert at head */
    else if (!pos || (pos&&(pos->id != nnode->id))) {
        *list = nnode;
        nnode->next = pos;
    }
}

/*  find the one or two elements that include two specific vertex nodes  */

int findTriangV(int count,int start,int end,int *lem, int ne, int *tri_1, int *tri_2)
{
   int i,j,jp1,cnt,found=0;
   cnt=count*2;
// printf("findTriangV: count = %i start = %i end = %i ne  = %i\n",count,start,end,ne);
   for (i=0;i<ne;i++) {
      for (j=0;j<3;j++){
         jp1=(j+1) % 3;
         if ((lem[i*6+j]==start)  && (lem[i*6+jp1]==end)){
   //       if(want==1)printf("fTnA count = %i start = %i end = %i ne = %i\n",count,start,end,ne);
            tri_1[cnt]=i;
            tri_1[cnt+1]=j;
            found++;
            if(found==2)return 2;
         }
         else if ((lem[i*6+j]==end) && (lem[i*6+jp1]==start)) {
   //       if(want==1)printf("fTnB count = %i start = %i end = %i ne = %i\n",count,start,end,ne);
            tri_2[cnt]=i;
            tri_2[cnt+1]=j;
            found++;
            if(found==2)return 2;
         }
      }
   }
   return found;
}

/*  find the neighbouring element that includes a specific midpoint node */

int findTriangM(int count, int midp, int nabr, int *lem, int ne, int *tri_1)
{
   int i,j,cnt,found=0;
   cnt=count*2;
   for (i=0;i<ne;i++) {
      for (j=3;j<6;j++){
         if ((lem[i*6+j]==midp)  && (i!=nabr)){
   //       if(want==1)printf("fTnA count = %i midp = %i nabr = %i ne = %i\n",count,midp,nabr,ne);
            tri_1[cnt]=i;
            tri_1[cnt+1]=j;
            found++;
            return found;
         }
      }
   }
   return found;
}


/*
 *  find all the elements which include both the start node
 *  and one of the nodes in the nbs list
 */
void findTriangles(int start,int *nbs,int cnt,int *lem,int ne,
                   int *tri_nbs, int *tricnt)
{
    int i, j, k;
    unsigned char found=0;

    for (i=0;i<ne;i++) {
        if (lem[i*6]-1==start || lem[i*6+1]-1==start ||
                        lem[i*6+2]-1==start) {
            found=0;
            for (j=0;j<3 && !found;j++) {
                for (k=0;k<cnt && !found;k++) {
                    if (lem[i*6+j]-1==nbs[k]) {
                        tri_nbs[*tricnt]=i;
                        (*tricnt)++;
                        found=1;
                    }
                }
            }
        }
    }
}

/*
 *  find all the neighbours of node, on one side of fault
 *  don't include the prev or next node
 */
void findNeighbourPts(int node, int prev, int next, 
                      struct triangulateio *out,
                      int *cnt, int *nbnodes, int maxnbs)
{
    int i, err=0, nb;
    REAL x[3],y[3];
    double angz=0.0, ref, ang, maxang, pi_x_2=M_PI+M_PI;;
    *cnt = 0;
    /*
     * find coords for prev and node and next
     * segment node-prev will be zero for angle calcs
     * segment node-next will be max for angle calcs
     */
    x[0]= out->pointlist[node*2];
    y[0]= out->pointlist[node*2+1];
    if (prev==node) {
        x[1]= x[0]-( out->pointlist[next*2]-x[0]);
        y[1]= y[0]-( out->pointlist[next*2+1]-y[0]);
        CartesianToPolar(x[1]-x[0],y[1]-y[0],0.0,&ref,&angz);
    }
    else if (next==node) {
        x[1]= out->pointlist[prev*2];
        y[1]= out->pointlist[prev*2+1];
        CartesianToPolar(x[1]-x[0],y[1]-y[0],0.0,&ref,&angz);
    }
    else {
        x[1]= out->pointlist[prev*2];
        y[1]= out->pointlist[prev*2+1];
        CartesianToPolar(x[1]-x[0],y[1]-y[0],0.0,&ref,&angz);
    }
    if (next==node) {
        x[1]= x[0]-( out->pointlist[prev*2]-x[0]);
        y[1]= y[0]-( out->pointlist[prev*2+1]-y[0]);
        CartesianToPolar(x[1]-x[0],y[1]-y[0],0.0,&maxang,&angz);
    }
    else {
        x[1]= out->pointlist[next*2];
        y[1]= out->pointlist[next*2+1];
        CartesianToPolar(x[1]-x[0],y[1]-y[0],0.0,&maxang,&angz);
    }
    maxang=maxang-ref;
    if (maxang<=0.0) maxang += pi_x_2;
    for (i=0; i<out->numberofedges && !err; i++) {
        nb = -1;
        if (out->edgelist[i*2]==node) nb = out->edgelist[i*2+1];
        else if (out->edgelist[i*2+1]==node) nb = out->edgelist[i*2];
        if (nb!=-1 && nb!=next && nb!=prev) {
            x[2]= out->pointlist[nb*2];
            y[2]= out->pointlist[nb*2+1];
            CartesianToPolar(x[2]-x[0],y[2]-y[0],0.0,&ang,&angz);
            ang=ang-ref;
            if (ang<=0.0) ang += pi_x_2;
            if (ang>pi_x_2) ang -= pi_x_2;
            if (ang>0.0 && ang<maxang && *cnt<maxnbs) {
                nbnodes[*cnt]=nb;
                (*cnt)++;
            }
        }
    }
}

int CartesianToPolar(REAL x, REAL y, REAL z, double *angxy, double *angz)
{
    *angz = asin((double)z);
    if (x==0.0) {
        if (y>=0.0) *angxy = M_PI*0.5;
        else        *angxy = M_PI * 1.5;
    }
    else {
        *angxy = atan((double)y/x);
        if (x<0.0) *angxy += M_PI;
        if (*angxy < 0.0) *angxy += 2.0*M_PI;
        else if (*angxy > M_PI) *angxy -= 2.0*M_PI;
    }
    return(0);
}

