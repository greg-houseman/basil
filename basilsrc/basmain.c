/*--------------------------------------------------------------------
*    Basil / Sybil:   basmain.c  1.1  1 October 1998
*
*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
*    See README file for copying and redistribution conditions.
*--------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "trimesh.h"
#include "poly.h"
#include "version.h"
     
#define MALLOC_ERR  -1
#define MAXNAME 33

/*
* Global pointers for all arrays used in basil
*/
int     *lem,*lemngh,*ibc,*ibpole,*ibngh,*ibctyp,*nor,*ielfix,
*ifbc,*ifbc2,*ifeqv,*jfbc1,*jfbc2, *iseg,
*iendp1,*iendp2,
*lgem,*lgibc,*lgibcf,
*ik,*ik2,*ik4,
*istelp,
*ielle, *ipolyn, *imat, *leno;
int IntegerVars[64];
float     *ex,*ey,*uvp,*qbnd,*ssq,*frot,*dens,*tempt,
          *thdiss,*thdint,*vhb,*vold,*taperf,
          *qload,*exo,*eyo,*exlg,*eylg,*uxlg,*uylg,
          *stelpx,*stelpy,*uvelpx,*uvelpy,
          *exts,*eyts,*ssqts,*frotts,*uvpts,*elad ;
double   *stk,*stk2,*stk3,*stk4;
float RealVars[64];
float area, quality;

int find_maxnbs(float *ex,float *ey,int *lem,int ne,int nn);
void InitPtrs(), FreeArrays();
int AllocIntArray(int cnt, int **addr, int dflt);
int AllocFlArray(int cnt, float **addr, float dflt);
int AllocDblArray(int cnt, double **addr, double dflt);
int AllocateArrays(int *iv, int nk, int nk2, int nk4);
void arraydim_(int *,int *,int *,int *,int *,int *,
       int *,int *,int *,int *,int *,int *,int *,int *);
void elementvar_(float *,float *,float *);
void flags_(int *,int *,int *,int *,int *,int *,int *,int *,
            int *,int *,int *,int *,int *,int *);
void getdblptr_(int *cxref, double *cxarray, int *nx);
int  gethostname(char *name, size_t len);
void getintptr_(int *cxref, int *cxarray, int *nx);
void getrealptr_(int *cxref, float *cxarray, int *nx);
void initialread_(char *, int *,int *,int *,float *,char *,int *,
          char *,int *,char *,int *,char *,int *,char *,int *,
          int *, int *,int *,int *,float *,float *,float *,float *,
          int *, int *,int *,int *,int *,int *,int *,int *,int *,
          int *,int *,int *,int *,int *,int *,int *,
          int *,char *,int *,int *,int *,int *,
          float *,int *,int *, int *, char *);
void initcomments_(char *,int *);
void initdimensions_(int *,float *,
         int *,int *,int *,int *,int *,int *,int *,int *,
         int *,int *,int *,int *,int *,int *,
         float *,float *,float *,float *,float *,float *,
         int *,int *,int *,int *,int *,
         int *,int *,int *,int *,int *,int *,
         int *,int *,int *,int *,int *,int *,float *,
         int *,int *,int *,int *,int *,int *);
void initremeshdim_(int *, float *, int *,int *,int *,int *,
     int *,int *,int *,int *,int *,int *,int *,int *, int *,
     float *, float *, float *, int *,int *,int *,int *,int *);
void initstkdim_(int *IntegerVars,int *imshp, int*maxnbs,
                 int *kbw,int *k2bw,int *k4bw,int *nk,int *nk2,int *nk4);
void inputexists_(char*, int *, int *);
void locnod_(float *exts,float *eyts, float *ex,float *ey,
     int *olem,int *nor, int *leno, float *elad, int *oldne,
     int *oldnup, int *nup);
void nextfilename_(char *, int *, int *, int *, int*);
void ntrpltd_(int *iden,float *olddens,float *newdens,float *bgd,
     int *ivv,float *oldvhb,float *newvhb, float *bgv,
     int *ithdi,float *othdint,float *thdint, float *bgdi,
     float *exts,float *eyts,float *ex,float *ey,
     int *olem,int *lem,int *nor,int *omat,int *imat,
     int *oldne,int *oldnup,int *ne,int *nup,int *imreg);
void ntrplt_(float *uvpts,float *uvp,float *bgu,
     int *olem,int *nor, int *leno, float *elad, 
     int *oldne,int *oldnup, int *nup, int *iq );
int  readelleinput( char *, char *, int *, int **);
int  ReallocateSTKArrays();
void runbasil_(int *,int *,float *, char*,int *,char *,int *,
       char *,int *, char *,int *, char *,int *, char *,int *,
       int *,int *,int *,int *,int *,int *,
       int *,int *,int *,int *,int *,
       int *,int *,float *,int *,int *,int *,
       int *,int *,int *,int *,int *,int *,int *,
       float *,float *,float *,float *,float *,float *,int *,
       float *,float *,int *,float *,float *,int *,int *);
void setmat_(float *exts,float *eyts,float *ex,float *ey,
     int *olem,int *lem, int *nor, int *omat,int *mat,
     int *oldne,int *oldnup, int *ne,int *nup);
void stkarraydim_(int *,int *,int *);

int main(int argc, char **argv)
{
    char polyfile[81], polycomments[81], basilin[]="basil.in";
    char namebc[MAXNAME]="", bindir[MAXNAME]="", outdir[MAXNAME]="";
    char namedb[MAXNAME]="", namew[MAXNAME]="", namer[MAXNAME]="";
    int  lnbd, lnod, lndb, lnwf, lnr;
    int i, arginp=0, argcnt=0, nmlen=0, baslen=8;
    int infile=0, finished=0, valid=0, err=0, iread=0, filein=0, verb=0;
    int nx,ny,imshp,ivis,inflt,ivv,ipoly,iden,itemp,icr,igrav,ivold,imreg;
    int nup,nn,nmp,ne,nbp,nnnof,nmpnof,nfp=0,nfpf3=0,nseg=0,nellep=0,ncomp=0,ithdi=0;
    int kbw,k2bw,k4bw,nk,nk2,nk4;
    int nxl,nyl,npmp,nrmp,nsm,ilag,inul,nlag;
    int maxnbs=0,pbxy=0,nbpext=0;
    float xlen,ylen,pdist;
    float visdflt=1.0, se=1.0;

 /* printf("\nBasil version: %s\n",version);  */

    if (argc>1) arginp=1;
    if (arginp && (argv[1][0]=='-')) {
        if ((argv[1][1]!='v') && (argv[1][1]!='h')) {
          printf("Command line option %s not understood\n",argv[1]);
          argv[1][1] = 'h';
        }
        if (argv[1][1]=='h') {
          printf("Command line options:\n");
          printf("-v: print the Basil version and exit\n");
          printf("-h: print this message and exit\n");
          printf("FILE: a valid Basil input file\n");
          printf("basil.in: a file listing valid Basil input files,\n");
          printf("          one per line\n");
          printf("          This is the default if there are no options\n");
        }
        return(err);
    }
    /*   filein is set non-zero if basil.in exists 
              but command line takes precedence        */

    inputexists_(basilin,&baslen,&filein);

    /*   while finished is zero, look for next input file to execute  */

    while (!finished) {

    /*   get the next basil input filename from command line 
         valid is set non-zero if a filename is obtained       
         finished flag is set if no more input files provided */

        infile++;
        /*printf("Starting execution of job number %i\n",infile);*/
        nmlen=0;
        for (i=0;i<MAXNAME-1;i++) namebc[i]=' ';
        if (arginp) {
            argcnt++;
            valid=0;
            if (argcnt < argc){
              nmlen=strlen(argv[argcnt]);
              for (i=0;i<nmlen;i++) namebc[i]=argv[argcnt][i];
              inputexists_(namebc,&nmlen,&valid);
              if(!valid){
                printf("File not found: %.*s\n",nmlen,namebc);
                err = 1;
              }
            }
            else {
              finished=1;
            }
        }
      /*    if no command line arguments, use basil.in   */

        else if (filein) {
            valid=0;
            nextfilename_(namebc,&nmlen,&infile,&finished,&valid);
        }
        else {
            printf("Input filenames should be specified on command line or listed in basil.in\n");
            err=1;
            finished=1;
        }
        if (!valid || err) return(err);
      /*printf("\nProceeding with execution for input file: %.*s\n",nmlen,namebc);*/

      /*  current input file now determined            */
      /*  Initialise pointers on first file            */

        if(infile==1)InitPtrs();
        FreeArrays();
        for (i=0;i<64;i++) IntegerVars[i]=0;
        for (i=0;i<64;i++) RealVars[i]=0.0;
        for (i=0;i<81;i++) polyfile[i]='\0';
        for (i=0;i<81;i++) polycomments[i]='\0';
        nx=0; ny=0; imshp=0; inflt=0; ivv=0; ipoly=0; ivis=0;
        nxl=16; nyl=32; npmp=80; nrmp=7; nsm=0; ilag=0;
        iden = icr = ithdi = ivold = imreg = itemp = igrav = 0;
        maxnbs = nellep = 0;
        valid = err = 0, iread=0;
        xlen=1.0,ylen=1.0;
        area=6.25e-4,quality=20;
        strcpy(polycomments,"");

      /* * The default values for the variables defining the
         * FE mesh and the Lagrangian mesh may be overwritten.
         * valid is set if the input file is found
         * err indicates an error was encountered reading input */

        initialread_(namebc,&nmlen,&infile,IntegerVars,RealVars,
               bindir,&lnbd,outdir,&lnod,namedb,&lndb,namew,&lnwf,
               namer,&lnr,&nx,&ny,&imshp,&inflt,&xlen,&ylen,&area,
               &quality,&nfp,&nfpf3,&nxl,&nyl,&npmp,&nrmp,&nsm,&ilag,&igrav,
               &iden,&itemp,&icr,&ivold,&imreg,&ivis,&ncomp,
               &ithdi,polyfile,&iread,&maxnbs,&infile,&pbxy,
               &pdist,&err,&valid,&verb,version);
   //   printf("after initialread: nfp = %i nfpf3 = %i\n",nfp,nfpf3);
        nlag=npmp*nsm;

      /*   if input file valid and not showing error  */

        if (valid && !err) {
          if (iread==0) {
            if (imshp==3) {
            /*
             * if imshp==3 then triangulate
             */
              if (trimesh(&inflt,&nup,&nn,&nmp,&nbp,&nnnof,&nmpnof,
                          &ne,&nfp,&nfpf3,&nseg,&ivis,&maxnbs,&nellep,&ielle,&ipolyn,
                          &lem,&lemngh,&ibc,&ibngh,&ibctyp,&ifbc,&ifbc2,&iseg,
                          &ex,&ey,&vhb,&xlen,&ylen,area,quality,
                          &visdflt,&se,pbxy,&nbpext,pdist,
                          polyfile,polycomments)) exit(1);
     //       printf("after trimesh pbxy = %i nbpext = %i nbp = %i\n",pbxy,nbpext,nbp);
              if (ipolyn!=NULL) ipoly=1;
            } /* mesh==3 */
            else if (strlen(polyfile)>0) {
              readelleinput(polyfile,polycomments,&nellep,&ielle);
            }
            i = strlen(polycomments);
            initcomments_(polycomments,&i);
            /*
             * copy the mesh values to the arrays of scalars
             * and calculate all the necessary array dimensions
             */
            initdimensions_(IntegerVars,RealVars,
                            &nup,&nn,&nmp,&nbp,&ne,&nfp,&nfpf3,&nseg,
                            &nx,&ny,&imshp,&inflt,&maxnbs,&nellep,
                            &xlen,&ylen,ex,ey,&visdflt,&se,
                            &nxl,&nyl,&npmp,&nrmp,&nsm,
                            &ithdi,&ilag,&inul,&ivis,&ipoly,
                            &iden,&itemp,&icr,&ivold,&imreg,&ncomp,
                            &pbxy,&pdist,
                            &kbw,&k2bw,&k4bw,&nk,&nk2,&nk4);
          } /* !iread */
          else {
            initstkdim_(IntegerVars,&imshp,&maxnbs,&kbw,&k2bw,&k4bw,
                                                   &nk, &nk2, &nk4);
          }
      /*
       * allocate memory required to run the basil program
       */
          if (AllocateArrays(IntegerVars,nk,nk2,nk4)!=0) {
            err=1;
            fprintf(stderr,"Malloc failed for input file %d\n",infile);
          }
          else {
            int dum=0;
/*    fprintf(stderr,"about to call runbasil nbp = %i\n",nbp);   */
            arraydim_(IntegerVars,&nup,&nbp,&ne,&nn,&dum,&nseg,&nfp,
                         &ny,&inul,&dum,&dum,&nrmp,&npmp);

            runbasil_(&err,IntegerVars,RealVars,namebc,&nmlen,namew,&lnwf,
                      namer,&lnr,bindir,&lnbd,outdir,&lnod,namedb,&lndb,
                      ielfix,ifbc,ifbc2,ifeqv,jfbc1,jfbc2,
                      iendp1,iendp2,lgem,lgibc,lgibcf,
                      ielle,ipolyn,tempt,&nup,&nfp,&nfpf3,
                      &kbw,&k2bw,&k4bw,&nk,&nk2,&nk4,&iread,
                      exo,eyo,exlg,eylg,uxlg,uylg,&inul,
                      stelpx,stelpy,istelp,uvelpx,uvelpy,&nlag,&verb);
          /* the following parameters are now passed by c_data-to_f:
                      ssq,frot,dens,
                      vhb,thdiss,thdint,imat,vold,
                      exts,eyts,ssqts,frotts,uvpts     */
          }
        }       /*  if input file valid and no error   */
    }         /*  while !finished                    */
    return(err);
}         /*  end of main                        */

void checkproblem_()
{
}
/*
* this routine accesses the global array pointers
* This routine should only be called from a Fortran subroutine
* It does not reset the array sizes just passes them through
* When it returns, the ptrs cxref, cyref etc should point to valid data
* which can be manipulated by the Fortran code
*/
void c_data_to_f_(int *cxref, int *cyref, int *cxtsref, int *cytsref,
      int *cssqref, int *cfrotref, int *cssqtsref, int *cfrottsref,
          int *cqbndref, int *cqloadref, int *cuvpref, int *cuvptsref,
      int *cdensref, int *cvhbref, int *cvoldref, int *cthdintref,
          int *cthdissref, int *ctaperfref,
      int *cstkref, int *cstk2ref, int *cstk3ref, int *cstk4ref,
      int *csegref,int *clemref,int *cimatref, 
      int *cibcref, int *cibpoleref, int *cibctypref,
      int *cibnghref, int *cnorref,
      int *cikref, int *cik2ref, int *cik4ref,
      int *intv)
{
    int NUP,NE,NN,NBP,NROWS,NFP,NSEG,NUL,NBL,NEL,NSM,NPM;
    int NY;
    int NK,NK2,NK4;
    int size;
    /*
    * get the flag values and array dimensions from the
    * Fortran routine
    * (This allows the program to maintain only one copy of
    * the indices definitions for the arrays of scalars)
    flags_(intv,&ILAG,&IFLT,&IMSH,&IGRAV,&ICR,&ITEMP,&IDEN,
          &IVIS,&IVOLD,&IMREG,&ITHDI,&NCOMP,&NFPF3);
    */
    arraydim_(intv,&NUP,&NBP,&NE,&NN,&NROWS,&NSEG,&NFP,
              &NY,&NUL,&NBL,&NEL,&NSM,&NPM);
    stkarraydim_(&NK,&NK2,&NK4);
/* fprintf(stdout,"in c_data_to_f, NK = %i NK2 = %i NK4 = %i\n",NK,NK2,NK4); */
    /* get fortran ref to C array */
    getrealptr_(cxref,ex,&NUP);
    getrealptr_(cyref,ey,&NUP);
    getrealptr_(cxtsref,exts,&NUP);
    getrealptr_(cytsref,eyts,&NUP);
    getrealptr_(cssqref,ssq,&NUP);
    getrealptr_(cssqtsref,ssqts,&NUP);
    getrealptr_(cfrotref,frot,&NUP);
    getrealptr_(cfrottsref,frotts,&NUP);
    getrealptr_(cuvpref,uvp,&NROWS);
    getrealptr_(cuvptsref,uvpts,&NROWS);
    getrealptr_(cqloadref,qload,&NROWS);
    size = NE * 7;
    getrealptr_(cdensref,dens,&size);
    getrealptr_(cthdintref,thdint,&size);
    getrealptr_(cthdissref,thdiss,&size);
    size = NE * 8;
    getrealptr_(cvhbref,vhb,&size);
    getrealptr_(cvoldref,vold,&size);
    size = NSEG * 3;
    getintptr_(csegref,iseg,&size);
    size = NE * 6;
    getintptr_(clemref,lem,&size);
    getintptr_(cimatref,imat,&NE);
    getintptr_(cibcref,ibc,&NBP);
    getrealptr_(ctaperfref,taperf,&NBP);
    size = NBP * 2;
    getrealptr_(cqbndref,qbnd,&size);
    getintptr_(cibpoleref,ibpole,&size);
    getintptr_(cibctypref,ibctyp,&size);
    getintptr_(cibnghref,ibngh,&size);
    getintptr_(cnorref,nor,&NUP);
    /*
    getrealptr_(cstkref,stk,&NK);
    getrealptr_(cstk2ref,stk2,&NK2);
    getrealptr_(cstk3ref,stk3,&NFP);
    getrealptr_(cstk4ref,stk4,&NK4); */
    getdblptr_(cstkref,stk,&NK);
    getdblptr_(cstk2ref,stk2,&NK2);
    getdblptr_(cstk3ref,stk3,&NFP);
    getdblptr_(cstk4ref,stk4,&NK4);
    getintptr_(cikref,ik,&NK);
    getintptr_(cik2ref,ik2,&NK2);
    getintptr_(cik4ref,ik4,&NK4);
    return;
}
/*
******************************************************************
*/
/*
* this routine accesses the global pointers
* and global vars area, quality
*/
void remesh_data_(int *ibp, int * newiseg, int *newnbp, 
                   int *newnseg, int *err)
{
    int i, j;
    int nx,ny,imshp,inflt,ivv,iden,itemp,icr,igrav,ivold,imreg;
    int nrows,nul,nbl,nel,npm;
    int maxnbs;
    int nup,nn,nmp,ne,nbp,nnnof,nmpnof,nfp=0,nfpf3=0;
    int nseg=0,nellep=0,ncomp=0,ithdi=0;
    int nxl,nyl,npmp,nrmp,nsm,ilag,iq;
    int oldne,oldnup,oldnfp,oldnn,oldnrows,*olem,*omat;
    float bgd=0,bgv=0,bgu=0,bgdi=0;
    float *odens=0,*ovhb=0,*othdint=0;
    float xlen, ylen;
    float exbnd[*newnbp], eybnd[*newnbp];
    float area,quality;

    int indxibc[*newnbp];
    int num_pt_attr=8;
    float pt_attributes[num_pt_attr * (*newnbp)];
    /*
    * get the flag values and array dimensions from the
    * Fortran routine
    * (This allows the program to maintain only one copy of
    * the indices definitions for the arrays of scalars)
    */
    flags_(IntegerVars,&ilag,&inflt,&imshp,&igrav,&icr,&itemp,&iden,
          &ivv,&ivold,&imreg,&ithdi,&ncomp,&nfpf3);
    arraydim_(IntegerVars,&nup,&nbp,&ne,&nn,&nrows,&nseg,&nfp,
              &ny,&nul,&nbl,&nel,&nsm,&npm);
    elementvar_(RealVars,&area,&quality);
    
/*fprintf(stdout,"in remesh_data, area = %f quality = %f\n",area,quality);*/
    /*
    * Find the index in the ibc for the points in ibp
    for (i=0;i<*newnbp;i++) {
    */
    int count=0;
    for (i=0;i<*newnseg;i++) {
        indxibc[i]=-1;
        for (j=0;j<nbp && indxibc[i]==-1;j++) {
            if (ibc[j]==ibp[i]) {
                indxibc[i]=j;
                count++;
            }
        }
    }
/* if (count!=*newnbp) {
printf("error finding ibc index ");
} */
    /*
     * tri_remesh will free the old arrays and allocate new
     * Save the old LEM array for ntrplt
     * Save the old SSQ and FROT arrays for ntrplt
     * Save the old DENS and/or VHB arrays for ntrpltd, if necessary
     * Save the old IMAT for setmat, if necessary
     */
    olem=0;
    omat=0;
    odens=0;
    othdint=0;
    ovhb=0;
    if ((*err=AllocIntArray(6*ne,&olem,0))!=0) return;
    for (i=0;i<ne*6;i++) olem[i] = lem[i];
    if (iden>0) {
      if ((*err=AllocFlArray(7*ne,&odens,0.0))!=0) return;
      for (i=0;i<ne*7;i++) odens[i] = dens[i];
    }
    if (ithdi>0) {
      if ((*err=AllocFlArray(7*ne,&othdint,0))!=0) return;
      for (i=0;i<ne*7;i++) othdint[i] = thdint[i];
    }
    if (ivv>0) {
      if ((*err=AllocFlArray(8*ne,&ovhb,0))!=0) return;
      for (i=0;i<ne*8;i++) ovhb[i] = vhb[i];
    }
    if (imreg>0) {
      if ((*err=AllocIntArray(ne,&omat,0))!=0) return;
      for (i=0;i<ne;i++) omat[i] = imat[i];
    }
    oldne = ne;
    oldnn = nn;
    oldnup = nup;
    oldnfp = nfp;
    oldnrows = nrows;
    /*
     * Pass through qbnd, ibctyp, ibpole, taperf as pt attributes
     */
    for (i=0;i<*newnbp;i++) {
        exbnd[i] = ex[ibp[i]-1];
        eybnd[i] = ey[ibp[i]-1];
        if (indxibc[i]==-1) {
            pt_attributes[i*num_pt_attr]=0;
            pt_attributes[i*num_pt_attr+1]=0;
            pt_attributes[i*num_pt_attr+2]=0;
            pt_attributes[i*num_pt_attr+3]=0;
            pt_attributes[i*num_pt_attr+4]=0;
            pt_attributes[i*num_pt_attr+5]=0;
            pt_attributes[i*num_pt_attr+6]=0;
            pt_attributes[i*num_pt_attr+7]=0;
        }
        else {
            pt_attributes[i*num_pt_attr]=qbnd[indxibc[i]];
            pt_attributes[i*num_pt_attr+1]=qbnd[indxibc[i]+nbp];
            pt_attributes[i*num_pt_attr+2]=ibctyp[indxibc[i]];
            pt_attributes[i*num_pt_attr+3]=ibctyp[indxibc[i]+nbp];
            pt_attributes[i*num_pt_attr+4]=ibpole[indxibc[i]];
            pt_attributes[i*num_pt_attr+5]=ibpole[indxibc[i]+nbp];
            pt_attributes[i*num_pt_attr+6]=taperf[indxibc[i]];
            pt_attributes[i*num_pt_attr+7]=1;
        }
    }
    *err =  tri_remesh( &inflt,&nup,&nn,&nmp,&nbp,&nnnof,&nmpnof,
                        &ne,&nfp,&nseg,&ivv,&maxnbs,
                        &nellep,&ielle,&ipolyn,
                        &lem,&lemngh,
                        &ibc,&ibngh,&ibctyp,&ifbc,&ifbc2,&ibpole,&iseg,
                        newiseg,*newnseg,ibp,*newnbp,num_pt_attr,
                        &imat, imreg,
                        &ex,&ey,&vhb,&qbnd,&taperf,
                        exbnd, eybnd, pt_attributes,
                        &xlen, &ylen,
                        area, quality
                        );
    imshp=3;
    if (*err==0) {
        initremeshdim_(IntegerVars,RealVars,
                        &nup,&nn,&nmp,&nbp,&ne,&nfp,&nseg,
                        &nx,&ny,&imshp,&inflt,&maxnbs,&nellep,
                        &xlen,&ylen,ey,
                        &nxl,&nyl,&npmp,&nrmp,&nsm
                        );
        arraydim_(IntegerVars,&nup,&nbp,&ne,&nn,&nrows,&nseg,&nfp,
                  &ny,&nul,&nbl,&nel,&nsm,&npm);

/*   set the region identifiers before interpolating fields

        if (imreg>0) {
          if (imat!=0) free(imat); imat=0;
          if ((*err=AllocIntArray(ne,&imat,0))!=0){
            printf("Problem allocating IMAT array, check SEGSQLIM\n");
            return; 
          }
          setmat_(exts,eyts,ex,ey,olem,lem,nor,
                  omat,imat,&oldne,&oldnup,&ne,&nup);
        } */

      /*   interpolate quantities defined on element nodes -
           first locate new nodes in old mesh     */

        if (leno!=0) free(leno); leno=0;
        if ((*err=AllocIntArray(nup,&leno,0))!=0) return;
        if (elad!=0) free(elad); elad=0;
        if ((*err=AllocFlArray(3*nup,&elad,0.0))!=0) return;
        locnod_(exts,eyts,ex,ey,olem,nor,leno,elad,&oldne,&oldnup,&nup);
   
      /*   successively interpolate U, V, P (all in UVP)  */

        printf("Interpolating nodal arrays: ncomp=%i \n",ncomp);
        if (uvp!=0) free(uvp); uvp=0;
        if ((*err=AllocFlArray(nrows,&uvp,0.0))!=0) return;
        bgu=0.; iq=2;
        ntrplt_(uvpts,uvp,&bgu,olem,nor,leno,elad,&oldne,&oldnup,&nup,&iq);
        ntrplt_(uvpts+oldnup,uvp+nup,&bgu,olem,nor,leno,elad,&oldne,
                &oldnup,&nup,&iq);
        if(ncomp>0){
          iq=1;
          ntrplt_(uvpts+2*oldnup,uvp+2*nup,&bgu,olem,nor,leno,elad,&oldne,
                  &oldnn,&nn,&iq);
          iq=2;
        }

      /*   and, as required interpolate SSQ and FROT, SSQTS and FROTS
           also re-set here for use in 2-step time loop   */

        printf("Interpolating nodal arrays: icr=%i \n",icr);
        if(icr>0 ){
          if(icr<3) {
            if (ssq!=0) free(ssq); ssq=0;
            if ((*err=AllocFlArray(nup,&ssq,0.0))!=0) return;
            ntrplt_(ssqts,ssq,&bgu,olem,nor,leno,elad,&oldne,&oldnup,&nup,&iq);
            if (ssqts!=0) free(ssqts); ssqts=0;
            if ((*err=AllocFlArray(nup,&ssqts,0.0))!=0) return;
            for (i=0;i<nup;i++) ssqts[i]=ssq[i];
            printf("interpolated SSQ ok (?)\n");
          }
          if(icr!=2) {
            if (frot!=0) free(frot); frot=0;
            if (frotts!=0) free(frotts); frotts=0;
            if ((*err=AllocFlArray(nup,&frot,0.0))!=0) return;
            if ((*err=AllocFlArray(nup,&frotts,0.0))!=0) return;
            ntrplt_(frotts,frot,&bgu,olem,nor,leno,elad,&oldne,&oldnup,&nup,&iq);
            printf("interpolated FROT ok (?)\n");
          }
        }
/*      ntrplt_(uvpts,uvp,&bgu,ossq,ssq,&bgs,ofro,frot,&bgf,
                exts,eyts,ex,ey,olem,nor,
                &oldne,&oldnup,&oldnfp,&oldnn,&oldnrows,
                &nup,&nfp,&nn,&nrows,&inflt,&ncomp,&icr);  */

/*   interpolate quantities defined on interpolation points (DENS,VHB)
     suggested values for bgd, bgv, outside range of expected values */

        if (iden>0) {
          if (dens!=0) free(dens); dens=0;
          if ((*err=AllocFlArray(ne*7,&dens,0.0))!=0) return;
          bgd=1.;
        }
        if (ivv>0) {
          if (vhb!=0) free(vhb); vhb=0;
          if ((*err=AllocFlArray(ne*8,&vhb,0.0))!=0) return;
          bgv=-10;
        }
        if (ithdi>0) {
          if (thdint!=0) free(thdint); thdint=0;
          if (thdiss!=0) free(thdiss); thdiss=0;
          if ((*err=AllocFlArray(ne*7,&thdint,0.0))!=0) return;
          if ((*err=AllocFlArray(ne*7,&thdiss,0.0))!=0) return;
          bgdi=0;
        }
        printf("Interpolating dens, vhb, thdint: ivv=%i iden=%i ithdi=%i\n"
                ,ivv,iden,ithdi);
        ntrpltd_(&iden,odens,dens,&bgd,&ivv,ovhb,vhb,&bgv,
                 &ithdi,othdint,thdint,&bgdi,
                 exts,eyts,ex,ey,olem,lem,nor,omat,imat,
                 &oldne,&oldnup,&ne,&nup,&imreg);
        
      /*   free up old arrays used in re-interpolation */

        if (odens!=0) free(odens); odens=0;
        if (othdint!=0) free(othdint); othdint=0;
        if (ovhb!=0) free(ovhb); ovhb=0;
        if (omat!=0) free(omat); omat=0;
        if (olem!=0) free(olem); olem=0;
        if (nor!=0) free(nor); nor=0;

      /* reallocate and re-initialize new arrays whose size is changed */

        if ((*err=AllocIntArray(nup,&nor,0))!=0) return;
        for (i=0;i<nup;i++) nor[i]=i+1;
        if (exts!=0) free(exts); exts=0;
        if ((*err=AllocFlArray(nup,&exts,0.0))!=0) return;
        for (i=0;i<nup;i++) exts[i]=ex[i];
        if (eyts!=0) free(eyts); eyts=0;
        if ((*err=AllocFlArray(nup,&eyts,0.0))!=0) return;
        for (i=0;i<nup;i++) eyts[i]=ey[i];
        if (uvpts!=0) free(uvpts); uvpts=0;
        if ((*err=AllocFlArray(nrows,&uvpts,0.0))!=0) return;
        for (i=0;i<nrows;i++) uvpts[i]=uvp[i];
        if (qload!=0) free(qload); qload=0;
        if ((*err=AllocFlArray(nrows,&qload,0.0))!=0) return;
        if (vold!=0) free(vold); vold=0;
        if ((*err=AllocFlArray(ne*8,&vold,0.0))!=0) return;

        if ((*err= ReallocateSTKArrays())!=0) return;
        
    }
    else return;
}
/*
 * this routine accesses the global pointers
 */
int AllocateArrays(int *iv, int NK, int NK2, int NK4)
{
    int err=0;
    int ILAG,IFLT,IMSH,ITHDI,NCOMP,NFPF3;
    int IGRAV,ICR,ITEMP,IDEN,IVIS,IVOLD,IMREG;
    int NUP,NBP,NE,NN,NROWS,NFP,NY;
    int NUL,NBL,NEL,NSM,NPM;
    int NSEG;

    /*
     * get the flag values and array dimensions from the
     * Fortran routine
     * (This allows the program to maintain only one copy of
     * the indices definitions for the arrays of scalars)
     */
    flags_(iv,&ILAG,&IFLT,&IMSH,&IGRAV,&ICR,&ITEMP,&IDEN,
              &IVIS,&IVOLD,&IMREG,&ITHDI,&NCOMP,&NFPF3);
    arraydim_(iv,&NUP,&NBP,&NE,&NN,&NROWS,&NSEG,&NFP,&NY,
                 &NUL,&NBL,&NEL,&NSM,&NPM);
 // stkarraydim_(&NK,&NK2,&NK4);
 // printf("before arrays allocated: NUP = %i NK2 = %i NFP = %i\n",NUP,NK2,NFP);
    /*
     * allocate the arrays required, depending on the flags set
     */
    if ((err=AllocIntArray(6*NE,&lem,0))!=0) return(err);
    if ((err=AllocFlArray(NUP,&ex,0.0))!=0) return(err);
    if ((err=AllocFlArray(NUP,&ey,0.0))!=0) return(err);
    if ((err=AllocIntArray(NBP,&ibc,0))!=0) return(err);
    if ((err=AllocIntArray(2*NBP,&ibngh,0))!=0) return(err);
    if (NSEG>0)
      if ((err=AllocIntArray(3*NSEG,&iseg,0))!=0) return(err);
    if ((err=AllocIntArray(2*NBP,&ibctyp,-1))!=0) return(err);
    if ((err=AllocIntArray(2*NBP,&ibpole,0))!=0) return(err);
    if (IVIS!=0 && vhb==NULL) {
      if ((err=AllocFlArray(8*NE,&vhb,1.0))!=0) return(err);
    }
    if (ITHDI!=0 && thdiss==NULL) {
      if ((err=AllocFlArray(7*NE,&thdiss,0.0))!=0) return(err);
    }
    if (ITHDI!=0 && thdint==NULL) {
      if ((err=AllocFlArray(7*NE,&thdint,0.0))!=0) return(err);
    }
    if (IVOLD) {
      if ((err=AllocFlArray(8*NE,&vold,0.0))!=0) return(err);
    }
    if (IMREG) {
      if ((err=AllocIntArray(NE,&imat,0))!=0) return(err);
    }
    if ((err=AllocIntArray(NUP,&nor,0))!=0) return(err);
    if (IGRAV==4 || ICR==1 || ICR==2) {
      if ((err=AllocIntArray(NUP,&ielfix,0))!=0) return(err);
      if ((err=AllocFlArray(NUP,&ssq,0.0))!=0) return(err);
      if ((err=AllocFlArray(NUP,&ssqts,0.0))!=0) return(err);
    }
    if (ICR==1 || ICR==3) {
      if ((err=AllocFlArray(NUP,&frot,0.0))!=0) return(err);
      if ((err=AllocFlArray(NUP,&frotts,0.0))!=0) return(err);
    }
    if ((err=AllocFlArray(NROWS,&uvp,0.0))!=0) return(err);
    if ((err=AllocFlArray(NROWS,&qload,0.0))!=0) return(err);
    if ((err=AllocFlArray(2*NBP,&qbnd,0.0))!=0) return(err);
    if ((err=AllocFlArray(NBP,&taperf,0.0))!=0) return(err);
    if (IDEN) {
      if ((err=AllocFlArray(7*NE,&dens,0.0))!=0) return(err);
    }
    if (ITEMP) {
      if ((err=AllocFlArray(NUP,&tempt,1.0))!=0) return(err);
    }
    if ((err=AllocFlArray(NUP,&exts,0.0))!=0) return(err);
    if ((err=AllocFlArray(NUP,&eyts,0.0))!=0) return(err);
    if ((err=AllocFlArray(NROWS,&uvpts,0.0))!=0) return(err);

    if ((err=AllocIntArray(NK,&ik,0))!=0) return(err);
    if ((err=AllocIntArray(NK2,&ik2,0))!=0) return(err);
    if ((err=AllocIntArray(NK4,&ik4,0))!=0) return(err);
    if ((err=AllocDblArray(NK,&stk,0.0))!=0) return(err);
    if ((err=AllocDblArray(NK2,&stk2,0.0))!=0) return(err);
    if ((err=AllocDblArray(NK4,&stk4,0.0))!=0) return(err);

    if (IFLT) {
      if ((err=AllocIntArray(NFP,&ifbc,0))!=0) return(err);
      if ((err=AllocIntArray(NFP,&ifbc2,0))!=0) return(err);
      if ((err=AllocIntArray(NFP,&ifeqv,0))!=0) return(err);
      if ((err=AllocIntArray(NFP,&jfbc1,0))!=0) return(err);
      if ((err=AllocIntArray(NFP,&jfbc2,0))!=0) return(err);
      if ((err=AllocDblArray(NFP,&stk3,0.0))!=0) return(err);
      if ((err=AllocIntArray(NN,&iendp1,0))!=0) return(err);
      if ((err=AllocIntArray(NN,&iendp2,0))!=0) return(err);
      if ((err=AllocFlArray(NUP,&exo,0.0))!=0) return(err);
      if ((err=AllocFlArray(NUP,&eyo,0.0))!=0) return(err);
    }
    if (ILAG) {
      if ((err=AllocIntArray(6*NEL,&lgem,0))!=0) return(err);
      if ((err=AllocFlArray(NUL,&exlg,0.0))!=0) return(err);
      if ((err=AllocFlArray(NUL,&eylg,0.0))!=0) return(err);
      if ((err=AllocFlArray(NUL,&uxlg,0.0))!=0) return(err);
      if ((err=AllocFlArray(NUL,&uylg,0.0))!=0) return(err);
      if ((err=AllocIntArray(NBL,&lgibc,0))!=0) return(err);
      if ((err=AllocIntArray(NBL,&lgibcf,0))!=0) return(err);
      if (ILAG==1 || ILAG==2) { /* markers */
        if ((err=AllocFlArray(NSM*NPM,&stelpx,0.0))!=0) return(err);
        if ((err=AllocFlArray(NSM*NPM,&stelpy,0.0))!=0) return(err);
        if ((err=AllocFlArray(NSM*NPM,&uvelpx,0.0))!=0) return(err);
        if ((err=AllocFlArray(NSM*NPM,&uvelpy,0.0))!=0) return(err);
        if ((err=AllocIntArray(NSM*NPM,&istelp,0))!=0) return(err);
      }
    }
    return(err);
}
/*
***************************************************************************
*/
int ReallocateSTKArrays()
{
    int err=0;
    int NK,NK2,NK4;
    /*
     * get the flag values and array dimensions from the
     * Fortran routine
     * (This allows the program to maintain only one copy of
     * the indices definitions for the arrays of scalars)
     */
    stkarraydim_(&NK,&NK2,&NK4);
    if (ik!=0) free(ik); ik=0;
    if (ik2!=0) free(ik2); ik2=0;
    if (ik4!=0) free(ik4); ik4=0;
    if (stk!=0) free(stk); stk=0;
    if (stk2!=0) free(stk2); stk2=0;
    if (stk3!=0) free(stk3); stk3=0;
    if (stk4!=0) free(stk4); stk4=0;
    if ((err=AllocIntArray(NK,&ik,0.0))!=0) return(err);
    if ((err=AllocIntArray(NK2,&ik2,0.0))!=0) return(err);
    if ((err=AllocIntArray(NK4,&ik4,0.0))!=0) return(err);
    if ((err=AllocDblArray(NK,&stk,0.0))!=0) return(err);
    if ((err=AllocDblArray(NK2,&stk2,0.0))!=0) return(err);
    if ((err=AllocDblArray(NK4,&stk4,0.0))!=0) return(err);
    return(err);
}
#if XY
but this will necessitate having fortran renew its pointer values
#endif
/*
***************************************************************************
*/
int AllocIntArray(int cnt, int **addr, int dflt)
{
    /*
     * malloc memory for cnt integers and initialise the array to 0
     */
    int i;
    if (*addr==NULL) {
        if (cnt<2 || (*addr=(int *)malloc(cnt*sizeof(int)))==NULL)
                                         return(MALLOC_ERR);
        for (i=0; i<cnt; i++) (*addr)[i] = dflt;
    }
    return(0);
}
/*
***************************************************************************
*/
int AllocFlArray(int cnt, float **addr, float dflt)
{
    /*
     * malloc memory for cnt float vals and initialise the array to 0
     */
    int i;
    if (*addr==NULL) {
        if (cnt<2 || (*addr=(float *)malloc(cnt*sizeof(float)))==NULL)
                                         return(MALLOC_ERR);
        for (i=0; i<cnt; i++) (*addr)[i] = dflt;
    }
    return(0);
}

int AllocDblArray(int cnt, double **addr, double dflt)
{
    /*
     * malloc memory for cnt doubles and initialise the array to 0.0
     */
    int i;
    if (*addr==NULL) {
        if (cnt<2 || (*addr=(double *)malloc(cnt*sizeof(double)))==NULL)
                                         return(MALLOC_ERR);
        for (i=0; i<cnt; i++) (*addr)[i] = dflt;
    }
    return(0);
}

/*
 * this routine accesses the global pointers
 */
void FreeArrays()
{
    if (lem!=0) free(lem); lem=0;
    if (lemngh!=0) free(lemngh); lemngh=0;
    if (ibc!=0) free(ibc);ibc=0;
    if (ibpole!=0) free(ibpole); ibpole=0;
    if (ibngh!=0) free(ibngh); ibngh=0;
    if (ibctyp!=0) free(ibctyp); ibctyp=0;
    if (iseg!=0) free(iseg); iseg=0;
    if (nor!=0) free(nor); nor=0;
    if (ielfix!=0) free(ielfix); ielfix=0;
    if (ifbc!=0) free(ifbc); ifbc=0;
    if (ifbc2!=0) free(ifbc2); ifbc2=0;
    if (ifeqv!=0) free(ifeqv); ifeqv=0;
    if (jfbc1!=0) free(jfbc1); jfbc1=0;
    if (jfbc2!=0) free(jfbc2); jfbc2=0;
    if (iendp1!=0) free(iendp1); iendp1=0;
    if (iendp2!=0) free(iendp2); iendp2=0;
    if (lgem!=0) free(lgem); lgem=0;
    if (lgibc!=0) free(lgibc); lgibc=0;
    if (lgibcf!=0) free(lgibcf); lgibcf=0;
    if (ik!=0) free(ik); ik=0;
    if (ik2!=0) free(ik2); ik2=0;
    if (ik4!=0) free(ik4); ik4=0;
/*
    if (nomov!=0) free(nomov); nomov=0;
    if (jnoal2!=0) free(jnoal2); jnoal2=0;
    if (ielf1!=0) free(ielf1); ielf1=0;
*/
    if (istelp!=0) free(istelp); istelp=0;
    if (ielle!=0) free(ielle); ielle=0;
    if (ipolyn!=0) free(ipolyn); ipolyn=0;

    if (ex!=0) free(ex); ex=0;
    if (ey!=0) free(ey); ey=0;
    if (uvp!=0) free(uvp); uvp=0;
    if (qbnd!=0) free(qbnd); qbnd=0;
    if (taperf!=0) free(taperf); taperf=0;
    if (ssq!=0) free(ssq); ssq=0;
    if (frot!=0) free(frot); frot=0;
    if (dens!=0) free(dens); dens=0;
    if (tempt!=0) free(tempt); tempt=0;
    if (vhb!=0) free(vhb);vhb=0;
    if (thdiss!=0) free(thdiss);thdiss=0;
    if (thdint!=0) free(thdint);thdint=0;
    if (imat!=0) free(imat);imat=0;
    if (vold!=0) free(vold);vold=0;
    if (qload!=0) free(qload); qload=0;
    if (exo!=0) free(exo); exo=0;
    if (eyo!=0) free(eyo); eyo=0;
    if (exlg!=0) free(exlg); exlg=0;
    if (eylg!=0) free(eylg); eylg=0;
    if (uxlg!=0) free(uxlg); uxlg=0;
    if (uylg!=0) free(uylg); uylg=0;
    if (stelpx!=0) free(stelpx); stelpx=0;
    if (stelpy!=0) free(stelpy); stelpy=0;
    if (uvelpx!=0) free(uvelpx); uvelpx=0;
    if (uvelpy!=0) free(uvelpy); uvelpy=0;
    if (stk!=0) free(stk); stk=0;
    if (stk2!=0) free(stk2); stk2=0;
    if (stk3!=0) free(stk3); stk3=0;
    if (stk4!=0) free(stk4); stk4=0;
    if (exts!=0) free(exts); exts=0;
    if (eyts!=0) free(eyts); eyts=0;
    if (ssqts!=0) free(ssqts); ssqts=0;
    if (frotts!=0) free(frotts); frotts=0;
    if (uvpts!=0) free(uvpts); uvpts=0;

}

/*
 * this routine accesses the global pointers
 */
void InitPtrs()
{
    lem=0;
    lemngh=0;
    ibc=0;
    ibngh=0;
    ibpole=0;
    ibctyp=0;
    iseg=0;
    nor=0;
    ielfix=0;
    ifbc=0;
    ifbc2=0;
    ifeqv=0;
    jfbc1=0;
    jfbc2=0;
    iendp1=0;
    iendp2=0;
    lgem=0;
    lgibc=0;
    lgibcf=0;
    ik=0;
    ik2=0;
    ik4=0;
/*
    nomov=0;
    jnoal2=0;
    ielf1=0;
*/
    istelp=0;
    ielle=0;
    ipolyn=0;

    ex=0;
    ey=0;
    uvp=0;
    qbnd=0;
    taperf=0;
    ssq=0;
    frot=0;
    dens=0;
    tempt=0;
    vhb=0;
    thdiss=0;
    thdint=0;
    imat=0;
    vold=0;
    qload=0;
    exo=0;
    eyo=0;
    exlg=0;
    eylg=0;
    uxlg=0;
    uylg=0;
    stelpx=0;
    stelpy=0;
    uvelpx=0;
    uvelpy=0;
    stk=0;
    stk2=0;
    stk3=0;
    stk4=0;
    exts=0;
    eyts=0;
    ssqts=0;
    frotts=0;
    uvpts=0;

}

/* int find_maxnbs(int *lem,int ne,int nn)
{
    int i,j,a,max=0;
    int numnbs[nn];

    for (i = 0; i < nn; i++) numnbs[i]=0;

    for (i = 0; i < ne; i++) {
        for (j = 0; j < 3; j++) {
          a=lem[i*6+j]-1;
          numnbs[a]++;
        }
    }
    for (i = 0; i < nn; i++) if (numnbs[i]>max) max=numnbs[i];
    return(max);
} */

int find_maxnbs(float *ex,float *ey,int *lem,int ne,int nn)
//
// this routine counts the number of elements attached to each 
// vertex node, numnbp is a kind of parity check used to find
// if the loop of adjoining elements is closed or interrupted
// by a boundary or fault.  If interrupted, numnbs is incremented
// by 1 to take account of additional 2 node connections.
{
    int i,j,j1,j2,ia,ib,ic,neb=0,max=0,inod=0;
    int numnbs[nn],numnbp[nn];

    for (i = 0; i < nn; i++) {numnbs[i]=numnbp[i]=0;}

    for (i = 0; i < ne; i++) {
        for (j = 0; j < 3; j++) {
          j1=(j+1)%3;
          j2=(j+2)%3;
          ia=lem[i*6+j]-1;
          ib=lem[i*6+j1]-1;
          ic=lem[i*6+j2]-1;
          numnbs[ia]++;
          numnbp[ia]=numnbp[ia]+ib-ic;
        }
    }
    for (i = 0; i < nn; i++) {
      if (numnbp[i] != 0) {
        numnbs[i]++;
//      neb++;
      }
      if (numnbs[i]>max) {
        max=numnbs[i];
        inod=i+1;
//      printf("find_maxnbs: max = %i inod = %i ex = %f ey = %f\n",
//                          max,inod,ex[inod-1],ey[inod-1]);
      }
    }
//  printf("find_maxnbs: neb = %i\n",neb);
    return(max);
}

