/*--------------------------------------------------------------------
 *    Basil / Sybil:   tripoly.c  1.1  13 May 1997
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tripoly.h"
#include "poly.h"
#ifndef SEEK_SET
#define SEEK_SET 0
#endif

#ifndef _STDLIB_H_
extern void *malloc();
extern void free();
#endif /* _STDLIB_H_ */
int WritePoly(struct triangulateio *out,char *name);

/*****************************************************************************

  tripoly()   Triangulate polygons

*****************************************************************************/
int tripoly( struct triangulateio *out, struct triangulateio *inp, char *infile,
             struct flagvals *flags, char *optionaltxt,
             int *maxindx, int **ielle)
{
  char str[32], buf[81];
  struct triangulateio mid, mid2;
  int err=0;
  int nbm=0,internal_bnd=0;
  void showvar();

  initio( &mid );
  initio( &mid2 );

  if ((err=readpolyinput(inp,infile,flags,optionaltxt,maxindx,ielle,
                         &nbm,&internal_bnd))
         !=0) return (err);
/* printf("Now Triangulate"); */

  /********************************/
  /******* Setup Structures *******/
  /********************************/
  /* Make necessary initializations so that Triangle can return a   */
  /*   triangulation in `mid' and a refined triangulation in 'out'. */

  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `out'.                                    */

  /**************************************/
  /******* Triangulate and Refine *******/
  /**************************************/
  /* Triangulate the points.            */
  /*   flags as defined in the triangle routine:
          -p  Triangulates a Planar Straight Line Graph (.poly file).
          -r  Refines a previously generated mesh.
          -q  Quality mesh generation.  A minimum angle may be specified[20].
          -a  Applies a maximum triangle area constraint.
          -A  Applies attributes to identify elements in certain regions.
          -c  Encloses the convex hull with segments.
          -e  Generates an edge list.
          -v  Generates a Voronoi diagram.
          -n  Generates a list of triangle neighbors.
          -g  Generates an .off file for Geomview.
          -B  Suppresses output of boundary information.
          -P  Suppresses output of .poly file.
          -N  Suppresses output of .node file.
          -E  Suppresses output of .ele file.
          -I  Suppresses mesh iteration numbers.
          -O  Ignores holes in .poly file.
          -X  Suppresses use of exact arithmetic.
          -z  Numbers all items starting from zero (rather than one).
          -o2 Generates second-order elements (6 nodes per triangle).
          -Y  Suppresses boundary segment splitting.
          -S  Specifies maximum number of added Steiner points.
          -i  Uses incremental method, rather than divide-and-conquer.
          -F  Uses Fortune's sweepline algorithm, rather than d-and-c.
          -l  Uses vertical cuts only, rather than alternating cuts.
          -s  Force segments into mesh by splitting (instead of using CDT).
          -C  Check consistency of final mesh.
          -Q  Quiet:  No terminal output except errors.
          -V  Verbose:  Detailed information on what I'm doing.
  */
  /*   Switches are chosen for the triangle routine as follows:  */
  /*   Read a PSLG (p), preserve the convex hull (c), number everything */
  /*   from  zero (z), assign a regional attribute to each element (A), */
  /*   and  produce an edge list (e), and a triangle neighbor list (n). */

#if XY 

  strcpy(str,"QpzAen");
  /*strcpy(str,"pczAen");*/
  if (flags->bndpts==0) {
    strcat(str,"Y");
  }
  printf("Triangulating 1st time, str =%s\n",str);
  triangulate(str, inp, &mid, (struct triangulateio *) NULL);

  /* Attach area constraints to the triangles in preparation for */
  /*   refining the triangulation.                               */

  /* Needed only if -r and -a switches used: */
/*
  if ((mid.trianglearealist = (REAL *)malloc(
                   mid.numberoftriangles * sizeof(REAL))==0) {
      fprintf(stderr,"Malloc failed in trimesh\n");
      return(1);
  }
  mid.trianglearealist[0] = 3.0;
  mid.trianglearealist[1] = 1.0;
*/

  /* Refine the triangulation according to the attached constraints.    */
  /*   Switches are chosen for the triangle routine as follows:         */
  /*   Read a PSLG (p), refine a previos mesh (r), quality mesh         */
  /*   generation (q), number everything from  zero (z), suppress output */
  /*   of poly file (P), produce an edge list (e), and a triangle neighbor */
  /*   list (n) and 6 nodes per element */

  strcpy(str,"QpzrAen");
  if (flags->bndpts==0) {
    strcat(str,"Y");
  }
  printf("Triangulating 2nd time, str =%s\n",str);
  triangulate(str, &mid, &mid2, (struct triangulateio *) NULL);

/*
  printf("Triangulating 3rd time\n");
  triangulate("prAeq10zn", &mid2, &mid, (struct triangulateio *) NULL);
*/

/* printf("Triangulating last time\n"); */
  strcpy(str,"QpzrAPn");
  /*strcpy(str,"przcPn");*/
  /*triangulate("przqPno2", &mid2, out, (struct triangulateio *) NULL);*/
  if (flags->area>0.0) {
    sprintf(buf,"a%f",flags->area);
    strcat(str,buf);
  }
  if (flags->quality>0) {
    sprintf(buf,"q%d",flags->quality);
    strcat(str,buf);
  }
  if (flags->bndpts==0) {
    strcat(str,"Y");
  }
  if (flags->midptnodes) strcat(str,"o2");
  triangulate(str, &mid2, out, (struct triangulateio *) NULL);
#endif

  if (nbm>0) strcpy(str,"QpzAen");
  else strcpy(str,"QpzAePn");
  /*strcpy(str,"QprzcPn");*/
  /*triangulate("przqPno2", &mid2, out, (struct triangulateio *) NULL);*/
  if (flags->area>0.0) {
    sprintf(buf,"a%f",flags->area);
    strcat(str,buf);
  }
  if (flags->quality>0) {
    sprintf(buf,"q%d",flags->quality);
    strcat(str,buf);
  }
  /*
   * if an area constraint is set but quality is not (ie q0)
   * can have problems with narrow polygons - gives large maxnbs
   * which causes malloc failure even though triangulation is
   * successful.  In this case, works better with a small q
   */
  else if (flags->area>0.0) {
    sprintf(buf,"q5");
    strcat(str,buf);
  }
  if (flags->bndpts==0) {
    strcat(str,"Y");
    if (internal_bnd) strcat(str,"Y");
  }
  if (flags->midptnodes) strcat(str,"o2");
/*printf("Triangulating now: options = %s\n",str);*/
  triangulate(str, inp, out, (struct triangulateio *) NULL);

  /************************************/
  /******* Output for debugging *******/
  /************************************/

  /*printf("Initial triangulation:\n\n");*/
  /*showvar(&mid, 1, 1, 1, 1, 1, 0);*/
  /*printf("Refined triangulation:\n\n");*/
  /*showvar(out, 0, 1, 0, 0, 0, 0);*/

  /***********************/
  /******* Cleanup *******/
  /***********************/
  /* Free all allocated arrays, including those allocated by Triangle. */

  /*cleanio( &inp );*/
  cleanio( &mid );
  cleanio( &mid2 );

  return 0;
}

int readpolyinput( struct triangulateio *inp, char *infile,
             struct flagvals *flags, char *optionaltxt,
             int *maxindx, int **ielle, int*nbm, int *internal_bnd)
{
  char str[32], buf[81];
  FILE *fp;
  long curr_pos;
  int i,j,k,m,ip,is,ir,p1,p2,bmark,istart=0;
  int np,ns,natt,ndim,nh,nr,npp;
  int rtnval=0;
  float xval, yval, attr, area;
  float tmp=0;
  float eps=1.0e-6, xtmp=0.0, ytmp=0.0;

  np=ns=natt=ndim=nh=nr=npp=0;
  str[0] = '\0';

  /**************************/
  /******* Read Input *******/
  /**************************/

  if((fp=fopen(infile,"r"))==NULL) {
      rtnval = fprintf(stderr,"cannot open file %s\n",infile);
      return(1);
  }

  /* Input File Format:
   First line:    <# of pts> <dim (=2)> <# of attr> <# of bndry markers (0 or 1)
   Next np lines: <point #> <x> <y> [attributes] [boundary marker]
   One line:      <# of segments> <# of boundary markers (0 or 1)
   Next ns lines: <segment #> <endpoint> <endpoint> [boundary marker]
   One line:      <# of holes>
   Next nh lines: <hole #> <x> <y>
   One line: <# of regional attributes and/or area constraints>
   Next nr lines:  <constraint #> <x> <y> <attrib> <max area>
   Optional line: Keyword followed by ascii data which will be interpreted
               according to the keyword
  */

  if (fgets(buf,80,fp)==0) {
      rtnval = fprintf(stderr,"Problem reading poly file \n");
      return(1);
  }
  if (sscanf(buf,"%d %d %d %d",&np, &ndim, &natt, nbm)!=4) {
      rtnval = fprintf(stderr,"Problem reading poly file %s\n",buf);
      return(1);
  }
/*  printf("%d %d %d %d",np, ndim, natt, *nbm);
    printf("sizeof REAL %d \n",sizeof(REAL));   */

  inp->numberofpoints = np;
  inp->numberofpointattributes = natt;
  if ((inp->pointlist = (REAL *) malloc(
             inp->numberofpoints * 2 * sizeof(REAL)))==0) {
      rtnval = fprintf(stderr,"Malloc failed in tripoly\n");
      return(1);
  }
  if (inp->numberofpointattributes != 0) {
    if ((inp->pointattributelist = (REAL *) malloc(inp->numberofpoints *
                                          inp->numberofpointattributes *
                                          sizeof(REAL)))==0) {
        rtnval = fprintf(stderr,"Malloc failed in trimesh\n");
        return(1);
    }
    for (i=0;i<np*natt;i++) inp->pointattributelist[i] = 0.0;
  }
  if ((inp->pointmarkerlist = (int *) malloc(
             inp->numberofpoints * sizeof(int)))==0) {
      rtnval = fprintf(stderr,"Malloc failed in trimesh\n");
      return(1);
  }

  for (i=0;i<np;i++){
     inp->pointmarkerlist[i]=0;
     if (fscanf(fp,"%d %f %f",&ip, &xval, &yval)!=3) {
         rtnval=fprintf(stderr,"Problem reading poly file point %d\n",ip);
         return(1);
     }
     if (i==0) istart=ip; /* numbering from zero or 1 in file */
     inp->pointlist[i*2] = xval;
     inp->pointlist[i*2+1] = yval;
     for (m=0;m<i;m++){
         xtmp = inp->pointlist[m*2];
         ytmp = inp->pointlist[m*2+1];
         if ((fabs(xval-xtmp)<eps) && (fabs(yval-ytmp)<eps)) {
             rtnval=fprintf(stderr,"Poly file point %d coincident with point %d\n",
                 i+istart,m+istart);
             return(1);
         }
     }
     if(natt>0){
        for (j=0;j<natt;j++){
           if (fscanf(fp,"%f",&attr)!=1) {
               rtnval=fprintf(stderr,"Problem reading point %d attribute %d\n",
                         i+istart,j);
               return(1);
           }
           inp->pointattributelist[i*natt+j] = attr;
        }
     }
     if(*nbm>0){
        if (fscanf(fp,"%d",&bmark)!=1) {
            rtnval=fprintf(stderr,"Problem reading point %d marker\n",i+istart);
            return(1);
        }
        inp->pointmarkerlist[i] = bmark;
        /*if (bmark==99) internal_bnd = 1;*/
     }
  }

  rtnval=fscanf(fp,"^\n");
  if (fscanf(fp,"\n")!=0) { // no conversions specified
       rtnval=fprintf(stderr,"Problem reading end of points\n");
       return(1);
  }
  if (fgets(buf,80,fp)==0) {
      rtnval=fprintf(stderr,"End of poly file \n");
      return(1);
  }
  if (sscanf(buf,"%d %d",&ns, nbm)!=2) {
      rtnval=fprintf(stderr,"Problem reading poly file %s\n",buf);
      return(1);
  }
/*printf("tripolyread: ns = %d nbm = %d\n",ns, *nbm); */
  inp->numberofsegments = ns;
  if ((inp->segmentlist = (int *) malloc(
                 inp->numberofsegments * 2 * sizeof(int)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in trimesh\n");
      return(1);
  }
  if ((inp->segmentmarkerlist = (int *) malloc(
                 inp->numberofsegments * sizeof(int)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in trimesh\n");
      return(1);
  }
  for (i=0;i<ns;i++){ 
     inp->segmentmarkerlist[i] = 0;
     if (fscanf(fp,"%d %d %d",&is, &p1, &p2)!=3) { 
         rtnval=fprintf(stderr,"Problem reading poly file segment %d\n",is);
         return(1);
     }
     if (i==0) istart=is; /* numbering from zero or 1 in file */
     inp->segmentlist[i*2] = p1-1;      
     inp->segmentlist[i*2+1] = p2-1;
         /* -1 in 2 lines above because point number is array address+1 */
     if(*nbm>0){
        if (fscanf(fp,"%d",&bmark)!=1) {
            rtnval=fprintf(stderr,"Problem reading poly file segment %d\n",i+istart);
            return(1);
        }
        inp->segmentmarkerlist[i] = bmark;
     }
  }

  if (fscanf(fp,"%d",&nh)!=1) {
      rtnval=fprintf(stderr,"Problem reading poly file - no. of holes \n");
      return(1);
  }
  inp->numberofholes = nh;
  if (nh > 0) {
    if ((inp->holelist = (REAL *) malloc(
                 inp->numberofholes * 2 * sizeof(REAL)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in trimesh\n");
      return(1);
    }
    for (i=0;i<nh;i++){ 
       if (fscanf(fp,"%d %f %f",&is, &xval, &yval)!=3) { 
            rtnval=fprintf(stderr,"Problem reading poly file hole %d\n",i);
            return(1);
        }
        inp->holelist[i*2] = xval;      
        inp->holelist[i*2+1] = yval;
    }
  }

  if (fscanf(fp,"%d",&nr)!=1) {
      rtnval=fprintf(stderr,"Problem reading poly file - no. of regions \n");
      return(1);
  }
  inp->numberofregions = nr;
  if (nr > 0) {
    if ((inp->regionlist = (REAL *) malloc(
                   inp->numberofregions * 4 * sizeof(REAL)))==0) {
        rtnval=fprintf(stderr,"Malloc failed in trimesh\n");
        return(1);
    }
    for (i=0;i<nr;i++){ 
        if (fscanf(fp,"%d %f %f %f %f",&ir, &xval, &yval, &attr, &area)!=5) { 
            rtnval=fprintf(stderr,"Problem reading poly file region %d\n",ir);
            return(1);
        }
        if (i==0) istart=ir; /* numbering from zero or 1 in file */
        inp->regionlist[i*4] = xval;      
        inp->regionlist[i*4+1] = yval;
        inp->regionlist[i*4+2] = attr;
        inp->regionlist[i*4+3] = area;
    }
  }
  rtnval=fscanf(fp,"^\n");
  if (fscanf(fp,"\n")!=0) { // no conversions specified
       rtnval=fprintf(stderr,"Problem reading end of regions\n");
       return(1);
  }
  /* optional keyword lines */
  while (fgets(buf,80,fp)!=0) {
    sscanf(buf,"%s",str);
    if (!strcmp(str,ELLE_KEY)) {
        if (optionaltxt) {
            strcpy(optionaltxt,buf);
            i=strlen(optionaltxt);
            if (optionaltxt[i-1]=='\n') optionaltxt[i-1] = '\0';
        }
        /* check if elle-poly numbers follow */
        curr_pos = ftell(fp);
        if (fgets(buf,80,fp)!=0) {
            if (sscanf(buf,"%d %d",&npp,&j)==2) {
                if (npp>0) {
                    *maxindx=j+1;
                    if ((*ielle = (int *)
                         malloc( (*maxindx) * sizeof(int)))==0) {
                      rtnval=fprintf(stderr,"Malloc failed in trimesh\n");
                      return(1);
                    }
                    for (i=0;i<*maxindx;i++) (*ielle)[i] = -1;
                    for (i=0;i<npp;i++){ 
                      if (fscanf(fp,"%d %d",&j, &k)!=2) {
                        rtnval=fprintf(stderr,"Problem reading elle pts\n");
                        free(*ielle);
                        *ielle=0;
                        *maxindx=0;
                        return(1);
                      }
                      if (j<*maxindx) (*ielle)[j] = k;
                      else rtnval=fprintf(stderr,"elle index out of range %d\n",j);
                    }
                }
            }
            else fseek( fp,curr_pos,SEEK_SET );
        }
    }
    else if (!strcmp(str,AREA_KEY)) {
        sscanf(buf,"%s %f",str,&tmp);
        if (tmp>0.0) flags->area=(float)tmp;
    }
    else if (!strcmp(str,QUALITY_KEY)) {
        sscanf(buf,"%s %d",str,&i);
		/* if flag field set <0 then don't use file value */
        if (flags->quality<0) flags->quality=0;
        else if (i>0) flags->quality=i;
    }
    rtnval=fscanf(fp,"^\n");
    if (fscanf(fp,"\n")!=0) { // no conversions specified
       rtnval=fprintf(stderr,"Problem reading end of regions\n");
       return(1);
    }
  }
  fclose(fp);
//return(0);
  return(rtnval);
}

int readelleinput( char *infile, char *optionaltxt, int *maxindx, int **ielle)
{
  char str[32], buf[81];
  FILE *fp;
  long curr_pos;
  int i,j,k;
  int npp;
  int rtnval=0;

  npp=0;
  str[0] = '\0';

  if((fp=fopen(infile,"r"))==NULL) {
      rtnval=fprintf(stderr,"cannot open file %s\n",infile);
      return(1);
  }
  /* optional keyword lines */
  while (fgets(buf,80,fp)!=0) {
    sscanf(buf,"%s",str);
    if (!strcmp(str,ELLE_KEY)) {
        if (optionaltxt) {
            strcpy(optionaltxt,buf);
            i=strlen(optionaltxt);
            if (optionaltxt[i-1]=='\n') optionaltxt[i-1] = '\0';
        }
        /* check if elle-poly numbers follow */
        curr_pos = ftell(fp);
        if (fgets(buf,80,fp)!=0) {
            if (sscanf(buf,"%d %d",&npp,&j)==2) {
                if (npp>0) {
                    *maxindx=j+1;
                    if ((*ielle = (int *)
                         malloc( (*maxindx) * sizeof(int)))==0) {
                      rtnval=fprintf(stderr,"Malloc failed in readelle\n");
                      return(1);
                    }
                    for (i=0;i<*maxindx;i++) (*ielle)[i] = -1;
                    for (i=0;i<npp;i++){ 
                      if (fscanf(fp,"%d %d",&j, &k)!=2) {
                        rtnval=fprintf(stderr,"Problem reading elle pts\n");
                        return(1);
                      }
                      if (j<*maxindx) (*ielle)[j] = k;
                      else rtnval=fprintf(stderr,"elle index out of range %d\n",j);
                    }
                }
            }
            else fseek( fp,curr_pos,SEEK_SET );
        }
    }
    rtnval=fscanf(fp,"^\n");
    if (fscanf(fp,"\n")!=1) {
       rtnval=fprintf(stderr,"Problem reading end of elle pts\n");
       return(1);
    }
  }
  fclose(fp);
//return(0);
  return(rtnval);
}
/*****************************************************

  tripolypts()   Triangulate polygons (from point list)

*****************************************************/

int tripolypts( struct triangulateio *out,
		struct flagvals *flags,
             int num, float *xvals, float *yvals )
{
  char str[30], buf[10];
  struct triangulateio inp;
  int np,ns,natt,nbm,ndim,nh,nr;
  int i;
  int rtnval=0;

  np=ns=natt=nbm=ndim=nh=nr=0;

  /**************************/
  /******* Setup Input *******/
  /**************************/

  initio( &inp );
  initio( out );

  inp.numberofpoints = np = num;
  if ((inp.pointlist = (REAL *) malloc(
             inp.numberofpoints * 2 * sizeof(REAL)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in tripoly\n");
      return(1);
  }
  if ((inp.pointmarkerlist = (int *) malloc(
             inp.numberofpoints * sizeof(int)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in trimesh\n");
      return(1);
  }
  for (i=0;i<np;i++){
      inp.pointlist[i*2] = xvals[i];
      inp.pointlist[i*2+1] = yvals[i];
      inp.pointmarkerlist[i] = 1;
  }

  inp.numberofsegments = ns = num;
  if ((inp.segmentlist = (int *) malloc(
                 inp.numberofsegments * 2 * sizeof(int)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in trimesh\n");
      return(1);
  }
  if ((inp.segmentmarkerlist = (int *) malloc(
                 inp.numberofsegments * sizeof(int)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in trimesh\n");
      return(1);
  }
  for (i=0;i<ns;i++){ 
      inp.segmentlist[i*2] = i;      
      inp.segmentlist[i*2+1] = (i+1)%ns;
      inp.segmentmarkerlist[i] = 1;
  }

  /*
   * no holes
   */
  /*   Switches are chosen for the triangle routine as follows:  */
  /*   Read a PSLG (p), number everything from  zero (z), */
  /*   and  produce an edge list (e), and a triangle neighbor list (n).*/

  strcpy(str,"pzAen");
  if (flags->bndpts==0) {
      strcat(str,"Y");
  }
  if (flags->area>0.0) {
      sprintf(buf,"a%f",flags->area);
      strcat(str,buf);
  }
  triangulate(str, &inp, out, (struct triangulateio *) NULL);


  /***********************/
  /******* Cleanup *******/
  /***********************/
  /* Free all allocated arrays, including those allocated by Triangle.
*/

  cleanio( &inp );

//return 0;
  return(rtnval);
}
/*****************************************************

  tripolydata()   Triangulate polygons (from point, seg data )

 *****************************************************/

int tripolydata( struct triangulateio *out,struct triangulateio *inp,
               struct flagvals *flags,
                 int numpts, int numseg, int num_pt_attr, int *segs,
                 float *xvals, float *yvals,
                 float *pt_attributes,
                 float *rgnattribs, int numrgns )

{
  char str[32], buf[81];
  int bndmark;
  int i,j,min,max;
  int rtnval=0;

  /**************************/
  /******* Setup Input *******/
  /**************************/

  initio( inp );
  initio( out );
  inp->numberofpoints = numpts;
  inp->numberofpointattributes = num_pt_attr;
  if ((inp->pointlist = (REAL *) malloc(
             inp->numberofpoints * 2 * sizeof(REAL)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in tripolydata\n");
      return(1);
  }
  if ((inp->pointmarkerlist = (int *) malloc(
             inp->numberofpoints * sizeof(int)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in tripolydata\n");
      return(1);
  }
  if ((inp->pointattributelist = (REAL *) malloc(
             inp->numberofpoints * num_pt_attr *sizeof(REAL)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in tripolydata\n");
      return(1);
  }
  for (i=0;i<numpts;i++){
        bndmark=1;
        inp->pointlist[i*2] = xvals[i];
        inp->pointlist[i*2+1] = yvals[i];
        if (inp->pointattributelist) {
          for (j=0;j<num_pt_attr;j++){
            inp->pointattributelist[num_pt_attr*i+j] =
                          pt_attributes[num_pt_attr*i+j];
          }
          bndmark = inp->pointattributelist[num_pt_attr*(i+1)-1];
        }
        inp->pointmarkerlist[i] = bndmark;
  }

    /*
     * determine if the seg endpts are numbered 0...numseg-1 or 1...numseg
     */
  for (i=0,min=numseg,max=0;i<numseg;i++){
      if (segs[i*3]<min) min=segs[i*3];
      if (segs[i*3+1]<min) min=segs[i*3+1];
      if (segs[i*3]>max) max=segs[i*3];
      if (segs[i*3+1]>max) max=segs[i*3+1];
  }
//printf("tripolydata, min and max of seg array: %i %i\n",min,max);
  inp->numberofsegments = numseg;
  if ((inp->segmentlist = (int *) malloc(
                 inp->numberofsegments * 2 * sizeof(int)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in tripolydata\n");
      return(1);
  }
  if ((inp->segmentmarkerlist = (int *) malloc(
                 inp->numberofsegments * sizeof(int)))==0) {
      rtnval=fprintf(stderr,"Malloc failed in tripolydata\n");
      return(1);
  }
  for (i=0;i<numseg;i++){
      inp->segmentlist[i*2] = segs[i*3]-min;
      inp->segmentlist[i*2+1] = segs[i*3+1]-min;
      inp->segmentmarkerlist[i] = segs[i*3+2];
  }
  if (numrgns>0) {
      inp->numberofregions=numrgns;
      if ((inp->regionlist = (REAL *) malloc(
                 inp->numberofregions * 4 * sizeof(REAL)))==0) {
        rtnval=fprintf(stderr,"Malloc failed in trimesh\n");
      return(1);
  }
  for (i=0;i<numrgns;i++){
      inp->regionlist[i*4]=rgnattribs[i*4];
      inp->regionlist[i*4+1]=rgnattribs[i*4+1];
      inp->regionlist[i*4+2]=rgnattribs[i*4+2];
      inp->regionlist[i*4+3]=rgnattribs[i*4+3];
  }
}

/*    write the current poly file to triin.inp   */
//WritePoly(inp,"triin.poly");

/*   Switches are chosen for the triangle routine as follows:  */
/*   Read a PSLG (p), number everything from  zero (z), */
/*   produce an edge list (e), and a triangle neighbor list (n).*/

/* use V in place of Q in options string to get detailed info on grid */

  strcpy(str,"QpzAen");
  if (flags->area>0.0) {
    sprintf(buf,"a%f",flags->area);
    strcat(str,buf);
  }
  if (flags->quality>0) {
    sprintf(buf,"q%d",flags->quality);
    strcat(str,buf);
  }
  /*
   * if an area constraint is set but quality is not (ie q0)
   * can have problems with narrow polygons - gives large maxnbs
   * which causes malloc failure even though triangulation is
   * successful.  In this case small q is better
   */
  else if (flags->area>0.0) {
    sprintf(buf,"q5");
    strcat(str,buf);
  }
//  Y flag specifies no adjustment of points on the boundary
  if (flags->bndpts==0) {
    strcat(str,"Y");
  }
  if (flags->midptnodes) strcat(str,"o2");
//printf("Triangulating now: options = %s\n",str);
  triangulate(str, inp, out, (struct triangulateio *) NULL);

  /*WritePoly(out,"out.poly");*/
  /*WriteBndPoly(out,"bnd.poly");*/

//return 0;
  return(rtnval);
}
/*
 * This routine writes the triangulateio in
 * poly file format - numbered 0...numberof-1
 */
int WritePoly(struct triangulateio *out,char *name)
{
    int i, j;
    int rtnval=0;
    FILE *fp;

    if ((fp=fopen(name,"w"))==0)
        return(1);
    rtnval=fprintf(fp,"%d 2 %d ", out->numberofpoints,
                                  out->numberofpointattributes);
    rtnval=fprintf(fp,"%d\n", (out->pointmarkerlist!=0?1:0));
    for (i=0;i<out->numberofpoints;i++) {
#ifdef SINGLE
        rtnval=fprintf(fp,"%d %f %f",i,out->pointlist[i*2],
#else
        rtnval=fprintf(fp,"%d %lf %lf",i,out->pointlist[i*2],
#endif
                                     out->pointlist[i*2+1]);
        for (j=0;j<out->numberofpointattributes;j++)
#ifdef SINGLE
          rtnval=fprintf(fp," %f",
#else
          rtnval=fprintf(fp," %lf",
#endif
            out->pointattributelist[i*out->numberofpointattributes+j]);
        if (out->pointmarkerlist) 
            rtnval=fprintf(fp," %d", out->pointmarkerlist[i]);
        rtnval= fprintf(fp,"\n");
    }
    if (out->numberofedges>0) {
        rtnval=fprintf(fp,"%d ",out->numberofedges);
        rtnval=fprintf(fp,"%d\n", (out->edgemarkerlist!=0?1:0));
        for (i=0;i<out->numberofedges;i++) {
            rtnval=fprintf(fp,"%d %d %d ",i,out->edgelist[i*2],
                                      out->edgelist[i*2+1]);
            if (out->edgemarkerlist)
                rtnval=fprintf(fp," %d", out->edgemarkerlist[i]);
            rtnval=fprintf(fp,"\n");
        }
    }
    else  if (out->numberofsegments>0) {
        rtnval=fprintf(fp,"%d ",out->numberofsegments);
        rtnval=fprintf(fp,"%d\n", (out->segmentmarkerlist!=0?1:0));
        for (i=0;i<out->numberofsegments;i++) {
            rtnval=fprintf(fp,"%d %d %d ",i,out->segmentlist[i*2],
                                      out->segmentlist[i*2+1]);
            if (out->segmentmarkerlist)
                rtnval=fprintf(fp," %d", out->segmentmarkerlist[i]);
            rtnval=fprintf(fp,"\n");
        }
    }
    rtnval=fprintf(fp,"0\n");
    if (out->numberofregions>0) {
        rtnval=fprintf(fp,"%d\n",out->numberofregions);
        for (i=0;i<out->numberofregions;i++) {
            rtnval=fprintf(fp,"%d %f %f %f %f\n",i,out->regionlist[i*4],
                           out->regionlist[i*4+1],out->regionlist[i*4+2],                                  out->regionlist[i*4+3]);
        }
     }
    fprintf(fp,"0\n");
    fprintf(fp,"0\n");
    fclose(fp);
//  return(0);
    return(rtnval);
}
/*
 * This routine assumes writes the triangulateio boundary
 * points and segments in
 * poly file format - numbered 0...numberof-1
 */
int WriteBndPoly(struct triangulateio *out,char *name)
{
    int i, j;
//  int newindex[out->numberofpoints];
    int rtnval=0;
    FILE *fp;

    if ((fp=fopen(name,"w"))==0)
        return(1);
    rtnval=fprintf(fp,"%d 2 %d ", out->numberofpoints,
                                  out->numberofpointattributes);
    rtnval=fprintf(fp,"%d\n", (out->pointmarkerlist!=0?1:0));
    int bndpts=0;
    for (i=0;i<out->numberofpoints;i++) {
//      newindex[i]=-1; 
//      if (out->pointmarkerlist[i]>0) newindex[i]=bndpts;
#ifdef SINGLE
        rtnval=fprintf(fp,"%d %f %f",i,out->pointlist[i*2],
#else
        rtnval=fprintf(fp,"%d %lf %lf",i,out->pointlist[i*2],
#endif
                                     out->pointlist[i*2+1]);
        for (j=0;j<out->numberofpointattributes;j++)
#ifdef SINGLE
          rtnval=fprintf(fp," %f",
#else
          rtnval=fprintf(fp," %lf",
#endif
            out->pointattributelist[i*out->numberofpointattributes+j]);
        rtnval=fprintf(fp," %d", out->pointmarkerlist[i]);
              rtnval=fprintf(fp,"\n");
        if (out->pointmarkerlist[i]>0)/* { */
          bndpts++;
    }
    if (out->numberofsegments>0) {
        rtnval=fprintf(fp,"%d ",out->numberofsegments);
        rtnval=fprintf(fp,"%d\n", (out->segmentmarkerlist!=0?1:0));
        for (i=0;i<out->numberofsegments;i++) {
            rtnval=fprintf(fp,"%d %d %d ",i,out->segmentlist[i*2],
                                      out->segmentlist[i*2+1]);
            if (out->segmentmarkerlist)
                rtnval=fprintf(fp," %d", out->segmentmarkerlist[i]);
            rtnval=fprintf(fp,"\n");
        }
    }
    rtnval=fprintf(fp,"0\n");
    rtnval=fprintf(fp,"0\n");
    fclose(fp);
 // return(0);
    return(rtnval);
}
void initio( io )
struct triangulateio *io;
{
  /* Not needed if -N switch used. */
  io->pointlist=(REAL *)NULL;
  /* Not needed if -N switch used or number of point attributes is zero: */
  io->pointattributelist=(REAL *)NULL;
  /* Not needed if -N or -B switch used. */
  io->pointmarkerlist=(int *)NULL;

  io->numberofpoints=0;
  io->numberofpointattributes=0;

  /* Not needed if -E switch used. */
  io->trianglelist=(int *)NULL;
  /* Not needed if -E switch used or number of triangle attributes is zero: */
  io->triangleattributelist=(REAL *)NULL;
  io->trianglearealist=(REAL *)NULL;

  io->neighborlist=(int *)NULL; /* Needed only if -n switch used. */
  io->numberoftriangles=0;
  io->numberofcorners=0;
  io->numberoftriangleattributes=0;

  /* Needed only if segments are output (-p or -c) and -P not used: */
  io->segmentlist=(int *)NULL;
  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  io->segmentmarkerlist=(int *)NULL;
  io->numberofsegments=0;

  io->holelist=(REAL *)NULL;
  io->numberofholes=0;
  io->regionlist=(REAL *)NULL;
  io->numberofregions=0;

  /* Needed only if -e switch used. */
  io->edgelist=(int *)NULL;
  /* Needed if -e used and -B not used. */
  io->edgemarkerlist=(int *)NULL;
  io->normlist=(REAL *)NULL;
  io->numberofedges=0;

}

void cleanio( struct triangulateio *io )
{
  if (io->pointlist!=NULL) free(io->pointlist);
  if (io->pointattributelist!=NULL) free(io->pointattributelist);
  if (io->pointmarkerlist!=NULL) free(io->pointmarkerlist);
  if (io->trianglelist!=NULL) free(io->trianglelist);
  if (io->triangleattributelist!=NULL) free(io->triangleattributelist);
  if (io->trianglearealist!=NULL) free(io->trianglearealist);
  if (io->neighborlist!=NULL) free(io->neighborlist);

  if (io->segmentlist!=NULL) free(io->segmentlist);
  if (io->segmentmarkerlist!=NULL) free(io->segmentmarkerlist);
  if (io->holelist!=NULL) free(io->holelist);
  if (io->regionlist!=NULL) free(io->regionlist);

  if (io->edgelist!=NULL) free(io->edgelist);
  if (io->edgemarkerlist!=NULL) free(io->edgemarkerlist);
  if (io->normlist!=NULL) free(io->normlist);

  initio(io);
}


/*****************************************************************************/
/*                                                                           */
/*  showvar()   Print the input or output.                                    */
/*                 for debugging                                             */
/*                                                                           */
/*****************************************************************************/

void showvar(io, markers, reporttriangles, reportneighbors, reportsegments,
            reportedges, reportnorms)

struct triangulateio *io;
int markers;
int reporttriangles;
int reportneighbors;
int reportsegments;
int reportedges;
int reportnorms;
{
  int i, j;

  for (i = 0; i < io->numberofpoints; i++) {
    printf("Point %4d:", i);
    for (j = 0; j < 2; j++) {
      printf("  %.6g", io->pointlist[i * 2 + j]);
    }
    if (io->numberofpointattributes > 0) {
      printf("   attributes");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      printf("  %.6g",
             io->pointattributelist[i * io->numberofpointattributes + j]);
    }
    if (markers) {
      printf("   marker %d\n", io->pointmarkerlist[i]);
    } else {
      printf("\n");
    }
  }
  printf("\n");

  if (reporttriangles || reportneighbors) {
    for (i = 0; i < io->numberoftriangles; i++) {
      if (reporttriangles) {
        printf("Triangle %4d points:", i);
        for (j = 0; j < io->numberofcorners; j++) {
          printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
        }
        if (io->numberoftriangleattributes > 0) {
          printf("   attributes");
        }
        for (j = 0; j < io->numberoftriangleattributes; j++) {
          printf("  %.6g", io->triangleattributelist[i *
                                         io->numberoftriangleattributes + j]);
        }
        printf("\n");
      }
      if (reportneighbors) {
        printf("Triangle %4d neighbors:", i);
        for (j = 0; j < 3; j++) {
          printf("  %4d", io->neighborlist[i * 3 + j]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportsegments) {
    for (i = 0; i < io->numberofsegments; i++) {
      printf("Segment %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->segmentlist[i * 2 + j]);
      }
      if (markers) {
        printf("   marker %d\n", io->segmentmarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportedges) {
    for (i = 0; i < io->numberofedges; i++) {
      printf("Edge %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->edgelist[i * 2 + j]);
      }
      if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
        for (j = 0; j < 2; j++) {
          printf("  %.6g", io->normlist[i * 2 + j]);
        }
      }
      if (markers) {
        printf("   marker %d\n", io->edgemarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }
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
