#ifndef _trimesh_h
#define _trimesh_h

/*#ifdef  ANSI_DECLARATORS*/

int trimesh( int *iflt,int *nup,int *nn,int *nmp,int *nbp,int *nnnof,int *nmpnof,
             int *ne,int *nfp,int *nfpf3,int *nseg,int *ivis,int *maxnbs,
             int *nellep,int **elle_adr,int **ipolyn,
             int **lem_adr,int **lemneighbor_adr,
             int **ibc_adr,int **ibcnb_adr,int **ibctyp,int **ifbc_adr,int **ifbc2_adr,
			 int **iseg_adr,
             float **ex_adr,float **ey_adr,float **vhb_adr,
             float *xlen,float *ylen,float area,float quality,
             float *visdflt, float *se,int pbxy,int *nbpext, float pdist,
             char *infile, char *comments);
int tri_remesh( int *iflt,int *nup,int *nn,int *nmp,int *nbp,int *nnnof,int *nmpnof,
             int *ne,int *nfp,int *nfpf3,int *nseg,int *ivv,int *maxnbs,
             int *nellep,int **ielle_adr,int **ipolyn_adr,
             int **lem_adr,int **lemneighbor_adr,
             int **ibc_adr,int **ibcnb_adr,int **ibctyp_adr,int **ifbc_adr,int **ifbc2_adr,
             int **ifeqv_adr,int **jfbc1_adr,int **jfbc2_adr,int **iseg_adr,
             int *isegnew, int newnseg, int *ibp, int newnbp, int num_pt_attr,
             int **imat, int imreg,
             float **ex_adr,float **ey_adr,float **vhb_adr,float **qbnd_adr,
             float *exbnd, float *eybnd, float *pt_attributes,
             float *xlen, float *ylen,
             float area, float quality, float *visdflt,float *se,int pbxy,int *nbpext,
             int npbpt);
unsigned char is_inex_bnd_pt(int id,int *iseg, int *label, int nseg);
void imatpp_(int *ibc, int *nbpfinal);
void smatpp_(float *qbnd, int *nbp2);

#endif

