/*--------------------------------------------------------------------
 *    Basil / Sybil:   sybfile.c  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#ifdef SUN
#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "sybfile.h"
#include "cmndefs.h"
#include "filedat.h"
#include "data.h"
#include "errnum.h"
#include "error.h"
#include "arrays.h"
int  ReadFlArray(int, float**, FILE*, unsigned char);
int  ReadIntArray(int, int**, FILE*, unsigned char);
int  read_data_blk(FILE*, char*, int*, float*, int**, float**);
int  skip_data_blk(FILE*);
void swapco_(float*, float*, float*, int*, int*, int*);

/*
 * byte swapping between little_endian and big_endian
 * 0->3, 1->2, 2->1, 3->0
 */
#define swap_4_bytes(a) {\
 char tmp;\
 bcopy(a+3,&tmp,1);\
 bcopy(a,a+3,1);\
 bcopy(&tmp,a,1);\
 bcopy(a+1,&tmp,1);\
 bcopy(a+2,a+1,1);\
 bcopy(&tmp,a+2,1);\
 }
/*
 * alternative - seems a bit slower
#define swap_4_bytes(tmp,a,b) {\
 unsigned long a,b;\
 bcopy(tmp,(char *)(&a),4);\
b = ( ((a) << 24) | \
              (((a) << 8) & 0x00ff0000) | \
              (((a) >> 8) & 0x0000ff00) | \
            ((a) >>24) );\
 bcopy((char *)(&b),tmp,4);\
}
 */
   

#ifdef SUN
void free_mem();
int init_file_data();
#else
void free_mem( int **arraysi,float **arraysf,int maxl,int maxf );
int init_file_data( file_data **data,char *name,int max);
#endif

int file_open(file_data** data, char* name, int max, char* mode)
{
    int err=0;

    if (!(err=init_file_data( data,name,max ))) {
        if (((*data)->fp = fopen( name,mode ))==NULL) err= OPEN_ERR;
    }
    return( err );
}

#ifdef SUN
int init_file_data( data,name,max )
file_data **data;
char *name;
int max;
#else
int init_file_data( file_data **data,char *name,int max)
#endif
{

    if (*data==NULL) *data=(file_data *)malloc(sizeof(file_data));
    if (*data!=NULL) {
        strncpy((*data)->fname,name,max);
        (*data)->fp = NULL;
        (*data)->rec_req = 1;
        (*data)->rec_curr = 0;
        (*data)->ref_req = 1;
        (*data)->ref_curr = 0;
        (*data)->rec_max = 0;
        return(0);
    }
    return(MALLOC_ERR);
}

int read_data(FILE *fp, char *string_vars, int *int_vars,
              float *fl_vars, int **arraysi, float **arraysf,
              float *xmin, float *xmax, float *ymin, float *ymax,
              int *current, int req, int *maxrec, int rotate )
{
    int i,err=0,rec_num=0,rotate_velocity=1;
    long curr_pos;
    struct stat buf;

    curr_pos = ftell( fp );

    if (*current>req) {
        rewind( fp );
        *current=0;
    }
    rec_num = *current;
    while (rec_num<req-1) {
        if ((err=skip_data_blk( fp ))) {
            if (err==EOF_ERR) *maxrec = rec_num;
            return( err );
        }
        else rec_num++;
    }
    if ((err=read_data_blk( fp,string_vars, int_vars,fl_vars,arraysi,arraysf))) {
        if (err==READ_ERR) {
            clear_arrays();
            rewind( fp );
            *current=0;
        }
        else if (err==EOF_ERR) {
            *maxrec = rec_num;
            fseek( fp,curr_pos,SEEK_SET );
        }
        return( err );
    }
    /*if (req==1) {*/
    if (*maxrec==0) {
        curr_pos = ftell( fp );
        err = fstat(fileno(fp),&buf);
        *maxrec = (float)buf.st_size/curr_pos + 0.5;
        *maxrec=0;
        while ((err=skip_data_blk( fp ))==0) {
            rec_num++;
        }
        if (err==EOF_ERR) {
            *maxrec=rec_num+1;
            fseek( fp,curr_pos,SEEK_SET );
            err=0;
        }
        else return( err );
        /*
         * if all records are the same length (ie IGRID!=4)
        curr_pos = ftell( fp );
        err = fstat(fileno(fp),&buf);
        *maxrec = (float)buf.st_size/curr_pos + 0.5;
         */
    }
    *current = req;
    if (rotate) {
         printf("swapco_ is called, rotate = %i\n",rotate);
         swapco_(arraysf[EX],arraysf[EY],arraysf[UVP],
                                &int_vars[NUP],&rotate,&rotate_velocity );
         if (int_vars[ILAG]) {
              rotate_velocity = 0;
              if (int_vars[ILAG]==1||int_vars[ILAG]==3)
                  swapco_(arraysf[EXLG],arraysf[EYLG],NULL,
                                &int_vars[NUL],&rotate,&rotate_velocity );
              i = int_vars[NSM] * int_vars[NPM];
              if (int_vars[ILAG]==1||int_vars[ILAG]==2)
                  swapco_(arraysf[STELPX],arraysf[STELPY],NULL,
                                &i,&rotate,&rotate_velocity );
         }
    }
    find_max_min(arraysf[EX],arraysf[EY],int_vars[NN],xmin,xmax,ymin,ymax);
//  printf("read_data: xmin, xmax = %f %f ymin, ymax = %f %f\n",*xmin,*xmax,*ymin,*ymax);
    return( 0 );
}

#define FUDGE_64BIT 2

int read_data_blk(FILE* fp, char* string_vars, int* int_vars, float* fl_vars,
                  int** arraysi, float** arraysf)
{
    register unsigned char byte_swap;
    register unsigned char array_resize;
    char *ptr;
    int num,i,err;
    int cnt;
    int prev_int_vars[REC_2_CNT];
    unsigned char longrecmark=0;
    long curr_pos;

    curr_pos = ftell( fp );

    byte_swap = 0;
    array_resize = 0;
        /* read the string vars */
    if ((num = fread( &cnt, 1, sizeof(int), fp ))!=sizeof(int)) {
        if (feof(fp)) return(EOF_ERR);
        else return(READ_ERR);
    }
    if (cnt!=REC_1_CNT*sizeof(char)) {
        /*
         * wrong sex or incompatible solution format
         */
        swap_4_bytes((char *)(&cnt));
        if (cnt!=REC_1_CNT*sizeof(char)) {
        /*
         * extra bytes at end of last record?
         * better to keep reading while num==1&&cnt==0 ??
         */
            for (i=0;i<FUDGE_64BIT;i++) num=fread( &cnt, 1, sizeof(int), fp );
            if (cnt!=REC_1_CNT*sizeof(char)) {
                swap_4_bytes((char *)(&cnt));
                byte_swap = 1;
            }
            if (cnt!=REC_1_CNT*sizeof(char)) {
                if (feof(fp)) return(EOF_ERR);
                else return(READ_ERR);
            }
        }
        else byte_swap = 1;
    }
    num = fread( string_vars, 1, REC_1_CNT*sizeof(char), fp );
    num = fread( &cnt, 1, sizeof(int), fp );
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (cnt!=REC_1_CNT*sizeof(char)) {
        longrecmark = 1;
        /*
         * read the first record again
         */
        fseek( fp,curr_pos,SEEK_SET );
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (cnt!=REC_1_CNT*sizeof(char)) {
        /*
         * extra bytes at end of last record?
         * better to keep reading while num==1&&cnt==0 ??
         */
            for (i=0;i<FUDGE_64BIT;i++) num=fread( &cnt, 1, sizeof(int), fp );
            if (byte_swap) swap_4_bytes((char *)(&cnt));
        }
        num = fread( string_vars, 1, REC_1_CNT*sizeof(char), fp );
        num = fread( &cnt, 1, sizeof(int), fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

        /* read the integer vars */
    for (i=0;i<REC_2_CNT;i++) prev_int_vars[i] = int_vars[i];
    num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (cnt!=REC_2_CNT*sizeof(int)) return(READ_ERR);
    num = fread( int_vars, 1, REC_2_CNT*sizeof(int), fp );
    if (byte_swap) {
        /* let the user know byte-swapping is necessary */
        fprintf(stdout,"The solution requires byte swapping\n");
        for (i=0,ptr=(char *)int_vars;i<REC_2_CNT;ptr+=4,i++)
            swap_4_bytes(ptr);
    }
    for (i=0;i<REC_2_CNT && !array_resize;i++) { 
        if (prev_int_vars[i] != int_vars[i]) array_resize=1;
        clear_alloc_arrays();
    }
    num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

        /* read the float vars */
    num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (cnt!=REC_3_CNT*sizeof(float)) return(READ_ERR);
    num = fread( fl_vars, 1, REC_3_CNT*sizeof(float), fp );
    if (byte_swap)
        for (i=0,ptr=(char *)fl_vars;i<REC_3_CNT;ptr+=4,i++)
            swap_4_bytes(ptr);
    num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

    /*
     * (b) mandatory mesh / solution / boundary condition blocks
     */
    num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (cnt!=(6*int_vars[NE]+int_vars[NUP]+int_vars[NBP]+
                int_vars[NBP2]+int_vars[NBP2])*sizeof(int))
        return(READ_ERR);
    if ((err = ReadIntArray(6*int_vars[NE],&arraysi[LEM],fp,byte_swap)))
                                                        return(err);
    if ((err = ReadIntArray(int_vars[NUP],&arraysi[NOR],fp,byte_swap)))
                                                        return(err);
    if ((err = ReadIntArray(int_vars[NBP],&arraysi[IBC],fp,byte_swap)))
                                                        return(err);
    if ((err = ReadIntArray(int_vars[NBP2],&arraysi[IBNGH],fp,byte_swap)))
                                                        return(err);
    if ((err = ReadIntArray(int_vars[NBP2],&arraysi[IBCTYP],fp,byte_swap)))
                                                        return(err);
    num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

    num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (cnt!=int_vars[NUP]*sizeof(float)*2+(int_vars[NUVP]+int_vars[NBP2])
                                            *sizeof(int))
        return(READ_ERR);
    if ((err = ReadFlArray(int_vars[NUP],&arraysf[EX],fp,byte_swap)))
                                                        return(err);
    if ((err = ReadFlArray(int_vars[NUP],&arraysf[EY],fp,byte_swap)))
                                                        return(err);
    if ((err = ReadFlArray(int_vars[NUVP],&arraysf[UVP],fp,byte_swap)))
                                                        return(err);
    if ((err = ReadFlArray(int_vars[NBP2],&arraysf[QBND],fp,byte_swap)))
                                                        return(err);
    num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

    /*
     * (c) optional crustal thickness / rotation block
     */
    if (int_vars[ICR]) {
        num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (int_vars[ICR]==1) {
          if (cnt!=int_vars[NUP]*sizeof(int)) return(READ_ERR);
          if ((err = ReadIntArray(int_vars[NUP],&arraysi[IELFIX],
                                          fp,byte_swap))) return(err);
          num = fread( &cnt, 1, sizeof(int), fp );
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

          num = fread( &cnt, 1, sizeof(int), fp );
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
          if (byte_swap) swap_4_bytes((char *)(&cnt));
          if (cnt!=int_vars[NUP]*2*sizeof(float)) return(READ_ERR);
          if ((err = ReadFlArray(int_vars[NUP],&arraysf[SSQ],
                                      fp,byte_swap))) return(err);
          if ((err = ReadFlArray(int_vars[NUP],&arraysf[FROT],
                                      fp,byte_swap))) return(err);
        }
        else if (int_vars[ICR]==2) {
          if (cnt!=int_vars[NUP]*sizeof(float)) return(READ_ERR);
          if ((err = ReadIntArray(int_vars[NUP],&arraysi[IELFIX],
                                         fp,byte_swap))) return(err);
          num = fread( &cnt, 1, sizeof(int), fp );
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

          num = fread( &cnt, 1, sizeof(int), fp );
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
          if (byte_swap) swap_4_bytes((char *)(&cnt));
          if (cnt!=int_vars[NUP]*sizeof(float)) return(READ_ERR);
          if ((err = ReadFlArray(int_vars[NUP],&arraysf[SSQ],
                                      fp,byte_swap))) return(err);
        }
        else if (int_vars[ICR]==3) {
          if (cnt!=int_vars[NUP]*sizeof(float)) return(READ_ERR);
          if ((err = ReadFlArray(int_vars[NUP],&arraysf[FROT],
                                      fp,byte_swap))) return(err);
        }
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (d) optional viscosity arrays
     */
    if (int_vars[IVIS]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (cnt!=int_vars[NE]*8*sizeof(float)) return(READ_ERR);
        if ((err = ReadFlArray(8*int_vars[NE],&arraysf[VHB],
                                      fp,byte_swap))) return(err);
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (e) optional density, temperature distribution arrays
     */
    if (int_vars[IDEN]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (cnt!=int_vars[NE]*7*sizeof(float)) return(READ_ERR);
        cnt = int_vars[NE]*7;
        if ((err = ReadFlArray(7*int_vars[NE],&arraysf[DENS],
                                      fp,byte_swap))) return(err);
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }
    if (int_vars[ITEMP]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (cnt!=int_vars[NUP]*sizeof(float)) return(READ_ERR);
        cnt = int_vars[NUP];
        if ((err = ReadFlArray(int_vars[NUP],&arraysf[TEMPT],
                                      fp,byte_swap))) return(err);
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (f) optional Lagrangian arrays
     */
    if (int_vars[ILAG]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
            /*
             * Read mesh arrays
             */
        if (int_vars[ILAG]==1||int_vars[ILAG]==3) {
          if (cnt!=(int_vars[NEL]*6+int_vars[NBL]*2)*sizeof(int))
            return(READ_ERR);
          if ((err = ReadIntArray(6*int_vars[NEL],&arraysi[LGEM],
                                         fp,byte_swap))) return(err);
          if ((err = ReadIntArray(int_vars[NBL],&arraysi[LGIBC],
                                         fp,byte_swap))) return(err);
          if ((err = ReadIntArray(int_vars[NBL],&arraysi[LGIBCF],
                                         fp,byte_swap))) return(err);
          num = fread( &cnt, 1, sizeof(int), fp );
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
                /* next record */
          num = fread( &cnt, 1, sizeof(int), fp );
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
          if (byte_swap) swap_4_bytes((char *)(&cnt));
          if (int_vars[ILAG]==1) {
            if (cnt!=(int_vars[NUL]*4+int_vars[NPM]*int_vars[NSM]*2)
                                     *sizeof(float)) return(READ_ERR);
          }
          else
            if (cnt!=int_vars[NUL]*4 *sizeof(float)) return(READ_ERR);
          if ((err = ReadFlArray(int_vars[NUL],&arraysf[EXLG],
                                      fp,byte_swap))) return(err);
          if ((err = ReadFlArray(int_vars[NUL],&arraysf[EYLG],
                                      fp,byte_swap))) return(err);
          if ((err = ReadFlArray(int_vars[NUL],&arraysf[UXLG],
                                      fp,byte_swap))) return(err);
          if ((err = ReadFlArray(int_vars[NUL],&arraysf[UYLG],
                                      fp,byte_swap))) return(err);
        }
            /*
             * Read marker arrays
             */
        if (int_vars[ILAG]==1||int_vars[ILAG]==2) {
          if (int_vars[ILAG]==2)
            if (cnt!=(int_vars[NPM]*int_vars[NSM]*2)
                                     *sizeof(float)) return(READ_ERR);
          cnt = int_vars[NSM]*int_vars[NPM];
          if ((err = ReadFlArray(cnt,&arraysf[STELPX],
                                      fp,byte_swap))) return(err);
          if ((err = ReadFlArray(cnt,&arraysf[STELPY],
                                      fp,byte_swap))) return(err);
        }
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (g) optional fault arrays block
     */
    if (int_vars[IFLT]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (cnt!=int_vars[NFP]*5*sizeof(int)) return(READ_ERR);
        if ((err = ReadIntArray(int_vars[NFP],&arraysi[IFBC],
                                         fp,byte_swap))) return(err);
        if ((err = ReadIntArray(int_vars[NFP],&arraysi[IFBC2],
                                         fp,byte_swap))) return(err);
        if ((err = ReadIntArray(int_vars[NFP],&arraysi[IFEQV],
                                         fp,byte_swap))) return(err);
        if ((err = ReadIntArray(int_vars[NFP],&arraysi[JFBC1],
                                         fp,byte_swap))) return(err);
        if ((err = ReadIntArray(int_vars[NFP],&arraysi[JFBC2],
                                         fp,byte_swap))) return(err);
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (h) optional series arrays block
     */
    if (int_vars[MSINDX]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (cnt!=int_vars[MSINDX]*2*sizeof(int)) return(READ_ERR);
        fseek(fp,cnt+sizeof(int),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (i) optional elle arrays block
     */
    if (int_vars[NELLEP] > 3) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (cnt!=int_vars[NELLEP]*sizeof(int)) return(READ_ERR);
        if ((err = ReadIntArray(int_vars[NELLEP],&arraysi[ELLENODE],
                                         fp,byte_swap))) return(err);
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (j) optional polygon numbers for triangles block
     */
    if (int_vars[IPOLY]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (cnt!=int_vars[NE]*sizeof(int)) return(READ_ERR);
        if ((err = ReadIntArray(int_vars[NE],&arraysi[POLYN],
                                         fp,byte_swap))) return(err);
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }
    /*
     * (k) optional starting viscosity arrays
     */

    if (int_vars[IVOLD]) {

      num = fread( &cnt, 1, sizeof(int), fp );
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      if (byte_swap) swap_4_bytes((char *)(&cnt));
      if (cnt!=int_vars[NE]*8*sizeof(float)) return(READ_ERR);
      if ((err = ReadFlArray(8*int_vars[NE],&arraysf[VOLD],
                fp,byte_swap))) return(err);
      num = fread( &cnt, 1, sizeof(int), fp );
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

    }

    /*
     * (l) optional material property arrays
     */
    if (int_vars[IMREG]) {
      num = fread( &cnt, 1, sizeof(int), fp );
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      if (byte_swap) swap_4_bytes((char *)(&cnt));
      if (cnt!=int_vars[NE]*sizeof(int)) return(READ_ERR);
      if ((err = ReadIntArray(int_vars[NE],&arraysi[IMAT],
                 fp,byte_swap))) return(err);
      num = fread( &cnt, 1, sizeof(int), fp );
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (m) optional segment array
     */
    if (int_vars[NSEG]>0) {
      num = fread( &cnt, 1, sizeof(int), fp );
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      if (byte_swap) swap_4_bytes((char *)(&cnt));
      if (cnt!=int_vars[NSEG]*3*sizeof(int)) return(READ_ERR);
      if ((err = ReadIntArray(int_vars[NSEG]*3,&arraysi[ISEG],
                 fp,byte_swap))) return(err);
      num = fread( &cnt, 1, sizeof(int), fp );
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (n) optional ibpole/taperf arrays (ibpole/taperf not utilised in sybil)
     */

    if (int_vars[IVRESET]>0) {

      num = fread( &cnt, 1, sizeof(int), fp );
      if (byte_swap) swap_4_bytes((char *)(&cnt));
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      if (cnt!=int_vars[NBP]*2*sizeof(int)) return(READ_ERR);
/*    if (err = ReadIntArray(cnt,&arraysi[IBPOLE],fp,byte_swap)) return(err); */
      fseek(fp,cnt+sizeof(int),SEEK_CUR);
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

      num = fread( &cnt, sizeof(int), 1, fp );
      if (byte_swap) swap_4_bytes((char *)(&cnt));
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      if (cnt!=int_vars[NBP]*sizeof(float)) return(READ_ERR);
/*    if (err = ReadFlArray(cnt,&arraysf[TAPERF],fp,byte_swap)) return(err); */
      fseek(fp,cnt+sizeof(float),SEEK_CUR);
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

      if (int_vars[IPOLE]>0) {
        num = fread( &cnt, sizeof(int), 1, fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (cnt!=int_vars[IPOLE]*3*sizeof(float)) return(READ_ERR);
/*    if (err = ReadFlArray(cnt,&arraysf[POLEP],fp,byte_swap)) return(err); */
        fseek(fp,cnt+sizeof(float),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      }
    }

    return( 0 );
}

int skip_data_blk(FILE* fp)
{
    register unsigned char byte_swap;
    char *ptr;
    int num,i;
    int cnt,tmp_vars[REC_2_CNT];
    unsigned char longrecmark=0;
    long curr_pos;

    curr_pos = ftell( fp );

    byte_swap = 0;
        /* skip the string vars */
    if ((num = fread( &cnt, 1, sizeof(int), fp ))!=sizeof(int)) {
        if (feof(fp)) return(EOF_ERR);
        else return(READ_ERR);
    }
    if (cnt!=REC_1_CNT*sizeof(char)) {
        /*
         * wrong sex or incompatible solution format
         */
        swap_4_bytes((char *)(&cnt));
        if (cnt!=REC_1_CNT*sizeof(char)) {
        /*
         * extra bytes at end of last record?
         */
            for (i=0;i<FUDGE_64BIT;i++) num=fread( &cnt, 1, sizeof(int), fp );
            if (cnt!=REC_1_CNT*sizeof(char)) {
                swap_4_bytes((char *)(&cnt));
                byte_swap = 1;
            }
            if (cnt!=REC_1_CNT*sizeof(char)) {
                if (feof(fp)) return(EOF_ERR);
                else return(READ_ERR);
            }
        }
        else byte_swap = 1;
    }
    fseek(fp,cnt,SEEK_CUR);
    num = fread( &cnt, 1, sizeof(int), fp );
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (cnt!=REC_1_CNT*sizeof(char)) {
        longrecmark = 1;
        /*
         * read the first record again
         */
        fseek( fp,curr_pos,SEEK_SET );
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (cnt!=REC_1_CNT*sizeof(char)) {
        /*
         * extra bytes at end of last record?
         * better to keep reading while num==1&&cnt==0 ??
         */
            for (i=0;i<FUDGE_64BIT;i++) num=fread( &cnt, 1, sizeof(int), fp );
            if (byte_swap) swap_4_bytes((char *)(&cnt));
        }
        fseek(fp,cnt,SEEK_CUR);
        num = fread( &cnt, 1, sizeof(int), fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

        /* read the integer vars */
    num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (cnt!=REC_2_CNT*sizeof(int)) return(READ_ERR);
    num = fread( tmp_vars, 1, REC_2_CNT*sizeof(int), fp );
    if (byte_swap)
        for (i=0,ptr=(char *)tmp_vars;i<REC_2_CNT;ptr+=4,i++)
            swap_4_bytes(ptr);
    num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

        /* skip the float vars */
    num = fread( &cnt, 1, sizeof(int), fp );
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (cnt!=REC_3_CNT*sizeof(float)) return(READ_ERR);
    fseek(fp,cnt+sizeof(int),SEEK_CUR);
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

    /*
     * (b) skip mandatory mesh / solution / boundary condition blocks
     */
    num = fread( &cnt, 1, sizeof(int), fp );
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (cnt!=(6*tmp_vars[NE]+tmp_vars[NUP]+tmp_vars[NBP]+tmp_vars[NBP2]
                +tmp_vars[NBP2])*sizeof(int))
        return(READ_ERR);
    fseek(fp,cnt+sizeof(int),SEEK_CUR);
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

    num = fread( &cnt, 1, sizeof(int), fp );
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (cnt!=tmp_vars[NUP]*sizeof(float)*2+(tmp_vars[NUVP]+tmp_vars[NBP2])
                                            *sizeof(int))
        return(READ_ERR);
    fseek(fp,cnt+sizeof(int),SEEK_CUR);
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);


    /*
     * (c) optional crustal thickness / rotation block
     */
    if (tmp_vars[ICR]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (cnt!=tmp_vars[NUP]*sizeof(int)) return(READ_ERR);
        fseek(fp,cnt+sizeof(int),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

        if (tmp_vars[ICR]==1) {
            num = fread( &cnt, 1, sizeof(int), fp );
            if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
            if (byte_swap) swap_4_bytes((char *)(&cnt));
            if (cnt!=tmp_vars[NUP]*2*sizeof(float)) return(READ_ERR);
            fseek(fp,cnt+sizeof(int),SEEK_CUR);
            if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        }
        else if (tmp_vars[ICR]==2) { 
            num = fread( &cnt, 1, sizeof(int), fp );
            if (byte_swap) swap_4_bytes((char *)(&cnt));
            if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
            if (cnt!=tmp_vars[NUP]*sizeof(float)) return(READ_ERR);
            fseek(fp,cnt+sizeof(int),SEEK_CUR);
            if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        }
    }

    /*
     * (d) optional viscosity arrays
     */
    if (tmp_vars[IVIS]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (cnt!=tmp_vars[NE]*8*sizeof(float)) return(READ_ERR);
        fseek(fp,cnt+sizeof(int),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (e) optional density, temperature distribution arrays
     */
    if (tmp_vars[IDEN]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (cnt!=tmp_vars[NE]*7*sizeof(float)) return(READ_ERR);
        fseek(fp,cnt+sizeof(int),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }
    if (tmp_vars[ITEMP]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (cnt!=tmp_vars[NUP]*sizeof(float)) return(READ_ERR);
        fseek(fp,cnt+sizeof(int),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (f) optional Lagrangian arrays
     */
    if (tmp_vars[ILAG]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (tmp_vars[ILAG]==1||tmp_vars[ILAG]==3) {
            if (cnt!=(tmp_vars[NEL]*6+tmp_vars[NBL]*2)*sizeof(int))
                return(READ_ERR);
            fseek(fp,cnt+sizeof(int),SEEK_CUR);
            if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
            num = fread( &cnt, 1, sizeof(int), fp );
            if (byte_swap) swap_4_bytes((char *)(&cnt));
            if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        }
        if (tmp_vars[ILAG]==1) {
            if (cnt!=(tmp_vars[NUL]*4+tmp_vars[NPM]*tmp_vars[NSM]*2)
                                     *sizeof(float)) return(READ_ERR);
        }
        else if (tmp_vars[ILAG]==2) {
            if (cnt!=(tmp_vars[NPM]*tmp_vars[NSM]*2)
                                     *sizeof(float)) return(READ_ERR);
        }
        else if (tmp_vars[ILAG]==3) {
            if (cnt!=tmp_vars[NUL]*4*sizeof(float)) return(READ_ERR);
        }
        fseek(fp,cnt+sizeof(int),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (g) optional fault arrays block
     */
    if (tmp_vars[IFLT]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (cnt!=tmp_vars[NFP]*5*sizeof(int)) return(READ_ERR);
        fseek(fp,cnt+sizeof(int),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (h) optional series arrays block
     */
    if (tmp_vars[MSINDX]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (cnt!=tmp_vars[MSINDX]*2*sizeof(int)) return(READ_ERR);
        fseek(fp,cnt+sizeof(int),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (i) optional elle arrays block
     */
    if (tmp_vars[NELLEP] > 3) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (cnt!=tmp_vars[NELLEP]*sizeof(int)) return(READ_ERR);
        fseek(fp,cnt+sizeof(int),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

    /*
     * (j) optional polygon numbers for triangles block
     */
    if (tmp_vars[IPOLY]) {
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (cnt!=tmp_vars[NE]*sizeof(int)) return(READ_ERR);
        fseek(fp,cnt+sizeof(int),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }
    /*
     * (k) optional starting viscosity arrays
     */

    if (tmp_vars[IVOLD]) {
      num = fread( &cnt, 1, sizeof(int), fp );
      if (byte_swap) swap_4_bytes((char *)(&cnt));
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      if (cnt!=tmp_vars[NE]*8*sizeof(float)) return(READ_ERR);
      fseek(fp,cnt+sizeof(int),SEEK_CUR);
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }
    /*
     * (l) optional material property arrays
     */

    if (tmp_vars[IMREG]) {
      num = fread( &cnt, 1, sizeof(int), fp );
      if (byte_swap) swap_4_bytes((char *)(&cnt));
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      if (cnt!=tmp_vars[NE]*sizeof(int)) return(READ_ERR);
      fseek(fp,cnt+sizeof(int),SEEK_CUR);
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }
    /*
     * (m) optional segment array
     */

    if (tmp_vars[NSEG]>0) {
      num = fread( &cnt, 1, sizeof(int), fp );
      if (byte_swap) swap_4_bytes((char *)(&cnt));
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      if (cnt!=tmp_vars[NSEG]*3*sizeof(int)) return(READ_ERR);
      fseek(fp,cnt+sizeof(int),SEEK_CUR);
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }
    /*
     * (n) optional ibpole/taperf arrays (ibpole/taperf not utilised in sybil)
     */

    if (tmp_vars[IVRESET]>0) {

      num = fread( &cnt, 1, sizeof(int), fp );
      if (byte_swap) swap_4_bytes((char *)(&cnt));
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      if (cnt!=tmp_vars[NBP]*2*sizeof(int)) return(READ_ERR);
/*    if (err = ReadIntArray(cnt,&arraysi[IBPOLE],fp,byte_swap)) return(err); */
      fseek(fp,cnt+sizeof(int),SEEK_CUR);
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

      num = fread( &cnt, sizeof(int), 1, fp );
      if (byte_swap) swap_4_bytes((char *)(&cnt));
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      if (cnt!=tmp_vars[NBP]*sizeof(float)) return(READ_ERR);
/*    if (err = ReadFlArray(cnt,&arraysf[TAPERF],fp,byte_swap)) return(err); */
      fseek(fp,cnt+sizeof(float),SEEK_CUR);
      if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

      if (tmp_vars[IPOLE]>0) {
        num = fread( &cnt, sizeof(int), 1, fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (cnt!=tmp_vars[IPOLE]*3*sizeof(float)) return(READ_ERR);
/*    if (err = ReadFlArray(cnt,&arraysf[POLEP],fp,byte_swap)) return(err); */
        fseek(fp,cnt+sizeof(float),SEEK_CUR);
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      }
    }



    return( 0 );
}

int read_ref_blk( fp,arraysf,nup )
FILE *fp;
float **arraysf;
int nup;
{
    register unsigned char byte_swap;
    char *ptr;
    int num,i,err;
    int cnt,tmp_vars[REC_2_CNT];
    unsigned char longrecmark=0;
    long curr_pos;

    curr_pos = ftell( fp );

    byte_swap = 0;
        /* skip the string vars */
    num = fread( &cnt, 1, sizeof(int), fp );
    if (cnt!=REC_1_CNT*sizeof(char)) {
        /*
         * wrong sex or incompatible solution format
         */
        swap_4_bytes((char *)(&cnt));
        if (cnt!=REC_1_CNT*sizeof(char)) {
        /*
         * extra bytes at end of last record?
         */
            for (i=0;i<FUDGE_64BIT;i++) num=fread( &cnt, 1, sizeof(int), fp );
            if (cnt!=REC_1_CNT*sizeof(char)) {
                swap_4_bytes((char *)(&cnt));
                byte_swap = 1;
            }
            if (cnt!=REC_1_CNT*sizeof(char)) {
                if (feof(fp)) return(EOF_ERR);
                else return(READ_ERR);
            }
        }
        else byte_swap = 1;
    }
    fseek(fp,cnt,SEEK_CUR);
    num = fread( &cnt, 1, sizeof(int), fp );
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (cnt!=REC_1_CNT*sizeof(char)) {
        longrecmark = 1;
        /*
         * read the first record again
         */
        fseek( fp,curr_pos,SEEK_SET );
        num = fread( &cnt, 1, sizeof(int), fp );
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (cnt!=REC_1_CNT*sizeof(char)) {
        /*
         * extra bytes at end of last record?
         * better to keep reading while num==1&&cnt==0 ??
         */
            for (i=0;i<FUDGE_64BIT;i++) num=fread( &cnt, 1, sizeof(int), fp );
            if (byte_swap) swap_4_bytes((char *)(&cnt));
        }
        fseek(fp,cnt,SEEK_CUR);
        num = fread( &cnt, 1, sizeof(int), fp );
        if (byte_swap) swap_4_bytes((char *)(&cnt));
        if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    }

        /* read the integer vars */
    num = fread( &cnt, 1, sizeof(int), fp );
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (cnt!=REC_2_CNT*sizeof(int)) return(READ_ERR);
    num = fread( tmp_vars, 1, REC_2_CNT*sizeof(int), fp );
    if (byte_swap)
        for (i=0,ptr=(char *)tmp_vars;i<REC_2_CNT;ptr+=4,i++)
            swap_4_bytes(ptr);
    if (tmp_vars[NUP]!=nup) return(REF_ERR);
    num = fread( &cnt, 1, sizeof(int), fp );
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

        /* skip the float vars */
    num = fread( &cnt, 1, sizeof(int), fp );
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (cnt!=REC_3_CNT*sizeof(float)) return(READ_ERR);
    fseek(fp,cnt+sizeof(int),SEEK_CUR);
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

    num = fread( &cnt, 1, sizeof(int), fp );
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (cnt!=(6*tmp_vars[NE]+tmp_vars[NUP]+tmp_vars[NBP]+tmp_vars[NBP2]
                +tmp_vars[NBP2])*sizeof(int))
        return(READ_ERR);
    fseek(fp,cnt+sizeof(int),SEEK_CUR);
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

    /*
     * read EXREF EYREF
     */
    num = fread( &cnt, 1, sizeof(int), fp );
    if (byte_swap) swap_4_bytes((char *)(&cnt));
    if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
    if (cnt!=tmp_vars[NUP]*sizeof(float)*2+(tmp_vars[NUVP]+tmp_vars[NBP2])
                                            *sizeof(int))
        return(READ_ERR);
    cnt = tmp_vars[NUP];
    if ((err = ReadFlArray(cnt,&arraysf[EXREF],
                                      fp,byte_swap))) {
        if (arraysf[EXREF]!=NULL) free(arraysf[EXREF]);
        return(err);
    }
    if ((err = ReadFlArray(cnt,&arraysf[EYREF],
                                      fp,byte_swap))) {
        if (arraysf[EXREF]!=NULL) free(arraysf[EXREF]);
        if (arraysf[EYREF]!=NULL) free(arraysf[EYREF]);
        return(err);
    }
    /*
     *  if Lagrangian arrays saved - read them
     */
    if (tmp_vars[ILAG]==1||tmp_vars[ILAG]==3) {
       cnt=(tmp_vars[NUVP]+tmp_vars[NBP2])*sizeof(int);
       fseek(fp,cnt+sizeof(int),SEEK_CUR);
       if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      /*
       * (c) optional crustal thickness / rotation block
       */
      if (tmp_vars[ICR]) {
          num = fread( &cnt, 1, sizeof(int), fp );
          if (byte_swap) swap_4_bytes((char *)(&cnt));
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
          if (cnt!=tmp_vars[NUP]*sizeof(int)) return(READ_ERR);
          fseek(fp,cnt+sizeof(int),SEEK_CUR);
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);

          if (tmp_vars[ICR]==1) {
              num = fread( &cnt, 1, sizeof(int), fp );
              if (byte_swap) swap_4_bytes((char *)(&cnt));
              if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
              if (cnt!=tmp_vars[NUP]*2*sizeof(float)) return(READ_ERR);
              fseek(fp,cnt+sizeof(int),SEEK_CUR);
              if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
          }
          else if (tmp_vars[ICR]==2) { 
              num = fread( &cnt, 1, sizeof(int), fp );
              if (byte_swap) swap_4_bytes((char *)(&cnt));
              if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
              if (cnt!=tmp_vars[NUP]*sizeof(float)) return(READ_ERR);
              fseek(fp,cnt+sizeof(int),SEEK_CUR);
              if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
          }
      }

      /*
       * (d) optional viscosity arrays
       */
      if (tmp_vars[IVIS]) {
          num = fread( &cnt, 1, sizeof(int), fp );
          if (byte_swap) swap_4_bytes((char *)(&cnt));
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
          if (cnt!=tmp_vars[NE]*8*sizeof(float)) return(READ_ERR);
          fseek(fp,cnt+sizeof(int),SEEK_CUR);
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      }

      /*
       * (e) optional density, temperature distribution arrays
       */
      if (tmp_vars[IDEN]) {
          num = fread( &cnt, 1, sizeof(int), fp );
          if (byte_swap) swap_4_bytes((char *)(&cnt));
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
          if (cnt!=tmp_vars[NE]*7*sizeof(float)) return(READ_ERR);
          fseek(fp,cnt+sizeof(int),SEEK_CUR);
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      }
      if (tmp_vars[ITEMP]) {
          num = fread( &cnt, 1, sizeof(int), fp );
          if (byte_swap) swap_4_bytes((char *)(&cnt));
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
          if (cnt!=tmp_vars[NUP]*sizeof(float)) return(READ_ERR);
          fseek(fp,cnt+sizeof(int),SEEK_CUR);
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
      }

      /*
       * (f) optional Lagrangian arrays
       */
      if (tmp_vars[ILAG]) {
          num = fread( &cnt, 1, sizeof(int), fp );
          if (byte_swap) swap_4_bytes((char *)(&cnt));
          if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
              /*
               * Read mesh arrays
               */
          if (tmp_vars[ILAG]==1||tmp_vars[ILAG]==3) {
            if (cnt!=(tmp_vars[NEL]*6+tmp_vars[NBL]*2)*sizeof(int))
              return(READ_ERR);
            fseek(fp,cnt+sizeof(int),SEEK_CUR);
            if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
                  /* next record */
            num = fread( &cnt, 1, sizeof(int), fp );
            if (byte_swap) swap_4_bytes((char *)(&cnt));
            if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
            if (tmp_vars[ILAG]==1) {
              if (cnt!=(tmp_vars[NUL]*4+tmp_vars[NPM]*tmp_vars[NSM]*2)
                                       *sizeof(float)) return(READ_ERR);
            }
            else
              if (cnt!=tmp_vars[NUL]*4 *sizeof(float)) return(READ_ERR);
            if ((err = ReadFlArray(tmp_vars[NUL],&arraysf[EXLREF],
                                        fp,byte_swap))) return(err);
            if ((err = ReadFlArray(tmp_vars[NUL],&arraysf[EYLREF],
                                        fp,byte_swap))) return(err);
            fseek(fp,2*cnt+sizeof(int),SEEK_CUR);
            if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
          }
              /*
               * Read marker arrays
               */
          if (tmp_vars[ILAG]==1||tmp_vars[ILAG]==2) {
            if (tmp_vars[ILAG]==2)
              if (cnt!=(tmp_vars[NPM]*tmp_vars[NSM]*2)
                                       *sizeof(float)) return(READ_ERR);
            cnt = tmp_vars[NSM]*tmp_vars[NPM];
            fseek(fp,2*cnt+sizeof(int),SEEK_CUR);
            if (longrecmark) fseek(fp,sizeof(int),SEEK_CUR);
          }
      /*
       * reference arrays not read from remainder of record
       */
      }
    }
    return( 0 );
}

int read_reference(FILE* fp, int* int_vars, float** arraysf, int current, 
                   int* ref, int req, int rotate)
{
    int err=0,rec_num,rotate_velocity=0;
    long curr_pos;

    curr_pos = ftell( fp );
    rec_num = current;
    if (rec_num>=req) {
        rewind( fp );
        rec_num=0;
    }
    while (rec_num<req-1) {
        if ((err=skip_data_blk( fp )))
            return( err );
        else rec_num++;
    }
    if ((err=read_ref_blk( fp,arraysf,int_vars[NUP] )))
        return( err );
    fseek( fp,curr_pos,SEEK_SET );
    *ref = req;
    if (rotate) {
        swapco_(arraysf[EXREF],arraysf[EYREF],NULL,
                                &int_vars[NUP],&rotate,&rotate_velocity );
        if (int_vars[ILAG]==1||int_vars[ILAG]==3)
             swapco_(arraysf[EXLREF],arraysf[EYLREF],NULL,
                                &int_vars[NUL],&rotate,&rotate_velocity );
    }
    return( 0 );
}

int ReadIntArray(int cnt, int** addr, FILE* fp, unsigned char byte_swap)
{
    int i;
    char *ptr;

    if (*addr==NULL)
        if ((*addr=(int *)malloc(cnt*sizeof(int)))==NULL)
                                         return(MALLOC_ERR);
    if (fread( *addr, sizeof(int), cnt, fp )!= cnt)
        return(READ_ERR);
    if (byte_swap)
        for (i=0,ptr=(char *)*addr;i<cnt;ptr+=4,i++)
            swap_4_bytes(ptr);
    return(0);
}

int ReadFlArray(int cnt, float** addr, FILE* fp, unsigned char byte_swap)
{
    int i;
    char *ptr;

    if (*addr==NULL)
        if ((*addr=(float *)malloc(cnt*sizeof(float)))==NULL)
                                         return(MALLOC_ERR);
    if (fread( *addr, sizeof(float), cnt, fp )!= cnt)
        return(READ_ERR);
    if (byte_swap)
        for (i=0,ptr=(char *)*addr;i<cnt;ptr+=4,i++)
            swap_4_bytes(ptr);
    return(0);
}

#define XY 0
#if XY
#endif

#ifdef SUN
void free_mem( arraysi,arraysf,maxl,maxf )
int **arraysi;
float **arraysf;
int maxl,maxf;
#else
void free_mem( int **arraysi,float **arraysf,int maxl,int maxf )
#endif
{
    int i;

    for (i=0;i<=(int)maxl;i++ ) free(arraysi[i]);
    for (i=0;i<=(int)maxf;i++ ) free(arraysf[i]);
}
