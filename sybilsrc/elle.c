#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef XSYB
#include <X11/Xlib.h>  /* for type XFontStruct in Plot_info - tidy */
#endif

#include "types.h"
#include "cmndefs.h"
#include "string_utils.h"
#include "elle.h"
#include "errnum.h"
#include "error.h"
#include "routines.h"

typedef struct {
    float x;
    float y;
} Coords;
int PlotXY(Coords*, Coords*);

typedef struct {
    Coords cellBBox[4];
    float xoffset;
    float yoffset;
    float xlength;
    float ylength;
    float cum_xoffset;
} CellData;

CellData Unitcell;

int PlotElleRegions(FILE*, float*, float*, int);
int ReadElleLocations(FILE*, float**, float**, int*);

void InitCellData();
extern void ElleFileChosen();
extern void SetClipToPlot(), UnsetClipToPlot();
void FindBBox();
void ElleCheckUnit(float *xpts,float *ypts,int num,
                   int *xflags,int *yflags,
                   Coords *bl,Coords *tr);

int plot_elle_regions(float *pwindo)
{
    int err=0,size;
    float *xpts, *ypts;
    float xmin,xmax;

    InitCellData(&Unitcell);
    rewind(Plot_info.elle_file->fp);
    err = ReadElleLocations(Plot_info.elle_file->fp,&xpts,&ypts,&size);
    if (!err) {
        if (Settings.rescale) {
            if (Unitcell.cellBBox[0].x < Unitcell.cellBBox[3].x)
                 xmin = Unitcell.cellBBox[0].x;
            else xmin = Unitcell.cellBBox[3].x;
            if (Unitcell.cellBBox[1].x > Unitcell.cellBBox[2].x)
                 xmax = Unitcell.cellBBox[1].x;
            else xmax = Unitcell.cellBBox[2].x;
            pwindo[XCMINREF]=pwindo[XCMIN] = xmin;
            pwindo[XCMAXREF]=pwindo[XCMAX] = xmax;
            pwindo[XCMINREF]=pwindo[YCMIN] = Unitcell.cellBBox[0].y;
            pwindo[XCMAXREF]=pwindo[YCMAX] = Unitcell.cellBBox[2].y;
#if XY
            rangex = pwindo[XCMAX]-pwindo[XCMIN];
            rangey = pwindo[YCMAX]-pwindo[YCMIN];
            if (fabs(Settings.zoom-Settings.prev_zoom)>1.0e-4) {
                centrex = pwindo[XCMIN] + rangex/2;
                centrey = pwindo[YCMIN] + rangey/2;
                rangex /= Settings.zoom;
                rangey /= Settings.zoom;
                pwindo[XCMIN] = centrex-rangex/2;
                pwindo[XCMAX] = pwindo[XCMIN]+rangex;
                pwindo[YCMIN] = centrey-rangey/2;
                pwindo[YCMAX] = pwindo[YCMIN]+rangey;
            }
#endif
            init_box( pwindo,&Settings );
            Settings.rescale = 0;
        }
        /*else*/
            /*init_box( pwindo,&Settings );*/
        rewind(Plot_info.elle_file->fp);
        SetClipToPlot();
        err = PlotElleRegions(Plot_info.elle_file->fp,xpts,ypts,size);
    }
    if (!err) Save_plot_data( NULL,NULL );
    if (Settings.clip_to_cell==0) UnsetClipToPlot();
    return(err);
}

int ReadElleLocations(FILE* fp, float** xpts, float** ypts, int* size)
{
    char line[SYB_FILENAME_MAX], str[256];
    int err=0, max = 1024;
    int i, num, nn, end_opts=0, finished=0, key;
    double tmp;

    if ((*xpts = (float *)malloc(max * sizeof(float)))==NULL)
        return(MALLOC_ERR);
    if ((*ypts = (float *)malloc(max * sizeof(float)))==NULL)
        return(MALLOC_ERR);

    if (fgets(line,SYB_FILENAME_MAX-1,fp)==NULL) return(READ_ERR);
    while (!feof(fp) && !finished) {
        /*
         * find keywords
         */
        sscanf(line,"%s",str);
        validate(str,&key,FileKeys);
        switch (key) {
        case E_OPTIONS :
            while (!feof(fp) && !err && !end_opts) {
                if ((num = fscanf(fp,"%s", str))!=1) return(READ_ERR);
                else if (str[0] == '#') dump_comments( fp );
                else {
                    validate(str,&i,run_option_terms);
                    switch(i) {
                    case RO_CELLBBOX :
                        if (fscanf(fp,"%f %f\n",&Unitcell.cellBBox[0].x,
                                                &Unitcell.cellBBox[0].y)!=2)
                                        err = READ_ERR;
                        else if (fscanf(fp,"%f %f\n",&Unitcell.cellBBox[1].x,
                                                &Unitcell.cellBBox[1].y)!=2)
                                        err = READ_ERR;
                        else if (fscanf(fp,"%f %f\n",&Unitcell.cellBBox[2].x,
                                                &Unitcell.cellBBox[2].y)!=2)
                                        err = READ_ERR;
                        else if (fscanf(fp,"%f %f\n",&Unitcell.cellBBox[3].x,
                                                &Unitcell.cellBBox[3].y)!=2)
                                        err = READ_ERR;
                        Unitcell.xlength = Unitcell.cellBBox[1].x -
                                              Unitcell.cellBBox[0].x;
                        Unitcell.ylength = Unitcell.cellBBox[3].y -
                                              Unitcell.cellBBox[0].y;
                        Unitcell.xoffset = Unitcell.cellBBox[3].x -
                                              Unitcell.cellBBox[0].x;
                        Unitcell.yoffset = Unitcell.cellBBox[0].y -
                                              Unitcell.cellBBox[1].y;
                        if (Unitcell.xoffset >= 1.0) {
                            Unitcell.xoffset =
                                (float)modf((double)Unitcell.xoffset,&tmp);
                            Unitcell.cum_xoffset += (float)tmp;
                        }
                        else Unitcell.cum_xoffset = Unitcell.xoffset;
                        Unitcell.yoffset = Unitcell.cellBBox[1].y -
                                            Unitcell.cellBBox[0].y;
                        break;
                    case RO_SSOFFSET:
                        if (fscanf(fp,"%f\n",&Unitcell.xoffset)!=1)
                        break;
                    case RO_CUMSSOFFSET:
                        if (fscanf(fp,"%f\n",&Unitcell.cum_xoffset)!=1)
                        break;
                    case RO_SWITCHDIST:
                    case RO_MAXNODESEP:
                    case RO_MINNODESEP:
                    case RO_SPEEDUP   :
                        dump_comments(fp);
                        break;
                    default: end_opts = 1;
                        break;
                    }
                }
            }
            break;
        case LOCATION :
            while (!feof(fp) && !err) {
                num = fscanf(fp,"%s", str);
                if (str[0] == '#') dump_comments( fp );
                else if (str[0]<'0' || str[0]>'9') return(0);
                else {
                    *size = nn = atoi(str);
                    if (nn>max) {
                        do { max*=1.5; } while (nn>max);
                        if ((*xpts = (float *)realloc(*xpts,max*sizeof(float)))
                                         ==NULL) return(MALLOC_ERR);
                        if ((*ypts = (float *)realloc(*ypts,max*sizeof(float)))
                                         ==NULL) return(MALLOC_ERR);
                    }
                    if ((num = fscanf(fp,"%f %f\n",&((*xpts)[nn]),
                                   &((*ypts)[nn]))) != 2) return(READ_ERR);
                }
                num=fscanf(fp,"\n");
            }
            finished=1;
            break;
        }
        finished=(fgets(line,SYB_FILENAME_MAX-1,fp)==NULL);
    }
    return(err);
}

int PlotElleRegions(FILE* fp, float* xpts, float* ypts, int size)
{
    char line[SYB_FILENAME_MAX], str[256];
    int err=0, i, j, k, s, u, pen_state;
    int num, nn, numnodes, finished=0, key;
 // int first;
    int col, outline;
    int xflags[2],yflags[2];
    float eps=1.5e-4;
    float xoffset, incrx, incry;
    float *xp=0, *yp=0, *xporig=0, *yporig=0;
    Coords xy, prev;
    Coords bl,tr,tmpbl,tmptr;

    col = Settings.fg;
    outline = 2;  /* no fill */
    if (fgets(line,SYB_FILENAME_MAX-1,fp)==NULL) return(READ_ERR);
    while (!feof(fp) && !finished) {
        /*
         * find keywords
         */
        num=sscanf(line,"%s",str);
        validate(str,&key,FileKeys);
        if (key==REGIONS || key==FLYNNS) {
            incrx = Unitcell.xlength;
            incry = Unitcell.ylength;
            xoffset = Unitcell.xoffset;
            if ((Unitcell.cum_xoffset<-eps || Unitcell.cum_xoffset>eps) &&
                    xoffset>Unitcell.xlength-eps &&
                        xoffset<Unitcell.xlength+eps)
                xoffset=0;
            while (!feof(fp) && !err) {
                num=fscanf(fp,"%s",str);
                if (str[0] == '#') dump_comments( fp );
                else if (str[0]<'0' || str[0]>'9') return(0);
                else {
                    if (key==REGIONS)
                        if ((num = fscanf(fp,"%d %d",&s,&u))!=2)
                            return(READ_ERR);
                    if ((num = fscanf(fp,"%d", &numnodes))!=1)
                        return(READ_ERR);
                    if ((num = fscanf(fp,"%d", &nn))!=1) return(READ_ERR);
                    if ((xp=(float *)malloc(numnodes*sizeof(float)))==0)
                        return(MALLOC_ERR);
                    if ((yp=(float *)malloc(numnodes*sizeof(float)))==0)
                        return(MALLOC_ERR);
                    if ((xporig=(float *)malloc(numnodes*sizeof(float)))==0)
                        return(MALLOC_ERR);
                    if ((yporig=(float *)malloc(numnodes*sizeof(float)))==0)
                        return(MALLOC_ERR);
     //             first = nn;
                    yflags[0] = yflags[1] = 0;
                    prev.x = xpts[nn]; prev.y = ypts[nn];
                    xp[0] = xporig[0] = xpts[nn];
                    yp[0] = yporig[0] = ypts[nn];
                    pen_state = PEN_UP;
                    plotu_(&xpts[nn],&ypts[nn],&pen_state);
                    for (i=1;i<numnodes;i++) {
                        if ((num = fscanf(fp,"%d", &nn))!=1)
                            return(READ_ERR);
                        xy.x = xpts[nn]; xy.y = ypts[nn];
                        PlotXY(&xy,&prev);
                        xp[i] = xy.x;  yp[i] = xy.y;
                        xporig[i] = xy.x; yporig[i] = xy.y;
                        prev = xy;
                    }
                    fillpoly_(xp,yp,&numnodes,&col,&outline);
                    if (Settings.plot_opts.dble) {
                        for (i=0;i<numnodes;i++) 
                            xp[i] += incrx; 
                        fillpoly_(xp,yp,&numnodes,&col,&outline);
                    }
                    ElleCheckUnit(xp,yp,numnodes,xflags,yflags,&bl,&tr);
    /*
     * if xoffset >= incrx/2 may need to draw offset back a unit
     */
                    int start=0;
                    if (yflags[0]||yflags[1]) start=-1;
                    for (i=start;i<xflags[0]+1;i++) {
                        for (j=0;j<numnodes;j++) {
                            xp[j] = xporig[j] + i*incrx;
                            yp[j] = yporig[j];
                        }
                        tmpbl = bl;
                        tmptr = tr;
                        tmpbl.x += i*incrx;
                        tmptr.x += i*incrx;
                        if (tmpbl.y<Unitcell.cellBBox[TOPRIGHT].y &&
                                tmptr.y>Unitcell.cellBBox[BASELEFT].y)
                            fillpoly_(xp,yp,&numnodes,&col,&outline);
                        for (k=0; k< yflags[0]; k++) {
                            for (j=0;j<numnodes;j++) xp[j] += (k+1)*xoffset;
                            for (j=0;j<numnodes;j++) yp[j] += (k+1)*incry;
                            tmpbl.x += (k+1)*xoffset;
                            tmpbl.y += (k+1)*incry;
                            tmptr.x += (k+1)*xoffset;
                            tmptr.y += (k+1)*incry;
                            if (tmpbl.x<Unitcell.cellBBox[TOPRIGHT].x &&
                                    tmptr.x>Unitcell.cellBBox[BASELEFT].x) 
                                fillpoly_(xp,yp,&numnodes,&col,&outline);
                        }
                        for (j=0;j<numnodes;j++) {
                            xp[j] = xporig[j] + i*incrx;
                            yp[j] = yporig[j];
                        }
                        tmpbl.x = bl.x + i*incrx;
                        tmpbl.y = bl.y;
                        tmptr.x = tr.x + i*incrx;
                        tmptr.y = tr.y;
                        for (k=0; k< yflags[1]; k++) {
                            for (j=0;j<numnodes;j++) xp[j] -= (k+1)*xoffset;
                            for (j=0;j<numnodes;j++) yp[j] -= (k+1)*incry;
                            tmpbl.x -= (k+1)*xoffset;
                            tmpbl.y -= (k+1)*incry;
                            tmptr.x -= (k+1)*xoffset;
                            tmptr.y -= (k+1)*incry;
                            if (tmpbl.x<Unitcell.cellBBox[TOPRIGHT].x &&
                                    tmptr.x>Unitcell.cellBBox[BASELEFT].x)
                                fillpoly_(xp,yp,&numnodes,&col,&outline);
                        }
                        if (start==0) start=1;
                    }
                    for (i=start;i<xflags[1]+1;i++) {
                        for (j=0;j<numnodes;j++) {
                            xp[j] = xporig[j] - i*incrx;
                            yp[j] = yporig[j];
                        }
                        tmpbl.x = bl.x - i*incrx;
                        tmpbl.y = bl.y;
                        tmptr.x = tr.x - i*incrx;
                        tmptr.y = tr.y;
                        if (tmpbl.y<Unitcell.cellBBox[TOPRIGHT].y &&
                                tmptr.y>Unitcell.cellBBox[BASELEFT].y)
                            fillpoly_(xp,yp,&numnodes,&col,&outline);
                        for (k=0; k< yflags[0]; k++) {
                            for (j=0;j<numnodes;j++) xp[j] += (k+1)*xoffset;
                            for (j=0;j<numnodes;j++) yp[j] += (k+1)*incry;
                            tmpbl.x += (k+1)*xoffset;
                            tmpbl.y += (k+1)*incry;
                            tmptr.x += (k+1)*xoffset;
                            tmptr.y += (k+1)*incry;
                            if (tmpbl.x<Unitcell.cellBBox[TOPRIGHT].x &&
                                    tmptr.x>Unitcell.cellBBox[BASELEFT].x)
                                fillpoly_(xp,yp,&numnodes,&col,&outline);
                        }
                        for (j=0;j<numnodes;j++) {
                            xp[j] = xporig[j] - i*incrx;
                            yp[j] = yporig[j];
                        }
                        tmpbl.x = bl.x - i*incrx;
                        tmpbl.y = bl.y;
                        tmptr.x = tr.x - i*incrx;
                        tmptr.y = tr.y;
                        for (k=0; k< yflags[1]; k++) {
                            for (j=0;j<numnodes;j++) xp[j] -= (k+1)*xoffset;
                            for (j=0;j<numnodes;j++) yp[j] -= (k+1)*incry;
                            tmpbl.x -= (k+1)*xoffset;
                            tmpbl.y -= (k+1)*incry;
                            tmptr.x -= (k+1)*xoffset;
                            tmptr.y -= (k+1)*incry;
                            if (tmpbl.x<Unitcell.cellBBox[TOPRIGHT].x &&
                                    tmptr.x>Unitcell.cellBBox[BASELEFT].x)
                                fillpoly_(xp,yp,&numnodes,&col,&outline);
                        }
                        if (start==0) start=1;
                    }
                    for (i=start;i<yflags[1]+1;i++) {
                        for (j=0;j<numnodes;j++) xp[j] = xporig[j]+i*xoffset;
                        for (j=0;j<numnodes;j++) yp[j] = yporig[j]+i*incry;
                        tmpbl.x = bl.x+i*xoffset;
                        tmpbl.y = bl.y+i*incry;
                        tmptr.x = tr.x+i*xoffset;
                        tmptr.y = tr.y+i*incry;
                        if (tmpbl.x<Unitcell.cellBBox[TOPRIGHT].x &&
                                tmptr.x>Unitcell.cellBBox[BASELEFT].x)
                            fillpoly_(xp,yp,&numnodes,&col,&outline);
                        if (start==0) start=1;
                    }
                    for (i=start;i<yflags[0]+1;i++) {
                        for (j=0;j<numnodes;j++) xp[j] = xporig[j]+i*xoffset;
                        for (j=0;j<numnodes;j++) yp[j] = yporig[j]+i*incry;
                        tmpbl.x = bl.x+i*xoffset;
                        tmpbl.y = bl.y+i*incry;
                        tmptr.x = tr.x+i*xoffset;
                        tmptr.y = tr.y+i*incry;
                        if (tmpbl.x<Unitcell.cellBBox[TOPRIGHT].x &&
                                tmptr.x>Unitcell.cellBBox[BASELEFT].x)
                            fillpoly_(xp,yp,&numnodes,&col,&outline);
                    }
                    if (xp){free(xp); xp=0;}
                    if (yp){free(yp); yp=0;}
                    if (xporig){free(xporig); xporig=0;}
                    if (yporig){free(yporig); yporig=0;}
                }
                num=fscanf(fp,"\n");
            }
            finished=1;
        }
        finished=(fgets(line,SYB_FILENAME_MAX-1,fp)==NULL);
    }
    return(err);
}

#define MAX_SIZE 0.75
int PlotXY(xy, prevxy)
Coords *xy, *prevxy;
{
    float unitsize_x,unitsize_y,xoffset,yoffset;

    /*
     * assumes that the unit cell remains a parallelogram
     */
    unitsize_x = Unitcell.xlength;
    unitsize_y = Unitcell.ylength;
    xoffset = Unitcell.xoffset;
    yoffset = Unitcell.yoffset;
    if ((xy->y - prevxy->y) >= unitsize_y*0.5) {
        xy->y -= unitsize_y;
        xy->x -= xoffset;
        while ((xy->y - prevxy->y) >= unitsize_y*MAX_SIZE) {
            xy->y -= unitsize_y;
            xy->x -= xoffset;
        }
    }
    else if ((xy->y - prevxy->y) < -unitsize_y*0.5) {
        xy->y += unitsize_y;
        xy->x += xoffset;
		while ((xy->y - prevxy->y) < -unitsize_y*MAX_SIZE) {
		    xy->y += unitsize_y;
		    xy->x += xoffset;
        }
    }
    if ((xy->x - prevxy->x) >= unitsize_x*0.5) {
		xy->x -= unitsize_x;
		xy->y -= yoffset;
        while ((xy->x - prevxy->x) >= unitsize_x*MAX_SIZE) {
            xy->x -= unitsize_x;
            xy->y -= yoffset;
        }
	}
    else if ((xy->x - prevxy->x) < -unitsize_x*0.5) {
		xy->x += unitsize_x;
		xy->y += yoffset;
        while ((xy->x - prevxy->x) < -unitsize_x*MAX_SIZE) {
            xy->x += unitsize_x;
            xy->y += yoffset;
        }
	}
    return(0);
}

void ElleCheckUnit(float *xpts,float *ypts,int num,
                   int *xflags,int *yflags,
                   Coords *bl,Coords *tr)
{
    int i;
    float minx,maxx,miny,maxy,xp,gamma;
    float eps;
    float xoffset;
    /*CellData unitcell;*/

    eps = 1.5e-6;
    eps = 5e-5;
    /*ElleCellBBox(&unitcell);*/
    maxy = (float)Unitcell.cellBBox[TOPLEFT].y;
    miny = (float)Unitcell.cellBBox[BASELEFT].y;
    maxx = (float)Unitcell.cellBBox[BASERIGHT].x;
    minx = (float)Unitcell.cellBBox[BASELEFT].x;
    xoffset = (float)(Unitcell.cellBBox[TOPLEFT].x-
                      Unitcell.cellBBox[BASELEFT].x);
    if ((Unitcell.cum_xoffset<-eps || Unitcell.cum_xoffset>eps) &&
                            xoffset>Unitcell.xlength-eps &&
                            xoffset<Unitcell.xlength+eps)
        xoffset=Unitcell.xlength;
    gamma = xoffset * (float)Unitcell.ylength;
    xflags[0] = xflags[1] = 0;
    yflags[0] = yflags[1] = 0;
    bl->x = tr->x = xpts[0];
    bl->y = tr->y = ypts[0];
    for (i=0;i<num;i++) {
        if (ypts[i]<miny) yflags[0]=1;
        if (ypts[i]>maxy) yflags[1]=1;
        xp = (ypts[i] - miny)*gamma + Unitcell.cellBBox[0].x;
        if (xpts[i]<(xp-eps)) xflags[0]=1;
        xp += Unitcell.xlength;
        if (xpts[i]>(xp+eps)) xflags[1]=1;
        if (ypts[i]<bl->y) bl->y = ypts[i];
        if (ypts[i]>tr->y) tr->y = ypts[i];
        if (xpts[i]<bl->x) bl->x = xpts[i];
        if (xpts[i]>tr->x) tr->x = xpts[i];
    }
    /* this hasn't been checked for -ve xoffset */
    while (bl->y<(miny-Unitcell.ylength*yflags[0])) yflags[0]++;
    while (tr->y>(maxy+Unitcell.ylength*yflags[1])) yflags[1]++;
    while (bl->x<(minx-(Unitcell.xlength-xoffset)*xflags[0]))
        xflags[0]++;
    while (tr->x>(maxx+(Unitcell.xlength-xoffset)*xflags[1]))
        xflags[1]++;
}



void InitCellData(CellData *cell)
{
    Coords DefaultTopLeft;
    Coords DefaultTopRight;
    Coords DefaultBaseLeft;
    Coords DefaultBaseRight;

    DefaultTopLeft.x = 0.0;
    DefaultTopLeft.y = 1.0;
    DefaultTopRight.x = 1.0;
    DefaultTopRight.y = 1.0;
    DefaultBaseLeft.x = 0.0;
    DefaultBaseLeft.y = 0.0;
    DefaultBaseRight.x = 1.0;
    DefaultBaseRight.y = 0.0;
    Unitcell.cellBBox[0] = DefaultBaseLeft;
    Unitcell.cellBBox[1] = DefaultBaseRight;
    Unitcell.cellBBox[2] = DefaultTopRight;
    Unitcell.cellBBox[3] = DefaultTopLeft;
    Unitcell.xlength = Unitcell.cellBBox[1].x -
                                  Unitcell.cellBBox[0].x;
    Unitcell.ylength = Unitcell.cellBBox[3].y -
                                  Unitcell.cellBBox[0].y;
    Unitcell.xoffset = Unitcell.cellBBox[3].x -
                                  Unitcell.cellBBox[0].x;
    Unitcell.cum_xoffset = 0;

}

