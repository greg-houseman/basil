#include <iostream>
#include <list>
#include <algorithm>
#include <utility>
#include <cstdio>
#include <cmath>

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */
#include "polyutils.h"

using std::find;
using std::cout;
using std::endl;

extern "C" {
int crossingstest_( double  *pgon, int *nverts, double  point[2],
                    int *res );
}
int CartesianToPolar(double x, double y, double z,
                     double *angxy, double *angz);
double polar_angle(double x,double y);
double polyArea(double *xy,int n);

int getpolydata(struct triangulateio *io,polydata *polys,int nr)
{
    bool found=0;
    int segs_left=0, *seg_count;
    int i, j, dum, nseg, end1, end2, num_pts, err=0;
    int *edgelist=0, *markerlist=0;
    double *ppts, pt[2];
    list<int> segments;
    list<int>::iterator it;
    polydata tmp;

    tmp.numpts = 0; tmp.ppts = 0;
    for (i=0;i<nr;i++) {
        polys[i].numpts=0;
        polys[i].ppts=0;
    }
    if ((nseg = io->numberofedges)>0) {
        edgelist =  io->edgelist;
    }
    else if ((nseg = io->numberofsegments)>0) {
        edgelist =  io->segmentlist;
        markerlist =  io->segmentmarkerlist;
    }
    else return(1);
        
    if ((seg_count = new int[nseg])==0)
        return(1);
    for (i=0;i<nseg;i++) seg_count[i]=0;
    segs_left = nseg;
    if (io->normlist!=0)
        for (i=0; i<nseg; i++) {
            if (edgelist[i*2]==-1 || edgelist[i*2+1]==-1)
                seg_count[i] = 1;
        }
    i=0;
    while (segs_left>0) {
        findPolygon(io->pointlist,io->normlist,
                    edgelist,nseg,
                    seg_count,&segs_left,segments);
                                                                                
        if ((num_pts = segments.size()) > 2) {
            tmp.numpts = num_pts;
            if ((tmp.ppts = (double *)malloc(sizeof(double)*num_pts*2))
                ==0) return(1);
            it=segments.begin();
            end1 = edgelist[*it *2 ];
            end2 = edgelist[*it *2 +1];
            it++;
            if (end2 != edgelist[*it *2] &&
                    end2 != edgelist[*it *2 +1])
            end2 = end1;
            tmp.ppts[0] = io->pointlist[end2*2];
            tmp.ppts[1] = io->pointlist[end2*2+1];
            /*cout << end2 << ' ' << tmp_x[0] << ' ' << tmp_y[0] << '
 * ';*/
            /*cout << endl;      */
            i=1;
            while (it!=segments.end()) {
                if (end2 == edgelist[*it *2])
                    end2 = edgelist[*it *2 +1];
                else end2 = edgelist[*it *2];
                /*cout << edgelist[*it *2 ] << ' ';*/
                /*cout << edgelist[*it *2 +1] << ' ';*/
                tmp.ppts[i*2] = io->pointlist[end2*2];
                tmp.ppts[i*2+1] = io->pointlist[end2*2+1];
                /*cout << end2 << ' ' << tmp_x[i] << ' ' << tmp_y[i] '
 * ';*/
                /*cout << endl;      */
                it++; i++;
            }
            // there must be a better way to stop it finding 
            // the cell boundary
            found=(polyArea(tmp.ppts,tmp.numpts)<0);
            for (j=0;j<nr&&!found;j++){
                if (polys[j].numpts==0) {
                   pt[0] = io->regionlist[j*4];
                   pt[1] = io->regionlist[j*4+1];
                   if ((found=crossingstest_(tmp.ppts,&tmp.numpts,pt,&dum))) {
                       polys[j].numpts=tmp.numpts;
                       polys[j].ppts=tmp.ppts;
                   }
                }
            }
        }
    }
    delete[] seg_count;
    return(err);
}

void findPolygon(REAL *pointlist, REAL *normlist, int *edgelist,
                 int nseg, int *seg_count,int *segs_left,
                 list<int> &processed)
{
    unsigned char finished = 0, valid=0;
    int startindx=0, i=0, index[4];
    int end1, end2, startnode, next, curr, tmp;
    list<int>::iterator it;
    std::pair<int,double> trial, nb;
    list<int> seg_index;
    list< std::pair<int,double> > ordered_list;
    list< std::pair<int,double> >::iterator ito;

    startindx=0;
    startnode=-1;
    while (!finished && startindx<nseg) {
        processed.clear();
        // find a valid starting segment
        while (startindx<nseg && (seg_count[startindx]==2 ||
               edgelist[startindx*2]==-1 || 
               edgelist[startindx*2+1]==-1)) startindx++;
        if (startindx<nseg) {
            // process the starting segment
            startnode = curr = end1 = edgelist[startindx*2];
            next = end2 =edgelist[startindx*2+1];
            index[0] = index[1] = startindx;
            processed.push_back(index[1]);
            seg_count[index[1]]++;
            if (seg_count[index[1]]==2) (*segs_left)--;
            /*
             * reverse order to prevent starting a
             * future traversal in the same direction
             */
            if (curr==edgelist[index[1]*2]) {
                edgelist[index[1]*2]=edgelist[index[1]*2+1];
                edgelist[index[1]*2+1]=curr;
            }
            valid = 1;
        }
        while (!finished &&startindx<nseg && valid) {
            /*
             * make a list of nbs in angle order,
             * including curr (angles are +ve anticlockwise
             * around node). The next node is the one
             * before curr in the list
             */
            i=0;
            nb.first = curr;
            nb.second =  polar_angle(
                         pointlist[curr*2]-pointlist[next*2],
                         pointlist[curr*2+1]-pointlist[next*2+1]);
            seg_index.clear();
            seg_index.push_back(index[1]);
            ordered_list.clear();
            ordered_list.push_back(nb);
            while (i<nseg) {
                if (seg_count[i]<2 && i!=index[1] &&
                     ( edgelist[i*2]==next ||
                         edgelist[i*2+1]==next) ) {
                       if (edgelist[i*2]==next)
                           trial.first = edgelist[i*2+1];
                       else trial.first = edgelist[i*2];
                       if (trial.first==-1)
                         trial.second =  polar_angle(
                           normlist[i*2],
                           normlist[i*2+1]);
                       else
                         trial.second =  polar_angle(
                           pointlist[trial.first*2]-
                                            pointlist[next*2],
                           pointlist[trial.first*2+1]-
                                            pointlist[next*2+1]);
                       ito=ordered_list.begin(); it = seg_index.begin();
                       while (ito!=ordered_list.end() &&
                               trial.second > (*ito).second) { ito++;it++; }
                       ordered_list.insert(ito,trial);
                       seg_index.insert(it,i);
                }
                i++;
            }
/*
                for(ito=ordered_list.begin();ito!=ordered_list.end();ito++)
                   cout<<ito->first<<' '<<ito->second<<"  ";
                cout<<endl;
*/
            /*
             * The list will always include curr
             * if size==1, reached dead-end
             */
            if (ordered_list.size()>1) {
#if XY
                if (ordered_list.back().first==curr) {
                    nb=ordered_list.front();
                    index[1] = seg_index.front();
                }
                else {
                    ito=ordered_list.begin(); it = seg_index.begin();
                    while (ito!=ordered_list.end() &&
                          ito->first!=curr) { ito++; it++; }
                    nb=*(++ito);
                    index[1] = *(++it);
                }
#endif
                if (ordered_list.front().first==curr) {
                    nb=ordered_list.back();
                    index[1] = seg_index.back();
                }
                else {
                    ito=ordered_list.begin(); it = seg_index.begin();
                    while (ito!=ordered_list.end() &&
                          ito->first!=curr) { ito++; it++; }
                    ito--; it--;
                    nb=*(ito);
                    index[1] = *(it);
                }
       
                processed.push_back(index[1]);
                seg_count[index[1]]++;
                if (seg_count[index[1]]==2) (*segs_left)--;
                index[0] = index[1];
                curr = next;
                next = nb.first;
                if (curr==edgelist[index[1]*2]) {
                    edgelist[index[1]*2]=edgelist[index[1]*2+1];
                    edgelist[index[1]*2+1]=curr;
                }
                finished = (next==startnode);
                if (next==-1) valid = 0;
            }
            else {
            // not valid starting segment, try again
                valid = 0;
            }
        }
    }
    if (!finished) {
        if (*segs_left) {
            printf("findPolygon-no start segment found\n");
            // and do what???
        }
    }
}

int CartesianToPolar(double x, double y, double z,
                     double *angxy, double *angz)
{
    double eps = 1e-9;
    *angz = asin(z);
    if (x==0.0) {
        if (y>=0.0) *angxy = M_PI*0.5;
        else        *angxy = M_PI * 1.5;
    }
    else {
        *angxy = atan(y/x);
        if (x<0.0) *angxy += M_PI;
        if (*angxy < -eps) *angxy += 2.0*M_PI;
        else if (*angxy > (2.0*M_PI-eps)) *angxy -= 2.0*M_PI;
    }
    return(0);
}

double polyArea(double *xy, int n)
{
    register int i, j;
    double ai, area = 0;
    if (n < 3) return 0;
    for (i = n-1, j = 0; j < n; i = j, j++) {
        ai = xy[i*2] * xy[j*2+1] - xy[j*2] * xy[i*2+1]; /*+ve counterclockwise*/
        area += ai;
    }
    return (area/2);
}

double polar_angle(double x,double y)
{
    int err = 0;
    double theta=0.0, angz=0.0, z=0.0;
                                                                                
    err=CartesianToPolar(x,y,z,&theta,&angz);
    if (theta<0) theta += 2.0*M_PI;
    return(theta);
}
