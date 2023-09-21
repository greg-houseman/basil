#ifndef _polyutils_h
#define _polyutils_h
#include <list>
#include "triangle.h"
#include "tripoly.h"
#include "polydata.h"

using std::list;

void findPolygon(REAL *pointlist, REAL *normlist,int *edgelist, int nseg,
                 int *seg_count,int *segs_left,
                 list<int> &processed);
#endif
