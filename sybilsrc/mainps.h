/*--------------------------------------------------------------------
 *    Basil / Sybil:   mainps.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

/* 
 * Structure to hold page defining variables -
 * equivalent info to that in plotWidget structure
 */

#ifndef _sybilps_h
#define _sybilps_h
typedef struct {
    int width_in_pts;
    int height_in_pts;
    int width_in_cells;
    int height_in_cells;
    int cell_width_in_pts;
    int cell_height_in_pts;
    float cell_width_in_cm;
    float cell_height_in_cm;
    cell_data *cell;
    float scale_factor;
} page_data;

#define DEFAULT_OUTPUT_FILE   "sybil.ps"
#define DEFAULT_OUTPUT_EXT   ".ps"

#endif /*_sybilps_h */
