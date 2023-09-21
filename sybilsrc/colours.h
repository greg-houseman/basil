/*--------------------------------------------------------------------
 *    Basil / Sybil:   colours.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#define STEPS        9
/*
 * No. of colours in the colourmap is PLOT_COLOURS + USEFUL_SIZE
 * which are defined in cmndefs.h
 * 8 basic colours are defined by the Red, Green, Blue arrays
 * and the plotting colours are defined in the Rgb arrays where
 * there are STEPS r g b entries. The program calculates rgb
 * values for the intermediate colours between each entry, giving
 * a total of PLOT_COLOURS.
 * The entries in the arrays should be in the range 0->255
 */
/*
 * Useful colours
 * background,foreground,red,grn,blue,cyn,yellow,magenta
 */
/*
 * greyscale
 */
#ifdef XSYB
static int       Red0[] = { 0,255,165, 135, 40,145,180,215};
static int     Green0[] = { 0,255,165, 135, 40,145,180,215};
static int      Blue0[] = { 0,255,165, 135, 40,145,180,215};
#else
static int       Red0[] = { 255,0,165, 135, 40,145,180,215};
static int     Green0[] = { 255,0,165, 135, 40,145,180,215};
static int      Blue0[] = { 255,0,165, 135, 40,145,180,215};
#endif

/*
 *  standard colourmap
 */
#ifdef XSYB
static int       Red1[] = {0,255,255,  0,  0,255,255,  0};
static int     Green1[] = {0,255,  0,255,  0,  0,255,255};
static int      Blue1[] = {0,255,  0,  0,255,255,  0,255};
#else
static int       Red1[] = {255,0,255,  0,  0,255,255,  0};
static int     Green1[] = {255,0,  0,255,  0,  0,255,255};
static int      Blue1[] = {255,0,  0,  0,255,255,  0,255};
#endif

/*
 * Plotting colours
 */
/*
 * greyscale
 * incr=16
 */
static int      Rgb0[STEPS*3]={    128,128,128,
                                   144,144,144,
                                   160,160,160,
                                   176,176,176,
                                   192,192,192,
                                   208,208,208,
                                   224,224,224,
                                   240,240,240,
                                   255,255,255
                                };
/*
 * greyscale
 * incr=24
static int      Rgb0[STEPS*3]={     43, 43, 43,
                                    67, 67, 67,
                                    91, 91, 91,
                                   115,115,115,
                                   139,139,139,
                                   163,163,163,
                                   187,187,187,
                                   211,211,211,
                                   235,235,235
                                };
 */
/*
 *  standard colourmap
 */
static int      Rgb1[STEPS*3]={       0,   0, 255,
                                    128, 128, 255,
                                    128, 255, 255,
                                      0, 255,   0,
                                    255, 255, 128,
                                    255, 255,   0,
                                    255, 128, 128,
                                    255,   0,   0,
                                    255, 128, 255
                                };
/*
 *  absolute value colourmap
static int      Rgb2[STEPS*3]={       0,   0, 255,
                                    128, 255, 255,
                                    255, 255, 128,
                                    255, 128, 128,
                                    255, 128, 255,
                                    255, 128, 128,
                                    255, 255, 128,
                                    128, 255, 255,
                                      0,   0, 255
                                };
 */

/*
 *  absolute value colourmap
 */
static int      Rgb2[STEPS*3]={       0,   0, 255,
                                    255,   0,   0,
                                    255,   0,   0,
                                    255,   0,   0,
                                    255,   0,   0,
                                    255,   0,   0,
                                    255,   0,   0,
                                    255,   0,   0,
                                      0,   0, 255
                                };
