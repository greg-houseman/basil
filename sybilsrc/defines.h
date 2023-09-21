/*--------------------------------------------------------------------
 *    Basil / Sybil:   defines.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#ifdef XSYB
#define fontdat XFontStruct *
#define colourdat XColor
#define BLACK_PEN  0
#else
#define fontdat int
#define colourdat int
#define BLACK_PEN  1
#endif
