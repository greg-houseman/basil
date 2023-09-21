/*--------------------------------------------------------------------
 *    Basil / Sybil:   data.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

/* 
 * characters
 */
#define REC_1_CNT    128
/* 
 * integer scalars
 */
#define REC_2_CNT     64
/* 
 * float scalars
 */
#define REC_3_CNT     64

/*
 * define indices into integer variable array
 * (mesh size,type,controls)
 */

#define NX          0
#define NY          1
#define IMSH        2
#define NCOMP       3
#define ICR         4
#define IVIS        5
#define IDEN        6
#define ILAG        7
#define IFLT        8
#define IGRAV       9
#define NE         10
#define NN         11
#define NMP        12
#define NUP        13
#define NBP        14
#define NBP2       15
#define NUP2       16
#define NUVP       17
#define NROWS      NUVP
#define IDATA      18
#define IVV        19
#define KSTEP      20
#define ISAVE      21
#define IDT0       22
#define KEXIT      23
#define ITSTOP     24
#define IPRINT     25
#define IWRITE     26
#define KSAVE      27
#define INDFIX     28
#define MPDEF      29
#define LFLAT      30
#define JLV        31
#define IDEFTYP    32
#define MROWS      33
#define NXMAX      34
#define NYMAX      35
#define NFP        36
#define IFCASE     37 /* not used */
#define NELLEP     38
#define IPOLY      39
#define IBCMOD     40
#define NEL        41
#define NUL        42
#define NSM        43
#define NPM        44
#define NBL        45
#define NRM        46
#define MSINDX     47
#define IMREG      48
#define IRHEOTYP   49
#define IVOLD      50
#define ITEMP      51
#define IDENS      52
#define ITOPO      53
#define NSEG       54
#define ITHDI	   55
#define IPOLE      56
#define IVRESET    57

/*
 * define indices into real variable array
 * (control params,physical params ) 
 */

#define XLEN        0
#define YLEN        1
#define WFIT        2 /* unused */
#define AC          3
#define TIME        4
#define BIG         5
#define SE          6
#define HLENSC      7
#define ARGAN       8
#define ARGANP      9
#define BRGAN      10
#define BRGANP     11
#define THRESH     12
#define AREA       13
#define VOLUM      14
#define VISCD      15
#define DVMX       16
#define AREFM      17
#define ACNL       18
#define TSAVE      19
#define TEXIT      20
#define BANGL      21
#define ERA        22
#define OMTOT      23
#define DEFV1      24
#define DEFV2      25
#define DEFV3      26
#define DEFV4      27
#define DEFV5      28
#define BCV1       29
#define BCV2       30
#define BCV3       31
#define BCV4       32
#define BCV5       33
#define RHOG1      34
#define RHOG2      35
#define RHOG3      36
#define RHOG4      37
#define RHOG5      38
#define VISP1      39
#define VISP2      40
#define VISP3      41
#define VISP4      42
#define VISP5      43
#define YLDSTR     44 /* unused */
#define VELXO      45
#define VELYO      46
#define TBXOFF     47
#define TBYOFF     48
#define STELPR     49
#define BDEPSC     50
#define TDIFF      51
#define ALPHA      52
#define RISOST     53
#define REFLEV     54
#define XREFM      55
#define YREFM      56
#define GAMMA      57
#define IELAREA    58
#define IELSEG     59
#define IELQUAL    60
#define IBETA      61
#define ITREF      62
#define IVC        63

/* define indices into vector of integer array addresses */
#define LEM         0
#define NOR         1
#define IBC         2
#define IBNGH       3
#define IBCTYP      4
#define IELFIX      5

#define LGEM        6
#define LGIBC       7
#define LGIBCF      8
#define LNOR        9
#define IFBC       10
#define IFBC2      11
#define IFEQV      12
#define JFBC1      13
#define JFBC2      14
#define ELLENODE   15
#define POLYN      16
#define IMAT       17
#define ISEG       18

#define IHELP      19 /* if changed, check PLOT_INT_ARRAYS */
#define IWORK      20
#define IMP        21

/* define indices into vector of float array addresses */
#define EX        0
#define EY        1
#define UVP       2
#define QBND      3
#define SSQ       4
#define FROT      5
#define VHB       6
#define DENS      7

#define BY        8
#define CX        9
#define DNDP     10
#define PNI      11
#define EXREF    12
#define EYREF    13

#define EXLG     14
#define EYLG     15
#define EXLREF   16
#define EYLREF   17
#define UXLG     18
#define UYLG     19
#define STELPX   20
#define STELPY   21
#define TEMPT    22
#define VOLD     23

#define AMESH    24 /* if changed, check PLOT_FL_ARRAYS */
#define BMESH    25
#define CMESH    26
#define DMESH    27
#define EMESH    28
#define SNTRP    29
#define BNTRP    30
#define ANTRP    31
#define CNTRP    32
#define DNTRP    33
#define BNDS     34

/* define sizes of vectors holding array addresses */
#define MAX_INT_ARRAYS    25
#define MAX_FL_ARRAYS     35

/* define index of first plot arrays (follow solution arrays) */
#define PLOT_INT_ARRAYS    IHELP
#define PLOT_FL_ARRAYS     AMESH
