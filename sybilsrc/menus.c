
/*--------------------------------------------------------------------
 *    Basil / Sybil:   menus.c  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

/* 
 *  menus.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
#ifdef SUN
#include <sys/param.h>
#endif
*/
#include <errno.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#include <X11/cursorfont.h>
#include <Xm/Xm.h>

/*
 * Public include files for widgets
 */
#include <Xm/MainW.h>
#include <Xm/RowColumn.h>
#include <Xm/PushBG.h>
#include <Xm/PushB.h>
#include <Xm/CascadeBG.h>
#include <Xm/CascadeB.h>
#include <Xm/Label.h>
#include <Xm/SeparatoG.h>
#include <Xm/SelectioB.h>
#include <Xm/MessageB.h>
#include <Xm/TextF.h>
#include <Xm/ToggleBG.h>
#include <Xm/FileSB.h>

#include "types.h"
#include "plot.h"
#include "menus.h"
#include "sybfile.h"
#include "globals.h"
#include "strain.h"
#include "deform.h"
#include "errnum.h"
#include "error.h"
#include "pref77.h"
#include "string_utils.h"
#include "log.h"
#include "help.h"
#include "routines.h"

#define XY 0
#define NO_VALUE        -10

#define OPEN              0
#define SAVE              1
#define SEP               2
#define EXIT              3

#define SOLN              0
#define LOG               1
#define STTGS             2
#define LOGPS             3

#define LABEL_INDX        0
#define COLOUR_INDX       1
#define LINE_INDX         2
#define DEL_INDX          3

#define EDIT_INDX         0
#define FONT_INDX         1
#define LRGE_FONT         5

#define CURR_CELL         11
#define ALL_CELLS         12

#define NEXT              0
#define PREV              1
#define SPEC              2
#define LAST              3
#define FIRST             4

/*
#define SOLID            -1
#define DASH              1
*/

#define VIEW              0
#define SOLNFILESDIR    "FD.sols"
#define ELLE_LABEL      "Elle"

#define MAX_FONTS         11

char Lastfilemask[SYB_FILENAME_MAX+1];
char Lastlogmask[SYB_FILENAME_MAX+1];
Cursor Waiting;
Widget Main;
Widget main_Menu, topLevel;
Widget LabelPrompt,RefPrompt,CurrPrompt;
/* Widget ZoomPrompt; */
Widget plotArea, currvalsArea;
XmFontList Curr_fontlist;
XmFontList Fontlists[MAX_FONTS];
XtAppContext App_context;
Dimension Minwidth;

/* These values are only used if the app_defaults file is not found */
String default_resources[] = {
        /* Plot resources*/
    "*cellWidthInPixels:      500",
    "*cellHeightInPixels:     200",
/*
    "*pixmapWidthInCells:       2",
    "*pixmapHeightInCells:      2",
*/
        /* Appearance Resources*/
    "*main_Menu*Foreground:     black",
    "*main_Menu*Background:     grey75",
    "*plotArea*Foreground:   white",
    "*plotArea*Background:   black",
    "current values*Foreground:   white",
    "currvalsArea*Background:   grey75",
    "*Help*Foreground:     white",
    "*Help*Background:   grey75",
    "*Xsybil*fontList: -*-helvetica-bold-r-*--*-120-*",
 /* "*Xsybil*fontList: -*-helvetica-r-*--*-120-*", */
        /* Fonts - Font5 used for title area text */
/*
    "*Helv 8*fontList: -*-helvetica-medium-r-*--*-80-*=TAG08H",
    "*Helv10*fontList: -*-helvetica-medium-r-*--*-100-*=TAG10H",
    "*Helv12*fontList: -*-helvetica-bold-r-*--*-120-*=TAG12H",
    "*Helv14*fontList: -*-helvetica-bold-r-*--*-140-*=TAG14H",
    "*Helv18*fontList: -*-helvetica-bold-r-*--*-180-*=TAG18H",
    "*Helv24*fontList: -*-helvetica-bold-r-*--*-240-*=TAG24H",
    "*Symbol10*fontList: -*-symbol-medium-r-*--*-100-*=TAG10S",
    "*Symbol12*fontList: -*-symbol-medium-r-*--*-120-*=TAG12S",
    "*Symbol14*fontList: -*-symbol-medium-r-*--*-140-*=TAG14S",
    "*Symbol18*fontList: -*-symbol-medium-r-*--*-180-*=TAG18S",
    "*Symbol24*fontList: -*-symbol-medium-r-*--*-240-*=TAG24S",
*/
    NULL
};
Widget BuildMenu();
void FileOpenChosen(),FileSaveChosen(),CurrChosen(),LabelChosen();
void RefChosen(),LabelDialog(),ExitDialog(),Exit_Proc();
void LocateChosen(),FontChosen(),ColourChosen(),LineChosen();
void DeleteChosen(),RescaleChosen(),HelpChosen();
void ZoomChosen();
void StrainArrowChosen(),StrainOptChosen();
void StressArrowChosen(),StressOptChosen();
void LayerOptChosen(),GravityOptChosen();
void RotationChosen(),RotationArrowChosen();
void VelocityArrowChosen(),VelocityOptChosen();
void VelocityChosen(),LayerChosen(),StrainChosen();
void DensityChosen(),GravityChosen();
void ElleVarChosen(),MeshVarChosen(),DeformVarChosen(),BndBoxChosen();
void ProfileMarkChosen(),StrnMrkChosen();
void ElleFileChosen();
void MeshParams();
void ReadLabel(),FindFont();
int  CreateColormap();
void CreateWorkingMessage(), ForceUpdate(), SetCursor();
void clear_menus();
void CreateFontlists();
void OpenDialog(), SaveDialog(), RefDialog(), CurrDialog();
void OpenElleDialog();
void PrintComments();
void Prfl2D_Dialog();
/* void ZoomDialog(); */
void  WorkingMsg(), ErrorMsg();
void warning_msg();
int Create_widgets(), Set_Colours();
/* int GetUserVals(); */
int GetDirSpec(), GetProfileVals(), UpdateMsgArea();
int Init_App(),Run_App();
void InitFont(),InitTitle(),Exit_App();
int Run_external_prog( char *command, char *msg, int len );
void read_ref_rec(),read_curr_rec(),sybflush_();
static void Syntax();
static void FileSaveSelectionOK(),FileOpenSelectionOK();
static void FileSelectionCancel();
/*static void ZoomPromptOK(); */
extern void clear_arrays(),CreateErrorMessage();
extern void LocateDialog();
extern void PrflptsDialog(),PrflvalsDialog(),CntrvalsDialog();
extern void PlotOptsDialog(),ZoomValsDialog();
extern void ChangeCell(), clear_X(), DeleteAllCells(), DeleteCell();
extern void ChangeCellConfig();
extern void toggled_uchar();

typedef struct {
    char *name;
    char *tag;
} FontData;

FontData labelfonts[] = {
    {"-*-helvetica-medium-r-*--*-80-*","TAG08H"},
    {"-*-helvetica-medium-r-*--*-100-*","TAG10H"},
    {"-*-helvetica-bold-r-*--*-120-*","TAG12H"},
    {"-*-helvetica-bold-r-*--*-140-*","TAG14H"},
    {"-*-helvetica-bold-r-*--*-180-*","TAG18H"},
    {"-*-helvetica-bold-r-*--*-240-*","TAG24H"},
    {"-*-symbol-medium-r-*--*-100-*","TAG10S"},
    {"-*-symbol-medium-r-*--*-120-*","TAG12S"},
    {"-*-symbol-medium-r-*--*-140-*","TAG14S"},
    {"-*-symbol-medium-r-*--*-180-*","TAG18S"},
    {"-*-symbol-medium-r-*--*-240-*","TAG24S"},
    {NULL}
};

        /* SUBMENUS */
/*
 * the menu labels should match the entries in log.h
 * (should be generated from log.h....)
 * Used when reading and writing log files
 */
MenuItemData fileopenopts[] = {
    { "Solution", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            FileOpenChosen, (XtPointer)SOLN, NULL, NO_VALUE },
    { "Log", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            FileOpenChosen, (XtPointer)LOG, NULL, NO_VALUE },
    {NULL},
};

MenuItemData filesaveopts[] = {
    { "Log", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            FileSaveChosen, (XtPointer)LOG, NULL, NO_VALUE },
    { "Log+Postscript", &xmPushButtonGadgetClass, NO_VALUE, 0,
       NO_VALUE, FileSaveChosen, (XtPointer)LOGPS, NULL, NO_VALUE },
/*
    { "Settings", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            FileSaveChosen, (XtPointer)STTGS, NULL, NO_VALUE },
*/
    {NULL},
};

MenuItemData strain_arrow_opts[] = {
    { "pstm", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainArrowChosen, (XtPointer)PSTR, NULL, NO_VALUE },
    { "pstd", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainArrowChosen, (XtPointer)SIGD, NULL, NO_VALUE },
    { "mssr", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainArrowChosen, (XtPointer)MSSR, NULL, NO_VALUE },
    { "ssft", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainArrowChosen, (XtPointer)SSFT, NULL, NO_VALUE },
    {NULL},
};

MenuItemData stress_arrow_opts[] = {
    { "taum", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressArrowChosen, (XtPointer)TAUM, NULL, NO_VALUE },
    { "taud", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressArrowChosen, (XtPointer)TAUD, NULL, NO_VALUE },
    { "sigd", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressArrowChosen, (XtPointer)SIGD, NULL, NO_VALUE },
    { "sigm", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressArrowChosen, (XtPointer)SIGM, NULL, NO_VALUE },
    {NULL},
};

MenuItemData strain_opts[] = {
    { "edxx", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)EDXX, NULL, NO_VALUE },
    { "edyy", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)EDYY, NULL, NO_VALUE },
    { "edzz", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)EDZZ, NULL, NO_VALUE },
    { "edxy", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)EDXY, NULL, NO_VALUE },
    { "psr1", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)PSR1, NULL, NO_VALUE },
    { "psr2", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)PSR2, NULL, NO_VALUE },
    { "msst", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)MSST, NULL, NO_VALUE },
    { "cang", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)CANG, NULL, NO_VALUE },
    { "tang", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)TANG, NULL, NO_VALUE },
    { "sang", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)SANG, NULL, NO_VALUE },
    { "dblc", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)DBLC, NULL, NO_VALUE },
    { "vort", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)VORT, NULL, NO_VALUE },
    { "ed2i", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)ED2I, NULL, NO_VALUE },
    { "vota", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrainOptChosen, (XtPointer)VOTA, NULL, NO_VALUE },
    /*{ "kvtn", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,*/
            /*StrainOptChosen, (XtPointer)KVTN, NULL, NO_VALUE },*/
    {NULL},
};

MenuItemData stress_opts[] = {
    { "taud", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)TAUD, NULL, NO_VALUE },
    { "taum", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)TAUM, NULL, NO_VALUE },
    { "taxx", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)TAXX, NULL, NO_VALUE },
    { "tayy", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)TAYY, NULL, NO_VALUE },
    { "tazz", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)TAZZ, NULL, NO_VALUE },
    { "taxy", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)TAXY, NULL, NO_VALUE },
    { "tau1", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)TAU1, NULL, NO_VALUE },
    { "tau2", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)TAU2, NULL, NO_VALUE },
    { "sixx", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)SIXX, NULL, NO_VALUE },
    { "siyy", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)SIYY, NULL, NO_VALUE },
    { "sizz", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)SIZZ, NULL, NO_VALUE },
    { "sig1", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)SIG1, NULL, NO_VALUE },
    { "sig2", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)SIG2, NULL, NO_VALUE },
    { "thdi", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)THDI, NULL, NO_VALUE },
    { "pres", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)PRES, NULL, NO_VALUE },
    { "brit", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)BRIT, NULL, NO_VALUE },
    { "bri2", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)BRI2, NULL, NO_VALUE },
    { "visc", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StressOptChosen, (XtPointer)VISC, NULL, NO_VALUE },
    {NULL},
};

MenuItemData velocity_opts[] = {
    { "Ux", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            VelocityOptChosen, (XtPointer)UX, NULL, NO_VALUE},
    { "Uy", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            VelocityOptChosen, (XtPointer)UY, NULL, NO_VALUE},
    { "Um", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            VelocityOptChosen, (XtPointer)UM, NULL, NO_VALUE},
    { "Ur", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            VelocityOptChosen, (XtPointer)UR, NULL, NO_VALUE},
    { "Uth", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            VelocityOptChosen, (XtPointer)UTH, NULL, NO_VALUE},
    {NULL},
};

MenuItemData deformopts[] = {
    { "ellp", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            DeformVarChosen, (XtPointer)ELLP, NULL, NO_VALUE },
    { "trel", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            DeformVarChosen, (XtPointer)TREL, NULL, NO_VALUE },
    { "tria", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            DeformVarChosen, (XtPointer)TRIA, NULL, NO_VALUE },
/*
    { "rota", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            DeformVarChosen, (XtPointer)ROTA, NULL, NO_VALUE },
    { "fdef", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            DeformVarChosen, (XtPointer)FDEF, NULL, NO_VALUE },
    { "tran", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            DeformVarChosen, (XtPointer)TRAN, NULL, NO_VALUE },
*/
    {NULL},
};

MenuItemData layer_opts[] = {
    { "thickness", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            LayerOptChosen, (XtPointer)THICKNESS, NULL, NO_VALUE},
    { "topog", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            LayerOptChosen, (XtPointer)TOPOG, NULL, NO_VALUE},
    { "gravPE", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            LayerOptChosen, (XtPointer)GRAVPE, NULL, NO_VALUE},
    {NULL},
};

MenuItemData gravityopts[] = {
    { "Bouguer", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            GravityOptChosen, (XtPointer)GRBOUG, NULL, NO_VALUE},
    { "free-air", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            GravityOptChosen, (XtPointer)GRFREE, NULL, NO_VALUE },
    { "topog", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
			GravityOptChosen, (XtPointer)TOPOG,  NULL, NO_VALUE },
    {NULL},
};

MenuItemData profile2Dopts[] = {
    { "Velocity", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, velocity_opts, PRFL2D },
    { "Strain", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, strain_opts, PRFL2D },
    { "Stress", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, stress_opts, PRFL2D },
    { "Layer", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, layer_opts, PRFL2D },
    { "Density", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            DensityChosen, (XtPointer)PRFL2D, NULL, NO_VALUE },
    {NULL},
};

MenuItemData profile1Dopts[] = {
    { "Velocity", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, velocity_opts, PRFL1D },
    { "Strain", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, strain_opts, PRFL1D },
    { "Stress", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, stress_opts, PRFL1D },
    { "Layer", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, layer_opts, PRFL1D },
    { "Density", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
             DensityChosen, (XtPointer)PRFL1D, NULL, NO_VALUE },
    { "Gravity(-Y)", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
             NULL, (XtPointer)0, gravityopts, PRFL1D },
    {NULL},
};

MenuItemData lgmeshopts[] = {
    { "elements", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            MeshVarChosen, (XtPointer)LGELMNT, NULL, NO_VALUE },
    { "boundary", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            MeshVarChosen, (XtPointer)LGBNDRY, NULL, NO_VALUE },
    {NULL},
};

MenuItemData meshopts[] = {
    { "elements", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            MeshVarChosen, (XtPointer)ELMNT, NULL, NO_VALUE },
    { "element+num", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            MeshVarChosen, (XtPointer)ELMNTNUM, NULL, NO_VALUE },
    { "boundary", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            MeshVarChosen, (XtPointer)BNDRY, NULL, NO_VALUE },
    { "internal_boundary", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            MeshVarChosen, (XtPointer)INTBND, NULL, NO_VALUE },
    { "internal_segments", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            MeshVarChosen, (XtPointer)INTLBNDRY, NULL, NO_VALUE },
    { "viscosity", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            MeshVarChosen, (XtPointer)VISCMSH, NULL, NO_VALUE },
    { "exponent", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            MeshVarChosen, (XtPointer)SEMSH, NULL, NO_VALUE },
    {NULL},
};

MenuItemData curropts[] = {
    { "Next", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            CurrChosen, (XtPointer)NEXT, NULL, NO_VALUE },
    { "Prev", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            CurrChosen, (XtPointer)PREV, NULL, NO_VALUE },
    { "Last", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            CurrChosen, (XtPointer)LAST, NULL, NO_VALUE },
    { "Specify", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            CurrChosen, (XtPointer)SPEC, NULL, NO_VALUE },
    {NULL},
};

MenuItemData refopts[] = {
    { "First", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            RefChosen, (XtPointer)FIRST, NULL, NO_VALUE },
    { "Specify", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            RefChosen, (XtPointer)SPEC, NULL, NO_VALUE },
    {NULL},
};

    /* number of entries matches USEFL_SIZE */
MenuItemData colouropts[] = {
    { "    ", &xmPushButtonWidgetClass, 0, 0, NO_VALUE,
            ColourChosen, (XtPointer)0, NULL, NO_VALUE },
    { "    ", &xmPushButtonWidgetClass, 1, 0, NO_VALUE,
            ColourChosen, (XtPointer)1, NULL, NO_VALUE },
    { "    ", &xmPushButtonWidgetClass, 2, 0, NO_VALUE,
            ColourChosen, (XtPointer)2, NULL, NO_VALUE },
    { "    ", &xmPushButtonWidgetClass, 3, 0, NO_VALUE,
            ColourChosen, (XtPointer)3, NULL, NO_VALUE },
    { "    ", &xmPushButtonWidgetClass, 4, 0, NO_VALUE,
            ColourChosen, (XtPointer)4, NULL, NO_VALUE },
    { "    ", &xmPushButtonWidgetClass, 5, 0, NO_VALUE,
            ColourChosen, (XtPointer)5, NULL, NO_VALUE },
    { "    ", &xmPushButtonWidgetClass, 6, 0, NO_VALUE,
            ColourChosen, (XtPointer)6, NULL, NO_VALUE },
    { "    ", &xmPushButtonWidgetClass, 7, 0, NO_VALUE,
            ColourChosen, (XtPointer)7, NULL, NO_VALUE },
    {NULL},
};

MenuItemData lineopts[] = {
    { "Solid", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            LineChosen, (XtPointer)SOLID, NULL, NO_VALUE },
    { "Dashed", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            LineChosen, (XtPointer)DASH, NULL, NO_VALUE },
    { "Dotted", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            LineChosen, (XtPointer)DOT, NULL, NO_VALUE },
    {NULL},
};

MenuItemData deleteopts[] = {
    { "Current Cell", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            DeleteChosen, (XtPointer)CURR_CELL, NULL, NO_VALUE },
    { "All Cells", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            DeleteChosen, (XtPointer)ALL_CELLS, NULL, NO_VALUE },
    {NULL},
};

MenuItemData profileopts[] = {
    { "1_D", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, profile1Dopts, PRFL1D},
    { "2_D", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)PRFL2D, profile2Dopts, NO_VALUE },
    { "Mark", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            ProfileMarkChosen, (XtPointer)PRFLMRK, NULL, NO_VALUE },
    {NULL},
};

        /* MENUS */
MenuItemData fileopts[] = {
    { "Open", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)OPEN, fileopenopts, NO_VALUE },
    { "Save", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)SAVE, filesaveopts, NO_VALUE },
    { "Sep", &xmSeparatorGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)SEP, NULL, NO_VALUE },
    { "Exit", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            ExitDialog, (XtPointer)EXIT, NULL, NO_VALUE },
    {NULL},
};

MenuItemData dataopts[] = {
    { "Current", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)NEXT, curropts, NO_VALUE },
    { "Reference", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)SPEC, refopts, NO_VALUE },
    {NULL},
};

MenuItemData xyplot_opts[] = {
    { "Mesh", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)MESH, meshopts, NO_VALUE },
    { "Deform", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)DEFORM, deformopts, NO_VALUE },
    { "BoundingBox", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            BndBoxChosen, (XtPointer)0, NULL, NO_VALUE },
    { "LGMesh", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)MESH, lgmeshopts, NO_VALUE },
    { "StrainMark", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            StrnMrkChosen, (XtPointer)0, NULL, NO_VALUE },
    { "Elle", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            ElleVarChosen, (XtPointer)ELLE, NULL, NO_VALUE },
    {NULL},
};

MenuItemData arrow_opts[] = {
    { "Velocity", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            VelocityArrowChosen, (XtPointer)ARROWS, NULL, NO_VALUE },
    { "Strain", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, strain_arrow_opts, ARROWS},
    { "Stress", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, stress_arrow_opts, ARROWS },
    { "Rotation", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            RotationArrowChosen, (XtPointer)ARROWS, NULL, NO_VALUE },
    {NULL},
};

MenuItemData contour_opts[] = {
    { "Velocity", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, velocity_opts, CNTRS },
    { "Strain", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, strain_opts, CNTRS },
    { "Stress", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, stress_opts, CNTRS },
    { "Layer", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, layer_opts, CNTRS },
    { "Density", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            DensityChosen, (XtPointer)CNTRS, NULL, NO_VALUE },
    {NULL},
};

    /* number of entries matches fontList data in default_resources */
MenuItemData fonts[] = {
    { "Helv 8", &xmPushButtonGadgetClass, NO_VALUE, 0, 0,
            FontChosen, (XtPointer)0, NULL, NO_VALUE },
    { "Helv10", &xmPushButtonGadgetClass, NO_VALUE, 0, 1,
            FontChosen, (XtPointer)1, NULL, NO_VALUE },
    { "Helv12", &xmPushButtonGadgetClass, NO_VALUE, 0, 2,
            FontChosen, (XtPointer)2, NULL, NO_VALUE },
    { "Helv14", &xmPushButtonGadgetClass, NO_VALUE, 0, 3,
            FontChosen, (XtPointer)3, NULL, NO_VALUE },
    { "Helv18", &xmPushButtonGadgetClass, NO_VALUE, 0, 4,
            FontChosen, (XtPointer)4, NULL, NO_VALUE },
    { "Helv24", &xmPushButtonGadgetClass, NO_VALUE, 0, 5,
            FontChosen, (XtPointer)5, NULL, NO_VALUE },
    { "Symbol10", &xmPushButtonGadgetClass, NO_VALUE, 0, 6,
            FontChosen, (XtPointer)6, NULL, NO_VALUE },
    { "Symbol12", &xmPushButtonGadgetClass, NO_VALUE, 0, 7,
            FontChosen, (XtPointer)7, NULL, NO_VALUE },
    { "Symbol14", &xmPushButtonGadgetClass, NO_VALUE, 0, 8,
            FontChosen, (XtPointer)8, NULL, NO_VALUE },
    { "Symbol18", &xmPushButtonGadgetClass, NO_VALUE, 0, 9,
            FontChosen, (XtPointer)9, NULL, NO_VALUE },
    { "Symbol24", &xmPushButtonGadgetClass, NO_VALUE, 0, 10,
            FontChosen, (XtPointer)10, NULL, NO_VALUE },
    {NULL},
};

MenuItemData labelopts[] = {
    { "Edit", &xmPushButtonGadgetClass, NO_VALUE, EDIT_INDX, NO_VALUE,
            LabelChosen, (XtPointer)0, NULL, NO_VALUE },
    { "Font", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer) 1, fonts, NO_VALUE },
    {NULL},
};

MenuItemData locateopts[] = {
    { "Next", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            LocateChosen, (XtPointer)NEXT, NULL, NO_VALUE },
    { "Prev", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            LocateChosen, (XtPointer)PREV, NULL, NO_VALUE },
    { "Specify", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            LocateChosen, (XtPointer)SPEC, NULL, NO_VALUE },
    {NULL},
};

MenuItemData optionopts[] = {
    { "Label", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, labelopts, NO_VALUE },
    { "Colour", &xmCascadeButtonGadgetClass, NO_VALUE, 10, NO_VALUE,
            NULL, (XtPointer)0, colouropts, NO_VALUE },
    { "Line", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, lineopts, NO_VALUE },
    { "Rescale", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            RescaleChosen, (XtPointer)0, NULL, NO_VALUE },
    { "Zoom", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            ZoomChosen, (XtPointer)&main_Menu, NULL, NO_VALUE },
    { "Delete", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, deleteopts, NO_VALUE },
    { "Plot", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            PlotOptsDialog, (XtPointer)&main_Menu, NULL, NO_VALUE },
    {NULL},
};

MenuItemData helpopts[] = {
    { "View", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            HelpChosen, (XtPointer)VIEW, NULL, NO_VALUE },
    {NULL},
};

MenuItemData mainopts[] = {
        /* Help button must be the last in the list */
    { "File", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, fileopts, NO_VALUE },
    { "Record", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, dataopts, NO_VALUE },
    { "XYPlot", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, xyplot_opts, NO_VALUE },
    { "Profile", &xmPushButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, profileopts, NO_VALUE },
    { "Arrow", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, arrow_opts, NO_VALUE },
    { "Contour", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, contour_opts, NO_VALUE },
    { "Locate", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, locateopts, NO_VALUE },
    { "Options", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, optionopts, NO_VALUE },
    { "Help", &xmCascadeButtonGadgetClass, NO_VALUE, 0, NO_VALUE,
            NULL, (XtPointer)0, helpopts, NO_VALUE },
    {NULL},
};

static void Syntax(argc, argv)
int argc;
char * argv[];
{

#if XY
    /* first argument is program name - skip that */
    for (i = 1; i < argc; i++) {
        if (!errs++) /* do first time through */
            fprintf(stderr, "sybil: command line option not understood:\n");
        fprintf(stderr, "option: %s\n", argv[i]);
    }
#endif

    fprintf(stderr, "Sybil understands Xt command line options as well as\n");
    fprintf(stderr, "these additional options:\n");
    fprintf(stderr, "Option   Valid Range\n");
    fprintf(stderr, "-pw      %d to %d    (specifies cells per row)\n",
        MINWIDTH, MAXWIDTH);
    fprintf(stderr, "-ph      %d to %d    (specifies cells per column)\n",
        MINHEIGHT,MAXHEIGHT);
    fprintf(stderr, "-i       name of log file \n");
    fprintf(stderr, "-h       print this message\n");
}
 
/*
 * Working Messages
 */
void WorkingMsg( display )
unsigned char display;
{
    char message[256];
    static Widget workingdialog;
    Arg args[6];
    int i=0;
    Dimension wdth;
    Position x,y;
    XmString xmstr;

    if (display==0 ) {
      if (workingdialog) XtUnmanageChild( workingdialog );
    }
    else  {
      CreateWorkingMessage(message);
#if XmVersion<1002
      xmstr=XmStringCreateSimple(message);
#else
      xmstr=XmStringCreateLocalized(message);
#endif
      if (!workingdialog) {
        XtSetArg(args[i],XmNmessageString,xmstr); i++;
        XtSetArg(args[i], XmNtitle, "STATUS"); i++;
        XtSetArg(args[i], XmNdefaultPosition, False); i++;
        XtVaGetValues(main_Menu, XmNx, &x, XmNwidth, &wdth, NULL);
        x += (Position)(wdth/2);
        XtSetArg(args[i], XmNx, x); i++;
        y = HeightOfScreen(XtScreen(main_Menu))/3;
        XtSetArg(args[i], XmNy, y); i++;
        XtSetArg(args[i], XmNdialogStyle, XmDIALOG_FULL_APPLICATION_MODAL); i++;
        workingdialog = XmCreateWorkingDialog(main_Menu,"working",args,i);
        XtUnmanageChild((Widget) XmMessageBoxGetChild(workingdialog,
                    XmDIALOG_OK_BUTTON));
        XtUnmanageChild((Widget) XmMessageBoxGetChild(workingdialog,
                    XmDIALOG_HELP_BUTTON));
        XtUnmanageChild((Widget) XmMessageBoxGetChild(workingdialog,
                    XmDIALOG_CANCEL_BUTTON));
      }
      else {
        XtVaSetValues(workingdialog,XmNmessageString,xmstr,NULL);
      }
      XtManageChild( workingdialog );
      XmStringFree(xmstr);
    }
}

void ForceUpdate(w)
Widget w;
{
    Widget shell, topshell;
    Window win, topwin;
    Display *dsply;
    XtAppContext cntxt = XtWidgetToApplicationContext(w);
    XWindowAttributes winattr;
    XEvent event;

    /* Find the shell */
    for (shell=w; !XtIsShell(shell); shell=XtParent(shell)) ;
    for (topshell=shell; !XtIsTopLevelShell(topshell);
                           topshell=XtParent(topshell)) ;
    /* Don't do anything unless both shells are realized */
    if (XtIsRealized(shell) && XtIsRealized(topshell)) {
        dsply = XtDisplay(topshell);
        win = XtWindow(shell);
        topwin = XtWindow(topshell);
        /* Wait for dialog to be mapped */
        while (XGetWindowAttributes(dsply,win,&winattr) &&
                    winattr.map_state != IsViewable ) {
            if (XGetWindowAttributes( dsply,topwin,&winattr) &&
                        winattr.map_state != IsViewable ) break;
            XtAppNextEvent( cntxt, &event );
            XtDispatchEvent(&event);
        }
    }
    /* next XSync() or XFlush() will get an expose event */
    XmUpdateDisplay(topshell);
}

void CreateWorkingMessage( message )
char *message;
{
    strcpy(message,"");
    switch(Plot_info.plot_type) {
    case ARROWS: strcpy(message,"Plotting arrows ");
                 strcat(message,Plot_info.variable);
                 break;
    case CNTRS:  strcpy(message,"Contouring ");
                 strcat(message,Plot_info.variable);
                 break;
    case PRFLMRK: strcpy(message,"Plotting profile line ");
                 break;
    case PRFL1D: strcpy(message,"Calculating 1_D profile ");
                 strcat(message,Plot_info.variable);
                 break;
    case PRFL2D: strcpy(message,"Calculating 2_D profile ");
                 strcat(message,Plot_info.variable);
                 break;
    case XYPLOT: strcpy(message,"Plotting ");
                 strcat(message,Plot_info.variable);
                 break;
    case PS    : strcpy(message,"Creating Postscript file ");
                 break;
    }
}

/*
 * Error Messages
 */
void ErrorMsg( msg, w, err_num )
char *msg;
Widget w;
int err_num;
{
    static Widget errdialog;
    Arg args[5];
    char message[SYB_FILENAME_MAX];
    int i=0;
    XmString xmstr;

    strcpy(message,"");
    strncpy(message,msg,sizeof(message)-21);
    CreateErrorMessage(message,err_num);
    xmstr=XmStringCreateSimple(message);
    if (w && XtIsRealized(topLevel)) {
        if (!errdialog) {
            XtSetArg(args[0],XmNmessageString,xmstr); i++;
               XtSetArg(args[i], XmNtitle, "ERROR MESSAGE"); i++;
                errdialog = XmCreateErrorDialog(w,"err",args,i);
            XtUnmanageChild((Widget) XmMessageBoxGetChild(errdialog,
                        XmDIALOG_HELP_BUTTON));
            XtUnmanageChild((Widget) XmMessageBoxGetChild(errdialog,
                        XmDIALOG_CANCEL_BUTTON));
        }
        else XtVaSetValues(errdialog,XmNmessageString,xmstr,NULL);
           XmStringFree(xmstr);
        XtManageChild( errdialog );
        XtPopup( XtParent( errdialog ), XtGrabNone );
    }
    else XtWarning(message);
}

/**********
 * FILE_OPT
 **********/

#if XY
/*
 * Save Log Edit callback function
 */
void Write_log(w, client_data, cbs)
Widget w;
XtPointer client_data;
XmSelectionBoxCallbackStruct *cbs;
{
    char *l, *logname;
    int err=0;
    FILE *log_fp;

    logname = (char *)client_data;
    XmStringGetLtoR(cbs->value, XmSTRING_DEFAULT_CHARSET, &l);
    strncpy(logname,l,MAX_LABEL_LEN-1);
    XtFree(l);

    if ((log_fp=fopen(logname,"w"))==NULL)
perror("sybil");
/*
                    ErrorMsg(logname,main_Menu,OPEN_ERR);
*/
    else {
        err = Write_plot_data(plotArea,log_fp);
        fclose( log_fp );
    }
}
#endif

static void FileSaveSelectionOK(w, client_data, cbs)
Widget w;
XtPointer client_data;
XmFileSelectionBoxCallbackStruct *cbs;/*always third param*/
{
    char *filename,*fname,command[SYB_FILENAME_MAX+256];
    char msg[MAX_LABEL_LEN+8];
    FILE *fp;
    int err=0,type=0,tmp;
    XtPointer ptr;

    if (!XmStringGetLtoR(cbs->value, XmFONTLIST_DEFAULT_TAG, &filename))
                return;
    XtVaGetValues(w, XmNuserData, &ptr, NULL);
    type = (long)ptr;

    if (filename[strlen(filename)-1]!='/') {
        if ( (Settings.filepath==0) && 
                (strstr(Plot_info.inp_file->fname,"FD.sols")==NULL) )
            ErrorMsg("Option filepath=0 but FD.sols not found ",
                       main_Menu,PATH_ERR);
            
        switch(type) {
        case LOG:   if ((fp=fopen(filename,"w"))==NULL)
                        ErrorMsg(filename,main_Menu,OPEN_ERR);
                    else {
                        err = Write_plot_data(plotArea,fp);
                        fclose( fp );
                        XtUnmanageChild(w);
                    }
                    break;
        case LOGPS: if ((fp=fopen(filename,"w"))==NULL)
                        ErrorMsg(filename,main_Menu,OPEN_ERR);
                    else {
                        SetCursor(main_Menu,Waiting);
                        err = Write_plot_data(plotArea,fp);
                        fclose( fp );
                        XtUnmanageChild(w);
                        /*if ((fname=strrchr(filename,'/'))==NULL)*/
                            fname=filename;
                        /*else fname++;*/
                        tmp = Plot_info.plot_type;
                        Plot_info.plot_type= PS;
                        snprintf(command,256+SYB_FILENAME_MAX,
                                 "%s -i %s &", "sybilps",fname);
                        strncat(msg,"sybilps: ",MAX_LABEL_LEN-1);
/*
                        err = Run_external_prog(command,msg,MAX_LABEL_LEN-1);
                        if (err!=0) {
                            ErrorMsg(msg,main_Menu,0);
                            err=0;
                        }
*/
                        int status = system(command);
                        int errnum = errno;
                        printf("%d %d \n",status, errnum);
                        if (status) {
                        if (errnum==EAGAIN || errnum==ENOENT)
                            ErrorMsg("Could not find sybilps",main_Menu,err);
                        else if (errnum)
                            ErrorMsg("Error running sybilps",main_Menu,err);
                        }
                        SetCursor(main_Menu,None);
                        Plot_info.plot_type= tmp;
                    }
                    break;
        default:    break;
        }
    }
    else ErrorMsg("The file selected is a directory",w,0);
    XtFree( filename );
}       /* FileSaveSelectionOK */

static void FileOpenSelectionOK(w, client_data, cbs)
Widget w;
XtPointer client_data;
XmFileSelectionBoxCallbackStruct *cbs;/*always third param*/
{
    unsigned char chng=0;
    char *filename,*fname,*buf;
    char input_str[SYB_FILENAME_MAX+1];
    int i,err=0,err_num=0,type=0;
    XtPointer ptr;
    XmString xmstr;
    input_options tmp;
    Dimension wdth,hgt,w1,h1;

    if (!XmStringGetLtoR(cbs->value, XmFONTLIST_DEFAULT_TAG, &filename))
                return;

 /* XtVaGetValues(w, XmNuserData, &type, NULL);  grief here */
    XtVaGetValues(w,XmNuserData,&ptr,NULL);
    type = (long)ptr;

    XtVaGetValues(w,XmNdirMask,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&buf);
    XmStringFree(xmstr);
    if (filename[strlen(filename)-1]!='/') {
        switch(type) {
        case SOLN:  if (Plot_info.inp_file!=NULL) {
                        if (Plot_info.inp_file->fp!=NULL) {
                            fclose(Plot_info.inp_file->fp);
                            clear_arrays();
                        }
                    }
                    fprintf(stdout,"\nOpening file: %s",filename);
                    if (!(err=file_open(&(Plot_info.inp_file),filename,
                                         SYB_FILENAME_MAX,"rb")))
                        err=CheckCurrentRec(Err_str,&err_num);
                    if (!err) {
                        if (Settings.plot_opts.label && Plot_info.title==NULL){
                            if ((fname=strrchr(Plot_info.inp_file->fname,'/'))
                                                                      ==NULL)
                                fname=Plot_info.inp_file->fname;
                            else fname++;
                            if (Settings.plot_opts.label) InitTitle( fname );
                        }
                        PrintComments();
                        UpdateMsgArea();
                        strncpy(Lastfilemask,buf,SYB_FILENAME_MAX);
                        XtUnmanageChild(w);
                    }
                    else ErrorMsg(filename,w,err_num);
                    break;
        case LOG:   if ((Plot_info.log_file_ptr=fopen( filename,"r" ))==NULL) 
                        err = OPEN_ERR;
                    if (!err) {
                        tmp = Settings;
                        err = initial_options(Plot_info.log_file_ptr,
                                                    &tmp,input_str);
                    }
                    if (!err) {
                        /*
                        if (tmp.rows!=Settings.rows||
                                          tmp.columns!=Settings.columns)
                         */
                            chng = 1;
                        if (tmp.orient!=Settings.orient) chng = 2;
                        Initial_Settings = Settings = tmp;
                        Settings.rescale = 1;
                        if (chng) {
                            if (chng>1) XtUnrealizeWidget(topLevel);
                            ChangeCellConfig(plotArea);
                            if (chng>1) {
                                XtVaGetValues(main_Menu,XtNheight,&h1,
                                                   NULL);
                                XtVaGetValues(plotArea,XtNwidth,&wdth,
                                                   XtNheight,&hgt,
                                                   NULL);
                                hgt += h1;
                                if (wdth<Minwidth) wdth = Minwidth;
                                XtMakeResizeRequest(topLevel,wdth,hgt,&w1,&h1);
                                XtRealizeWidget(topLevel);
                            }
                        }
                        /*
                         * these aren't set when initial opts read
                         */
                        setpencolor_(&Settings.fg);
                        setlinewidth_(&Settings.linewidth);
                        i=1; dashln_(&Settings.linestyle,&i);
                    }
                    if (!err) {
                        strncpy(Lastlogmask,buf,SYB_FILENAME_MAX);
                        XtUnmanageChild(w);
                        SetCursor(main_Menu,Waiting);
                        Plot_info.display_log = 1;
                        process_log_file(Plot_info.log_file_ptr,&Settings);
                        Plot_info.display_log = 0;
                        UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
                        UpdateMsgArea();
                        SetCursor(main_Menu,None);
                    }
                    else ErrorMsg(filename,w,err);
                    break;
        default:    break;
        }
    }
    else ErrorMsg("The file selected is a directory",w,0);
    XtFree( buf );
    XtFree( filename );
}       /* FileOpenSelectionOK */

/*
 * FileSelectionCancel() simply unmaps the file selection box.
 */
static void FileSelectionCancel(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    XtUnmanageChild(w);
}       /* FileSelectionCancel */

/*
 * Sets the string for the file mask to all files in
 * current directory/FD.sols if it exists else
 * current directory
 */
int GetDirSpec( path,file_type )
char *path;
int file_type;
{
    char ext[30];
    char *patherr;

    if (!(patherr = getcwd(path,SYB_FILENAME_MAX))) {
     /* fprintf(stderr,path); */
        fprintf(stderr,"%s\n",path);
        return(0);
    }
    switch(file_type) {
    case SOLN: strcpy(ext,"/");
               strcat(ext,SOLNFILESDIR);
               strcat(path,ext);
               if (access(path,R_OK)!=0) 
                 path[strlen(path)-strlen(SOLNFILESDIR)-1]='\0';
               strcpy(ext,"/*");
               strcat(path,ext);
               break;
    case LOGPS: /* fall through */
    case LOG : strcpy(ext,"/");
               strcat(ext,"*.log");
               strcat(path,ext);
               break;
    case ELLE: strcpy(ext,"/");
               strcat(ext,"*.elle");
               strcat(path,ext);
               break;
    default  : break;
    }
    return(1);
}

/* 
 * maps a file selection box
 */
void OpenDialog(title,file_type)
char *title;
int file_type;
{
    static Widget fileSelection;

    char *buf;
    XmString tmp1, tmp2=NULL;
    Arg args[5];
    int i=0;
    Cardinal argcount;

 /* fprintf(stdout,"in OpenDialog: %s %d \n",title,file_type); */

#if XmVersion<1002
    tmp1 = XmStringCreateSimple(title);
#else
    tmp1 = XmStringCreateLocalized(title);
#endif
    XtSetArg(args[i], XmNdialogTitle, tmp1); i++;
    XtSetArg(args[i], XmNuserData, file_type); i++;
    if (file_type==SOLN) buf = Lastfilemask;
    else buf = Lastlogmask;
    if (!fileSelection) {
      strcpy(Lastfilemask,"");
      strcpy(Lastlogmask,"");
      if (GetDirSpec( buf,file_type )) {
#if XmVersion<1002
          tmp2 = XmStringCreateSimple(buf);
#else
          tmp2 = XmStringCreateLocalized(buf);
#endif
          XtSetArg(args[i], XmNdirMask, tmp2); i++;
      }
      argcount=i;
      fileSelection= (Widget) XmCreateFileSelectionDialog(
            main_Menu,  /* parent widget */
            "fileSelection",  /* widget name */
            args,   /* argument list*/
            argcount   /* arglist size */
            );
      XtUnmanageChild((Widget) XmFileSelectionBoxGetChild(fileSelection,
                    XmDIALOG_HELP_BUTTON));
      XtAddCallback(fileSelection, XmNcancelCallback,
                    FileSelectionCancel, NULL);
      XtAddCallback(fileSelection, XmNokCallback, FileOpenSelectionOK, NULL);
    }
    else {
      if (strlen(buf)==0) GetDirSpec( buf,file_type );
#if XmVersion<1002
      tmp2 = XmStringCreateSimple(buf);
#else
      tmp2 = XmStringCreateLocalized(buf);
#endif
      XtSetArg(args[i], XmNdirMask, tmp2); i++;
      XtSetValues(fileSelection,args,i);
    }
    XmStringFree(tmp1); XmStringFree(tmp2);
    XtManageChild(fileSelection);
    XtPopup(XtParent(fileSelection), XtGrabNone);
}

/* 
 * maps a file selection box
 */
void SaveDialog(title,file_type)
char *title;
int file_type;
{
    static Widget fileSave;

    char buf[SYB_FILENAME_MAX+1];
    XmString tmp1, tmp2=NULL;
    Arg args[5];
    int i=0;

    tmp1 = XmStringCreateSimple(title);
    XtSetArg(args[i], XmNdialogTitle, tmp1); i++;
    if (GetDirSpec( buf,file_type )) {
#if XmVersion<1002
        tmp2 = XmStringCreateSimple(buf);
#else
        tmp2 = XmStringCreateLocalized(buf);
#endif
        XtSetArg(args[i], XmNdirMask, tmp2); i++;
    }
    XtSetArg(args[i], XmNuserData, file_type); i++;
    if (!fileSave) {
    fileSave= (Widget) XmCreateFileSelectionDialog(
            main_Menu,  /* parent widget */
            "fileSave",  /* widget name */
            args,   /* argument list*/
            i   /* arglist size */
            );
    XtUnmanageChild((Widget) XmFileSelectionBoxGetChild(fileSave,
                    XmDIALOG_HELP_BUTTON));
    XtAddCallback(fileSave, XmNcancelCallback,
                    FileSelectionCancel, NULL);
    XtAddCallback(fileSave, XmNokCallback, FileSaveSelectionOK, NULL);
    }
    else XtSetValues(fileSave,args,i);
    XmStringFree(tmp1); XmStringFree(tmp2);
    XtManageChild(fileSave);
    XtPopup(XtParent(fileSave), XtGrabNone);
}

void warning_msg(err_num,str)
int err_num;
char *str;
{
    ErrorMsg(str,main_Menu,err_num);
}

/*
 * function called when user confirms exit
 */
void Exit_Proc()
{
    Exit_App(0);
}

void Exit_App(err)
int err;
{
    if (Plot_info.inp_file!=NULL) {
        if (Plot_info.inp_file->fp!=NULL) fclose(Plot_info.inp_file->fp);
        free(Plot_info.inp_file);
    }
    if (Plot_info.elle_file!=NULL) {
        if (Plot_info.elle_file->fp!=NULL) fclose(Plot_info.elle_file->fp);
        free(Plot_info.elle_file);
    }
    if (Plot_info.title) {
        free(Plot_info.title);
        Plot_info.title = NULL;
    }
    clear_arrays();
    if (plotArea) clear_X(plotArea);
    exit(err); 
}

/*
 * exit button callback function
 */
void ExitDialog()
{ 
    Widget warndialog;
    Arg args[5];
    int i=0;
    XmString msg,label1,label2;

    label1 = XmStringCreateSimple("Yes");
    label2 = XmStringCreateSimple("No");
       msg=XmStringCreateSimple("Do you really want to quit?");
    XtSetArg(args[i],XmNmessageString,msg); i++;
    XtSetArg(args[i],XmNokLabelString,label1); i++;
    XtSetArg(args[i],XmNcancelLabelString,label2); i++;
       XtSetArg(args[i], XmNtitle, "SYBIL WARNING"); i++;
    warndialog = XmCreateWarningDialog(main_Menu,"warn",args,i);
    XtAddCallback( warndialog, XmNokCallback, Exit_Proc, NULL );
    XtUnmanageChild((Widget) XmMessageBoxGetChild(warndialog,
                    XmDIALOG_HELP_BUTTON));
       XmStringFree(msg);
       XmStringFree(label1);
       XmStringFree(label2);
    XtManageChild( warndialog );
    XtPopup( XtParent( warndialog ), XtGrabNone );
}

/*
 * File open button callback function
 */
void FileOpenChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    switch((long)client_data) {
    case SOLN: OpenDialog("Solution files",SOLN);
               break;
    case LOG:  OpenDialog("Log files",LOG);
               break;
    case STTGS:  
               break;
    default:   break;
    }
}

/*
 * File save button callback function
 */
void FileSaveChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{

    switch((long)client_data) {
    case LOG:  SaveDialog("Log files",LOG);
               break;
    case LOGPS:SaveDialog("Log files",LOGPS);
               break;
    case STTGS:  
               break;
    default:   break;
    }
}

/**********
 * DATA_IN
 **********/

/*
 * Reference record Edit callback function
 */
void read_ref_rec(w, p, cbs)
Widget w;
Widget p;
XmSelectionBoxCallbackStruct *cbs;
{
    char *s;
    int err=0, err_num=0;

    XmStringGetLtoR(cbs->value, XmSTRING_DEFAULT_CHARSET, &s);
    Plot_info.inp_file->ref_req = atoi(s);
    XtFree(s);
    /* try to read the requested record */
    SetCursor(main_Menu,Waiting);
    if ((err=CheckReferenceRec(Err_str,&err_num))) {
        SetCursor(main_Menu,None);
        ErrorMsg(Err_str,main_Menu,err_num);
    }
    UpdateMsgArea();
    SetCursor(main_Menu,None);
}

/*
 * maps a prompt dialog for reference record input
 */
void RefDialog()
{
    char s[5];

    XmString p = XmStringCreateSimple("Enter reference no.:");
    XmString d;
    Arg args[8];

    int i=0;
       XtSetArg(args[i], XmNselectionLabelString, p); i++;
       XtSetArg(args[i], XmNtitle, "REFERENCE"); i++;
    sprintf(s,"%2d",Plot_info.inp_file->ref_curr);
    d = XmStringCreateSimple(s);
       XtSetArg(args[i], XmNtextString, d); i++;
    if (!RefPrompt) {
        RefPrompt = XmCreatePromptDialog(main_Menu, "prompt", args, i);
        XtAddCallback(RefPrompt, XmNokCallback, read_ref_rec, 0);

        /* No help is available... */
        XtUnmanageChild((Widget) XmSelectionBoxGetChild(RefPrompt,
                    XmDIALOG_HELP_BUTTON));
    }
    else XtSetValues(RefPrompt,args,i);
        XmStringFree(p);
        XmStringFree(d); /* always destroy compound strings when done */


    XtManageChild(RefPrompt);
    XtPopup(XtParent(RefPrompt), XtGrabNone);
        /*
         * if using Motif1.2 - add code here
         * use XmNinitialFocus to specify that XmDIALOG_PROMPT
         * will have keyboard focus when dialog popped up
         */

}

/*
 * Reference button callback function
 */
void RefChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    int err=0, err_num=0;

    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL)
        ErrorMsg("No solution file open",main_Menu,0);
    else {
        switch((long)client_data) {
        case FIRST:Plot_info.inp_file->ref_req = 1;
                   /* try to read the requested record */
                   SetCursor(main_Menu,Waiting);
                   if ((err=CheckReferenceRec(Err_str,&err_num))) {
                       SetCursor(main_Menu,None);
                       ErrorMsg(Err_str,main_Menu,err_num);
                   }
                   UpdateMsgArea();
                   SetCursor(main_Menu,None);
                   break;
        case SPEC: RefDialog();
                   break;
        default:   break;
        }
    }
}


/*
 * Current record edit callback function
 */
void read_curr_rec(w, p, cbs)
Widget w;
Widget p;
XmSelectionBoxCallbackStruct *cbs;
{
    char *s;
    int err=0, err_num=0;

    XmStringGetLtoR(cbs->value, XmSTRING_DEFAULT_CHARSET, &s);
    Plot_info.inp_file->rec_req = atoi(s);
    XtFree(s);
    /* try to read the requested record */
    SetCursor(main_Menu,Waiting);
    if ((err=CheckCurrentRec(Err_str,&err_num))) {
        SetCursor(main_Menu,None);
        ErrorMsg(Err_str,main_Menu,err_num);
    }
    SetCursor(main_Menu,None);
    UpdateMsgArea();
}

/*
 * maps a prompt dialog for current record input
 */
void CurrDialog()
{
    char s[5];
    XmString d;

    XmString p = XmStringCreateSimple("Enter record no.:");
    Arg args[8];

    int i=0;
       XtSetArg(args[i], XmNselectionLabelString, p); i++;
       XtSetArg(args[i], XmNtitle, "CURRENT"); i++;
    sprintf(s,"%2d",Plot_info.inp_file->rec_curr);
    d = XmStringCreateSimple(s);
       XtSetArg(args[i], XmNtextString, d); i++;
    if (!CurrPrompt) {
           CurrPrompt = XmCreatePromptDialog(main_Menu, "prompt", args, i);
           XtAddCallback(CurrPrompt, XmNokCallback, read_curr_rec, 0);

        /* No help is available... */
        XtUnmanageChild((Widget) XmSelectionBoxGetChild(CurrPrompt,
                    XmDIALOG_HELP_BUTTON));
    }
    else XtSetValues(CurrPrompt,args,i);
       XmStringFree(p);
       XmStringFree(d); /* always destroy compound strings when done */


    XtManageChild(CurrPrompt);
    XtPopup(XtParent(CurrPrompt), XtGrabNone);
        /*
         * if using Motif1.2 - add code here
         * use XmNinitialFocus to specify that XmDIALOG_PROMPT
         * will have keyboard focus when dialog popped up
         */

}

/*
 * Current button callback function
 */
void CurrChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    int err=0, err_num=0;

    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL)
        ErrorMsg("No solution file open",main_Menu,0);
    else {
        switch((long)client_data) {
        case NEXT: Plot_info.inp_file->rec_req++;
                   /* try to read the requested record */
                   SetCursor(main_Menu,Waiting);
                   if ((err=CheckCurrentRec(Err_str,&err_num))) {
                       SetCursor(main_Menu,None);
                       ErrorMsg(Err_str,main_Menu,err_num);
                   }
                   UpdateMsgArea();
                   SetCursor(main_Menu,None);
                   break;
        case PREV: if (Plot_info.inp_file->rec_req-1>0) {
                       Plot_info.inp_file->rec_req--;
                       /* try to read the requested record */
                       SetCursor(main_Menu,Waiting);
                       if ((err=CheckCurrentRec(Err_str,&err_num))) {
                           SetCursor(main_Menu,None);
                           ErrorMsg(Err_str,main_Menu,err_num);
                       }
                       UpdateMsgArea();
                       SetCursor(main_Menu,None);
                   }
                   else ErrorMsg("Already reading first record",main_Menu,0);
                   break;
        case LAST: Plot_info.inp_file->rec_req=Plot_info.inp_file->rec_max;
                   /* try to read the requested record */
                   SetCursor(main_Menu,Waiting);
                   if ((err=CheckCurrentRec(Err_str,&err_num))) {
                      SetCursor(main_Menu,None);
                      ErrorMsg(Err_str,main_Menu,err_num);
                   }
                   UpdateMsgArea();
                   SetCursor(main_Menu,None);
                   break;
        case SPEC: CurrDialog();
                   break;
        default:   break;
        }
    }
}

/*************************
 * XY_PLOT ARROW  CONTOUR
 *************************/

/*
 * Arrow submenu callback function
 */
void RotationArrowChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    sprintf(Plot_info.variable,"%s","Rotation");
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.var_num= ROTATION;
    Plot_info.plot_type= ARROWS;
    Plot_info.plot_description= 0;
    RotationChosen(Plot_info.dflt_label1,Plot_info.var_num,
                   Plot_info.plot_type,Plot_info.plot_description);
}

void RotationChosen(name, option, plot, descrip)
char *name;
int option,plot, descrip;
{
    int err=0;

    SetCursor(main_Menu,Waiting);
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
        SetCursor(main_Menu,None);
        ErrorMsg("No solution file open",main_Menu,0);
    }
    else {
        UpdateMsgArea();
        if (plot==PRFL1D||plot==PRFL2D) {
            if (GetProfileVals(&descrip,Data_vars_int[NCOMP])==USER_CANCEL) {
                SetCursor(main_Menu,None);
                return;
            }
        }
        WorkingMsg(1);
        err = plot_rotation(Data_vars_int,Data_vars_fl,
                    Data_arrays_int,Data_arrays_fl,
                    name,option,plot,descrip,Contour_vals,Profile_vals,
                    Pwindo,Mesh,&Settings.fg,&Settings.linestyle
                    );
        SetCursor(main_Menu,None);
        WorkingMsg(0);
        if (err) ErrorMsg("Rotation ",main_Menu,err);
        else UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
    }
}

/*
 * Velocity contour submenu callback function
 */
void VelocityOptChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int len, plot_type, cancel=0;
    int choice = (long) client_data;
    XtPointer ptr;
    
/*
 * choice is the number used by the f77 routines - not button no.
 */
    XtVaGetValues(w,XmNuserData,&ptr,NULL);
    plot_type = (long)ptr;
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    len = sizeof(Plot_info.variable) - 1 - strlen(var_name);
    sprintf(Plot_info.variable,"%s%s%s","Velocity",".",var_name);
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.var_num= choice;
    Plot_info.plot_type= plot_type;
    switch(plot_type) {
    case CNTRS: Plot_info.plot_description= Settings.plot_opts.contour_plot;
                break;
    case PRFL1D:Plot_info.plot_description= DIM1;
                break;
    case PRFL2D:if (Plot_info.plot_description!=DIM2_X)
                    Plot_info.plot_description= DIM2_Y;
                break;
        /* catch garbage val  - FIX*/
    default:    break;
    }
    if (!cancel) VelocityChosen(Plot_info.dflt_label1,Plot_info.var_num,
                   Plot_info.plot_type,Plot_info.plot_description);
}

/*
 * Arrow submenu callback function
 */
void VelocityArrowChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int len=0;
    
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name); /* var_name alloc if ok */
    XmStringFree(xmstr);
    len = strlen(var_name);
    if (len > MAXNAME-1) len = MAXNAME-1;
/*
    while (strlen(var_name)<4) {
        strcat(var_name," ");
    }
    sprintf(Plot_info.variable,"%s%s%s","Velocity",".",var_name);
*/
    sprintf(Plot_info.variable,"%s","Velocity");
    strncpy(Plot_info.dflt_label1,var_name, len);
    XtFree(var_name);
    Plot_info.var_num= SYB_VELOCITY;
    Plot_info.plot_type= ARROWS;
    Plot_info.plot_description= 0;
    VelocityChosen(Plot_info.dflt_label1,Plot_info.var_num,
                   Plot_info.plot_type,Plot_info.plot_description);
}

void VelocityChosen(name, option, plot, descrip)
char *name;
int option,plot,descrip;
{
    int err=0;

    SetCursor(main_Menu,Waiting);
/*  fprintf(stdout," entered VelocityChosen, option = %i\n",option); */
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
        SetCursor(main_Menu,None);
        ErrorMsg("No solution file open",main_Menu,0);
    }
    else {
        UpdateMsgArea();
        if (plot==PRFL1D||plot==PRFL2D) {
            if (GetProfileVals(&descrip,Data_vars_int[NCOMP])==USER_CANCEL) {
                SetCursor(main_Menu,None);
                return;
            }
        }
        WorkingMsg(1);
        err = plot_velocity(Data_vars_int,Data_vars_fl,
                    Data_arrays_int,Data_arrays_fl,
                    name,option,plot,descrip,Contour_vals,Profile_vals,
                    Pwindo,Mesh,&Settings.fg,&Settings.linestyle
                    );
        SetCursor(main_Menu,None);
        WorkingMsg(0);
        if (err) ErrorMsg("Velocity ",main_Menu,err);
        else UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
    }
}

/*
 * Strain submenu callback function
 */
void StrainArrowChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int len;
    int choice = (long) client_data;
    
/*
 * choice is the number used by the f77 routines - not button no.
 */
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    len = sizeof(Plot_info.variable) - 1 - strlen(var_name);
    sprintf(Plot_info.variable,"%s%s%s","Strain",".",var_name);
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.var_num= choice;
    Plot_info.plot_type= ARROWS;
    Plot_info.plot_description= 0;
    StrainChosen(Plot_info.dflt_label1,Plot_info.var_num,
                 Plot_info.plot_type,&Plot_info.plot_description);
}

/*
 * Stress submenu callback function
 */
void StressArrowChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int len;
    int choice = (long) client_data;
    
/*
 * choice is the number used by the f77 routines - not button no.
 */
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    len = sizeof(Plot_info.variable) - 1 - strlen(var_name);
    sprintf(Plot_info.variable,"%s%s%s","Stress",".",var_name);
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.var_num= choice;
    Plot_info.plot_type= ARROWS;
    Plot_info.plot_description= 0;
    StrainChosen(Plot_info.dflt_label1,Plot_info.var_num,
                 Plot_info.plot_type,&Plot_info.plot_description);
}

/*
 * Stress submenu callback function
 */
void StressOptChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int len, plot_type, cancel=0;
    int choice = (long) client_data;
    XtPointer ptr;
/*
 * choice is the number used by the f77 routines - not button no.
 */
    XtVaGetValues(w,XmNuserData,&ptr,NULL);
    plot_type = (long)ptr;
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    len = sizeof(Plot_info.variable) - 1 - strlen(var_name);
    sprintf(Plot_info.variable,"%s%s%s","Stress",".",var_name);
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.var_num= choice;
    Plot_info.plot_type= plot_type;
    switch(plot_type) {
    case CNTRS: Plot_info.plot_description= Settings.plot_opts.contour_plot;
                break;
    case PRFL1D:Plot_info.plot_description= DIM1;
                break;
    case PRFL2D:if (Plot_info.plot_description!=DIM2_X)
                    Plot_info.plot_description= DIM2_Y;
                break;
        /* catch garbage val  - FIX*/
    default:    break;
    }
    if (!cancel) StrainChosen(Plot_info.dflt_label1,Plot_info.var_num,
                   Plot_info.plot_type,&Plot_info.plot_description);
}

/*
 * Strain submenu callback function
 */
void StrainOptChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int len, plot_type, cancel=0;
    int choice = (long) client_data;
    XtPointer ptr;
    
/*
 * choice is the number used by the f77 routines - not button no.
 */
    XtVaGetValues(w,XmNuserData,&ptr,NULL);
    plot_type = (long)ptr;
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    len = sizeof(Plot_info.variable) - 1 - strlen(var_name);
    sprintf(Plot_info.variable,"%s%s%s","Strain",".",var_name);
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.var_num= choice;
    Plot_info.plot_type= plot_type;
    switch(plot_type) {
    case CNTRS: Plot_info.plot_description= Settings.plot_opts.contour_plot;
                break;
    case PRFL1D:Plot_info.plot_description= DIM1;
                break;
    case PRFL2D:if (Plot_info.plot_description!=DIM2_X)
                    Plot_info.plot_description= DIM2_Y;
                break;
        /* catch garbage val  - FIX*/
    default:    break;
    }
    if (!cancel) StrainChosen(Plot_info.dflt_label1,Plot_info.var_num,
                   Plot_info.plot_type,&Plot_info.plot_description);
}

/*
 * Strain submenu callback function
 */
void StrainChosen(name, option, plot, descrip)
char *name;
int option,plot,*descrip;
{ 
    int err=0;

    SetCursor(main_Menu,Waiting);
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
        SetCursor(main_Menu,None);
        ErrorMsg("No solution file open",main_Menu,0);
    }
    else {
        UpdateMsgArea();
        if (plot==PRFL1D||plot==PRFL2D) {
            if (GetProfileVals(descrip,Data_vars_int[NCOMP])==USER_CANCEL) {
                SetCursor(main_Menu,None);
                return;
            }
        }
        WorkingMsg(1);
        if ((err=plot_strain(Data_vars_int,Data_vars_fl,
                            Data_arrays_int,Data_arrays_fl,
                            Pwindo,Mesh,option,Contour_vals,Profile_vals,name,
                            plot,*descrip,&Settings.fg,&Settings.linestyle,
                            Arrow_colours
                            ))) {
            SetCursor(main_Menu,None);
            if (err) ErrorMsg("strain",main_Menu,err);
        }
        else {
            UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
            SetCursor(main_Menu,None);
        }
        WorkingMsg(0);
    }
}

/*
 * LGMesh and Mesh submenu callback function
 */
void MeshVarChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int err=0;
    int len, hgt=8;
    int choice = (long) client_data;
    fontdat curr_font;
    
    SetCursor(main_Menu,Waiting);
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    len = sizeof(Plot_info.variable) - 1 - strlen(var_name);
     /*
      * these should pull out parent button's label and be tied 
      * to structures in log.h
      */
    if (choice==LGELMNT||choice==LGBNDRY)
        sprintf(Plot_info.variable,"%s%s%s","LGMesh",".",var_name);
    else sprintf(Plot_info.variable,"%s%s%s","Mesh",".",var_name);
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    switch(choice) {
    case ELMNT  : if (Settings.plot_opts.dble) choice = ELMNTDBL;
                  break;
    case BNDRY  : if (Settings.plot_opts.dble) choice = BNDRYDBL;
                  break;
    case INTBND : if (Settings.plot_opts.dble) choice = INTBNDDBL;
                  break;
    case VISCMSH: if (Settings.plot_opts.dble) choice = VISCDBL;
                  break;
    case SEMSH: if (Settings.plot_opts.dble) choice = SEDBL;
                  break;
    default:      break;
    }
    Plot_info.var_num= choice;
    Plot_info.plot_type= XYPLOT;
    Plot_info.plot_description= 0;
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
        SetCursor(main_Menu,None);
        ErrorMsg("No solution file open",main_Menu,0);
    }
    else if ((choice==INTBND || choice==INTBNDDBL) &&
					Data_vars_int[IMREG]==0) {
		SetCursor(main_Menu,None);
        ErrorMsg("Region numbers not defined in file,\n internal boundaries cannot be plotted",main_Menu,0);
    }
    else {
        UpdateMsgArea();
        if (choice==ELMNTNUM) {
            /* set font to small Helvetica for element labels */
            curr_font = Plot_info.font;
            len = strlen("Helvetica");
            setfont_("Helvetica",&len,&hgt);
        }
        /*WorkingMsg(1);*/
        err=plot_mesh(choice,Data_vars_int,Data_vars_fl,
                        Data_arrays_int,Data_arrays_fl,
                        Pwindo,Mesh,&Settings.fg
                        );
        if (choice==ELMNTNUM) Plot_info.font = curr_font;
        SetCursor(main_Menu,None);
        if (err) ErrorMsg(Err_str,main_Menu,err);
        else UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
        /*WorkingMsg(0);*/
    }
}

/*
 * Elle submenu callback function
 */
void ElleVarChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int len;
    int choice = (long) client_data;
    
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    len = sizeof(Plot_info.variable) - 1 - strlen(var_name);
     /*
      * these should pull out parent button's label and be tied 
      * to structures in log.h
      */
    sprintf(Plot_info.variable,"%s","Elle");
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.var_num= choice;
    Plot_info.plot_type= XYPLOT;
    Plot_info.plot_description= 0;
    UpdateMsgArea();
    OpenElleDialog(main_Menu,"Elle files",ELLE);
}

void ElleFileChosen()
{
    int err=0;

    if (Plot_info.elle_file->fp!=NULL) {
        SetCursor(main_Menu,Waiting);
        /*WorkingMsg(1);*/
        err=plot_elle_regions(Pwindo);
        if (err) {
            SetCursor(main_Menu,None);
            ErrorMsg(Err_str,main_Menu,err);
        }
        else {
            UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
            SetCursor(main_Menu,None);
        }
        /*WorkingMsg(0);*/
    }
}
/*
 * Strain markers submenu callback function
 */
void StrnMrkChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    sprintf(Plot_info.variable,"%s","StrainMark");
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.plot_type= XYPLOT;
    Plot_info.plot_description= 0;
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
        ErrorMsg("No solution file open",main_Menu,0);
    }
    else {
        WorkingMsg(1);
        plot_strain_markers(Data_vars_int,Data_vars_fl,
                                        Data_arrays_fl,Pwindo);
        UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
        WorkingMsg(0);
    }
}

/*
 * Bounding box submenu callback function
 */
void BndBoxChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
/*
    sprintf(Plot_info.variable,"%s%s%s","Velocity",".",var_name);
*/
    sprintf(Plot_info.variable,"%s","BndBox");
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.plot_type= XYPLOT;
    Plot_info.plot_description= 0;
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
        ErrorMsg("No solution file open",main_Menu,0);
    }
    else {
        WorkingMsg(1);
        DrawBndBox(Pwindo,Profile_vals);
        UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
        WorkingMsg(0);
    }
}


/*
 * Deform submenu callback function
 */
void DeformVarChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int err=0;
    int len;
    int choice = (long) client_data;
    
    SetCursor(main_Menu,Waiting);
/*
 * choice is the number used by the f77 routines - not button no.
 */
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    len = sizeof(Plot_info.variable) - 1 - strlen(var_name);
    sprintf(Plot_info.variable,"%s%s%s","Deform",".",var_name);
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.var_num= choice;
    Plot_info.plot_type= XYPLOT;
    Plot_info.plot_description = 0;
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
        err=1;
        SetCursor(main_Menu,None);
        ErrorMsg("No solution file open",main_Menu,0);
    }
    if (!err) {
        if (Plot_info.inp_file->ref_curr==0) { 
            err=1;
            SetCursor(main_Menu,None);
            ErrorMsg("No reference record selected",main_Menu,0);
        }
    }
    if (!err) {
        UpdateMsgArea();
        WorkingMsg(1);
        if ((err=plot_deformation(Data_vars_int,Data_vars_fl,
                Data_arrays_int,Data_arrays_fl,
                Pwindo,Mesh,choice,Plot_info.plot_description,
                Contour_vals,Profile_vals,Plot_info.dflt_label1,
                &Settings.fg
                ))) {
            SetCursor(main_Menu,None);
            ErrorMsg("deform",main_Menu,err);
        }
        else {
            UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
        }
        WorkingMsg(0);
    }
    SetCursor(main_Menu,None);
}

/*
 * Layer submenu button callback function
 */
void LayerOptChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int len, plot_type, cancel=0;
    int choice = (long) client_data;
    XtPointer ptr;
    
    XtVaGetValues(w,XmNuserData,&ptr,NULL);
    plot_type = (long)ptr;
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    len = sizeof(Plot_info.variable) - 1 - strlen(var_name);
    sprintf(Plot_info.variable,"%s%s%s","Layer",".",var_name);
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.var_num= choice;
    Plot_info.plot_type= plot_type;
    switch(plot_type) {
    case CNTRS: Plot_info.plot_description= Settings.plot_opts.contour_plot;
                break;
    case PRFL1D:Plot_info.plot_description= DIM1;
                break;
    case PRFL2D:if (Plot_info.plot_description!=DIM2_X)
                    Plot_info.plot_description= DIM2_Y;
                break;
        /* catch garbage val  - FIX*/
    default:    break;
    }
    if (!cancel) LayerChosen(Plot_info.dflt_label1,Plot_info.var_num,
                   Plot_info.plot_type,&Plot_info.plot_description);
}

/*
 * Layer submenu callback function
 */
void LayerChosen(name, option, plot, descrip)
char *name;
int option,plot,*descrip;
{
    int err=0;

    SetCursor(main_Menu,Waiting);
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
        SetCursor(main_Menu,None);
        ErrorMsg("No solution file open",main_Menu,0);
    }
    else {
        if (plot==PRFL1D||plot==PRFL2D) {
            if (GetProfileVals(descrip,Data_vars_int[NCOMP])==USER_CANCEL) {
            SetCursor(main_Menu,None);
            return;
            }
        }
        /*if (Contour_vals[CNTR_LVL]==0.0)  Contour_vals[CNTR_LVL]=0.35;*/
        WorkingMsg(1);
        if ((err = plot_layer(Data_vars_int,Data_vars_fl,
                    Data_arrays_int,Data_arrays_fl,
                    Pwindo,Mesh,name,plot,*descrip,
                    Contour_vals,Profile_vals,option,
                    &Settings.fg,&Settings.linestyle
                    ))) {
            SetCursor(main_Menu,None);
            if (err) ErrorMsg("Layer",main_Menu,err);
        }
        else {
            UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
        }
        SetCursor(main_Menu,None);
        WorkingMsg(0);
    }
}

#if XY
#endif
/*
 * Density submenu button callback function
 */
void DensityChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int choice = (long) client_data;
    int err=0;

    
    if (XtIsRealized(main_Menu)) {
        XDefineCursor(XtDisplay(topLevel),
                            XtWindow(XtParent(main_Menu)),Waiting);
        XFlush( XtDisplay(topLevel));
    }
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    sprintf(Plot_info.variable,"%s","Density");
    XmStringFree(xmstr);
    XtFree(var_name);
    Plot_info.plot_type= choice;
    switch(choice) {
    case CNTRS: Plot_info.plot_description= Settings.plot_opts.contour_plot;
                break;
    case PRFL1D:Plot_info.plot_description= DIM1;
                break;
    case PRFL2D:if (Plot_info.plot_description!=DIM2_X)
                    Plot_info.plot_description= DIM2_Y;
                break;
        /* catch garbage val  - FIX*/
    default:    break;
    }
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
        SetCursor(main_Menu,None);
        ErrorMsg("No solution file open",main_Menu,0);
    }
    else {
        if (choice==PRFL1D||choice==PRFL2D) {
            if (GetProfileVals(&Plot_info.plot_description,
                                Data_vars_int[NCOMP])==USER_CANCEL) {
                SetCursor(main_Menu,None);
                return;
            }
        }
        WorkingMsg(1);
        if ((err = plot_density(Data_vars_int,Data_vars_fl,
                    Data_arrays_int,Data_arrays_fl,
                    Pwindo,Mesh,Plot_info.variable,Plot_info.plot_type,
                    Plot_info.plot_description,Contour_vals,Profile_vals,
                    &Settings.fg,&Settings.linestyle
                    ))) {
            SetCursor(main_Menu,None);
            if (err) ErrorMsg("Density",main_Menu,err);
        }
        else {
            UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
        }
        SetCursor(main_Menu,None);
        WorkingMsg(0);
    }
}

/*
 * Gravity(-Y) submenu button callback function
 */
void GravityOptChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int len, plot_type, cancel=0;
    int choice = (long) client_data;
    XtPointer ptr;
    
    XtVaGetValues(w,XmNuserData,&ptr,NULL);
    plot_type = (long)ptr;
    XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
    XmStringGetLtoR(xmstr,"",&var_name);
    XmStringFree(xmstr);
    len = sizeof(Plot_info.variable) - 1 - strlen(var_name);
    sprintf(Plot_info.variable,"%s%s%s","Gravity(-Y)",".",var_name);
    strncpy(Plot_info.dflt_label1,var_name, MAXNAME-1);
    XtFree(var_name);
    Plot_info.var_num= choice;
    Plot_info.plot_type= plot_type;
    switch(plot_type) {
    case CNTRS: Plot_info.plot_description= Settings.plot_opts.contour_plot;
                break;
    case PRFL1D:Plot_info.plot_description= DIM1;
                break;
    case PRFL2D:if (Plot_info.plot_description!=DIM2_X)
                    Plot_info.plot_description= DIM2_Y;
                break;
        /* catch garbage val  - FIX*/
    default:    break;
    }
    if (!cancel) GravityChosen(Plot_info.dflt_label1,Plot_info.var_num,
                   Plot_info.plot_type,&Plot_info.plot_description);
}

/*
 * Gravity(-Y) submenu callback function
 */
void GravityChosen(name, option, plot, descrip)
char *name;
int option,plot,*descrip;
{
    int err=0;
    float zmeas;

    SetCursor(main_Menu,Waiting);
    if (Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) {
        SetCursor(main_Menu,None);
        ErrorMsg("No solution file open",main_Menu,0);
    }
    else if (Data_vars_int[IMREG]==0) {
		SetCursor(main_Menu,None);
        ErrorMsg("Region numbers not defined in file,\n gravity calculation not available",main_Menu,0);
    }
    else {
        if (plot==PRFL1D||plot==PRFL2D) {
		    if(option==TOPOG)zmeas=1.0;
		    if(option==GRBOUG||option==GRFREE)zmeas=1.05;
		    if(Profile_vals[PRFL_X1]==0)Profile_vals[PRFL_X1]=Pwindo[XCMIN];
		    if(Profile_vals[PRFL_X2]==0)Profile_vals[PRFL_X2]=Pwindo[XCMAX];
		    Profile_vals[PRFL_Y1]=zmeas;
		    Profile_vals[PRFL_Y2]=zmeas;
            if (GetProfileVals(descrip,Data_vars_int[NCOMP])==USER_CANCEL) {
            SetCursor(main_Menu,None);
            return;
            }
        }
        WorkingMsg(1);
        if ((err = plot_gravity(Data_vars_int,Data_vars_fl,
                    Data_arrays_int,Data_arrays_fl,
                    Pwindo,Mesh,plot,*descrip,
                    Contour_vals,Profile_vals,
                    Plot_info.dflt_label1,Plot_info.var_num,
                    &Settings.fg,&Settings.linestyle
                    ))) {
            SetCursor(main_Menu,None);
            if (err) ErrorMsg("Gravity(-Y)",main_Menu,err);
        }
        else {
            UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
        }
        SetCursor(main_Menu,None);
        WorkingMsg(0);
    }
}
/**********
 * PROFILE
 **********/

/*
 * Profile Mark callback function
 */
void ProfileMarkChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{ 
    char *var_name;
    XmString xmstr;
    int choice = (long) client_data;

    SetCursor(main_Menu,Waiting);
    if ((Plot_info.inp_file==NULL||Plot_info.inp_file->fp==NULL) &&
       (Plot_info.elle_file==NULL||Plot_info.elle_file->fp==NULL)) {
        SetCursor(main_Menu,None);
        ErrorMsg("No solution file open",main_Menu,0);
    }
    else {
        XtVaGetValues(w,XmNlabelString,&xmstr,NULL);
        XmStringGetLtoR(xmstr,"",&var_name);
        sprintf(Plot_info.variable,"%s","Mark");
        XmStringFree(xmstr);
        XtFree(var_name);
        Plot_info.plot_type= choice;

        if (GetProfileVals(&Plot_info.plot_type,Data_vars_int[NCOMP])!=
             USER_CANCEL) {
            DrawProfileLine(Profile_vals,Pwindo,Data_vars_int[NCOMP]);
        }
    }
    SetCursor(main_Menu,None);
}

/**********
 * LABEL
 **********/

/*
 * Label Edit callback function
void ReadLabel(w, client_data, cbs)
Widget w;
XtPointer client_data;
XmSelectionBoxCallbackStruct *cbs;
{
    char *l, *str;

    str = (char *)client_data;
    XmStringGetLtoR(cbs->value, XmSTRING_DEFAULT_CHARSET, &l);
    strncpy(str,l,MAX_LABEL_LEN-1);
    XtFree(l);
}
 */
void ReadLabel(w,dest_str)
Widget w;
char *dest_str;
{
    char *str;

    str = XmTextFieldGetString(XmSelectionBoxGetChild(w,XmDIALOG_TEXT));
    strncpy(dest_str,str,MAX_LABEL_LEN-1);
    XtFree(str);
}

/*
 * maps a prompt dialog for label text input
 */
void LabelDialog(title,prompt,dest,callbk)
char *title,*prompt,*dest;
void (*callbk);
{
    Widget toggle,toggle_box;
    XmString p = XmStringCreateLtoR(prompt,XmFONTLIST_DEFAULT_TAG);
#if XmVersion<1002
    XmString d = XmStringCreateSimple(dest);
    XmString ok_str = XmStringCreateSimple("Done");
#else
    XmString d = XmStringCreateLocalized(dest);
    XmString ok_str = XmStringCreateLocalized("Done");
#endif
    Arg args[8];

    int i=0;
    XtSetArg(args[i], XmNselectionLabelString, p); i++;
    XtSetArg(args[i], XmNtitle, title); i++;
    XtSetArg(args[i], XmNtextString, d); i++;
    XtSetArg(args[i], XmNokLabelString, ok_str); i++;
    if (!LabelPrompt) {
        LabelPrompt = XmCreatePromptDialog(main_Menu, "prompt", args, i);

        if (callbk!=NULL) XtAddCallback(LabelPrompt,XmNokCallback,
                    (XtCallbackProc)callbk,dest);
        /* No help is available... */
        XtUnmanageChild((Widget) XmSelectionBoxGetChild(LabelPrompt,
                    XmDIALOG_HELP_BUTTON));
        /* Don't need cancel button... */
        XtUnmanageChild((Widget) XmSelectionBoxGetChild(LabelPrompt,
                    XmDIALOG_CANCEL_BUTTON));

        XtVaSetValues((Widget) XmSelectionBoxGetChild(LabelPrompt,
                        XmDIALOG_TEXT),XmNfontList,Curr_fontlist,NULL);
        /*
         * create ToggleButton gadget for setting bg option
         */
        toggle_box=XtVaCreateManagedWidget("toggle",
            xmRowColumnWidgetClass, LabelPrompt, NULL );
            toggle = XtVaCreateManagedWidget("Fill background",
                xmToggleButtonGadgetClass, toggle_box,
                XmNuserData, (XtPointer)&Settings.text_bg,
                NULL );

            XmToggleButtonGadgetSetState(toggle,
                           (Settings.text_bg==1), 0);
            XtAddCallback( toggle,XmNvalueChangedCallback,
                            toggled_uchar, (XtPointer)&Settings.text_bg );
    }
    /*else XtSetValues(LabelPrompt,args,i);*/
    XmStringFree(p);
    XmStringFree(ok_str);
    XmStringFree(d); /* always destroy compound strings when done */

    XtManageChild(LabelPrompt);
    XtPopup(XtParent(LabelPrompt), XtGrabNone);
}

/*
 * maps a prompt dialog for zoom input
 */
/*   This routine replaced by call to GetUserVals
void ZoomDialog(w,title,prompt,curr,retval)
Widget w;
char *title,*prompt;
int curr;
float *retval;
{
    XmString p = XmStringCreateLtoR(prompt,XmFONTLIST_DEFAULT_TAG);
    Arg args[4];

    int i=0;
    XtSetArg(args[i], XmNselectionLabelString, p); i++;
    XtSetArg(args[i], XmNtitle, title); i++;
    if (!ZoomPrompt) {
        ZoomPrompt = XmCreatePromptDialog(w, "zoom", args, i);

        XtUnmanageChild((Widget) XmSelectionBoxGetChild(ZoomPrompt,
                    XmDIALOG_HELP_BUTTON));
        XtAddCallback(ZoomPrompt, XmNokCallback, ZoomPromptOK,
NULL);

    }
    XmStringFree(p);

    XtManageChild(ZoomPrompt);
    XtPopup(XtParent(ZoomPrompt), XtGrabNone);
} */

/*
 * Font menu button callback function
 * UPDATE for Motif 1.2
 */
void FontChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    char *tag;
    XmFontContext cntxt;
    XmFontType type_return;
    XmFontListEntry entry;
    XtPointer tmp;

    XtVaGetValues( w, XmNfontList, &Curr_fontlist, NULL );
    if (!XmFontListInitFontContext(&cntxt,Curr_fontlist)) {
        XtWarning("Could not allocate font context");
    }
    else {
        entry = XmFontListNextEntry( cntxt );
        tag = XmFontListEntryGetTag( entry );
        tmp = XmFontListEntryGetFont(entry,&type_return);
        if (type_return==XmFONT_IS_FONT)
            Plot_info.font = (XFontStruct *) tmp;
        if (tag[strlen(tag)-1]=='H')
            strcpy(Plot_info.fontname,"Helvetica");
        else if (tag[strlen(tag)-1]=='S')
            strcpy(Plot_info.fontname,"Symbol");
        XtFree(tag);
        XmFontListFreeFontContext(cntxt);
        
        /* set the label widget's font */
        if (LabelPrompt) 
            XtVaSetValues((Widget) XmSelectionBoxGetChild(LabelPrompt,
                        XmDIALOG_TEXT),XmNfontList,Curr_fontlist,NULL);
    }
}

/*
 * Label button callback function
 */
void LabelChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    switch((long)client_data) {
    case EDIT_INDX: LabelDialog("Label",
                    "Enter label\nIndicate position using the mouse cursor\nClick left mouse button to display text",
                                Plot_info.curr_label,NULL);
               break;
    default:   break;
    }
}

/**********
 * LOCATE
 **********/
/*
 * Locate button callback function
 */
void LocateChosen(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{ 
    int max,cell_number;

    max = (Plot_info.max_row+1)*(Plot_info.max_col+1)-1;
    cell_number = Plot_info.curr_cell;
    switch((long)client_data) {
    case NEXT: if (Plot_info.curr_cell<max) {
                    Plot_info.curr_cell++;
               }
               else {
                   Plot_info.curr_cell=0;
                   ErrorMsg("moved to first cell",main_Menu,0);
               }
               break;
    case PREV: if (Plot_info.curr_cell-1>=0) {
                    Plot_info.curr_cell--;
               }
               else {
                   Plot_info.curr_cell=max;
                   ErrorMsg("moved to last cell",main_Menu,0);
               }
               break;
    case SPEC: LocateDialog( main_Menu,plotArea );
               break;
    default:   break;
    }
    if (cell_number!=Plot_info.curr_cell) 
        ChangeCell( cell_number,Plot_info.curr_cell );
}

/**********
 * COLOUR
 **********/

/*
 * Colour callback function
 */
void ColourChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    Settings.fg = (long) client_data;
    setpencolor_(&Settings.fg);
    Save_option((XtPointer)&Settings.fg,O_FOREGRND,SYB_INT);
}

/**********
 * LINE
 **********/

/*
 * Line callback function
 */
void LineChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    int width=1;

    Settings.linestyle = (long) client_data;
    dashln_(&Settings.linestyle,&width);
    Save_option((XtPointer)&Settings.linestyle,O_LINESTYLE,SYB_INT);
}

/**********
 * RESCALE
 **********/

/*
 * Rescale callback function
 */
void RescaleChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    /* XtUngrabPointer(w,time()); XtUngrabKeyboard(w,time()); */
    Settings.rescale = 1;
    Save_option((XtPointer)&Settings.rescale,O_RESCALE,SYB_BOOL);
/*  following parameters will be reset when new solution is read */
/*  Pwindo[XCMINREF]=Pwindo[XCMIN];
    Pwindo[XCMAXREF]=Pwindo[XCMAX];
    Pwindo[YCMINREF]=Pwindo[YCMIN];
    Pwindo[YCMAXREF]=Pwindo[YCMAX]; */
/*  it is questionable whether init-box should be called at this point */
    init_box( Pwindo,&Settings );
    MeshParams( Pwindo[YCMAX]-Pwindo[YCMIN],Pwindo[XCMAX]-Pwindo[XCMIN],
                                                             Mesh );
    /*err=CheckCurrentRec(Err_str,&err_num);*/
    Settings.rescale = 0;
}

/**********
 * ZOOM
 **********/
/*
 * Zoom button callback function
 */
/* void ZoomChosen(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{ 
    float zoomfactor=Settings.zoom;
    Widget *parent = (Widget *)client_data;

    ZoomDialog(*parent,"Zoom Dialog","Zoom factor",Settings.zoom,&zoomfactor );
} */
/*
 * zoom ok button callback function  
   This routine now obsolete with new ZoomChosen
 */
/*static void ZoomPromptOK(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    char *str;
    float zoomfactor;
    XmSelectionBoxCallbackStruct *selection;
                                                                                
    selection=(XmSelectionBoxCallbackStruct *) call_data;
	XmStringGetLtoR(selection->value, 
				XmSTRING_DEFAULT_CHARSET, &str);
    zoomfactor = atof( str );
    XtFree(str);
    if (fabs(zoomfactor-Settings.zoom)>0.0001) {
        Settings.prev_zoom = Settings.zoom;
        Settings.zoom = zoomfactor;
        Settings.rescale = 1;
        init_box(Pwindo,&Settings);
        Settings.rescale = 0;
        Save_option((XtPointer)&Settings.zoom,O_ZOOM,SYB_FLOAT);
    }
} */

/**********
 * DELETE
 **********/

/*
 * Delete callback function
 */
void DeleteChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    switch((long)client_data) {
    case CURR_CELL: DeleteCell(Plot_info.curr_cell);
                    break;
    case ALL_CELLS: DeleteAllCells();
                    break;
    default:        break;
    }
}

/**********
 * HELP
 **********/

/*
 * Help callback function
 */
void HelpChosen(w, client_data, call_data)
Widget w;
XtPointer client_data;
XtPointer call_data;
{
    Pixel fg,bg;

    switch((long)client_data) {
    case VIEW: XtVaGetValues( main_Menu,
                    XmNforeground,&fg,
                    XmNbackground,&bg,
                    NULL);
               HelpView(w,fg,bg);
               break;
    default:   break;
    }
}

/*
 *  Initialize font to default size or closest menu font
 */
void InitFont(w)
Widget w;
{
    Widget wdgt,pulldown;
    XmFontContext cntxt;
    XtPointer tmp;
    XmFontType type_return;
    XmFontListEntry entry;

    /*
     * Get the size of the last font in the menu -
     * assumed to be the largest
     * Set size of title area accordingly
     */
    wdgt = XtNameToWidget(main_Menu,mainopts[OPTIONS].label);
    if (wdgt) {
        XtVaGetValues( wdgt,XmNsubMenuId,&pulldown,NULL );
        wdgt = XtNameToWidget(pulldown,optionopts[LABEL_INDX].label);
        XtVaGetValues( wdgt,XmNsubMenuId,&pulldown,NULL );
        wdgt = XtNameToWidget(pulldown,labelopts[FONT_INDX].label);
    }
    if (wdgt) {
        XtVaGetValues( wdgt,XmNsubMenuId,&pulldown,NULL );
        wdgt = XtNameToWidget(pulldown,fonts[LRGE_FONT].label);
    }
    if (wdgt) {
        XtVaGetValues( wdgt, XmNfontList, &Curr_fontlist, NULL );
        if (!XmFontListInitFontContext(&cntxt,Curr_fontlist)) {
            XtWarning("Could not allocate font context");
        }
        else {
            entry = XmFontListNextEntry( cntxt );
            tmp = XmFontListEntryGetFont(entry,&type_return);
            if (type_return==XmFONT_IS_FONT)
                Plot_info.large_font = (XFontStruct *) tmp;
            Plot_info.title_offset = (Plot_info.large_font->ascent+
                                    Plot_info.large_font->descent)*1.5;
        }
        XmFontListFreeFontContext(cntxt);
    }
    /*
     * Initialize the font to the default font in .options
     */
    FindFont("Helvetica",Settings.dfltfonthgt);
    Plot_info.dflt_font = Plot_info.font;
}

/*
 * Find closest menu font to hgt
 */
void FindFont(char* name, int hgt)
{
    unsigned char found=0;
    char testtag[10],*tag;
    int i, last;
    Widget wdgt,pulldown;
    XtPointer tmp;
    XmFontType type_return;
    XmFontContext cntxt;
    XmFontListEntry entry;

    /*
     * Check font menu for a matching font size
     * Set the current font if successful
     */
    if (!strcmp(name,"Symbol"))sprintf(testtag,"%s%02dS","TAG",hgt);
    else sprintf(testtag,"%s%02dH","TAG",hgt);
    last = strlen(testtag) -1;
    wdgt = XtNameToWidget(main_Menu,mainopts[OPTIONS].label);
    if (wdgt) {
        XtVaGetValues( wdgt,XmNsubMenuId,&pulldown,NULL );
        wdgt = XtNameToWidget(pulldown,optionopts[LABEL_INDX].label);
        XtVaGetValues( wdgt,XmNsubMenuId,&pulldown,NULL );
        wdgt = XtNameToWidget(pulldown,labelopts[FONT_INDX].label);
    }
    if (wdgt) {
        XtVaGetValues( wdgt,XmNsubMenuId,&pulldown,NULL );
        for (i=0;i<(XtNumber(fonts)-1)&&!found;i++) {
            wdgt = XtNameToWidget(pulldown,fonts[i].label);
            if (wdgt) {
                XtVaGetValues( wdgt, XmNfontList, &Curr_fontlist, NULL );
                if (!XmFontListInitFontContext(&cntxt,Curr_fontlist)) {
                    XtWarning("Could not allocate font context");
                }
                else {
                    while ((entry = XmFontListNextEntry( cntxt)) && !found ) {
                    /*entry = XmFontListNextEntry( cntxt );*/
                    tag = XmFontListEntryGetTag( entry );
                    if (strcmp(testtag,tag)<=0 && testtag[last]==tag[last]) {
                        tmp = XmFontListEntryGetFont(entry,&type_return);
                        if (type_return==XmFONT_IS_FONT)
                            Plot_info.font = (XFontStruct *) tmp;
                        found = 1;
                    }
                    XtFree(tag);
                    }
                }
            }
        }
        XmFontListFreeFontContext(cntxt);
    }
}

/*
 * These fns display the appropriate dialog to
 * request input from the user
 * Events for the main window are not processed until
 * the dialog is completed
 */
int GetUserVals(float *user_vals, int dialogflag, unsigned char prfl_flag, int *plot, char *labels )
{
/*
 * answer indicates that the user has acknowledged
 * the window and whether ok or cancel
 */
    int answer=0;

    if (!XtIsRealized(topLevel)||Plot_info.display_log) return(0);

    switch(dialogflag) {
    case CNTRVALS: CntrvalsDialog( main_Menu,user_vals,&answer );
                   break;
    case PRFLVALS: PrflvalsDialog( main_Menu,user_vals,&answer );
                   break;
    case PRFLPTS:  PrflptsDialog( main_Menu,user_vals,&answer,prfl_flag );
                   break;
    case PRFLDIR:  Prfl2D_Dialog( main_Menu,user_vals,plot,&answer );
                   break;
    case ZOOMVALS: ZoomValsDialog( main_Menu,user_vals,&answer );
                   break;
    }
    while (!answer||XtAppPending(App_context))
        XtAppProcessEvent(App_context, XtIMAll);
    return(answer);
}

int GetProfileVals(plot_descrip,ncomp)
int *plot_descrip;
int ncomp;
{
    char dfltlabels[2][MAXNAME];
    int ret_val=0;
    int i,len=1,err=0;
    float tmpmax, tmpmin, xref, yref;
    float current_vals[MAXPRFLVALS];
/*  fprintf(stdout,"GetProfileVals: *plot descrip =%i\n",*plot_descrip); */

/*    XREFM and YREFM usage depends on context */
    xref=Data_vars_fl[XREFM];
    yref=Data_vars_fl[YREFM];
    if (Data_vars_int[IDEFTYP]>110) {
       xref=0.;
       yref=0.;
    }

/*  copy profile end points from previous plot */
    for (i=0 ;i<MAXPRFLVALS;i++)
            current_vals[i]=Profile_vals[i];

/*  if no previous profile set default profile endpoints from loaded mesh */
    if (current_vals[PRFL_X1]==0 && current_vals[PRFL_Y1]==0 &&
        current_vals[PRFL_X2]==0 && current_vals[PRFL_Y2]==0) {
        if (*plot_descrip==DIM1) {
           current_vals[PRFL_X1]=current_vals[PRFL_X2]=
                             (Pwindo[XCMAX]+Pwindo[XCMIN])/2;
           current_vals[PRFL_Y1]=Pwindo[YCMIN];
           current_vals[PRFL_Y2]=Pwindo[YCMAX];
        }
        else if (*plot_descrip==DIM2_X) {
           current_vals[PRFL_Y1]=current_vals[PRFL_Y2]=
                             (Pwindo[YCMAX]+Pwindo[YCMIN])/2;
           current_vals[PRFL_X1]=current_vals[LWR_LIM]=Pwindo[XCMIN];
           current_vals[PRFL_X2]=current_vals[UPR_LIM]=Pwindo[XCMAX];
        }
        else {
           current_vals[PRFL_X1]=current_vals[PRFL_X2]=
                             (Pwindo[XCMAX]+Pwindo[XCMIN])/2;
           current_vals[PRFL_Y1]=current_vals[LWR_LIM]=Pwindo[YCMIN];
           current_vals[PRFL_Y2]=current_vals[UPR_LIM]=Pwindo[YCMAX];
        }

/*  project (x,y) to long,lat if needed for initial values */

        if (ncomp==-1) {
              projectdeg_(&current_vals[PRFL_X1],&current_vals[PRFL_Y1],
                       &xref,&yref,&len,&ncomp,&err);
              projectdeg_(&current_vals[PRFL_X2],&current_vals[PRFL_Y2],
                       &xref,&yref,&len,&ncomp,&err);
        }

    }
    if (*plot_descrip==DIM1){
        /* allow user to adjust values or cancel */
        ret_val = GetUserVals( &current_vals[PRFL_X1],PRFLPTS,0,
                                    plot_descrip,(char*)dfltlabels );
                                          
    }
    else if (*plot_descrip==PRFLMRK) {
        ret_val = GetUserVals( &current_vals[PRFL_X1],PRFLPTS,1,NULL,NULL );
    }
    else {     /* applies to 2D profile  */
        ret_val = GetUserVals( &current_vals[PRFL_MIN_U],PRFLDIR,0,
                                                  plot_descrip, NULL );
    }
    if (ret_val!=USER_CANCEL) {
        for (i=0;i<MAXPRFLVALS;i++) 
                       Profile_vals[i]=current_vals[i];
 /*   fprintf(stdout,"Profile_vals = current_vals in GetProfileVals\n"); */

 /*  project coordinates from long, lat to (x-y) if required
       (this conditional changed from IDEFTYP = 110) */

        if (ncomp==-1) {
            projectxy_(&Profile_vals[PRFL_X1],&Profile_vals[PRFL_Y1],
                       &xref,&yref,&len,&ncomp,&err);
            projectxy_(&Profile_vals[PRFL_X2],&Profile_vals[PRFL_Y2],
                       &xref,&yref,&len,&ncomp,&err);
        /*  fprintf(stdout,"converting to (x-y) in GetProfileVals\n"); */
            if (*plot_descrip==DIM2_X) {
              tmpmin = tmpmax = yref;
              projectxy_(&Profile_vals[LWR_LIM],&tmpmin,
                         &xref,&yref,&len,&ncomp,&err);
              projectxy_(&Profile_vals[UPR_LIM],&tmpmax,
                         &xref,&yref,&len,&ncomp,&err);
            }
            else if (*plot_descrip==DIM2_Y) {
              tmpmin = tmpmax = xref;
              projectxy_(&tmpmin,&Profile_vals[LWR_LIM],
                         &xref,&yref,&len,&ncomp,&err);
              projectxy_(&tmpmax,&Profile_vals[UPR_LIM],
                         &xref,&yref,&len,&ncomp,&err);
            }
        }
    }
    return(ret_val);
}
/*  ZoomChosen upgraded here, replaces old version */
void ZoomChosen(w, client_data, call_data)
Widget w;
XtPointer client_data, call_data;
{
    int ret_val=0;
    float zoom_vals[3];

/*  set previously stored values of zoom factor, x centre, ycentre
     stored values should be in degrees if spherical option  */
    zoom_vals[0]=Settings.zoom;
    zoom_vals[1]=Settings.xcentre;
    zoom_vals[2]=Settings.ycentre;

/* allow user to adjust values or cancel */
    ret_val = GetUserVals( &zoom_vals[0],ZOOMVALS,1,NULL,NULL );

/* reset the stored values */
    if (ret_val!=USER_CANCEL) {
        Settings.prev_zoom = Settings.zoom;
        Settings.zoom=zoom_vals[0];
        Settings.xcentre=zoom_vals[1];
        Settings.ycentre=zoom_vals[2];
        Save_option((XtPointer)&Settings.zoom,O_ZOOM,SYB_FLOAT);
        Save_option((XtPointer)&Settings.xcentre,O_XCENTRE,SYB_FLOAT);
        Save_option((XtPointer)&Settings.ycentre,O_YCENTRE,SYB_FLOAT);
/* clipping is automatically switched on if zoom > 1.1 */
        if(Settings.zoom > 1.1){
          Settings.clip_to_cell = 1;
          Save_option((XtPointer)&Settings.clip_to_cell,
                                             O_CLIPTOCELL,SYB_BOOL);
        }
/*  call init_box to invoke modified zoom settings */
        Settings.rescale = 0;
        init_box(Pwindo,&Settings);
    }
}

void InitTitle( title )
char *title;
{
    XFontStruct *oldfont;

    if ((Plot_info.title=(char *)malloc((strlen(title)+1)*sizeof(char)))
                                                            ==NULL)
        ErrorMsg("Title",main_Menu,MALLOC_ERR);
    else {
        strcpy(Plot_info.title,title);
        oldfont = Plot_info.font;
        Plot_info.font = Plot_info.large_font;
        DrawTitle( Plot_info.title );
        Plot_info.font = oldfont;
        UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
    }
}

/*
 * update the data displayed below the drawing area
 */
int UpdateMsgArea()
{
    XmString str;
    char buf[1024],*fname,blank[10];

    strcpy(blank,"         ");
    strcpy(buf,blank);
    if (Plot_info.inp_file==NULL) fname=blank;
    else {
        if ((fname=strrchr(Plot_info.inp_file->fname,'/'))==NULL)
            fname=Plot_info.inp_file->fname;
        else fname++;
        sprintf(buf, "%s curr=%3d time=%5.3f req=%3d ref=%3d max=%3d", fname,
                                        Plot_info.inp_file->rec_curr,
                                        Data_vars_fl[TIME],
                                        Plot_info.inp_file->rec_req,
                                        Plot_info.inp_file->ref_curr,
                                        Plot_info.inp_file->rec_max);
    }
#if XmVersion<1002
    str = XmStringCreateSimple(buf);
#else
    str = XmStringCreateLocalized(buf);
#endif
    XtVaSetValues(currvalsArea, XmNlabelString, str, NULL);
    XmStringFree(str);
    return(0);
}

/*
 * Create the menu, the cascade button
 * that owns the menu, and then the submenu items.
 */
Widget BuildMenu(parent, menu_title, spacing, data)
Widget parent;
char *menu_title;
int spacing;
MenuItemData *data;
{
    Widget pulldown, cascade, widget;
    static int owner;
    int i;
    XmString str;
    Arg args[5];
    Pixel bg, fg,top_shdw, bottom_shdw,select_col;
    extern XColor Colours[];
    extern Colormap cmap;

    i=0;
    if (spacing>0) {
        XtSetArg(args[i], XmNmarginWidth, spacing); i++;
        XtSetArg(args[i], XmNmarginHeight, spacing); i++;
        XtSetArg(args[i], XmNspacing, spacing); i++;
    }
    pulldown = XmCreatePulldownMenu(parent, "plldwn", args, i);

    str = XmStringCreateSimple(menu_title);
    cascade = XtVaCreateManagedWidget(menu_title,
        xmCascadeButtonGadgetClass, parent,
        XmNsubMenuId,   pulldown,
        XmNlabelString, str,
        NULL);

    XmStringFree(str);

    /* add the menu items */
    for (i = 0; data[i].label != NULL; i++) {
        widget = 0;
        /*
         * If submenu exists, create the pull-right menu
         * If owner is set let submenu items store their origin
         */
        if (data[i].submenu) {
         // if (strcmp(data[i].label,ELLE_LABEL) || Settings.elle) {
            if (data[i].submenu_data!=NO_VALUE) owner=data[i].submenu_data;
            else owner=NO_VALUE;
            widget = BuildMenu(pulldown,
                            data[i].label, data[i].spacing, data[i].submenu);
          //}
        }
        else
          //if (strcmp(data[i].label,ELLE_LABEL) || Settings.elle) 
            widget = XtVaCreateManagedWidget(data[i].label,
                     *data[i].class, pulldown,
                     NULL);
        if (widget) {
          if (data[i].bg!=NO_VALUE) {
              bg = Colours[data[i].bg].pixel;
              XmGetColors( XtScreen(pulldown),cmap,bg,
                              &fg, &top_shdw, &bottom_shdw, &select_col);
              XtVaSetValues(widget,
                      XmNbackground, bg,
                      XmNborderColor, fg,
                      XmNtopShadowColor, top_shdw,
                      XmNbottomShadowColor, bottom_shdw,
                      XmNarmColor, select_col,
                      NULL);
          }
          /*
           * note - activateCallback
           */
          if (data[i].callback)
              XtAddCallback(widget, XmNactivateCallback,
                  data[i].callback, data[i].callback_data);
          if (owner!=NO_VALUE)
              XtVaSetValues(widget,
                      XmNuserData, owner, NULL);
          if (data[i].fontindex!=NO_VALUE)
              XtVaSetValues(widget,
                      XmNfontList, Fontlists[data[i].fontindex],
                      NULL);
      }
    }
    return cascade;
}

int Create_widgets( top )
Widget top;
{
    Widget mainWindow, button;
    Pixel fg,bg;
    int i;
    Dimension wdth;

    /* create main window */
    Main = mainWindow = XtVaCreateManagedWidget( "mainWindow",
            xmMainWindowWidgetClass, top,
			XmNresizable, True,
            NULL);

    main_Menu = XmCreateMenuBar(
        mainWindow,   /* parent widget*/
        "main_Menu",  /* widget name */
        NULL,
        0
        );
    XtManageChild( main_Menu );

    Minwidth=0;
    for ( i=0;mainopts[i].label!=NULL;i++ ) {
        button = BuildMenu( main_Menu,mainopts[i].label,
                                mainopts[i].spacing,
                                mainopts[i].submenu );
        XtVaGetValues( button,XmNwidth,&wdth,NULL);
        Minwidth += wdth;
    }
    Minwidth += (Minwidth/(XtNumber(mainopts)-1));
        /*
         * help button must be the last one created.
         * Set the resource so that it is placed on
         * the far right of the menu bar
         */
    XtVaSetValues( main_Menu,XmNmenuHelpWidget,button,NULL );

    InitFont(top);

    plotArea = XtVaCreateManagedWidget( "plotArea",
            plotWidgetClass, mainWindow,
            NULL);

    currvalsArea = XtVaCreateManagedWidget("current values",
    xmLabelWidgetClass, mainWindow,
    XmNalignment,   XmALIGNMENT_BEGINNING,
    NULL);
  
    XmMainWindowSetAreas(
        mainWindow,   /* widget */
        main_Menu,      /* menu */
        NULL,          /* command */
        NULL,          /* h scroll */
        NULL,          /* v scroll */
        plotArea      /* region */
        );

    XtVaSetValues(mainWindow,
        XmNmessageWindow, currvalsArea,
        NULL);

    XtVaGetValues(main_Menu,
        XmNforeground, &fg,
        XmNbackground, &bg,
        NULL);

    XtVaSetValues(currvalsArea,
        XmNforeground, fg,
        XmNbackground, bg,
        NULL);
        
    return(0);
}

void clear_menus()
{
    int i;
    Widget button,pulldown;

        /* reset sensitivity for xyplot menu */
    button=XtNameToWidget(main_Menu,mainopts[XY_PLOT].label);
    XtVaGetValues( button,XmNsubMenuId,&pulldown,NULL );
    button=XtNameToWidget(pulldown,xyplot_opts[0].label);
    XtVaGetValues( button,XmNsubMenuId,&pulldown,NULL );
    i=0;
    while (meshopts[i].label!=NULL) {
        button=XtNameToWidget(pulldown, meshopts[i].label);
        XtVaSetValues(button, XmNsensitive, True, NULL);
        i++;
    }
}

#if XY
int systemcall(s)
char *s;
{
    int done=0;
    int status,pid,wval,w;

    if ((pid = fork()) == 0) {
        execlp("sh", "sh", "-c", s, (char *)0);
    while (!done) {
        wval=wait(&status);
        fprintf(stderr,"%d %d   ",wval,status);
    }
        exit(127);
    }
    /*while ((wval=wait(&status)) != pid && wval != -1) ;*/
                        /* what to do if wval==error?? */
    if (wval==-1) status = -1;
    return( status );
}
#endif

#if XY
 /usr/include/asm-generic/errno-base.h

/* SPDX-License-Identifier: GPL-2.0 WITH Linux-syscall-note */
#define	EPERM		 1	/* Operation not permitted */
#define	ENOENT		 2	/* No such file or directory */
#define	ESRCH		 3	/* No such process */
#define	EINTR		 4	/* Interrupted system call */
#define	EIO		 5	/* I/O error */
#define	ENXIO		 6	/* No such device or address */
#define	E2BIG		 7	/* Argument list too long */
#define	ENOEXEC		 8	/* Exec format error */
#define	EBADF		 9	/* Bad file number */
#define	ECHILD		10	/* No child processes */
#define	EAGAIN		11	/* Try again */
#define	ENOMEM		12	/* Out of memory */
#define	EACCES		13	/* Permission denied */
#define	EFAULT		14	/* Bad address */
#define	ENOTBLK		15	/* Block device required */
#define	EBUSY		16	/* Device or resource busy */
#define	EEXIST		17	/* File exists */
#define	EXDEV		18	/* Cross-device link */
#define	ENODEV		19	/* No such device */
#define	ENOTDIR		20	/* Not a directory */
#define	EISDIR		21	/* Is a directory */
#define	EINVAL		22	/* Invalid argument */
#define	ENFILE		23	/* File table overflow */
#define	EMFILE		24	/* Too many open files */
#define	ENOTTY		25	/* Not a typewriter */
#define	ETXTBSY		26	/* Text file busy */
#define	EFBIG		27	/* File too large */
#define	ENOSPC		28	/* No space left on device */
#define	ESPIPE		29	/* Illegal seek */
#define	EROFS		30	/* Read-only file system */
#define	EMLINK		31	/* Too many links */
#define	EPIPE		32	/* Broken pipe */
#define	EDOM		33	/* Math argument out of domain of func */
#define	ERANGE		34	/* Math result not representable */

#endif
int Run_external_prog( char *command, char *msg, int len )
{
    int status=0;
    int err=0;

    pid_t pid = fork();
    if ( pid == 0 ) { //child
    /*
       run command.
       -b enables asynchronous notification of completion, else EAGAIN or
          EINTR commonly returned if command not found etc
     */
        execlp("sh", "sh", "-b", command, (char *)0);
    // Child will only return if error
        strncat(msg,strerror(errno),len);
        exit(127);
    }
    if ( pid >0 ) { //parent
        int wstatus=0;
        status = wait(&wstatus);
        if (WIFEXITED(wstatus) && (status = WEXITSTATUS(wstatus))) {
            status = WEXITSTATUS(wstatus);
            err = errno;
            if (err) {
                strncat(msg,strerror(errno),len);
                status = -1;
            }
            else status = 0;
        }
        if (WIFSIGNALED(wstatus))
            status = WTERMSIG(wstatus);
    }
    return( status );
}

void SetCursor(w,cs)
Widget w;
Cursor cs;
{
    if (XtIsRealized(w)) {
        XDefineCursor(XtDisplay(topLevel),
                            XtWindow(XtParent(w)),cs);
        XFlush( XtDisplay(topLevel));
    }
}

void sybflush_()
{
    UpdateCell( plotArea,Plot_info.curr_cell,DRAW );
}

int Init_App(int argc, char** argv, FILE** optsfp)
{
    int err=0;
    char optsfile[256],LogFile[256],PsFile[256],*dir;
    char input_str[SYB_FILENAME_MAX+1];
    extern Display *Dsply;

    static XrmOptionDescRec table[] = {
        {"-pw",            "*pixmapWidthInCells", XrmoptionSepArg, NULL},
        {"-ph",            "*pixmapHeightInCells", XrmoptionSepArg, NULL},
               /* just to make motif ignore -i option */
        {"-ilbl",            "*labelList", XrmoptionSepArg, NULL},
    };
    
    XtToolkitInitialize();
    App_context = XtCreateApplicationContext();
    XtAppSetFallbackResources(App_context,default_resources);
    if ((Dsply = XtOpenDisplay(App_context,NULL,NULL,"Xsybil",
                           table,XtNumber(table),&argc,argv))==NULL)
        error_msg(err,"Unable to open display");
#if XY
    topLevel = XtVaAppInitialize(
        &App_context,       /* Application context */
        "Xsybil",          /* application class name */
        table,              /* command line option list */
        XtNumber(table),
        &argc, argv,        /* command line args */
        default_resources,  /* for missing app-defaults file */
        NULL);              /* terminate varargs list */
#endif

    strcpy(optsfile,DEFAULT_OPTIONS_FILE);
    strcpy(LogFile,"");
    strcpy(PsFile,"");
    if (argc > 1) {
/*
        Syntax( argc,argv );
        return(1);
*/
       if ((err=ParseOptions(argc,argv,LogFile,PsFile))) {
           Syntax( argc,argv );
           if (err==HELP_ERR) return(1);
       }
    }
    if (LogFile[0] != '\0') {
        if ((*optsfp = fopen(LogFile,"r"))!=NULL) {
            if ((err = initial_options(*optsfp,&Settings,input_str)))
                error_msg(err,LogFile);
        }
        else error_msg(OPEN_ERR,LogFile);
    }
    else {
        /* check in current dir for options file */
        if ((*optsfp = fopen(optsfile,"r"))==NULL) {
            /* check in user's home dir for options file */
            if ((dir = getenv("HOME"))!=NULL) {
                strcpy(optsfile,dir);
                strcat(optsfile,"/");
                strcat(optsfile,DEFAULT_OPTIONS_FILE);
                *optsfp = fopen(optsfile,"r");
            }
        }
        if (*optsfp!=NULL) {
            if ((err = initial_options(*optsfp,&Settings,input_str)))
                error_msg(err,input_str);
            fclose(*optsfp);
        }
        *optsfp = NULL;
    }
    return( 0 );
}

int Set_Colours(int* red, int* green, int* blue, int* irgbv, int steps)
{
    extern Display *Dsply;

    if ((Plot_info.max_colours=CreateColormap( Dsply,
                                               XDefaultScreenOfDisplay(Dsply),
                                               red,green,blue,irgbv,
                                               steps)) == 0) return( 1 );
/*
    if ((Plot_info.max_colours=CreateColormap( XtDisplay(topLevel),
                                               XtScreen(topLevel),
                                               red,green,blue,irgbv,
                                               steps)) == 0) return( 1 );
*/
    Colour_range[0] = USEFL_SIZE;
    Colour_range[1] = Plot_info.max_colours;
    Arrow_colours[0] = 7; /* cyan */
    Arrow_colours[1] = 5; /* magenta */
    return(0);
}

void CreateFontlists(fontlists)
XmFontList fontlists[];
{
    int i;
    XFontStruct *font;
    XmFontListEntry entry;

i = XtNumber(labelfonts);
    for (i=0; i<MAX_FONTS && i<(XtNumber(labelfonts)); i++) {
        font = XLoadQueryFont(XtDisplay(topLevel),labelfonts[i].name);
        entry = XmFontListEntryCreate(labelfonts[i].tag,
                      XmFONT_IS_FONT,font);
        fontlists[i] = XmFontListAppendEntry(NULL,entry);
        XtFree((char *)entry);
    }
}

int Run_App(FILE* logfp)
{
    int dum=1;
    extern Colormap cmap;
    extern Display *Dsply;
    topLevel = XtVaAppCreateShell(NULL,"Xsybil",
                                applicationShellWidgetClass,
                                Dsply,
                                XtNcolormap,cmap,
                                NULL);
    CreateFontlists(Fontlists);
    if (Create_widgets( topLevel )) return(1);

    setpencolor_(&Settings.fg);
    setlinewidth_(&Settings.linewidth);
    dashln_(&Settings.linestyle,&dum);

    Waiting = XCreateFontCursor( XtDisplay(topLevel),XC_watch);
    /* for plotting area with black bg */
    XRecolorCursor(  XtDisplay(topLevel),Waiting,&Colours[1],
                     &Colours[0]);

    XtRealizeWidget(topLevel);

    if (logfp) {
        SetCursor(main_Menu,Waiting);
        Plot_info.display_log = 1;
        process_log_file(logfp,&Settings);
        Plot_info.display_log = 0;
        UpdateMsgArea();
        SetCursor(main_Menu,None);
    }
    else {
        ChangeCell( 0,Plot_info.curr_cell );
    }
    
    XtAppMainLoop(App_context);
    return(0);
}
