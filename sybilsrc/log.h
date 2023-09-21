/*--------------------------------------------------------------------*
 *    Basil / Sybil:   log.h  1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

#ifndef _log_h
#define _log_h

#define LABEL_MRK  '"'
#define COMMENT_MRK  #
#define MAX_OPTION_NAME 20
#define CMNTS_FILE "logcomments.txt"
/*
 * allowed option types
 */
#define SYB_BOOL   201
#define SYB_INT    202
#define SYB_FLOAT  203

/*
 * options - used when reading and writing options
 */
#define O_ORIENT      101
#define O_MARKCELL    102
#define O_ROWS        103
#define O_COLS        104
#define O_FONTHGT     105
#define O_XMARG       106
#define O_YMARG       107
#define O_ZOOM        108
#define O_XCENTRE     109
#define O_YCENTRE     110
#define O_PGSCL       111
#define O_TL0         112
#define O_RHOC        113
#define O_RHOM        114
#define O_Z0          115
#define O_TMAX        116
#define O_ELEV0       117
#define O_BGAM0       118
#define O_BGAM1       119
#define O_NDIVR       120
#define O_S0          121
#define O_DELS        122
#define O_DELOM       123
#define O_NX3         124
#define O_MP          125
#define O_MPE         126
#define O_NP          127
#define O_CNTRPLT     128
#define O_PRFLPTS     129
#define O_FLIP        130
#define O_DBLE        131
#define O_LABEL       132
#define O_STIP        133
#define O_SOLNROT     134
#define O_FOREGRND    135
#define O_LINESTYLE   136
#define O_LINEWDTH    137
#define O_COLMAP      138
#define O_VERBOSE     139
#define O_PG_XMARG    140
#define O_PG_YMARG    141
#define O_PAPERSZ     142
#define O_COLBAR      143
#define O_TIC         144
#define O_VISCMIN     145
#define O_VISCMAX     146
#define O_SEMIN       147
#define O_SEMAX       148
#define O_CLIPTOCELL  149
#define O_FILEPATH    150
#define O_RESCALE     151 /* options written only if < O_RESCALE */
#define O_ELLE        152
#define O_TEXTBG      153
#define O_AXISXMIN    154
#define O_AXISXMAX    155
#define O_AXISYMIN    156
#define O_AXISYMAX    157

static valid_terms keywords[] = {
                           {"File", FILE_OPT},
                           {"Data_In", DATA_IN},
                           {"XYPlot", XY_PLOT},
                           {"Profile", PRFL},
                           {"Arrow", ARROWS},
                           {"Contour", CNTRS},
                           {"Label", LABEL},
                           {"Locate", LOCATE},
                           {"Title", TITLE},
                           {"Options", OPTIONS},
                           {NULL}
                          };


static valid_terms xyplots[] = {
                           {"Mesh", MESH},
                           {"Deform", DEFORM},
                           {"BndBox",BBOX},
                           {"LGMesh",LGMESH},
                           {"StrainMark",STRN_MRKR},
                           {"Elle",ELLE},
                           {NULL}
                          };
static valid_terms meshplots[] = {
                           {"elements", ELMNT},
                           {"element+num", ELMNTNUM},
                           {"boundary", BNDRY},
                           {"internal_boundary", INTBND},
                           {"internal_segments", INTLBNDRY},
                           {"viscosity", VISCMSH},
                           {"exponent", SEMSH},
                           {NULL}
                          };
static valid_terms deformplots[] = {
                           {"ellp", ELLP},
                           {"trel", TREL},
                           {"tria", TRIA},
                           {NULL}
                          };
static valid_terms lgmeshplots[] = {
                           {"elements", LGELMNT},
                           {"boundary", LGBNDRY},
                           {NULL}
                          };
static valid_terms arrowplots[] = {
                           {"Velocity", SYB_VELOCITY},
                           {"Strain", SYB_STRAIN},
                           {"Stress", SYB_STRESS},
                           {"Rotation", ROTATION},
                           {NULL}
                          };
static valid_terms profileplots[] = {
                           {"1_D", DIM1},
                           {"2_D_X", DIM2_X},
                           {"2_D_Y", DIM2_Y},
                           {"Mark",PRFLMRK},
                           {NULL}
                          };
static valid_terms contourplots[] = {
                           {"Velocity", SYB_VELOCITY},
                           {"Strain", SYB_STRAIN},
                           {"Stress", SYB_STRESS},
                           {"Thickness", THICKNESS},
                           {"Layer", LAYER},
                           {"Density", DENSITY},
                           {"Gravity(-Y)", GRAVITY},
                           {NULL}
                          };
static valid_terms layerplots[] = {
                           {"thickness", THICKNESS},
                           {"topog", TOPOG},
                           {"gravPE", GRAVPE},
                           {NULL}
                          };
static valid_terms gravityplots[] = {
                           {"Bouguer", GRBOUG},
                           {"topog", TOPOG},
                           {"free-air", GRFREE},
                           {NULL}
                          };
static valid_terms velocityplots[] = {
                           {"Ux", UX},
                           {"Uy", UY},
                           {"Um", UM},
                           {"Ur", UR},
                           {"Uth", UTH},
                           {NULL}
                          };
static valid_terms strain_arrow_plots[] = {
                           {"pstr", PSTR},
                           {"pstm", PSTR},
                           {"pstd", SIGD},
                           {"mssr", MSSR},
                           {"ssft", SSFT},
                           {"taum", TAUM},
                           {"taud", TAUD},
                           {"sigd", SIGD},
                           {"sigm", SIGM},
                           {NULL}
                          };
static valid_terms strain_cntr_plots[] = {
                           {"edxx", EDXX},
                           {"edyy", EDYY},
                           {"edzz", EDZZ},
                           {"edxy", EDXY},
                           {"psr1", PSR1},
                           {"psr2", PSR2},
                           {"msst", MSST},
                           {"cang", CANG},
                           {"tang", TANG},
                           {"sang", SANG},
                           {"dblc", DBLC},
                           {"vort", VORT},
                           {"ed2i", ED2I},
                           {"vota", VOTA},
                           {"taud", TAUD},
                           {"TAUD", TAUD},
                           {"taum", TAUM},
                           {"taxx", TAXX},
                           {"tayy", TAYY},
                           {"tazz", TAZZ},
                           {"taxy", TAXY},
                           {"tau1", TAU1},
                           {"tau2", TAU2},
                           {"sixx", SIXX},
                           {"siyy", SIYY},
                           {"sizz", SIZZ},
                           {"sig1", SIG1},
                           {"sig2", SIG2},
                           {"thdi", THDI},
                           {"pres", PRES},
                           {"brit", BRIT},
                           {"bri2", BRI2},
                           {"visc", VISC},
                           {NULL}
                          };
static valid_terms plot_types[] = { 
			      {"shade",SHADE},
                              {"lines",LINES},
                              {NULL}
                            };
static valid_terms option_terms[] = { 
				{"orientation",O_ORIENT},
                                {"paper_size",O_PAPERSZ},
                                {"mark_cell",O_MARKCELL},
                                {"filepath",O_FILEPATH},
                                {"rows",O_ROWS},
                                {"columns",O_COLS},
                                {"font_height",O_FONTHGT},
                                {"xmargin",O_XMARG},
                                {"ymargin",O_YMARG},
                                {"page_xmargin",O_PG_XMARG},
                                {"page_ymargin",O_PG_YMARG},
                                {"page_scale",O_PGSCL},
                                {"zoom",O_ZOOM},
                                {"xcentre",O_XCENTRE},
                                {"ycentre",O_YCENTRE},
                                {"tl0",O_TL0},
                                {"rhoc",O_RHOC},
                                {"rhom",O_RHOM},
                                {"z0",O_Z0},
                                {"tmax",O_TMAX},
                                {"elev0",O_ELEV0},
                                {"bgam0",O_BGAM0},
                                {"bgam1",O_BGAM1},
                                {"ndivr",O_NDIVR},
                                {"s0",O_S0},
                                {"dels",O_DELS},
                                {"delom",O_DELOM},
                                {"nx3",O_NX3},
                                {"mp",O_MP},
                                {"mpe",O_MPE},
                                {"np",O_NP},
                                {"contour_plot",O_CNTRPLT},
                                {"profile_pts",O_PRFLPTS},
                                {"flip",O_FLIP},
                                {"dble",O_DBLE},
                                {"label",O_LABEL},
                                {"stipple",O_STIP},
                                {"solution_rot",O_SOLNROT},
                                {"foreground",O_FOREGRND},
                                {"linestyle",O_LINESTYLE},
                                {"linewidth",O_LINEWDTH},
                                {"colour",O_COLMAP},
                                {"bar",O_COLBAR},
                                {"visc_min",O_VISCMIN},
                                {"visc_max",O_VISCMAX},
                                {"se_min",O_SEMIN},
                                {"se_max",O_SEMAX},
                                {"clip_to_cell",O_CLIPTOCELL},
                                {"rescale",O_RESCALE},
                                {"verbose",O_VERBOSE},
                                {"elle",O_ELLE},
                                {"ticmarks",O_TIC},
                                {"text_bg",O_TEXTBG},
                                {"x_axis_min",O_AXISXMIN},
                                {"x_axis_max",O_AXISXMAX},
                                {"y_axis_min",O_AXISYMIN},
                                {"y_axis_max", O_AXISYMAX},
                                { NULL }};
static valid_terms paper_terms[] = {   
				{"A4",A4_PAPER},
                                {"USLetter", US_PAPER},
                                { NULL} };
static valid_terms colmap_terms[] = { 
				{"greyscale",GREY_MAP},
                                {"standard", STD_MAP},
                                {"absolute", ABS_MAP},
                                { NULL} };
static valid_terms linestyle_terms[] = { 
				{"solid",SOLID},
                                {"dash", DASH},
                                {"dot", DOT},
                                { NULL} };
static valid_terms colourbar_terms[] = { 
				{"vertical",CB_VERT},
                                {"horizontal", CB_HORIZ},
                                {"none", CB_NONE},
                                { NULL }};
static valid_terms ticmark_terms[] = { 
				{"centred", J_CENTRE},
                                {"internal",J_BASE},
                                {"external", J_TOP},
                                {"none", TIC_NONE},
                                { NULL }};
static char *cntr_terms[] = { "colour","line","limits","stipple", NULL };
static char *colour_vals[] = { "min", "max", "scale", "bar", NULL };
static char *line_terms[] = { "min", "max", "scale", "shift", NULL };
static char *line_vals[] = { "level", "step", "num", NULL };
static char *stipple_vals[] = { "type", NULL };
static char *prfl_terms[] = { "line","pts","limits", NULL };
static char *prfl_pts[] = { "x1","y1","x2","y2", NULL };
static char *prfl_limits[] = { "lower","upper", NULL };
#endif
