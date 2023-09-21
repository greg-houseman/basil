/* id values
#define  REGIONS    20
#define  NODES      21
#define  COMMENT    22
#define  EXPAND     23
#define  COLOUR     24
#define  LOCATION   25
#define  E_OPTIONS  26 */

/* define unitcell corners */
#define BASELEFT  0
#define BASERIGHT 1
#define TOPRIGHT  2
#define TOPLEFT   3

#define  E_OPTIONS  20
#define  REGIONS    21
#define  COMMENT    22
#define  FLYNNS     23
#define  PARENTS    24
#define  EXPAND     25
#define  COLOUR     26

#define  LOCATION   30
#define  VELOCITY   31
#define  VEL_X      32
#define  VEL_Y      33
/* the next 2 defs are also in attrib.h - fix */
#define  CONC_A     38
#define  NUM_NB     39

#define  STRESS     40
#define  TAU_XX     41
#define  TAU_YY     42
#define  TAU_ZZ     43
#define  TAU_XY     44
#define  TAU_1      45
#define  PRESSURE   46

#define  NUM_STRESS_VALS 6


#define  CAXIS      50
#define  CAXIS_X    51
#define  CAXIS_Y    52
#define  CAXIS_Z    53

#define  ENERGY     60
#define  GBE_LUT    61

#define VISCOSITY   70

#define  MINERAL    80
#define  QUARTZ     81
#define  FELDSPAR   82

#define  STRAIN     90
#define  INCR_S     91
#define  BULK_S     92

#define  NUM_STRAIN_VALS 2

/*
 *general dummy flynn attributes
 */
#define  SPLIT     101

#define  GRAIN     102

static valid_terms FileKeys[] = {
                           {"OPTIONS", E_OPTIONS},
                           {"REGIONS", REGIONS},
                           {"FLYNNS", FLYNNS},
                           {"PARENTS", PARENTS},
                           {"EXPAND", EXPAND},
                           {"COLOUR", COLOUR},
                           {"LOCATION", LOCATION},
                           {"VELOCITY", VELOCITY},
                           {"STRESS", STRESS},
                           {"STRAIN", STRAIN},
                           {"CAXIS", CAXIS},
                           {"ENERGY", ENERGY},
                           {"VISCOSITY", VISCOSITY},
                           {"GBENERGY_LUT", GBE_LUT},
                           {"MINERAL", MINERAL},
                           {"GRAIN", GRAIN},
                           {"SPLIT", SPLIT},
                           {NULL}
                          };
/*typedef struct { char *name; int id; } valid_terms;*/

/*
static valid_terms FileKeys[] = {
                           "OPTIONS", E_OPTIONS,
                           "REGIONS", REGIONS,
                           "NODES", NODES,
                           "EXPAND", EXPAND,
                           "COLOUR", COLOUR,
                           "LOCATION", LOCATION,
                           NULL
                          }; */
/*
 * definitions used when accessing option fields
 */
#define RO_DISPLAY     201
#define RO_COLOURMAP   202
#define RO_STAGES      203
#define RO_VERBOSE     204
#define RO_INITFUNC    205
#define RO_RUNFUNC     206
#define RO_FILE        207
#define RO_EXTRAFILE   208
#define RO_SWITCHDIST  209
#define RO_MAXNODESEP  210
#define RO_MINNODESEP  211
#define RO_SPEEDUP     212
#define RO_SVFREQ      213
#define RO_SVFILEROOT  214
#define RO_CELLBBOX    215
#define RO_SSOFFSET    216
#define RO_CUMSSOFFSET 217
#define RO_PRESSURE    218
#define RO_TEMPERATURE 219

static valid_terms run_option_terms[] = {
                               {"Stages",RO_STAGES},
                               {"Display", RO_DISPLAY},
                               {"Verbose", RO_VERBOSE},
                               {"Init Function", RO_INITFUNC},
                               {"Run Function", RO_RUNFUNC},
                               {"SaveInterval", RO_SVFREQ},
                               {"Save_File Root",RO_SVFILEROOT},
                               {"SwitchDistance",RO_SWITCHDIST},
                               {"MaxNodeSeparation",RO_MAXNODESEP},
                               {"MinNodeSeparation",RO_MINNODESEP},
                               {"SpeedUp",RO_SPEEDUP},
                               {"CellBoundingBox",RO_CELLBBOX},
                               {"SimpleShearOffset",RO_SSOFFSET},
                               {"CumulativeSimpleShear",RO_CUMSSOFFSET},
                               {"Temperature",RO_TEMPERATURE},
                               {"Pressure",RO_PRESSURE},
                               {NULL}
                              };
