#define ELLE_KEY "ELLE"
#define AREA_KEY "AREA"
#define QUALITY_KEY "QUALITY"
#define VISCOSITY_KEY "VISCOSITY"
#define SE_KEY "S_EXPONENT"
#define FLYNN_REAL_ATTRIB_KEY "F_ATTRIB_A"
#define UNODES_KEY "UNODES"
#define U_VISCOSITY_KEY "U_VISCOSITY"
#define DEFAULT "Default"

/*
 * arbitrary minimum area for triangulation
 * smaller values can cause basil to fail array allocation
 */
#define AREA_MIN 0.00008

/*
 * this func is in c_funcs.c
 */
int lines_intersect( double x1, double y1,   /* First line segment */
                     double x2, double y2,

                     double x3, double y3,   /* Second line segment */
                     double x4, double y4,

                     double *x,
                     double *y         /* Output value:
                                * point of intersection */
               );
