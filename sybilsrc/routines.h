void arrow_(float *, float *, float *, int *, int *, int *, int *, int *, int *,
            int *, int *, int *, float *, float *, float *, float *, float *,
            float *, int *, int *);
void arwtwo_(float *, float *, float *, float *, float *, float *, int *, int *, 
             int *, int *, int *, int *, int *, int *, int *, int *, int *, float *,
             float *, float *, float *, float *, float *, int *, char *, int *);
void calcprofil_(float*, float*, int*, int*, float*, float*, float*, 
                 float*, float*, float*, float*, float*, int*);
int  CheckCurrentRec(char*, int*);
int  CheckReferenceRec(char*, int*);
void clipregionu_(float*, float*, float*, float*);
void cntour_(float*, float*, int*, int*, int*, int*, int*, int*, int*,
             int*, int*, int*, float*, float*, float*, float*, float*,
             float*, float*, float*, float*);
int  contour(int*, float*, int**, float**, float*, int*, int, float*,
             float*, char*, int*, int*);
void createcolormap_(int*, int*, int*, int*, int*, int*, int*, int*);
int  create_xy_plot(FILE*, char*);
int  create_arrow_plot(FILE*, char*);
int  create_contour_plot(FILE*, char*);
int  create_optional_label(char*, FILE*);
int  create_profile_plot(FILE*, char*);
void c3code_(float*, int*, int*, int*, int*, int*, int*, int*, float*, float*, 
             float*, float*, float*, float*, float*, float*);
void dashln_(int*, int*);
void drawautolabels_(char*, int*, float*, int*, float*, int*);
int  DrawBndBox(float*, float*);
void drawcolourbar_(int*, int*, float*, float*, float*, float*, float*,
                    float*, float*, float*, int*);
void drawlabel_(float*, float*, int*, int*, int*, char*, int*);
void drawlabelcm_(float*, float*, int*, int*, int*, char*, int*);
/*void DrawPreviewLine(Widget, int, int, int, int); */
int  DrawProfileLine(float*, float*, int);
void drawrectangleu_(float *, float *, float *, float *);
void drawticmarks_(float*, float*, float*, float*, float*, int*, int*);
void DrawTitle(char*);
int  dump_comments(FILE*);
void edges_(float*, float*, float*, float*, float*, float*, int*,
            int*, int*, int*, int*, int*, int*);
void fdput_(int *, float *, float *, float *, float *, int *, int *, int *, int *,
            int *, int *, int *, int *, int *, int *);
void fileopen_(char*, int*, int*);
int  file_open(file_data**, char*, int, char*);
void fillpoly_(float *, float *, int *, int *, int *);
void FindFont(char*, int);
void formatnumber_(char *, float *, int *);
int  GetUserVals(float *, int, unsigned char, int *, char *);
void gravcp_(int *, int *, float *, float *, float *, float *, int *, float *, float *,
             int *, int *, float *, float *, float *, float *, float *, int *);
void hplots_(int*, int*, int*, int*);
void intbnd_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *,
             int *, float *, float *, float *, float *, float *, int *, int *,
             int *, float *, float *, int *);
int  initial_options(FILE*, input_options*, char*);
int  Init_App(int, char**, FILE**);
void init_box(float*, input_options*);
void integrattrpz_(float *, int *, int *, float *, float *, char *, int *);
void integrat_(float *, int *, int *, float *, float *, char *, int *);
void labelcolourbar_(float *, float *, float *, float *, float *, float *, int *);
void linecntour_(float *, float *, int *, int *, int *, int *, int *, int *, 
                 float *, float *, float *, float *, float *, float *, int *, 
                 int *, float *, float *);
void lprint_(float *, float *, float *, int *, int *, int *, int *, float *);
void markcolourbar_(float *, float *, int *, float *, float *, int *, float *,
                    float *, float *, float *);
void meshdf_(int *, int *, char *, float *, float *, float *, float *, int *, int *,
             float *, int *, int *, int *, float *, float *, int *);
void NewMesh();
void ntrpld_(float *, float *, float *, float *, float *, float *, float *, float *,
             float *, float *, float *, float *, int *, int *, int *, int *, int *,
             int *, int *, int *, int *, int *);
void ntrpln_(float *, float *, int *, float *, float *, float *, float *, float *,
             float *, float *, float *, float *, float *, float *, float *, int *,
             int *, int *, int *, int *, int *, int *, int *, char *, int *);
void ntrplnd_(float *, float *, int *, float *, float *, float *, float *, float *,
              float *, float *, float *, float *, float *, int *, int *, int *, int *,
              int *, int *, int *);
void ntrplt_(float*, float*, float*, float*, float*, float*, float*, float*,
             float*, float*, float*, float*, int*, int*, int*, int*, int*, int*,
             int*, int*, int*, int*);
int  ParseOptions(int, char**, char*, char*);
void plmesh_(int*, int*, int*, int*, int*, float*, float*, int*, int*, int*, int*,
             float*, float*, float*, int*, int*, int*, int*, float*, float*, float*,
             float*, float*, float*, float*, int*, float*, float*, int*);
void plotu_(float*, float*, int*);
int  plot_arrows(int, char*, float**, int**, float*, int*, float*, int*, int*, int*);
int  plot_deformation(int*, float*, int**, float**, float*, int*, int, int, float*, 
                      float*, char*, int*);
int  plot_density(int*, float*, int**, float**, float*, int*, char*, int, 
                 int, float*, float*, int*, int*);
int  plot_elle_regions(float*);
int  plot_gravity(int*, float*, int**, float**, float*, int*, int, int, float*, float*, 
                 char*, int, int*, int*);
int  plot_layer(int*, float*, int**, float**, float*, int*, char*, int, int, float*,
                float*, int, int*, int*);
int  plot_mesh(int, int*, float*, int**, float**, float*, int*, int*);
int  plot_rotation(int*, float*, int**, float**, char*, int, int, int, float*, float*, 
                   float*, int*, int*, int*);
int  plot_strain(int*, float*, int**, float**, float*, int*, int, float*, float*, char*,
                 int, int, int*, int*, int*);
int  plot_strain_markers(int*, float*, float**, float*);
int  plot_velocity(int*, float*, int**, float**, char*, int, int, int, float*, float*, 
                   float*, int*, int*, int*);
void plseg_(int*, int*, int*, int*, float*, float*, int*, int*, int*, int*, int*,
            float*, float*, float*, float*, int*, float*, float*, int*);
int  process_log_file(FILE*, input_options*);
int  Profile(float*, int*, int, float*, int, char*, float*, int*, int*); 
int  profile_2D(char*, int, float**, int**, float*, float*, int*, int,
                int*, float*, float*, float*);
void projectxy_(float*, float*, float*, float*, int*, int*, int*);
void projectdeg_(float*, float*, float*, float*, int*, int*, int*);
void rangxy_(float*, int*, int*, int*, int*, int*, int*, int*, int*, float*, int*,
             int*, float*, int*, int*, int*, float*);
void rotor_(float*, float*, float*, float*, int*, int*, int*, int*, int*,
            int*, int*, int*, float*, float*, float*, float*);
int  read_arrow_data(FILE*, char*, int*);
int  read_contour_data(FILE*, char*, int*);
int  read_contour_values(FILE*, char*);
int  read_data(FILE*, char*, int*, float*, int**, float**, float*, float*, float*, float*,
              int*, int, int*, int);
int  read_filename(FILE*, char*);
int  read_label(char*, char*, int*, int*, float*, float*, FILE*);
int  read_options(FILE*, char*, input_options*, unsigned char);
int  read_profile_data(FILE*, char*, int*);
int  read_profile_values(FILE*, char*, int);
int  read_reference(FILE*, int*, float**, int, int*, int, int);
int  read_title(FILE*, char*, char**);
int  read_xyplot_data(FILE*, char*, int*);
int  Run_App(FILE*);
int  Save_plot_data(float*, float*);
void scale_(float*, float*, float*, float*, float*, float*, float*, float*);
int  set_cell(FILE*);
void setfont_(char*, int*, int*);
void setlinewidth_(float*);
void setpencolor_(int*);
int  Set_Colours(int*, int*, int*, int*, int);
int  set_max_cntrs(float*, float);
int  set_records(FILE*, file_data*);
void shadrt_(float*, float*, float*, float*);
void snorm_(float*, int*, int*, float*, float*);
void stmark_(float*, float*, int*, int*, int*, float*, int*, float*, float*, float*, float*);
void stmesh_(int*,int*, float*, float*, float*, float*, 
/*           float*, float*, float*, float*, float*, */
             float*, float*, float*, float*, float*, float*,
             float*, float*, int*, int*, int*, int*, int*, int*, float*,
             int*, int*, float*, int*, int*, int*, int*, int*, int*, int*, 
             int*, int*, int *,int*, int*, float*, float*, float*, float*,
             float*, float*, float*, float*, int*, char*, int *);
void sybflush_();
int  ticstep(float, float, float *, float *, float *);
void unsetclipregionu_();
void vprint_(float *, float *, float *, int *, int *, int *, int *, float *);
int DataPrintLine(char *,float *,float *, float, float, float, float, int, int,
                  int );
int DataPrintMesh(char *label,float *data,int *flags,float *pwindo,
                  int nx3,int ny3,float xref,float yref,int ncomp,
                  int verbose);
int DataPrintArrow(char *label,float *c2t,float *s2t, float *wx,float *wy,
                   int m1,int m2,int mp,int n1,int n2,int np,int nx3,
                   float *pwindo, float xref, float yref, int ncomp,
                   int verbose);
int DataPrintMeshVelo(char *label, float *data1, float *data2,
                      int m1, int m2, int mp, int n1, int n2, int np,
                      int nx3, int ny3,
                      float *pwindo, float xref, float yref, int ncomp,
                      int verbose);
