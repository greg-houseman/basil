.TH basil l
.SH NAME
basil \- 2-D Finite Deformation Simulations
.SH SYNTAX
.BR basil [input.file]
.SH DESCRIPTION
.PP
The program basil uses the finite element method to
calculate quantities that describe stress and strain in non-linear viscous
materials, up to strains of order 100%.  The calculations describe
very viscous Earth materials which undergo irreversible large-strain
deformation at high temperature and over long time periods, under the
influence of body forces and surface tractions.
.PP
The program permits a spatially variable Newtonian or non-Newtonian
viscosity in a 2-D geometry with boundary conditions on traction 
and/or velocity.  It is also possible to include a single fault or
discontinuity in the problem in a dynamically self consistent way.  The
2-D deformation field represents either plane-strain deformation, or
it permits a specified distribution of normal stress in the third
direction.  The latter is referred to as the thin viscous sheet 
formulation when the normal force is due to gravity acting on 
variations of the layer thickness.  Plane-stress calculations are
a specific case of the thin viscous sheet formulation.
.PP
This document provides specific instructions and examples of how to use
basil.  Basil calculations are often quite lengthy and so are normally
run in background, saving a solution to disc at various stages during
the deformation.  Graphical output of these solutions is obtained using the
program sybil (refer man pages on sybil).
.PP
The programs basil and sybil have been developed mainly at 
Monash University since 1988, and before that at ANU and Harvard.
The present set of programs has been developed mainly by Greg Houseman,
Terence Barr and Lynn Evans.  Important contributions by others,
particularly by Philip England during the early stages of development
are gratefully acknowledged.  Comments on, or questions about, basil
and sybil should be sent to basil@earth.monash.edu.au
.PP
A User Manual which gives background theory and worked examples 
of finite deformation problems solved using basil is under development.
.SH PROGRAM SUMMARY
.PP
.IP \(bu 5
.BR BASIL
sets up the finite element mesh, obtains a time-zero solution and
proceeds in a time-stepping mode to compute successive solutions as
finite deformation proceeds.  A subset of the solutions at selected
times are saved to a disc file of specified binary format.
.IP \(bu 5
.BR SYBIL
reads the solution files, and enables the user to
plot selected quantities in the form of contour plots of strain, stress
or displacement quantities, arrow plots of vectors or principal
stresses, mesh plots or other derived quantities.  The program is run
interactively or from a pre-saved set of instructions referred to
as the log file.  The output is directed to the screen in an X-window.
.IP \(bu 5
.BR SYBILPS
produces a postscript output from a set of instructions created by
SYBIL.  The postscript files can be viewed on the screen using
ghostview or printed using lpr.

.SH CAVEAT 
.PP
The purpose of this document is to provide guidance and reference
material to new users.  These programs are under current development,
however, so the specifications set out below apply to the current state
of the programs and probably will change with future development.


.SH SETTING UP A NEW PROJECT AND RUNNING BASIL

.PP
1. Assuming that the executables for basil and sybil have been
set up under /usr/local/basil, define the following environment
variable: setenv BASILPATH /usr/local/basil
and add $BASILPATH/bin to your execution path.
.PP
2. Set up a project directory for files related to the current
basil project, e.g.: mkdir Proj1; cd Proj1
.PP
3. Under the project directory, establish subdirectories to hold
(a) binary solution files: mkdir FD.sols, and (b) ascii output
files: mkdir FD.out.  It is up to the user to decide on names 
and locations for these directories, but our examples will show
the above names used.
.PP
4. Copy in to Proj1 an example of a basil input file.  The example
that you start with would typically be chosen from 
$BASILPATH/examples, firstly to test basil by rerunning
a previously solved problem, and then to be modified as the
input file for a new problem, e.g.: 
cp $BASILPATH/examples/compression/compinput .
.PP
5. Edit or list the input file (assumed now to be called compinput)
in order to examine the structure of the file.  All lines are comments
to the user unless the line commences with a recognised command word.
The first command line in this file is the OUTPUT command.
You will require an input file similar to this for each job that you
wish to run, and you will need to become familiar with the commands
and command arguments described below in order to direct the 
calculations as you require.
.PP
6. Commence execution in background with standard output redirected 
to a file (BAS.scr) and a command line argument of
the input filename, by: nice basil compinput > BAS.scr &
Or the program may be run in foreground with standard output direct 
to the screen by: basil compinput
.PP
7. Multiple input files may be queued for sequential operation by
setting up separate input files for each job, listing the names of
all input files (1 per line) in a file called basil.in, and then
invoking basil without a command line input parameter.
.PP
8. Standard output (directed to BAS.scr in the above example) is
designed so that the user can monitor the progress of the calculations.
It is useful for debugging problems but is not intended as a permanent
record of the calculation.
If the calculation stops, examine the last few lines of standard output
to examine whether the calculation stopped normally or because of
some unexpected problem (such as an input parameter set in error or loss
of convergence).  If running in background, monitor BAS.scr to see
that the calculation is proceeding as expected.
.PP
9. After a solution has been written to disc, graphical output can be
examined using sybil (refer man pages on sybil).  Even if the basil
calculation is not yet finished, solutions that are already written
may be examined to check that the calculation is proceeding as expected.
.PP
10. Each input file that is acted on by basil creates 3 or 4 output
files, with name determined by the BIN parameter in the OUTPUT command.
The primary output files are the binary solution file saved in the 
directory BINDIR and the .out file saved in the directory OUTDIR.
These two files provide a relatively complete record of a given
task.  The input file is copied unchanged to the beginning of the .out
file, so it is unnecessary to also archive the input file separately.

.SH SETTING UP AN INPUT FILE
.PP
A large number of different parameters are required to define a
basil finite element calculation.  The values of these
parameters are controlled entirely by the input file for a given
calculation.  Although all parameters in basil have default values,
the defaults in almost all cases define the simplest action.
An example of an input file is shown below.  Parameter
values are organised in groups defined by a Command Keyword.
Each parameter that requires a non-default value is redefined here.
The simplest way to set up a new input file is to edit an
existing file from a previous job.
.PP
Basil will execute in sequence the set of input files that are listed 
(one per line) in the file 'basil.in'.  This file may be modified 
during execution, so that jobs may be added to the queue, or removed 
before execution.  Input files may also be edited while in the queue
before execution.  If only one input file is to be run, basil may be
called with the name of the input file as a parameter.  The .out file
will indicate unambiguously which parameter values were used in a given
calculation.

.SH EXAMPLE OF AN INPUT FILE
.PP
.EX
 #1. Setup commands - 
 #---------------------------------------------------------
 OUTPUT   BIN=indent BINDIR=FD.sols OUTDIR=FD.out DUMP=debug
 LABEL    Continental collision, rigid eastern boundary
 MESH     NX=16 NY=40 TYPE=0 FAULT=0
 #READ     BIN=output, RECORD=2, RESET
 GEOMETRY XLEN = 1.4, YLEN = 2.0, NCOMP = 0, IGRAV = 4
 VISDENS  SE = 3.0
 FORCE    RHOG= 1.0 ;
 LAYER    THICKNESS, ROTATION &
          ARGAN=1.0, BRGAN=0.0, THRESH=0.0314285714
 BCOND
 SOLVE    AC = 5.0E-5, ACNL = 5.0E-3, WFIT = 1.0 &
          NITER = 50000, ITSTOP = 100
 LAGRANGE MARKERS
 SERIES
 #DEFORM   TYPE=0, BANGL=90.0, ERA=0.125, DEFV=1.0 ;
 
 #2. Timestep commands
 #---------------------------------------------------------
 STEPSIZE TYPE=EXPLICIT, IDT0=20, MPDEF=5
 SAVE     KSAVE=50, TSAVE=0.1
 STOP     KEXIT=100, TEXIT=0.5
 REMESH   IGRID=3

 #3. Boundary condition commands
 #---------------------------------------------------------
 #The following lines state Boundary Conditions
 #for the continental collision problem, as used in the
 #Rubey paper (rigid eastern boundary).

 #North, East and West boundaries are made rigid
 ON X = 1.4 : UX = 0.0
 ON X = 1.4 : UY = 0.0
 ON Y = 2.0 : UX = 0.0
 ON Y = 2.0 : UY = 0.0
 ON Y = 0.0 : UX = 0.0
 ON Y = 0.0 : UY = 0.0

 #The X-component velocity is tapered either side of
 #a central indenting block.
 #Taper function:
 #TP = 2 is cos**2 (zero gradient at both ends)
 
 ON X = 0.0 FOR Y = 1.5 TO 2.0   : UX = 0.0
 ON X = 0.0 FOR Y = 1.25 TO 1.5  : UX = 1.0 TO 0.0 : TP = 2
 ON X = 0.0 FOR Y = 0.75 TO 1.25 : UX = 1.0
 ON X = 0.0 FOR Y = 0.5 TO 0.75  : UX = 0.0 TO 1.0 : TP = 2
 ON X = 0.0 FOR Y = 0.0 TO 0.5   : UX = 0.0
 ON X = 0.0 : UY = 0.0
 
 #4. Physical Parameter commands-
 #---------------------------------------------------------
 #Define a special region representing the Tarim Basin
 REG E { 461,501,541,422,462,502,542,383,423,463,&
         503,344,384,424,464,345,385,425 } : &
         VR = 0 ; VC = 10.0 ;
 
 #5. Strain marker commands
 #---------------------------------------------------------
 #define a set of elliptical strain markers 
 MARKERS  ROWS=5,COLS=10,R=0.02,XMIN=0.1,XMAX=0.6, &
          YMIN=0.5,YMAX=1.5
.EE

.SH INPUT FILE SPECIFICATION
.PP
The input file consists of a sequence of commands separated by optional comments.
.PP
Each command is identified to the program by a command keyword.
Any line not commencing with a command keyword is treated as a comment.  A command line 
may be commented out by inserting '#' in front of the keyword.  Valid Command
keywords are listed below.  Except for the OUTPUT command (which is required)
each command is optional.
.PP
Each command has its own set of valid parameter keywords.  Parameter
keywords identify real, integer, or string  variables.  
Some parameter keywords identify an array of real variables (terminated by ';').  
Default values are set for all parameters, so each parameter is optional.
.PP
Parameters are separated by blank space or comma.  Avoid using tabs in
the input file.  Arrays must be terminated with a semicolon preceded by
a blank space ' ;'.  All valid keywords are in uppercase.
A continuation character '&' is used to continue commands across
multiple lines.
.PP
In the event of an unrecognised parameter, or a syntax error while reading 
a command, the command line will be echoed to Standard Out with an indicator 
as to where the error occurred, and processing of the current input file stops.
.PP
Only the first command (OUTPUT) is essential - the program will use
defaults for all other parameters, unless overridden by command.  
An input file is disabled by commenting out the OUTPUT command.
.PP
Most of the basil commands are used to set parameter values and/or to
activate subroutines.  
These actions occur in the order in which they are given in the input file.  
The relative sequence of statements such as BCOND and DEFORM is significant.  
In this sequence, the boundary conditions will be applied to a rectangular region 
which will then be deformed.
.PP
Each command should only appear once in the file, with the exception of 
boundary condition (ON), physical parameter (REG) and strain marker (MARKER)
commands for which multiple statements are expected.  The effect of these
statements is cumulative.
The input file is opened separately to read the boundary conditions (ON statements), 
physical parameters (REG statements) and strain marker data (MARKERS statements).  
.PP
There is a recommended sequence for the input instructions as follows 
(although in some cases the order does not affect the results of the calculation):
.EX
 1. setup commands
 2. timestep commands
 3. boundary condition commands
 4. physical parameter commands
 5. strain marker commands
.EE

.SH SETUP COMMANDS
.PP
Valid Setup commands are:
.PP
OUTPUT   MESH     POLY      LABEL     READ      GEOMETRY   FORCE
LAYER    BCOND    DEFORM    VISDENS   LAGRANGE  SOLVE      SERIES
.PP
These commands define both the physical problem under investigation
and control parameters.
Refer to the above input file example for illustration of how some
of the commands are used.  Each command is specified as follows, with
default values of parameter values given in brackets:
.PP
OUTPUT   BIN = filename  OUTDIR = directory  & 
         BINDIR = directory  DUMP = filename
.PP
The OUTPUT command directs output to nominated files and directories.
This command is essential; the remainder of the input file is skipped
if the first command is not the OUTPUT command.  BIN is the name of
the binary solution file; it is filed in directory BINDIR [FD.sols].
DUMP [debug] (also in BINDIR) is the name of a binary solution file containing 
only the most recent solution record.  Its primary use is for debugging and/or monitoring progress of a long job.  Two ascii files BIN.out and BIN.dat are also created and filed in OUTDIR [FD.out].  BIN.out is the primary 
record of the calculation.  BIN.dat contains timeseries data, only if 
defined by the SERIES command (refer below).  If these files already exist they will be
overwritten without further checking, so it is prudent to disable an input
file for a long job that has completed successfully.
.PP
MESH   NX = #, NY = #, TYPE = #, FAULT = #
.PP
The MESH command defines the finite element mesh on which the problem
is to be solved.  If TYPE = 0,1,2 it defines the number of finite element
triangles used in the calculations and therefore directly determines the
size of the matrix equation that is solved by basil, and hence the 
computational resource requirements (memory and execution time) of the
current job.  The mesh is generated within a rectangular region which may
be later deformed into a different shape (refer DEFORM below).  NX and NY
define the approximate number of triangles spanning the rectangle in each of
X and Y directions.  Either NX [0] or NY [0] must be set. If you leave
NX (NY) unset, then it will be set by the program using the chosen value of
NY (NX) to give triangles that are most nearly equilateral.  NX and NY
should both be even numbers > 4.  TYPE [0] defines the topology of the
triangular mesh:
    0 for triangles with hexagonal symmetries,
    1 for a rectangular mesh subdivided by diagonals in one direction,
    2 for a rectangular mesh subdivided by diagonals in alternating directions.
    3 for a mesh defined by polygons in a poly file (see POLY) 
FAULT [0] defines whether the region is cut by a fault (with duplication of nodes
on either side of the fault).  FAULT = 0 means no fault. FAULT = 1 defines a fault on both of
the X boundaries and a solution that is periodic in the X direction. FAULT = 2  defines
a fault on X=XLEN and the mesh is duplicated on the other side of the fault.
.PP
POLY   FILE = filename
.PP
If the MESH TYPE is 3 then the program will generate a mesh by
triangulating the polygons defined in this file. The format must
conform to the .poly file format defined for Triangle (J.R.  Shewchuk,1996)
First line: <# of points> <dimension (2)> <# of attributes> <# of boundary markers (0/1)>
Following lines: <point #> <x> <y> [attributes] [boundary marker]
One line: <# of segments> <# of boundary markers (0/1)>
Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
One line: <# of holes>
Following lines: <hole #> <x> <y>
Optional line: <# of regional attributes and/or area constraints>
Optional following lines: <constrain #> <x> <y> <attrib> <max area>
.PP
LABEL  This comment is stored in the solution file header
.PP
The LABEL  command is used to input a 1-line comment which is saved in
the binary solution file header.  The comment is up to 80 characters in length, and does not require any delimiting marks.
.PP
READ  BIN = filename  RECORD = n  RESET/CONTINUE
.PP
The READ command is used to read in an old solution stored in the
file identified by BIN.  BIN must be specified and is assumed to be in BINDIR,
as specified in the OUTPUT statement. A given solution consists of a number of complete records, each saved at a different timelevel.  The parameter RECORD [1] 
is used to select the nth complete record in the solution.
The keyword RESET (default) causes time and step level parameters to be reset to zero.
The keyword CONTINUE assumes that the calculation is to be
continued without change to the controlling parameters.
.PP
GEOMETRY  XLEN = #, YLEN=#, NCOMP = #,  IGRAV= #
.PP
The GEOMETRY command sets parameters that define the shape of the solution region, the 2D flow assumption and the direction of gravity.
XLEN [1.0] and YLEN [1.0] are the X and Y-lengths respectively of the initial
rectangular solution region.
The length unit may be chosen arbitrarily to represent any real length.
NCOMP [0] specifies whether plane-strain (NCOMP=1) or thin viscous 
sheet (NCOMP=0) formulation is assumed.  IGRAV [0] defines the direction 
of gravity: +/-1 = +/-X, +/-2 = +/-Y, 3=arbitrary direction in the (X-Y) 
plane (not implemented yet), or 4 = - Z direction).  
.PP
FORCE  RHOG = #[#[#[#[#]]]] ;
.PP
The FORCE command defines the values of variables which are used
in the body force calculation.  RHOG is an array of between 1 and 5
numbers (note: must be terminated with ' ;').   RHOG(1) [1.0] 
is a scaling constant for the body force term.  The other 4
parameters in RHOG are not yet in use.
.PP
LAYER  THICKNESS, ROTATION, ARGAN = #, BRGAN = #, &
       THRESH = #, HLENSC = #, BDEPSC = #
.PP
The LAYER command is used to set parameters that are only used
in the thin viscous sheet formulation (when NCOMP=0).  THICKNESS
and ROTATION are two flags which, if set, indicate that the
layer thickness and rotation arrays (respectively) are to be 
calculated and saved in the binary solution file.
HLENSC [50.0] is the ratio of length unit to layer thickness
(refer GEOMETRY statement for choice of length unit).
BDEPSC [0.35] is the ratio of buoyant layer thickness to strength
bearing layer thickness (crustal thickness to lithospheric thickness
in geological terms).
ARGAN [0.0] is the Argand number (initially assumed uniform) 
that determines the magnitude of the buoyancy forces that act 
on variations in thickness of the buoyant layer. If the layer 
thickness exceeds a threshold THRESH[1.0], the Argand number 
is locally increased to ARGAN + BRGAN [0.0].
The latter two parameter are used to represent the sudden 
local increase in Argand number associated with the phenomenon 
of convective thinning of the lithosphere.
.PP
BCOND  IOFF=#, VELXO=#, VELYO=#, YLDSTR=#
.PP
The BCOND command tells the program to set boundary conditions at this
point in execution.  The program will then immediately read and interpret
all of the ON statements entered in the 3rd section of the input file (see
below).  In addition, periodic boundaries may apply in the X-
direction if IOFF[0] is set to 1.  On the periodic boundaries
there is an offset of VELXO[0.0] in the x velocity and VELYO[0.0] 
(note spelling: 'O', not zero) in the y velocity. YLDSTR[0.0] is 
the tangential stress on this boundary if it represents an unlocked
fault.
.PP
DEFORM TYPE=#, DEFV=#[,#[,#[,#[,#]]]] ; BANGL=#, &
       ERA=#, TILTDEG=#, NDFORM=#
.PP
The DEFORM command instructs the program to distort the initial
geometry of the solution region, by distorting the external boundary
of the initially rectangular 2D solution region, or perturbing the
thickness of the thin viscous layer.  These distortions do not disturb
the topology of the finite element mesh, so they have little impact on the 
formulation of the problem and so represent an easy way to address
problems with somewhat more complex geometry than the basic rectangle.
Each of the distortions is specific to a particular problem that has
been studied in the past.  To add a new kind of distortion, add another
IF THEN block to subroutine CFORM in cform.f, triggered by a new value of
the variable TYPE [0], and use the set of control parameters DEFV [0] to
control the effect of the distortion.
.PP
TYPE = 0 - no distortion will occur.  
TYPE = 1 - the middle section of the X = 0 boundary is stretched 
onto a circular arc, with a radius of curvature given by the variable 
ERA [0.0].  
TYPE = 2 - combines the effects of TYPES 1 and 3.
TYPE = 3 - causes the X = 0 boundary to be rotated clockwise through 
an angle in degrees BANGL [0.0].  Side boundaries and internal mesh
are stretched or squeezed as necessary.
Types 1-3 were designed for setting the initial geometry in the
indenter problem.
TYPE = 4 - takes the initial rectangular region and rotates part of it
clockwise through an angle in degrees BANGL [0.0].  The rotation is zero for
X < DEFV(1) and increases linearly to BANGL for X > DEFV(2).
Type 4 was designed to provide an initial condition for a subducted
slab problem.
TYPE = 5 - distorts the Y = 0 boundary of the rectangle by deflecting it
in the form of a harmonic perturbation with amplitude DEFV(1) and relative 
wavenumber DEFV(2).  The entire mesh may be translated in the -Y direction
by distance DEFV(3) and that part of the mesh with Y < 0 may be stretched
by factor DEFV(4) in the -Y direction.
Type 5 was designed to provide initial conditions for the Rayleigh-Taylor
problem, either with no lower layer (DEFV(3)=DEFV(4)=0), or with a lower
layer covered by a relatively coarse mesh (DEFV(3)=-1, DEFV(4)=4).
TYPE = 6 - distorts the initial layer thickness distribution, without
changing the external geometry of the mesh.  The perturbation to
crustal thickness is in the form of an X-direction wave, amplitude
DEFV(1), +/- a Y-direction wave, amplitude DEFV(2), or a hexagonal planform,
amplitude DEFV(3).  
TYPE = 7 - stretches the mesh within an unchanged external boundary, in
order to cause element boundaries to coincide with a region with
anomalous strength.  For X < DEFV(1), the Y-direction displacement of
the mesh nodes is maximum DEFV(4) on X = 0, between Y = DEFV(2) and
Y = DEFV(3).  The displacement is linearly tapered to zero towards
X = DEFV(1), Y = 0, and Y = YLEN.
TILTDEG > 0 distorts the rectangle into a parallelogram, rotating
the vertical sides of the rectangle through angle TILTDEG [0.0]
NDFORM[0] is used to specify compression of the mesh around fault tips.
.PP
VISDENS SE = #, RHO = #
.PP
The VISDENS command firstly is used to set background values of
the strain-rate vs stress exponent SE [1.0], and the density variable
RHO [1.0] which is used to determine the body force.  The VISDENS command
also cause the program to immediately read and interpret the region
statements (see REG below) which define more complicated distributions
of the physical variables defining viscosity coefficient and density.
.PP
LAGRANGE MARKERS, LGMESH, NXL = #, NYL = #, & 
         NPMP = #, NRMP = # FILE = filename
.PP
The LAGRANGE command initiates generation of Lagrangian markers
which provide a simple means of recording the finite deformation
of the medium by tracking the movement of tracers that are passively
advected by the deforming medium.
LAGRANGE has two optional arguments.  The argument LGMESH causes
an array of tracers on a uniform triangular mesh to be initiated.
NXL [16] and NYL [32] determine the number of tracers in each of X and Y directions.
The argument MARKERS causes the program to immediately read 
and interpret all the MARKER statements (see below)
which initiate arrays of circular strain markers that may be later
analysed to determine finite strain ellipsoids at specific locations.
NPMP [80] and NRMP [7] define the number of tracer points in each circular strain marker
and the number of points along the radius of the marker respectively.
These parameters are used to allocate memory for the markers.
FILE specifies a file where the current mesh will be saved at each
timestep. If the file exists, the Lagrange mesh will be read from it
rather than generating a new mesh. At each timestep.
.PP
SERIES  TIME, XX=#,# ; YY=#,# ; UX=#,# ; UY=#,# ; &
        XMIN, XMAX, YMIN, YMAX, UMIN, UMAX, VMIN, VMAX &
        TH=#,# ; THMIN, THMAX, PRES=#,# ; PRESMIN, PRESMAX
.PP
The SERIES command specifies that the BIN.dat file will be created.
This file is ascii format, containing one line of data values per
timestep.  The arguments to the SERIES command are all optional.
The selection and order of these arguments defines the data values
to be written on each line of the BIN.dat file.  The file does not
include a header; the selection of data values is recorded in the BIN.out
file under the SERIES command.  It is strongly recommended that the
first parameter selected is TIME.  XMIN, XMAX, YMIN, YMAX designate
the minimum and maximum values of the finite element node coordinates
in X and Y directions respectively.  UMIN, VMIN, UMAX, VMAX designate
the minimum and maximum values of the finite element node velocity
components respectively.  THMIN and THMAX designate minimum and maximum layer 
thickness.  PRESMIN and PRESMAX designate minimum and maximum
pressure values (and will be ignored if NCOMP=0). For all of the MIN and
MAX values, sign is significant, i.e. the MIN is the most negative and
the MAX is the most positive.
The keywords XX, YY, UX, UY, TH, or PRES indicate that the X coordinate, Y
coordinate, X velocity, Y velocity, thickness, or pressure at a particular
node is to be recorded. These keywords are therefore followed by the X 
and Y values which specify the node (each pair terminated by ' ;'). 
Note that it is simplest to identify such nodes before the mesh 
has been deformed by either the DEFORM command, or by timestepping.
If there are no parameters to the SERIES command, the default parameters are:
TIME, YY=0.0,0.0 ; UY = 0.0,0.0 ; YY = xlen,0.0 ; UY = xlen,0.0 ;
where xlen is the value defined by XLEN in the GEOMETRY statement.

.PP
TIMESTEP COMMANDS
.PP
The distinction between Setup commands and Timestep command is
a little blurred, but timestep commands are intended to affect
the definition of the problem during timestepping.  In that sense
they define operations that modify physical parameters or boundary
conditions at each timestep.  Valid Timestep commands are:
.PP
SOLVE STEPSIZE  STOP  SAVE  REMESH  BCONDMOD  DENSMOD  RHEOMOD
.PP
These commands control aspects of the timestepping.  Refer to the
example input file above for illustration of how the commands are
used.
.PP
SOLVE AC=#, ACNL=#, WFIT=#, NITER=#, ITSTOP=#
.PP
The SOLVE command specifies the parameters which control convergence of
the solution at a given time level.  AC [0.001] sets the convergence
criterion for the linear conjugate gradient iterations; ACNL [0.1]
sets the convergence criterion within the non-linear iterations
(only relevant if SE not 1.0).  In each case, a smaller value implies a more
accurate solution, but it is possible to set these numbers so small
that convergence is impossible in the presence of round-off error.
The optimum values are found by convergence experiments.
WFIT[1.0] may be used to damp instability in the non-linear 
iterations by combining a fraction of the new solution (WFIT) 
with a fraction of the old solution (1-WFIT). 
NITER[NROWSP*3] sets an upper limit on the number of iterations 
allowed in the conjugate gradient solver, ITSTOP[30] sets 
an upper limit on the number of iterations permitted in the 
non-linear iterative loop (only relevant if SE not 1.0).
.PP
STEPSIZE  TYPE = RK/EXPLICIT,  IDT0 = #,  MPDEF = #
.PP
The STEPSIZE command controls the type of timestep algorithm used and
sets parameters which govern the size of the timestep.  The TYPE
[EXPLICIT] variable specifies the type of timestep (RK = 2nd order
Runge-Kutta or EXPLICIT = 1st order forward-time difference method).  
The maximum size of the timestep is 1/IDT0 [25], and the maximum 
dimensionless strain which may occur in one timestep is 1/MPDEF [5].
.PP
STOP    KEXIT = #, TEXIT = #, IWRITE = #
.PP
The STOP command allows the user to specify a number of independent
criteria which will cause execution of the current input file to terminate
normally.  In general the program will then proceed to execute the next
input file in the queue.  Execution completes if 
(a) the number of timesteps >= KEXIT [100], or 
(b) the dimensionless time level >= TEXIT [0.4], or 
(c) the number of records written to a file exceeds IWRITE [10].
The latter criterion is designed to prevent a runaway job filling
up the output disc.  The job will also complete normally if the
magnitude of the timestep required to continue the computation 
becomes less than 0.002.  This parameter is not externally adjustable.  
If KEXIT=-1 the job will terminate normally after the mesh has been
created but no solution calculated. This is useful to check that the mesh
generation is satifactory.
Other completion criteria may be implemented later.
.PP
SAVE    KSAVE=#, TSAVE = #
.PP
The SAVE comand specifies how frequently a solution record is written 
to the binary solution file.  Each binary solution file consists of a 
sequence of solution records at different timelevels.  Each solution 
record is complete and contains all of the information necessary to 
reconstruct all aspects of the solution at that time level.  A solution
record is written if either of the following criteria is satisfied 
(a) (current time level - time level of solution last saved) >= TSAVE [0.08] 
(b) (current step number - step number of solution last saved) >= KSAVE [4]
.PP
REMESH  IGRID=# SEGMIN=# SHAPEMIN=#
.PP
The REMESH command controls whether the mesh regridding (or conditioning)
routines are called at each timestep.  IGRID ([0] if  FAULT=0, else [1])
specifies the type of regridding algorithm.   IGRID = 1 indicates
regridding of only the nodes on a fault.  IGRID = 2 indicates
regridding of the entire mesh (keeping the external boundary fixed).  
IGRID = 3 indicates a partial reconditioning of the mesh, in which the
most deformed triangles (least equilateral) are adjusted by moving
vertex nodes so as to improve the shape of the triangles.
This mesh conditioning allows a calculation to be extended under conditions
of extreme mesh deformation, but since it does not affect the topology of
the finite element mesh, the benefit is limited.  Properties such as 
layer thickness are re-interpolated onto the nodes that are moved.
Caution is required in using these options since they have only been
tested on a small number of example calculations. IGRID = 4 indicates complete remeshing within the limits of the external boundaries and preserving any internal boundaries.  The applied boundary conditions are preserved and the properties are interpolated onto the new mesh. SEGMIN and SHAPEMIN control the frequency of remeshing when IGRID is set to 4. SEGMIN sets a segment length threshold below which a node is removed from the boundary before remeshing.  SHAPEMIN is a limit for triangular element distortion based on the sum of the angle cosines (1.5 for equilateral elements and decreasing to 1 as the shape degenerates).
.PP
BCONDMOD  TYPE=#, OMTOT=#, BCV=#, [#, [#, [#, [#]]]];
.PP
The BCONDMOD command is used to instruct the program to modify the boundary 
conditions at every time level for problems in which time-dependent boundary 
conditions are required.  Initial boundary conditions are set by the BCOND 
statement (refer above).  If the velocities and or tractions on all boundary 
segments are constant with time, BCONDMOD should not be used.
As for the mesh deformation operation (DEFORM) the following operations have 
all been implemented for specific problems defined by TYPE [0].  Other types of boundary 
condition modification will require alterations to a set of source routines.
BCV [0] is an array of up to 5 variables used to specify the details of any of these operations.  
TYPE = 1 causes the lithostatic boundary condition (assumes NCOMP=0) to be
applied to any boundary on which both components of traction are specified.  
For this type of boundary, the constant normal traction must be resolved
into components whose relative magnitude depends on the current orientation
of the boundary segment.  These boundary traction components are re-evaluated
and reset at every timestep.
TYPE = 2 causes part of the boundary to be rotated as it is moved horizontally
(rotating indenter).  The rotation in degrees of the indenter (per 0.5
time units) is specified by the variable OMTOT [0.0].  For a rotating boundary
segment, the velocity components on this boundary element change continuously 
with time so they are reset after every timestep.  This option will only work 
under specific conditions.
TYPE = 3 causes both rotating indenter and the lithostatic boundary conditions
to be applied (to different boundary segments).  
TYPE = 4 causes specified boundary nodes to be changed from fixed velocity to 
fixed traction.  Nodes whose X-coordinate exceed BCV(2) and whose Y-coordinate
exceed -2 have the Y-component of traction set to BCV(3).
TYPE = 5 causes the velocity components on all boundary segments to be 
set to zero (used to switch off the indenter).
TYPE = 6 is similar to TYPE = 4, except that the change from fixed velocity
to fixed traction occurs when the X-coordinate exceeds BCV(2)+BCV(4)*TIME.
TYPES 4 and 6 are used to represent the effect of subduction on a dense slab.
.PP
DENSMOD
.PP
The DENSMOD command causes the density distribution to be modified at
every time level.  The density distribution is initially set using the
VISDENS command and REG statements and the density distribution is
advected with the deforming medium.  If the density is constant with
time at all points the DENSMOD command should not be used.  The DENSMOD
command is used to reset the density at a subset of mesh points, for
example to represent the effect of a depth-dependent phase transition 
as the mesh passes a given depth level. This command is not yet fully
implemented.  At present DENSMOD has no arguments, but invoking it
causes the subroutine DENSIT to be called at every timestep, and this
subroutine causes density values to be reset from RHOG(2) to RHOG(3)
when the Y-coordinate of a point becomes less than RHOG(4).  At present
the RHOG parameters are set using the FORCE command (refer above).
.PP
RHEOMOD
.PP
The RHEOMOD command is not yet implemented.  It is intended that this
command will be used to modify the rheological parameters during the
course of a calculation in a manner that is analogous to the way that
DENSMOD changes the density distribution.  Rheological parameters are
set using the VISDENS command and the REG commands, and at present
these parameters remain fixed to the deforming material throughout
the calculation.
.PP
BCONDMOD, DENSMOD and RHEOMOD assume that the boundary conditions, 
density, and rheological parameters are defined in the setup phase, 
and only require minor amendment during the timestepping loop.  These 
three commands work in a similar manner to the DEFORM command.  They 
are governed by an option number which specifies the particular
modification required, together with a number of parameter values that
are adjustable.  Specific code modifications will probably be required
for each new kind of modification operation that is required.
.PP


.SH BOUNDARY CONDITION COMMANDS
.PP
The third section of the input file contains the specification of
boundary conditions, all recognised as statements that commence with
the keyword ON.  A simple boundary condition statement looks like
.PP
ON X = 0.0 : UY = 1.0
.PP
which simply translates as the statement that on the boundary X = 0,
the Y-component of velocity is required to have the value 1.0.  The
routine that interprets these commands (VSBCON) determines the
bounding box for the rectangular solution region so that ON X = XMIN,
ON X = XMAX, ON Y = YMIN and ON Y = YMAX are valid instructions.
The rectangular region may be distorted into a non-rectangle
after boundary conditions have been set (using the DEFORM command, or
by timestepping); in that case the boundary conditions applying to a
node after it is moved are unchanged by the movement of the node.
.PP
For the problem to be completely specified in an analytical sense, we
require that every node on the external boundary is associated with 2
boundary conditions, one applying to the Y-component of velocity or
traction, one applying to the X-component of velocity or traction.  If,
after all boundary condition statements have been read, the program has
not detected the necessary number of boundary conditions, it will
finish execution of the current input file with a message which specifies
those boundary nodes on which boundary conditions have not yet been
specified.  If that occurs, check that the required boundary conditions
are correctly set, and restart the job.  Any statement commencing with
the keyword ON is recognised as a boundary condition specification.  
You can use as many ON statements as you require, but in general, at
least 8 are required (2 for each of the four boundaries).  If a
particular node(s), e.g. at a corner, are specified by two or more ON
statements, the latter specification takes effect.
.PP
The first part of an ON statement specifies the boundary segment to which the condition 
applies; it is terminated by a colon ":", and takes the form
.PP
ON X = # FOR Y = # TO # :    or   ON Y = # FOR X = # TO # :
.PP
the phrase commencing with the keyword FOR is optional, but is used to
restrict the application of the boundary condition to the line segment
between the 2 variables, or to specify the endpoints of a tapered
boundary condition (see below).  Note punctuation precisely; include 1
or more spaces between every syntax element.
.PP
Following the domain specification in the first part of the ON
statement is the 2nd or action component of the statement, of which an
example is:
.PP
TX = # TO # : TP = #
.PP
where TX in the above example specifes the X-component of traction.  In
place of TX, you can set TY (y-component of traction), 
UX (x- component of velocity), UY (y-component of velocity), 
FX (x-coefficient of friction), FY (y-coefficient of friction).  
In the latter cases (FX, FY) the component of traction T is 
linearly related to the component of velocity u, by the coefficient 
of friction m:
.PP
			T = - mu
.PP
In general you must specify, one of UX, TX or FX on each boundary node,
together with one of UY, TY, or FY on each boundary node.
.PP
If the boundary condition is constant along the boundary, the remaining
part of the ON statement beginning with the keyword TO 
is not required.  This last part of the ON statement is used to specify a
boundary condition which varies along the boundary segment.  
The key symbol TP defines a type of taper
function (TP = 1 for linear; TP = 2 for cos squared; TP = 3 for cos; TP = 4
for sin) that defines a smooth variation from the first value at the 
first end of the boundary segment to the 2nd value at the other end.
Application of the taper requires end points specified by the FOR phrase 
as described above.  For each of the trigonometric functions, the first
1/4 of a wavelength is used and TP = 2, 3, or 4 is chosen to obtain 
zero gradient at both end points, at the first end point, or at the 
second end point of the boundary segment, respectively.  
.PP
Another type of boundary is possible.  The fault boundary is designated
by the keywords FL and FU, which specify that the boundary segment is
part of a fault, either locked or unlocked respectively.  A locked
fault (FL) is specified by the condition that both components of
traction and both components of velocity are continuous from one side
of the fault to the other side.  An unlocked fault (FU) is specified by
the condition that the normal components of velocity and traction are
continuous, but the tangential component of traction is zero and the
tangential component of velocity may be discontinuous.  These keywords
may only be applied to a boundary which has been specified as a fault.
Periodic boundaries in one direction may be implemented using an external 
fault which is specified to be locked.

.SH PHYSICAL PARAMETER COMMANDS
.PP
This section of the input file is used to specify a priori variation of
the physical variables that describe the distribution of the
rheological coefficient B(x,y), the stress exponent n, and (where
required) the density distribution rho(x,y).  Any part of the solution 
region may be defined to have a particular distribution of physical
properties by means of the REG statements.  Multiple statements may be
used to define complicated distributions, with those statements coming
last in the sequence taking precedence.  These commands are only
required if one or more of the following statements applies to the problem 
under investigation:

(a) the stress exponent n  is spatially varying
(b) the dimensionless rheological coefficient B' differs from B' = 1, anywhere 
(c) the dimensionless density rho'  differs from rho' = 1, anywhere.
The REG statement is used to reset values of strain-rate vs stress exponent,
SE, rheological coefficient VC, or density RO, anywhere in the solution
region.  The background values of these parameters are set using the 
VISDENS command (refer above).
.PP
The physical parameters SE, VC and RHO are assumed to be piecewise continuous, 
and are set by means of a REG statement (stands for REGION).  An example is:
.PP
REG A : VR = 2,  VC = 0.1, 0.0, 2.302585 ;
.PP
The first part of the REG statement specifies the domain of applicability,
terminated by the colon ":"  Only points within this domain of 
applicability are affected by the subsequent statement.  There are 3 ways to
specify the geometric region:
.PP
REG A :  specifies the entire solution region
.PP
REG E {n1, n2, n3, .....} : specifies a list of element (E) numbers
You need to make a map of element numbers using sybil, in order for
this mode to be useful, and the element numbers change unpredictably
if you change the number of node points used in the calculation, so this
is not the recommended mode.
.PP
REG P {x1, y1, x2, y2, x3, y3, ......xn, yn} : specifies a polygonal
(P) region defined by the n vertices in the 2-D solution region
When reading this statement, points are read until the closing brace
('}') so the continuation mark ('&') should not be appended to these lines.
.PP
REG statements in all cases apply to entire triangular elements.  A
particular REG  P statement applies to a particular triangular element
if the centre of gravity of the element is contained within the
polygonal area.
.PP
Following the specification of the domain of applicability, the following
key symbols are recognised:
.PP
SE = #,  sets the value of the strain-rate vs stress exponent within the domain.
VR = #,  defines the viscosity rule to be applied in setting the viscosity
coefficient VC.
DR = #,  defines the density rule to be applied in setting the density RO.
In each case VR and DR may take values of 0 for constant VC or RO, 1 for linear
variation with respect to the spatial coordinates, e.g. 
RO = RO(1)+X*RO(2)+Y*RO(3), or 2 for exponential variation,
e.g. VC = VC(1)*EXP(X*VC(2)+Y*VC(3)).  VR = 3 cause the current viscosity
coefficient to be multiplied by VC(1).
The arrays VC and RO (terminated with ' ;') specify the coefficents to
be used in setting the linear or exponential variation of density or
rheological coefficient.

.SH STRAIN MARKER COMMANDS
.PP
The final section of the input file is used to specify strain-markers
if they are required, as specified by the LAGRANGE command (refer above).  
A MARKER command is used to specify an array of circular strain markers,
each composed of a set of points that are advected passively by the
deforming material.  Multiple MARKERS commands may be used, and their
effect is additive, but the total number of markers that may be set
is limited by the NRMP parameter in the LAGRANGE statement.  
An example of a strain marker command is:
.PP
MARKERS  ROWS=5, COLS=10, R=0.02, XMIN=-0.07, XMAX=0.41, &
         YMIN=0.545, YMAX=0.765
.PP
Each MARKERS command causes an array of circular strain-markers to be
inscribed in the deforming medium.  Each strain-marker consists of a
set of points in a circular distribution of radius R.  With finite
deformation the circles are deformed, into ellipses if the strain
field is homogeneous, or into more complex shapes if the strain field
is more complex.  The array of strain markers is distributed on an
initially rectangular grid, with COLS markers across the page between
XMIN and XMAX, and ROWS markers down the page between YMIN and YMAX.
For a more complex distribution of strain markers, use multiple MARKERS
statements. If R<=0 a single line marker with endpoints XMIN, YMIN and
XMAX, YMAX is created.

.SH RUN-TIME PROBLEMS
.PP
.SH a) failure to converge
.PP
The program may give a message that it cannot converge on a solution.
If it occurs at time zero this symptom may indicate an inconsistency in
the specification of the boundary conditions.  If it occurs after large
deformation it probably means that the mesh has undergone extreme
deformation in some part of the solution region.  In that case the
calculation is completed as far as BASIL is concerned.  You may try 
using REMESH  IGRID = 3 and rerunning the job.  This flag permits a limited form of
mesh regridding to occur during the timestepping and therby permits the
calculation to proceed further than it otherwise would.  At present
there is no simple way in basil to create a new well-conditioned mesh
inside an extremely deformed boundary of arbitrary geometry.
.PP
The first step in debugging an unsuccessful run is to determine
whether a time-zero solution was obtained.  If the time zero solution
is available open it in sybil, plot the mesh to check that the geometry
is as expected, and plot the velocity field to check that the 
deformation is consistent with the boundary conditions you set.
.PP
If the time-zero solution looks ok, but the job failed further down the
track, then examine the other solutions that were saved, particularly
the mesh, to determine if it has been excessively deformed.  Note that
the last solution to be obtained is always present in the debug file.
.PP
If the time-zero solution is not available, check whether the debug
file contains a record, and in particular whether the mesh as recorded
in this file has been properly set up.  Failure to converge at the
first timestep has occurred if the rheology is extremely
non-linear (SE > 3) and the initial guess at the solution is poor.
The standard strategy in such cases is to attempt to obtain a converged
time-zero solution for a lesser value of SE, e.g. SE = 1.0, and then
if necessary, use this solution as an initial condition for a calculation
with a higher SE value.
.PP
If the job didn't get this far, it is likely that there was a syntax
error in the input file, or conceivably a parameter set to a bad
value.  Check the Standard Out file and the BIN.out file for messages
indicating the source of such problems.
.PP
.SH b) excessively large velocities
.PP
If the velocities calculated at any time in the calculation become
excessively large, a number of problems may appear, including loss
of ascii output data as the format of output variables is overflowed.
If this occurs, you can scale back the flow rates in the solution by
setting a larger viscosity coefficient throughout the solution region
(using the REG statement).
.PP
.SH c) excessively small velocities
.PP
If you define a problem in which the driving forces are zero, or
nearly zero, the program will produce a solution with velocity
field everywhere, almost equal to zero.  It will then rapidly
integrate forward in time without changing the solution, but
perhaps producing lots of output.  Check the specification of
the problem and that boundary conditions and driving forces are
set properly.
