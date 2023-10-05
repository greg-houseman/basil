.TH sybil l
.SH NAME
sybil \- Graphical display of 2-D finite element solutions 
.SH SYNTAX
sybil [-i file.log]
.br
sybilps [-i file.log] [-o file.ps]
.SH DESCRIPTION
.PP
The programs sybil and sybilps are used for graphical output from
solution files written by the 2-D finite element program basil (refer man pages on
basil).  Basil solution files contain multiple records, each of 
which defines the current geometry, velocity field, strain-rate 
and stress distribution.  The sequential records in a file describe 
a solution at different times during the course of a finite deformation
calculation.
.PP
Sybil uses a menu-driven graphical user interface using an X-window.
Sybil is designed for interactive use, but also has the capacity to 
execute in background using a set of instructions previously saved 
(the so-called log file).  Sybilps uses the log file to construct 
a postscript image of the drawing for printing.
.PP
This document provides specific instructions and examples of how to use
sybil.  There is also a help window which can be opened from the top
right corner of the sybil window in order to obtain reference information.
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
of finite deformation problems solved using basil and sybil is 
under development.
.SH PROGRAM SUMMARY
.PP
.IP \(bu 5
.BR SYBIL
reads the solution files, and enables the user to interactively
plot selected quantities in the form of contour plots of strain, stress
or displacement quantities, arrow plots of vectors or principal
stresses, mesh plots or other derived quantities.  The program can also
be directed by a pre-saved set of instructions referred to
as the log file.  Output is directed to the screen in an X-window.
.IP \(bu 5
.BR SYBILPS
produces postscript output from the log file created by
SYBIL.  The postscript files can be viewed on the screen using
ghostview or printed using lpr.

.SH CAVEAT 
.PP
The purpose of this document is to provide guidance and reference
material to new users.  These programs are under current development,
however, so the specifications set out below apply to the current state
of the programs and probably will change with future development.

.SH RUNNING SYBIL
.PP
Sybil can be run from any directory, but on initial opening
it first looks for solutions in the sub-directory FD.sols
of the current directory.  FD.sols was specified by the BINDIR
variable in the basil input file, as the location for basil
solution files.  Sybil uses an input file in the current directory
called sybil.in, in which a number of format variables are set to 
default values as specified further below.  If sybil.in does not exist, 
sybil and sybilps will run with a set of global defaults.
.PP
On opening the sybil window, the following menu headers appear:
.PP
File, Record, XYPlot, Profile, Arrow, Contour, Locate, Options,
and, in the top right hand corner: Help.
.PP
Each of these headers describes a class of actions which are 
initiated by the user selecting an appropriate menu option.
Refer to SYBIL MENU SPECIFICATION below for details of these actions.
Unless otherwise specified, all menu selections require only
the left mouse button.
.PP
Most sybil drawings will contain multiple plots, either e.g., of
different quantities at the same time-level, or perhaps of the same
quantity at different time levels.  The user constructs the drawing
to use a set of rectangular cells on an A4 or USLetter page.  Each cell
may hold one or more separate drawings.  The two parameters called 'rows'
and 'columns' in the sybil.in file determine the number of cells to be used
in the drawing.  Use a number between 1 and 3 (or greater if required)
for each.  Also important here is the parameter called 'orientation'
which may take the value of PORTRAIT or LANDSCAPE.  These variables
should be set to the desired values before execution of sybil.
Refer to: INPUT FILE SPECIFICATION below for complete set of input variables.
.PP
Sybilps is very similar to sybil, but it is linked to a different
set of graphics libraries in order to produce postscript output
instead of screen output.  Sybilps does not run interactively; it
executes a set of instructions contained in a log file previously
constructed in an interactive sybil session.  Refer to: LOG FILE
SPECIFICATION below, for log file structure and content.
.PP
In using sybil, the axes convention is that x is horizontal and y
is vertical as plotted on the screen, regardless of the orientation 
of the page, or whether the solution has been rotated on input.

.SH EXAMPLE SYBIL SESSION
.PP
1. After Sybil is started up, the first action must be to open
a solution file.  From the 'File' menu, choose: 'Open -> Solution'.
Find a solution file using the file chooser that is then produced.  If
sybil was initiated from the project directory, a list of available
solutions will appear in the right hand column.  If not, type the
required directory in the Filter line at the top to find the right
directory, and then open a solution by highlighting it and hitting OK
(or double clicking on the solution).
.PP
2. The solution is automatically opened at the first record, usually
the time-zero solution and scale is automatically set so that the entire
mesh will appear within the plotting region.  It is now possible to
make a plot using most of the other options.  For example, from the 'XYPlot'
menu, choose 'Mesh -> elements', and a line drawing of the finite
element mesh is produced in the top-left cell.
.PP
3.  You can superpose multiple plots in the same cell, but the next step
usually is to direct plotting to another cell.  From the 'Locate'
menu, choose 'Next', 'Previous', or 'Select' to move the active
plotting area to another cell.  The default sequence is to fill
across the page in rows, but you can override this at any time.
Note also that while plotting is centred on a cell, drawing outside
the cell boundaries is not prevented.  Scale is preserved from
plot to plot, only changing the origin from one cell to another.
Scale can however be reset if required, as described below
under 'Options -> Rescale'.
.PP
4. Suppose you now want a velocity plot at the same time-level.
The next step is to choose 'Velocity' from the 'Arrow' menu, and
a vector plot of velocity arrows will now appear in the current
cell.
.PP
5.  If a previous plot in the current drawing cell was a mistake, 
you can erase it by moving the mouse to that cell and pressing the
middle button or by choosing 'Options->Delete'.
.PP
6.  Suppose you now want to plot contours of the x-component of
velocity. From the 'Contour' menu, choose 'Velocity -> Ux'
For contour plots, a secondary menu appears at this stage to confirm
selection of colour scale and / or contour levels.  Sybil uses a
standard colour scale with red indicating high values and blue
indicating low values.  Superimposed on the colour shading is a
set of contours if required.  Usually the colour scale is stretched
to fill the available range of values, but you can override the
this colour mapping using the first 2 slots to define the values
of the Min and Max colour.  Contours are drawn at constant increments
of 'Step', above and below 'Level'.
.PP
7. If you have now completed the required plots for
the first solution record, and wish to make some plots for the
later solution records,  choose 'Next' 'Previous', 'Last', 
or 'Specify' from the 'Record' menu, to select any of the available
records in the current solution.  An info line at the bottom of the
window shows the total number of records in this file (under max)
and the currently selected solution (under curr).  The model 'time'
parameter is also shown in this line. Any subsequent
plots will be for the newly selected solution record.  You
can also open a different solution file and choose any record
from that file in the course of constructing a sybil drawing.
.PP
8. Automatic labelling in sybil is simple and minimal and it
is often switched off (under 'Options->Plot->Label') in order 
to exercise precise control over labels.  
To add your own labels to the drawing, choose 'Label->Edit' from
the 'Options' menu.  Enter a string in the box that appears,
then use the left mouse button to position the lower-left corner
of that string on the drawing.  If the location is not quite right
pick a new location (but unwanted labels must be edited out of
the log file - see below). Choose 'Label->Font' to select from
Helvetica or Symbol fonts at different sizes.
.PP
9. Having constructed a drawing interactively, perhaps including
other kinds of plots not discussed in this example, we now save the
set of instructions that enables this plot to be easily recreated.
Choose 'Save->Log' from the 'File' menu in order to obtain a
directory / file specification window.  We suggest that the file
be saved with the suffix .log, so that the set of log files is
easily identified.
.PP
10.  At this point execution may finish by choosing 'Exit' from 
the 'File' menu.
.PP
11.  The log file may be edited in order to change any of the
specification of the drawings.  Take care with such changes if
you are not familiar with the structure of the log file. Some
types of edit are often necessary and easy to do: e.g., deleting
unwanted labels - just locate the string and delete that instruction
line.  Modifying a label content, or its address on the page is
also straightforward.  All labels are grouped at the end of the 
drawing in the order in which they were added to the drawing.
Another simple change is to change the solution filename.  In this
way a single log file can be used as a template to construct
similar log files for a set of drawings using different files.
Detailed structure of the log file is specified below.
.PP
11.  A drawing may be easily reconstructed from its logfile in 
2 ways.  The command 'sybil -i log.file' will recreate the
screen drawing without user intervention.  Alternatively, first
start sybil, then use 'Open->Log' from the 'File' menu to choose
a log file from the directory.  The chosen log file will be
automatically executed.  Further additions to the drawing can
be made interactively at this point, and a new log file then
saved.
.PP
12.  To obtain a postscript version of the drawing, use the
command 'sybilps -i log.file', or, from within sybil, 
choose 'Save->Log+Postscript' from the 'File' menu.  The latter action
saves the log file and then calls sybilps.  Note that there are minor
differences in the placement and size of labels with sybil and
sybilps, owing to differences in the available fonts.

.SH SYBIL MENU SPECIFICATION
.PP
.SH File
Open: opens a file, either 
.br
.in 6P
Solution: to open a solution file, or
.br
Log: to open and execute a previously saved log file
.br 
.in 3P
Save: writes a file, either
.br
.in 6P
Log: to write a log file that recreates the current drawing, or
.br
Log+Postscript: to write the log file, and then cause sybilps to 
execute and write a postscript file of the current drawing.
.br
.in 3P
Exit: terminates execution of sybil, after checking.
.SH Record
Current: selects the solution record for subsequent drawing,
.br
.in 6P
Next: selects the next record in the current file
.br
Previous: selects the previous record in the current file
.br
Last: selects the last record in the current file
.br
Specify: type input number of solution record (1 is the
first record and the value of max in the info line at the bottom
of the window is the last record)
.br
.in 3P
Reference: selects the solution record to define an initial geometry
relative to which subsequent deformations are measured,
.br
.in 6P
First: uses the first record of the current file
.br
Specify: uses the record number you identify, as above.
.br
.SH XYPlot
.PP
Mesh: creates an X-Y plot of the finite element mesh, either
.br
.in 6P
elements: outline of all elements in the mesh
.br
element+num: outline as above plus element numbers of some elements
.br
boundary: external mesh boundary only.
.br
viscosity: shows regions with anomalous viscosity
.br
.in 3P
Deform: draws deformation indicators relative to a previously saved
reference record (see Record->Reference above),
.br
.in 6P
ellp: draws finite deformation ellipses for a subset of triangles
.br
tria: draws the subset of triangles used for the ellipse calculations
.br
trel: draws ellipses and triangles, as specified by two preceding options.
.br
.in 3P
Bounding box: draws a rectangular box around the current plotting area,
excluding margins.
.br
LGMesh: draws the Lagrangian mesh that may be used by basil to track deformation,
.br
.in 6P
elements: draws the entire Lagrangian mesh
.br
boundary: draws the boundary only
.br
.in 3P
StrainMark: draws any strainmarkers defined by basil to track internal deformation
.SH Profile
.PP
1_D: constructs a 1-D profile of selected parameters between 2 points on the x-y plane.
After selecting a quantity to profile from the lists below, input of the x-y
coordinates of the two profile endpoints is requested.  Min and Max values on the
profile are then shown and, the user must input values for the Min and Max on the
vertical axis of the profile plot.  Default values are the min and max values
of the function on the first call, and scale is preserved for subsequent profile
plots.
.br
.in 6P
Velocity: select from the following parameters:
.br
.in 9P
Ux: x-component of velocity
.br
Uy: y component of velocity
.br
Ur: radial component of velocity (relative to local origin)
.br
Uth: tangential component of velocity
.br
.in 6P
Strain: select from the following quantities related to strain-rate
.br
.in 9P
edxx: the xx-component of strain-rate
.br
edyy: the yy-component of strain-rate
.br
edzz: the zz-component of strain-rate
.br
edxy: the xy-component of strain-rate
.br
psr1: maximum principal strain-rate in x-y plane
.br
psr2: minimum principal strain-rate in x-y plane
.br
msst: maximum shear strain-rate in the x-y plane
.br
cang: orientation of principal compressional axis
.br
tang: orientation of principal extensional axis
.br
sang: orientation of direction of maximum shear
.br
dblc: parameter that indicates type of faulting
.br
vort: vertical component of vorticity
.br
ed2i: 2nd invariant of the strain-rate tensor
.br
vota: ratio of vorticity to shear strain-rate
.br
.in 6P
Stress: select from the following quantities related to stress
.br
.in 9P
taud: direction of maximum deviatoric stress
.br
taum: maximum deviatoric shear stress
.br
taxx: the xx component of deviatoric stress
.br
tayy: the yy-component of deviatoric stress
.br
tazz: the zz-component of deviatoric stress
.br
taxy: the xy-component of deviatoric stress
.br
tau1: maximum principal deviatoric stress
.br
tau2: minimum principal deviatoric stress
.br
sixx: the xx-component of total stress
.br
siyy: the yy-component of total stress
.br
sizz: the zz-component of total stress
.br
sig1: maximum principal total stress
.br
sig2: minimum principal total stress
.br
thdi: thermal dissipation function 
.br
	(stress*strain-rate)
.br
pres: the pressure field
.br
brit: parameter indicating type of faulting
.br
bri2: parameter indicating type of faulting
.br
visc: effective viscosity parameter
.br
.in 6P
Thickness: layer thickness in the z-direction
.br
Density: density parameter
.br
.in 3P
2_D: constructs a horizontal or vertical profile of selected 
quantities integrated in the orthogonal direction.  The same
set of parameters can be profiled, as listed above for 1_D.
After choosing the quantity to be profiled, choose whether
an x-direction, or a y-direction profile is required, and
input the upper and lower limits of the orthogonal integration
variable.  E.g.  To show how the integral of Ux.dy varies
as a function of x, choose the x-direction profile and the
two limits are the extreme values of y for the integration.
Avoid integrating outside the solution region if meaningful
results are required.  Input of vertical axis limits is then
requested, as for 1_D profile.
.br
Mark: is used to draw a line in the current cell showing the
path of the 1-D profile in the x-y plane.  
The profile path can be plotted over a previously drawn contour
plot of the variable that is to be plotted, in order to
show where the profile has come from, and to check correlation
with features in the 2D plane.  After entering the coordinates
of the 2 endpoints, use preview to check position of the profile
which can then be amended as necessary.  The Mark command can
be used also for drawing any additional lines required on a
sybil drawing, regardless of whether a profile plot is required.

.SH Arrow
Velocity:  constructs an arrow plot, consisting of equally spaced
arrow heads, whose direction is parallel to the velocity field and
whose length is proportional to the magnitude of the velocity.
The spacing of arrows can be set under `Options->Plot->mp and np'
(refer 'Options' below).
.br
.in 6P
Strain: constructs an arrow plot, consisting of equally spaced
crosses, whose orientation indicates the directions of the
principal axes of strain-rate.  For direction-only plots all
crosses are the same size, for magnitude plots, the length of
each arm is proportional to the magnitude of the principal axis.
Options are:
.br
.in 9P
pstm: Principal strain-rate magnitudes
.br
pstd: Principal strain-rate directions only
.br
mssr: planes and magnitudes of maximum shear strain-rate
.br
msft: direction of strike-slip faulting planes
.br
.in 6P
Stress: constructs an arrow plot as for Strain (see above),
but with lengths of the arrows either constant or proportional 
to the relevant stress magnitudes.  Options are:
.in 9P
taum: directions and magnitudes of maximum shear stress
.br
taud: directions of maximum shear stress
.br
sigd: principal deviatoric stress, directions only
.br
sigm: principal deviatoric stress, magnitudes and directions.

.SH Contour
.PP
The Contour command is used to display 2-D quantities using
a graphical display in which areas are coloured according to the value
of the physical parameter being displayed, or equally spaced contours
of this physical parameter are drawn, or both.  The same set of
physical parameters listed above under 'Profile->1_D' is
available for contouring.  After choosing the parameter to
be contoured, the user is asked to confirm or enter a set of
parameters that affect the display.  Shown at the top left
are the actual minimum and maximum values of the quantity to be
contoured.  The two values at the top right, entered by the user,
determine the mapping of colour to number.  With the default
colour scale, blue is mapped to the minimum, and magenta to the 
maximum.  Contours, if present, are chosen in equal increments
of 'Step' above and below 'Level'.  A maximum number of contours
is also available.  The number labelled scale is only used to
multiply the labels that appear on the colour scale.  The settings
are preserved for successive calls to 'Contour' so that the
default setting usually is to keep the same scales as for the
preceding plot.  Various default parameters affecting the
display may be set using the 'Options' command (see below).

.SH Locate
.PP
The Locate command is used to specify which of the current
drawing cells is active, as follows:
.br
.in 6P
Next: Move to the next cell, either the cell to the right,
or if already at the right limit of the page, the left most
cell of the next row.
.br
Prev: Move to the previous cell.
.br
Specify: Increment or decrement the row or column number
of the current cell address using the buttons provided
in order to move directly to a particular cell.

.SH Options
.PP
Label: is used to add labels to the plot as follows:
.br
.in 6P
Edit: enter a string to be drawn on the plot.  Then click
with the left mouse button at the point on the drawing
where the lower left corner of the string is to be
located.  The same string can be added to the plot as
many times as required and the string can be edited
at any time in the labelling step.
.br
Font: choose one of the 2 fonts (Helvetica, Symbol) and 5
or 6 font sizes (8, 10, 12, 14, 18, 24) that are provided 
as options.
.br
.sp
.in 3P
Colour: choose one of the 8 standard colours to be used
in subsequent drawing and labelling operations.
.br
.sp
Line: choose solid, dashed, or dotted line for subsequent
drawing operations.
.br
.sp
Rescale: On making the first plot, the scale is set automatically
so that the entire solution region will appear within the plotting
cell, including a small margin on each side.  Subsequent plots are
drawn at the same scale, with the origin shifted to the new
drawing cell.  The scale will be automatically reset for the
next plot, only if this option is first selected.  Judicious
choice of cell numbers and margin widths may be necessary to
obtain a particular scale.
.br
.sp
Delete: deletes previous drawings from the screen.  Note that the
first option is also available at any time by locating the cursor 
within the current cell and pressing the middle button of the mouse.
.br
.in 6P
Current cell: only deletes contents of one cell
.br
All Cells: deletes contents of all cells
.br
.sp
.in 3P
Plot: requests input of certain options and parameters that
affect the visual display, as follows:
.in 6P
nx3: determines the number of interpolation points to be
used in contouring and other operations.  A number of order 100
is recommended. Greater numbers will improve the resolution of
the plot, but will require more time to plot.
.br
.sp
mp, np: in constructing arrow plots (see above) an arrow is
situated at every (mp)th interpolation point in the x-direction
and every (np)th interpolation point in the y-direction.  The
size of the maximum arrow is equal to the distance between
arrows.  Numbers of 4 or 8 are recommended.  Note that if nx3
is doubled, mp and np should be doubled to preserve the same
arrow spacing.
.br
.sp
profile_pts: the number of sample points to take along a given profile
.br
.sp
stipple: if > 0, stippling is added to contour plots in specified
regions, with reference to the Level parameter (see Contour above):
.br
stipple = 1 -> stipple if below Level
.br
stipple = 2 -> stipple if above Level
.br
stipple = 3 -> stipple if abs value above Level
.br
stipple = 4 -> stipple if abs value below level
.br
solution_rot: This parameter should be a number between 1 and 3, to
rotate the solution by quanta of 90 degrees anticlockwise.  The
rotation is done immediately that the solution is input, and any
reference to x and y-coordinates then applies to the solution
as displayed in its new reference frame.
.br
.sp
dble: Causes causes most area plots to be displayed twice,
with a horizontal offset equal to the width of the solution
region.  This option is only used for solutions with 
periodicity in the x-direction.
.br
.sp
label: causes all automatic labelling to be switched on or off.
This option does not affect any labels added manually, as above.
.br
.sp
flip: causes the x-axis to be reversed on the display.
.br
.sp
Contour Options: either or both of the following may be selected.
.br
.in 6P
Lines: contour lines are added to all contour plots
.br
Shading: shading by colour is added to all contour plots
.br
Bar: the colour bar for contour plots may be vertical, 
to the right of the plot, or horizontal, underneath the 
plot, or may be omitted.
.br
Tics: ticks on profile plots may be internal to the plot,
external, centered, or omitted.
.br
.sp
.in 3P
Verbose: outputs numbers (lots of them) to standard out, from 
profiles, contour plots, arrow plots and strain-marker plots, 
in order that other external processing may be applied to these data.
This option should only be turned on temporarily to extract a
particular data set and is not yet implemented uniformly.
.br
.sp
Mark Cell: shows the current plotting cell by a line drawn around
its perimeter.  This boundary line does not appear on the postscript
plot.

.SH Help
.PP
The Help command displays reference information about menu
options, parameter usage and structure of the Log file.  Information
appears in a separate window.
.br
.SH INPUT FILE SPECIFICATION
.PP
The input file, sybil.in, is normally read upon commencing execution and
is used to set key parameters for the drawing. In the case where
a log file is used, sybil.in is ignored and these parameters are
set by the log file.
.PP
In either case, the file commences with the word 'Options' followed by
a set of optional parameter assignments as listed below. The
file may also contain optional comment lines commencing with '#'.
These parameter values are read when sybil is opened, but the 
file is only found if sybil is opened in the directory containing
the file.  A subset of the parameter values may be reset during
execution of sybil, using the 'Options->Plot' command, as
described above.
.PP
To construct an up-to-date sybil.in file, open sybil and, before making
any plots, write a log file using 'File->Save->Log', saving it as
sybil.in. This file includes a set of default parameter values, together
with a short description and a list of possible values for each parameter.
.br
.SH LOG FILE SPECIFICATION
.PP
To be completed.
