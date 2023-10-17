# basil
2D viscous flow code; Lagrangian frame advection:
This package is used to solve general 2D viscous flow problems in a finite domain, allowing for spatially varying and/or non-Newtonian viscosity.
The program uses the finite element method and depends on the triangle package of Shewchuk to provide solutions for complex closed domains
with flexible boundary conditions on velocity or traction.
The program can provide solutions for plane-strain, thin viscous sheet (plane-stress), thin viscous shell, or axiymmetric cases.
The program can be installed on a linux system from the zip package on the GitHub website as follows:\
unzip basil-main.zip ; cd basil-main ; ./install.sh ;\
If successful a set of executables will be left in basil-main/bin.
You can access man pages for basil and sybil by adding an ammendment to the MANPATH environment variable,e.g. 
for bash, edit the .bashrc file to include lines as follows:\
export BASILPATH=/home/yourpath/basil-main\
export PATH=$PATH:$BASILPATH/bin\
export MANPATH=$MANPATH:$BASILPATH\
The linux compilation uses xmkmf (apt-get install xutils-dev) and the sybil program GUI uses OpenMotif (apt-get install libmotif-dev) which both need to be installed on the system.
An alternative compilation script ./InstallS does not require xmkmf and 
has been used on MacOS systems. Successful compilation using InstallS
usually requires editing of the files basilsrc/MakeSimple and sybilsrc/MakeSimple in order to specify correct location of libraries.
Use of sybil on a MacOS system also requires XCode and XQuartz.
Basil is a complex research tool for which documentation is under construction.  Effective use may require advice from the developers.  You can contact us using the "issues" button on the GitHub page.
