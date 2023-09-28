# basil
2D viscous flow code; Lagrangian frame
This package used to solve general 2D viscous flow problems in finite domain allowing for sptaially varying and /or non-Newtonian viscosity.
The program uses the finite element method and depends on the triangle package of Shewchuk to provide solutions for complex closed domains
with flexible boundary conditions on velocity or traction.
The program can provide solutions for plane-strain, thin viscous sheet (plane-stress), thin viscous shell, or axiymmetric cases.
The program can be installed on a linux system from the zip package on the GitHub website as follows:
"unzip basil-main.zip ; cd basil-main ; ./install.sh ;"
If successful a set of executables will be left in basil-main/bin.
The linux compilation  uses xmkmf and the GUI in the sybil program uses OpenMotif which both need to be installed on the linux system.
An alternative compilation script which has been used on MacOS is ./InstallS  This script does not require xmkmf but OpenMotif
still needs to be installed.  
Successful operation of InstallS will probably require editing of the files
basilsrc/MakeSimple and sybilsrc/MakeSimple in order to specify correct location of libraries on your system.
Basil is a complex research tool for which documentation is under construction.  Effective use may require advice from the devlopers.
