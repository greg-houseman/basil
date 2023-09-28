# basil
2D viscous flow code; Lagrangian frame advection:
This package is used to solve general 2D viscous flow problems in a finite domain, allowing for spatially varying and/or non-Newtonian viscosity.
The program uses the finite element method and depends on the triangle package of Shewchuk to provide solutions for complex closed domains
with flexible boundary conditions on velocity or traction.
The program can provide solutions for plane-strain, thin viscous sheet (plane-stress), thin viscous shell, or axiymmetric cases.
The program can be installed on a linux system from the zip package on the GitHub website as follows:
"unzip basil-main.zip ; cd basil-main ; ./install.sh ;"
If successful a set of executables will be left in basil-main/bin.
The linux compilation  uses xmkmf and the sybil program GUI uses OpenMotif which both need to be installed on the system.
An alternative compilation script ./InstallS does not require xmkmf.
InstallS has been used on MacOS systems but successful compilation using InstallS
usually requires editing of the files basilsrc/MakeSimple and sybilsrc/MakeSimple in order to specify correct location of libraries 
and use of sybil on a MacOS system requires XCode and XQuartz.
Basil is a complex research tool for which documentation is under construction.  Effective use may require advice from the developers.
