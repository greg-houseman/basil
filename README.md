# basil
2D viscous flow code; Lagrangian frame
This package used to solve general 2D viscous flow problems in finite domain allowing for sptaially varying and /or non-Newtonian viscosity.
The program uses the finite element method and depends on the triangle package of Shewchuk to provide solutions for complex closed domains
with flexible boundary conditions on velocity or traction.
The program can provide solutions for plane-strain, thin viscous sheet (plane-stress), thin viscous shell, or axiymmetric cases.
The program can be installed on a linux system from the zip package on the GitHub website as follows:
unzip basil-main.zip
cd basil-main
./install.sh
If successful a set of executables will be left in basil-main/bin.
The ilinux compilation  uses xmkmf which should be installed.  The GUI in the sybil program uses OpenMotif which will need to be installed.
An alternative compilation script used on MacOS is ./InstallS  This does not require xmkmf ,but will probably require editing of the files
basilsrc/MakeSimple and sybilsrc/MakeSimple to specify location of libraries.
Basil is a complex research tool for which documentation is under construction.  Effective use may require advice from the devlopers.
