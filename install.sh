#!/bin/sh
#
# provide version number to source directories
#
read vname < version
echo 'char version[8] = "' $vname '" ;' > basilsrc/version.h
echo 'char version[8] = "' $vname '" ;' > sybilsrc/version.h
#
# executables will be moved into bin directory
#
if ! [ -e bin ]
then
  mkdir bin
fi
#
# build the makefiles
#
xmkmf
if [ $? != 0 ]
then
  exit 1
fi
make Makefiles
if [ $? != 0 ]
then
  exit 1
fi

#
# build and install the binaries
#
make all
if [ $? != 0 ]
then
  exit 1
fi
make install
if [ $? != 0 ]
then
  exit 1
fi

make depend

make clean

make install.man
