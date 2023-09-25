#!/bin/bash
# Script to be executed in the topas source code directory

# Change paths as appropriate
FLAVOR="meshgeom"
BUILDDIR="build_$FLAVOR"
EXTENSIONS_PATH="/Users/isaacmeyer/topas_extensions/TOPAS-MeshGeom/"
GEANT4PATH="/Users/isaacmeyer/geant4.10.07.p03/build/source/externals/ptl/include"
export Geant4_DIR=/Applications/geant4.10.07.p03-install
export GDCM_DIR=/Applications/gdcm-install
# export CMAKE_PREFIX_PATH=/optional/Qt/path # (e.g. /Qt/5.15.2/clang_64/bin)

rm -rf $BUILDDIR
mkdir $BUILDDIR
cd $BUILDDIR

echo $(which cmake)
cmake -DTOPAS_TYPE="public"\
      -DCMAKE_PREFIX_PATH=$GEANT4PATH\
      -DTOPAS_EXTENSIONS_DIR=$EXTENSIONS_PATH\
      -DCMAKE_INSTALL_PREFIX="~/topas/topas-install-$FLAVOR"\
      ../
make -j12
