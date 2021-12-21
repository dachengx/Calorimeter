#!/bin/bash

source /opt/geant4.10.07-install/bin/geant4.sh
source /opt/root-install/bin/thisroot.sh
unset JUPYTER_PATH
unset JUPYTER_CONFIG_DIR

mkdir build
cd build

cmake -DGeant4_DIR=$G4DIR -DCMAKE_BUILD_TYPE=Debug ..

make -j

./exampleB4c -m run.mac > log.log

hadd gamma.root gamma_1.root gamma_2.root gamma_5.root gamma_10.root gamma_50.root gamma_100.root

cp *.root /tmp/.
