source /cvmfs/sft.cern.ch/lcg/views/LCG_106b/x86_64-el9-gcc11-opt/setup.sh
 export CXX=`which g++`
 export CC=`which gcc`
 cmake3 -DGeant4_DIR=/cvmfs/geant4.cern.ch/geant4/11.3/x86_64-el9-gcc11-optdeb-MT/lib64/Geant4-11.3.0/ ../
 make
