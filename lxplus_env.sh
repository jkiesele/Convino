
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "${SCRIPT}")

source /cvmfs/sft.cern.ch/lcg/views/LCG_92/x86_64-slc6-gcc62-opt/setup.sh

export PATH=$BASEDIR:$PATH
