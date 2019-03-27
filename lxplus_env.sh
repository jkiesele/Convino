
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "${SCRIPT}")

source /cvmfs/sft.cern.ch/lcg/views/ROOT-latest/x86_64-slc6-gcc49-opt/setup.sh

export PATH=$BASEDIR:$PATH
