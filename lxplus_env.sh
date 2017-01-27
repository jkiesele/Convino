
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "${SCRIPT}")

. /afs/cern.ch/sw/lcg/external/gcc/5.3/x86_64-slc6-gcc53-opt/setup.sh
cd /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.06/x86_64-slc6-gcc49-opt/root/
. bin/thisroot.sh
cd -


export LD_LIBRARY_PATH=$BASEDIR:$BASEDIR/armadillo/usr/local/lib64:$LD_LIBRARY_PATH
export PATH=$BASEDIR:$PATH
