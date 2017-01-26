
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "${SCRIPT}")

echo "using CMSSW environment to setup right root and gcc version (no real dependence on CMSSW)"

cd /afs/cern.ch/user/j/jkiesele/public/CMSSW_8_0_4
eval `scramv1 runtime -sh`
cd - 
export LD_LIBRARY_PATH=$BASEDIR:$BASEDIR/armadillo/usr/local/lib64:$LD_LIBRARY_PATH
export PATH=$BASEDIR:$PATH
