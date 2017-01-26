if [[ "${PWD##*/}" != "Convino" ]]; 
then 
echo Please source the script from the Convino base directory

else

BASEDIR=`pwd`

export DYLD_LIBRARY_PATH=$BASEDIR:$DYLD_LIBRARY_PATH
export PATH=$BASEDIR:$PATH

fi
