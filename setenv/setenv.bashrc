#!/bin/bash 

NAME=$(scutil --get ComputerName)
if [ $NAME = 'vicuna' ]; then
	echo '...setting for machine of "vicuna"'
elif [ $NAME = 'nebula' ]; then
	echo '...setting for machine of "nebula"'
	export HOPPY_PATH="/Users/enoto/work/drbv1/soft/git/hoppy"	
elif [ $NAME = 'llama' ]; then
	echo '...setting for machine of "llama"'
	export HOPPY_PATH="/Users/enoto/work/soft/hoppy"		
else
	echo 'no corresponding computer setup.'
fi

export PATH=$HOPPY_PATH/hoppy/script:$PATH
export PYTHONPATH=$HOPPY_PATH:$PYTHONPA

echo ComputerName    = $NAME
echo HOPPY_PATH = $HOPPY_PATH

alias hoppyinit="source $HOPPY_PATH/setenv/setenv.bashrc"

################
#     NICER    # 
################  
export PATH=$HOPPY_PATH/hoppy/nicer:$PATH
export NICER_BGD_MODEL_VERSION="BGMod_3C50"
export NICER_BGD_LIBRARY_PATH="/Users/enoto/Dropbox/enoto/research/nicer/data/mitbgd/"$NICER_BGD_MODEL_VERSION
export NICER_RESP_PATH="/Users/enoto/work/niresp"

echo NICER_BGD_MODEL_VERSION=$NICER_BGD_MODEL_VERSION
echo NICER_BGD_LIBRARY_PATH=$NICER_BGD_LIBRARY_PATH
echo NICER_RESP_PATH=$NICER_RESP_PATH