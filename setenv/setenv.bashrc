#!/bin/bash 

echo "#################"
echo "#     HOPPY     #"
echo "#################"

NAME=$(scutil --get ComputerName)
if [ $NAME = 'vicuna' ]; then
	echo '...setting for machine of "vicuna"'
	export HOPPY_PATH="/Users/enoto/work/soft/hoppy"	
	export HEASARC_REPOSITORY="/Users/enoto/work/repository/heasarc"		
elif [ $NAME = 'nebula' ]; then
	echo '...setting for machine of "nebula"'
	export HOPPY_PATH="/Users/enoto/work/drbv1/soft/git/hoppy"	
	export HEASARC_REPOSITORY="/Users/enoto/work/drbv1/reporitory/heasarc"	
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

export HEASARC_REMOTE_FTPDIR="https://heasarc.gsfc.nasa.gov/FTP"

echo "#################"
echo "#     NICER     #"
echo "#################"
export NICER_SOFT_PATH=$HOPPY_PATH/hoppy/nicer
export NICER_BGD_MODEL_VERSION=BGMod_3C50
export NICER_BGD_LIBRARY_PATH=/Users/enoto/Dropbox/enoto/research/nicer/data/mitbgd/$NICER_BGD_MODEL_VERSION
export NICER_RESP_PATH=/Users/enoto/work/niresp
export NICER_DATA_REPOSITORY=$HEASARC_REPOSITORY/data/nicer/data/obs
export NICER_DEFAULT_FPARAM_FILE=$NICER_SOFT_PATH/nipipeline/pipeline_setup_default.yaml
export NICER_SOURCE_LIST_FOR_DOWNLOAD=$NICER_SOFT_PATH/download_nicer_source.lst

export PATH=$NICER_SOFT_PATH:$PATH

echo NICER_SOFT_PATH=$NICER_SOFT_PATH
echo NICER_BGD_MODEL_VERSION=$NICER_BGD_MODEL_VERSION
echo NICER_BGD_LIBRARY_PATH=$NICER_BGD_LIBRARY_PATH
echo NICER_RESP_PATH=$NICER_RESP_PATH
echo NICER_DATA_REPOSITORY=$NICER_DATA_REPOSITORY
echo NICER_DEFAULT_FPARAM_FILE=$NICER_DEFAULT_FPARAM_FILE

echo "#################"
echo "#     RXTE      #"
echo "#################"
export RXTE_SOFT_PATH=$HOPPY_PATH/hoppy/rxte
export PATH=$RXTE_SOFT_PATH:$PATH

echo "#################"
echo "#     MAXI      #"
echo "#################"
export MAXI_SOFT_PATH=$HOPPY_PATH/hoppy/maxi
export PATH=$MAXI_SOFT_PATH:$PATH

echo "##################"
echo "#     TIMING     #"
echo "##################"
export TIMING_SOFT_PATH=$HOPPY_PATH/hoppy/timing
export PATH=$TIMING_SOFT_PATH:$PATH
echo TIMING_SOFT_PATH=$TIMING_SOFT_PATH
