#!/bin/bash 

NAME=$(scutil --get ComputerName)
if [ $NAME = 'vicuna' ]; then
	echo '...setting for machine of "vicuna"'
elif [ $NAME = 'nebula' ]; then
	echo '...setting for machine of "nebula"'
	export HOPPY_PATH="/Users/enoto/work/drbv1/soft/git/hoppy"	
elif [ $NAME = 'llama' ]; then
	echo '...setting for machine of "fireant"'
else
	echo 'no corresponding computer setup.'
fi

export PATH=$HOPPY_PATH/hoppy/script:$PATH

echo ComputerName    = $NAME
echo HOPPY_PATH = $HOPPY_PATH
export PYTHONPATH=$HOPPY_PATH:$PYTHONPATH