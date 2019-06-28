#!/bin/bash 

echo "#################"
echo "#     HOPPY     #"
echo "#################"

#export HOPPY_PATH=$(pwd)
export PYTHONPATH=$HOPPY_PATH:$PYTHONPATH

echo "----- FTOOLS (general) -----"
export PATH=$HOPPY_PATH/hoppy/ftools/cli:$PATH

echo "----- NICER -----"
export NICERPAGE_USERNAME="xxxxx"
export NICERPAGE_PASSWORD="xxxxx"
export NICERDATA_DECRYPT_PASSPHRASE="xxxxx" 
source /Users/enoto/Dropbox/01_enoto/library/hoppy/setenv/nicer_password.bashrc

export NICER_SEG_TEAM_WEBSITE="https://heasarc.gsfc.nasa.gov/docs/nicer/team_schedule/nicer_seg_team.html"
export NICER_DOWNLOAD_SOURCE_LIST=$HOPPY_PATH/data/nicer/nicer_download_source.lst

export PATH=$HOPPY_PATH/hoppy/nicer/cli:$PATH


