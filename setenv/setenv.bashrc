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
export NICER_VIS_TEAM_WEBSITE="https://heasarc.gsfc.nasa.gov/docs/nicer/team_schedule/nicer_vis_team.html"
#export NICER_DOWNLOAD_SOURCE_LIST=$HOPPY_PATH/data/nicer/nicer_download_source.lst
export NICER_DOWNLOAD_SOURCE_LIST=$HOPPY_PATH/data/nicer/nicer_download_source_gsheet.lst

export NICER_BKGDIR="/Users/enoto/Dropbox/01_enoto/research/nicer/data/mitbgd/nibackgen3C50_alpha/bg_model_3C50_RGv5"
export NICER_RMF_DEFAULT="/Users/enoto/Dropbox/01_enoto/research/nicer/response/190628_craig_v1/nicer-rmf-upd52-teamonly/nicer_upd_d52.rmf"
export NICER_ARF_DEFAULT="/Users/enoto/Dropbox/01_enoto/research/nicer/response/190628_craig_v1/nicer-arf-consim135-teamonly/nicer-arf-consim135-teamonly-array52.arf"

export PATH=$HOPPY_PATH/hoppy/nicer/cli:$PATH


