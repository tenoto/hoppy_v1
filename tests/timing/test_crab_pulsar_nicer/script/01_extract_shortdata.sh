#!/bin/sh -f 

# NICER Crab X-ray observation with radio
# Date : 2017-08-09 (MJD 57974, doy 221)
# NICER ObsID : 1013010104 (Exposure 10.1 ks, UTC 17:13-23:53)
# Usuda radio observatory: UTC 16:45-23:59 UTC
# 
# Usuda 2nd GTI from 2017221_UsdSch0123_GTItbl.txt
# 				Start 				End
# Usuda-MJD-UTC? 	57974.750000000000  57974.791655080973 (duraiton 3599.286457)
# NICER-MET 		113767202.000		113770800.999

echo $NICER_DATA_REPOSITORY

rm -rf data/ni1013010104_0mpu7_cl_crab_30sec.evt

fselect $NICER_DATA_REPOSITORY/2017_08/1013010104/xti/event_cl/ni1013010104_0mpu7_cl.evt.gz \
	data/ni1013010104_0mpu7_cl_crab_30sec.evt \
	"(TIME >= 113769700.0) && (TIME < 113769730.0)"

cp  $NICER_DATA_REPOSITORY/2017_08/1013010104/auxil/ni1013010104.orb.gz data

gunzip data/ni1013010104.orb.gz

