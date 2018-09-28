#!/bin/sh -f

export INPUT_FILE_DIR=/Volumes/NICER_BackupHDD_20180725/heasarc/data/nicer/data/obs/2017_08/1013010104/
echo "INPUT_FILE_DIR="$INPUT_FILE_DIR

nibarytime.py \
	$INPUT_FILE_DIR/xti/event_cl/ni1013010104_0mpu7_cl.evt.gz \
	83.633218 22.014464 $INPUT_FILE_DIR/auxil/ni1013010104.orb.gz \
	--refframe ICRS --ephem JPLEPH.430 \
	--outfits ni1013010104_0mpu7_cl_bary_DE430.evt
