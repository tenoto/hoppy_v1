#!/bin/sh -f 

nibarytime.py \
	data/ni1013010104_0mpu7_cl_crab_30sec.evt \
	83.633218 22.014464 data/ni1013010104.orb \
	--refframe FK5 --ephem JPLEPH.200 \
	--outfits out/ni1013010104_JBOephem_nibarytime.evt

