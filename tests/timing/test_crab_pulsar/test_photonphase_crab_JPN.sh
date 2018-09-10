#!/bin/sh -f

rm -f ni1013010104_0mpu7_cl_500kcnt_JPN_phase.png
rm -f ni1013010104_0mpu7_cl_500kcnt_JPN_phase.evt

# --addphase 
# --barytime 
# --absphase 

photonphase --ephem DE430 --orb data/ni1013010104.orb --addphase \
	--outfile ni1013010104_0mpu7_cl_500kcnt_JPN_phase.evt \
	--plot --plotfile ni1013010104_0mpu7_cl_500kcnt_JPN_phase.png \
	data/ni1013010104_0mpu7_cl_500kcnt.evt crab_JPNradio_ephemeris_MJD57974.par
