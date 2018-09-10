#!/bin/sh -f

rm -f ni1013010104_0mpu7_cl_500kcnt_JBO_phase.png
rm -f ni1013010104_0mpu7_cl_500kcnt_JBO_phase.evt

# --addphase 
# --barytime 
# --absphase 

photonphase --ephem DE200 --orb data/ni1013010104.orb --absphase \
	--outfile ni1013010104_0mpu7_cl_500kcnt_JBO_phase.evt \
	--plot --plotfile ni1013010104_0mpu7_cl_500kcnt_JBO_phase.png \
	data/ni1013010104_0mpu7_cl_500kcnt.evt crab_JBO_ephemeris_MJD57980.par
