#!/bin/sh -f 

rm -f ni1013010104_0mpu7_cl_500kcnt_bary.evt
barycorr infile=data/ni1013010104_0mpu7_cl_500kcnt.evt \
	outfile=ni1013010104_0mpu7_cl_500kcnt_bary.evt \
	ra=83.633218 dec=22.014464 \
	orbitfiles=data/ni1013010104.orb \
	refframe=ICRS ephem=JPLEPH.430 

rm -f ni1013010104_0mpu7_cl_500kcnt_bary_phase.evt 
faddphase_nu.py \
	-i ni1013010104_0mpu7_cl_500kcnt_bary.evt \
	-o ni1013010104_0mpu7_cl_500kcnt_bary_phase.evt \
	-n 29.639601201518 \
	-d -3.687105e-10 \
	-e 113702332.821496

rm -f ni1013010104_0mpu7_cl_500kcnt_bary_pulse.fits 
fplot_pulseprofile.py \
 	-i ni1013010104_0mpu7_cl_500kcnt_bary_phase.evt \
 	-o ni1013010104_0mpu7_cl_500kcnt_bary_pulse.fits \
 	-n 100 -c PHASE	