#!/bin/sh -f

#fplot_hist2d_pulseprofile.py \
#	data/ni1020130117_0mpu7_cl_o1u20cut_photonphase_sample.evt \
#	ni1020130117_0mpu7_cl_o1u20cut_photonphase_sample.pdf \
#	--emin 0.0 --emax 1.0  --ebin 20 --pbin 40

fplot_hist2d_pulseprofile.py \
	/Users/enoto/work/nicer/project/psrb0656/180918_phootnphase/nipulse/psrb0656p14_pulsephase.evt \
	ni1020130117_0mpu7_cl_o1u20cut_energy_vs_phsae.fits \
	--emin 0.0 --emax 1.4 --ebin 60 --pbin 30

	