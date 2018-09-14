#!/bin/sh -f

mkdir -p out 

photonphase --absphase --barytime --ephem DE200 \
	--orb data/ni1013010104.orb  \
	--outfile out/ni1013010104_JBOephem_photonphase.evt \
	--plotfile out/ni1013010104_JBOephem_photonphase.png \
	data/ni1013010104_0mpu7_cl_crab_30sec.evt \
	data/crab_JBO_ephemeris_MJD57980.par
