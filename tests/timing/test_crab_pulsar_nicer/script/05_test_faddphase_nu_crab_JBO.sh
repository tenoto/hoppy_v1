#!/bin/sh -f 

rm -f ni1013010104_JBOephem_nibarytime_phase_mjd.evt
faddphase_nu.py \
	out/ni1013010104_JBOephem_nibarytime.evt \
	57980.000000184372766 \
	29.639410077200001 \
	--nudot=-3.6870211e-10 \
	--nu2dot=9.17300618091751e-21 \
	--outfits out/ni1013010104_JBOephem_nibarytime_phase_mjd.evt \
	--flag_mjd  

rm -f out/ni1013010104_JBOephem_nibarytime_phase_mjd.fht 	
fplot_pulseprofile.py \
	out/ni1013010104_JBOephem_nibarytime_phase_mjd.evt \
	--outfits out/ni1013010104_JBOephem_nibarytime_phase_mjd.fht \
	--nbin 100 \
	--colname PULSE_PHASE \
	--title "Crab pulse profile (MJD)"    

# 57980 --> 114220732.816 
# 57980.000000184372766 --> 114220732.816 + 0.0159298069 
rm -f ni1013010104_JBOephem_nibarytime_phase.evt
faddphase_nu.py \
	out/ni1013010104_JBOephem_nibarytime.evt \
	114220732.831929807 \
	29.639410077200001 \
	--nudot=-3.6870211e-10 \
	--nu2dot=9.17300618091751e-21 \
	--outfits out/ni1013010104_JBOephem_nibarytime_phase.evt 

rm -f out/ni1013010104_JBOephem_nibarytime_phase.fht 	
fplot_pulseprofile.py \
	out/ni1013010104_JBOephem_nibarytime_phase.evt \
	--outfits out/ni1013010104_JBOephem_nibarytime_phase.fht \
	--nbin 100 \
	--colname PULSE_PHASE \
	--title "Crab pulse profile (MET)"    	