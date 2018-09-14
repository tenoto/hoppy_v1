#!/usr/bin/env python

import matplotlib.pyplot as plt
import astropy.io.fits as pyfits 
import numpy as np 

hdu_photonphase = pyfits.open('out/ni1013010104_JBOephem_photonphase.evt')
hdu_barycorrmet = pyfits.open('out/ni1013010104_JBOephem_nibarytime_phase.evt')
hdu_barycorrmjd = pyfits.open('out/ni1013010104_JBOephem_nibarytime_phase_mjd.evt')

print(len(hdu_photonphase[1].data['ABS_PHASE']))
flag_diff = (hdu_photonphase[1].data['ABS_PHASE']-hdu_barycorrmjd[1].data['PULSE_NUMBER']==0)
print(list(flag_diff).count(True))
#print(hdu_photonphase[1].data['TIME'][flag_diff==False])

flag_diff = (hdu_photonphase[1].data['ABS_PHASE']-hdu_barycorrmet[1].data['PULSE_NUMBER']==0)
print(list(flag_diff).count(True))
#print(hdu_photonphase[1].data['TIME'][flag_diff==False])

pulse_phase_diff_mjd = hdu_photonphase[1].data['PULSE_PHASE']-hdu_barycorrmjd[1].data['PULSE_PHASE']
print(pulse_phase_diff_mjd)

fig, ax = plt.subplots(figsize=(5,4))
#ax.bar(center, hist, align='center', width=width)
#ax.set_xticks(bins)
plt.hist(pulse_phase_diff_mjd/1e-5,bins=100,range=(-6,-2))
plt.xlabel('Pulse Phase Difference (1e-5)')
#plt.ylabel('Events')
plt.title('MJD')
fig.savefig("out/ni1013010104_JBOephem_nibarytime_phase_mjd_diff.pdf")


pulse_phase_diff_met = hdu_photonphase[1].data['PULSE_PHASE']-hdu_barycorrmet[1].data['PULSE_PHASE']
print(pulse_phase_diff_met)

fig, ax = plt.subplots(figsize=(5,4))
#ax.bar(center, hist, align='center', width=width)
#ax.set_xticks(bins)
plt.hist(pulse_phase_diff_met/1e-5,bins=100,range=(-6,-2))
plt.xlabel('Pulse Phase Difference (1e-5)')
#plt.ylabel('Events')
plt.title('MET')
fig.savefig("out/ni1013010104_JBOephem_nibarytime_phase_met_diff.pdf")
