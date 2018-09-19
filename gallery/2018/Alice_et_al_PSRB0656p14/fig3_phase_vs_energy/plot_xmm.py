#!/usr/bin/env python

PI_TO_KEV = 0.001
EVENT_FILE = '/Users/enoto/Dropbox/enoto/research/nicer/analysis/psrb0656/data/180621_kobayashi/phase_ana/pn_sr_bacor_phaseadded.evt'
ebin = 60
pbin = 30
emin = 0.0
emax = 1.4 

import os 
import sys 
import numpy as np 
import astropy.io.fits as pyfits

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rc
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['font.size'] = 18
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['axes.xmargin'] = '0' #'.05'
plt.rcParams['axes.ymargin'] = '.2'

cmd  = 'rm -f psrb0656p14_xmmpn_hist2d.pdf psrb0656p14_xmmpn_hist2d_pulsed.pdf;'
cmd += 'rm -f psrb0656p14_xmmpn_hist2d.eps psrb0656p14_xmmpn_hist2d_pulsed.eps'
print(cmd);os.system(cmd)

hdu = pyfits.open(EVENT_FILE)
nevents = len(hdu['EVENTS'].data)
phase  = hdu[1].data['PULSE_PHASE']
energy = (hdu[1].data['PI'] + np.random.rand(nevents)) * PI_TO_KEV
x = np.append(phase, phase+1.0)
y = np.append(energy, energy)

# https://matplotlib.org/api/_as_gen/matplotlib.pyplot.hist2d.html
fig = plt.figure(figsize=(10,8))
#cm = generate_cmap(['#00008B', '#aaaaab', '#FFFFFF', '#F4D793', '#F4A460'])
hist2d, xedges, yedges, image = plt.hist2d(x,y,
	#norm=mpl.colors.PowerNorm(gamma=1./2.),
	#norm=mpl.colors.LogNorm(),
	cmin=1.0,#vmin=4000,
	#cmap = cm,
	cmap='inferno',
	bins=[pbin,ebin],range=[[0.0,2.0],[emin,emax]])
plt.xlabel('Pulse Phase')
plt.ylabel('Energy (keV)')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Photon Counts')
fig.savefig('psrb0656p14_xmmpn_hist2d.pdf')
fig.savefig('psrb0656p14_xmmpn_hist2d.eps')

hist2d_pulsed_list = []
for profile in hist2d.T:
	if np.isnan(profile[0]):
		hist2d_pulsed_list.append(np.array([np.nan for i in range(len(profile))]))
		continue
	profile_max = np.max(profile)
	profile_min = np.min(profile)
	profile_mean = np.mean(profile)	
	#print(profile_max,profile_mean,profile_min)
	#print(profile)
	pulsed_profile = (profile - profile_mean)/profile_mean*100.0
	hist2d_pulsed_list.append(pulsed_profile)
	#hist2d_pulsed = np.append([hist2d_pulsed],[pulsed_profile])
print(hist2d)	
hist2d_pulsed = np.array(hist2d_pulsed_list)

plt.gcf().clear()
# https://matplotlib.org/api/_as_gen/matplotlib.pyplot.imshow.html
#fig = plt.figure(figsize=(14,0))
#ax  = fig.add_subplot(111)
img = plt.imshow(hist2d_pulsed,interpolation="nearest",aspect='auto',
	vmin=-35,vmax=35,
	origin="lower",cmap='inferno',extent=[0.0,2.0,emin,emax])
plt.xlabel('Pulse Phase')
plt.ylabel('Energy (keV)')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Pulsed Fluctuation (%)')
fig.savefig('psrb0656p14_xmmpn_hist2d_pulsed.pdf',dpi=200)
fig.savefig('psrb0656p14_xmmpn_hist2d_pulsed.eps',dpi=200)




