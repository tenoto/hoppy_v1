#!/usr/bin/env python

import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt 

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = '13'
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['xtick.top'] = 'True'
mpl.rcParams['ytick.right'] = 'True'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['axes.grid'] = 'True'
mpl.rcParams['axes.xmargin'] = '.05' #'.05'
mpl.rcParams['axes.ymargin'] = '.05'
mpl.rcParams['savefig.facecolor'] = 'None'
mpl.rcParams['savefig.edgecolor'] = 'None'
mpl.rcParams['savefig.bbox'] = 'tight'

df = pd.read_csv('data/xspec_multifit_fit_v180912a.csv')

fig, axes = plt.subplots(1,1,figsize=(9.6,3.6,))
axes.errorbar(
	df['MJD_DATEOBS'],df['rate_0p2to2p0keV'],
	yerr=[df['rate_0p2to2p0keV_err'],df['rate_0p2to2p0keV_err']],
	label='NICER (0.2-2.0 keV)',
	fmt='o',
	color='#D35400',
	markersize=2.8,
	markerfacecolor='#D35400',
	markeredgecolor='#D35400')
axes.legend(loc='upper left',shadow=False)
#axes.set_title()
axes.set_xlim(58030,58260)
axes.set_ylim(0.0,15.0)
axes.set_xlabel('MJD')
axes.set_ylabel('Count Rate (cnt s$^{-1}$)')
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('fig_psrb0656p14_nicer_lightcurve.eps')