#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pylab as plt 
from astropy.io import fits

path_to_file = '/Users/enoto/drop/01_enoto/research/nicer/analysis/proximacen/data/200718_Notsu_TESS/add_nicer/'
hdu1 = fits.open(path_to_file + 'tess2019112060037-s0011-0000000388857263-0143-s_lc_nicer.fits')
hdu2 = fits.open(path_to_file + 'tess2019140104343-s0012-0000000388857263-0144-s_lc_nicer.fits')

df1 = pd.DataFrame(hdu1['TRBLIST'].data)
#df1.drop(index='cow', columns='small')
mask_tess1 = (df1['PDCSAP_FLUX'] > 0)
PDCSAP_FLUX_mean1 = np.mean(np.array(df1['PDCSAP_FLUX'][mask_tess1]))
df1['PDCSAP_FLUX_normalized'] = df1['PDCSAP_FLUX'] / PDCSAP_FLUX_mean1

df1_sel = df1.loc[:,['MJD','cps_0.3_1.0keV','exp_1.0_3.0keV','PDCSAP_FLUX','PDCSAP_FLUX_normalized']]

df2 = pd.DataFrame(hdu2['TRBLIST'].data)
mask_tess2 = (df2['PDCSAP_FLUX'] > 0)
PDCSAP_FLUX_mean2 = np.mean(np.array(df2['PDCSAP_FLUX'][mask_tess2]))
df2['PDCSAP_FLUX_normalized'] = df2['PDCSAP_FLUX'] / PDCSAP_FLUX_mean2
df2_sel = df2.loc[:,['MJD','cps_0.3_1.0keV','exp_1.0_3.0keV','PDCSAP_FLUX','PDCSAP_FLUX_normalized']]

df = pd.concat([df1_sel, df2_sel])
mask = (df['exp_1.0_3.0keV'] > 0)

fontsize = 14
fig, axs = plt.subplots(2,1, figsize=(11.69,8.27), 
	sharex=True, gridspec_kw={'hspace': 0})
axs[0].step(df['MJD'][mask],df['cps_0.3_1.0keV'][mask],
	'o',mec='k',markersize=4,where='mid',color='#FF5733')
axs[0].set_ylabel('NICER rate (cps) 0.3-1.0 keV',fontsize=fontsize)
axs[0].set_yscale('log')
axs[1].step(df['MJD'],df['PDCSAP_FLUX_normalized'],
	'o',mec='k',markersize=4,where='mid',color='#1565C0')
#axs[1].set_yscale('log')
axs[1].set_ylabel('TESS optical intensity (arbitary unit)',fontsize=fontsize)
axs[1].set_ylim(0.995,1.08)
axs[1].set_xlabel('MJD (day)',fontsize=fontsize)
for ax in axs:
	ax.label_outer()	
	ax.minorticks_on()
	ax.xaxis.grid(True)
	ax.xaxis.grid(which='major', linestyle='--', color='#000000')
	ax.xaxis.grid(which='minor', linestyle='-.')	
	#ax.xaxis.set_minor_locator(dates.HourLocator())
	ax.tick_params(axis="both", which='major', direction='in', length=5)
	ax.tick_params(axis="both", which='minor', direction='in', length=3)
fig.align_ylabels(axs)
plt.tight_layout(pad=2)
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.savefig('curve.pdf')

for i in range(len(axs)):
	axs[i].set_xlim(58622,58624)
plt.savefig("tmp1.pdf")	

fig, axs = plt.subplots(1,1, figsize=(8.27,8.27))
plt.plot(df['PDCSAP_FLUX_normalized'][mask],df['cps_0.3_1.0keV'][mask],'o')
#	'o',mec='k',markersize=4,where='mid',color='#FF5733')
plt.xlabel('TESS optical intensity (arbitary)',fontsize=fontsize)
#axs[0].set_yscale('log')
plt.ylabel('NICER (cps) 0.3-1.0 keV',fontsize=fontsize)
plt.tight_layout(pad=2)
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.savefig('scatter.pdf')
