#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pylab as plt 

df = pd.read_csv('out/at2019wey_v200901/fit/xspec_multifit_fit.csv')
mjd = 0.5 * (df['MJD_DATEOBS']+df['MJD_DATEEND'])
mjd_err = 0.5 * (-df['MJD_DATEOBS']+df['MJD_DATEEND'])

# par_2 : kT
# par_3 : disk_bb 
# par_4 : Photon index
# par_7 : norm 

fig, axs = plt.subplots(7,1,figsize=(8.27,11.69),sharex=True, gridspec_kw={'hspace': 0})

axs[0].plot(mjd,df['rate_0p4to1p0keV'],'o',mec='black',label='0.4-1.0 keV')
axs[0].plot(mjd,df['rate_1p0to10p0keV'],'o',mec='black',label='1-10 keV')
axs[0].set_ylabel('Rate\n(counts/sec)')
axs[0].set_xlim(59068,59092)
axs[0].set_ylim(0,400)
axs[0].legend(loc='upper left')

axs[1].plot(mjd,df['rate_1p0to10p0keV']/df['rate_0p4to1p0keV'],'o',mec='black',mfc='gray')
axs[1].set_ylabel('Hardness\n(1-10/0.4-1.0)')
axs[1].set_ylim(1.8,3.8)

axs[2].plot(mjd,df['flux_0p4to1p0keV']/1e-10,'o',mec='black',label='0.4-1.0 keV')
axs[2].plot(mjd,df['flux_1p0to10p0keV']/1e-10,'o',mec='black',label='1-10 keV')
axs[2].set_ylabel('Absorbed flux\n($10^{-10}$ erg s$^{-1}$ cm$^{-2}$)')
axs[2].set_yscale('log')
axs[2].set_ylim(0.1,80.)
axs[2].legend(loc='upper left')

axs[3].plot(mjd,df['par_2'],'o',mec='black')
axs[3].set_ylabel('Disk blackbody\nkT (keV)')
axs[3].set_ylim(0.15,0.18)

axs[4].plot(mjd,np.sqrt(df['par_3']),'o',mec='black')
axs[4].set_ylabel('Inner disk radius\n'r'R$_{in}d_{10 kpc}/\sqrt{\cos {\theta}}$ (km)')
axs[4].set_ylim(380,900)

axs[5].plot(mjd,df['par_4'],'o',mec='black',mfc='orange')
axs[5].set_ylabel('Power-law \nPhoton index')

axs[6].plot(mjd,df['par_7']/100.0,'o',mec='black',mfc='orange')
axs[6].set_ylabel('Unabsorb PL 2-10 keV flux \n($10^{-10}$ erg s$^{-1}$ cm$^{-2}$) ')
axs[-1].set_xlabel('MJD')
for ax in axs:
	ax.label_outer()	
	#ax.minorticks_on()
	ax.xaxis.grid(True)
	ax.xaxis.grid(which='major', linestyle='--', color='#000000')
	ax.xaxis.grid(which='minor', linestyle='-.')	
	ax.tick_params(axis="both", which='major', direction='in', length=5)
	ax.tick_params(axis="both", which='minor', direction='in', length=3)			
outpdf = 'test.pdf'
fig.align_ylabels(axs)
plt.tight_layout(pad=2)
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"		
plt.savefig(outpdf)
