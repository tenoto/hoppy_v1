#!/usr/bin/env python

import pandas as pd 
import matplotlib as mpl
import matplotlib.pyplot as plt 
import datetime

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = '14'
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

df_maxi = pd.read_csv('data/J1520-571_00055058g_lc_1day_all_crab.dat')
df_asm  = pd.read_csv('data/rxte_asm_cirx1_1dayave_crab.dat')

fig, axes = plt.subplots(1,1,figsize=(9.6,3.6))
print(df_asm['MJD'])
print(df_asm['crab_intensity'])
axes.errorbar(df_asm['MJD'],df_asm['crab_intensity'],
	yerr=[df_asm['crab_intensity_uncertainty'],df_asm['crab_intensity_uncertainty']],
	label='RXTE/ASM (2-10 keV)',fmt='o',
	color="#35a16b",markersize=2.8,
	markerfacecolor="#35a16b",markeredgecolor="#35a16b")
axes.errorbar(df_maxi['MJDcenter'],df_maxi['crab_intensity_2-20keV'],
	yerr=[df_maxi['crab_intensity_2-20keV_err'],df_maxi['crab_intensity_2-20keV_err']],
	label='MAXI/GSC (2-20 keV)',fmt='o',
	color="#0041ff",markersize=2.8,
	markerfacecolor="#0041ff",markeredgecolor="#0041ff")
axes.legend(loc='upper left',shadow=False)
axes.set_xlim(50133,58390)
# MJD 50133 --> 1996-02-20 00:00:00.000 UTC
# MJD 58390 --> 2018-09-29 00:00:00.000 UTC
axes.set_ylim(0.0,4.0)	
axes.set_xlabel('MJD (day)')			
axes.set_ylabel('X-ray intensity (Crab unit)')
axes2 = axes.twiny()
date1 = datetime.datetime(1996,2,20)
date2 = datetime.datetime(2018,9,29)
axes2.set_xlim(date1,date2)
axes.vlines(58315,1.8,2.7,colors='k')
axes.vlines(54530,1.8,2.7,colors='k')
#axes2.set_xticks([0.0,0.5])
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('./out/cirx1_asm_maxi_longlc_crab_intensity.pdf',dpi=300)




