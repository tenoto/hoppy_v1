#!/usr/bin/env python

# http://134.160.243.88/home/asai/maxibatHR/ 

import pandas as pd 
import matplotlib as mpl
import matplotlib.pyplot as plt 

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
axes.set_ylim(0.0,4.0)	
axes.set_xlabel('MJD (day)')			
axes.set_ylabel('X-ray intensity (Crab unit)')
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('./out/cirx1_asm_maxi_longlc_crab_intensity.pdf',dpi=300)



