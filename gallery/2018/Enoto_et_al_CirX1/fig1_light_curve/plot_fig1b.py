#!/usr/bin/env python

# http://134.160.243.88/home/asai/maxibatHR/ 

import numpy as np
import pandas as pd 
import matplotlib as mpl
import matplotlib.pyplot as plt 
import datetime
import hoppy.timing.xtime as xtime 
import matplotlib.dates as mdates

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

plot_mjd_start = 58305.0
plot_mjd_stop  = 58345.0
plot_ymin = 0.0
plot_ymax = 22.0

xtime_start = xtime.XTime(plot_mjd_start,format='mjd',scale='utc')
yyyymmdd, hhmmss = xtime_start.get_utc('isot').split('T')
start_yyyy, start_mm, start_dd = yyyymmdd.split('-')
start_hour, start_minute, start_sec = hhmmss.split(':')
datetime_start = datetime.datetime(int(start_yyyy),int(start_mm),int(start_dd),
	int(start_hour),int(start_minute),0)

xtime_stop  = xtime.XTime(plot_mjd_stop,format='mjd',scale='utc')
yyyymmdd, hhmmss = xtime_stop.get_utc('isot').split('T')
start_yyyy, start_mm, start_dd = yyyymmdd.split('-')
start_hour, start_minute, start_sec = hhmmss.split(':')
datetime_stop = datetime.datetime(int(start_yyyy),int(start_mm),int(start_dd),
	int(start_hour),int(start_minute),0)

df = pd.read_csv('data/J1520-571-MAXI-BAT-HR_1orbit_lum.qdp',skiprows=1,skipfooter=13,
	delim_whitespace=True,names=['MJD','L36_MAXI','L36_MAXI_ERR','L36_BAT','L36_BAT_ERR','HARDNESS','HARDNESS_ERR'])
#print(df)

df_nicer = pd.read_csv('data/cirx1_nicer_pcfabs_diskbb_bbodyrad_Kabs_fit.csv')
Lx = 1.2e+32*9.4**2*(df_nicer['flux_1p0to10p0keV']/1e-12)/1e+37
Lx_err_min = 1.2e+32*9.4**2*(df_nicer['flux_1p0to10p0keV_err_min']/1e-12)/1e+37
Lx_err_max = 1.2e+32*9.4**2*(df_nicer['flux_1p0to10p0keV_err_max']/1e-12)/1e+37



np.where(df['MJD']>plot_mjd_start)

fig, axes = plt.subplots(1,1,figsize=(9.6,3.6))
axes.errorbar(df['MJD'],df['L36_MAXI']/10.0,
	yerr=[df['L36_MAXI_ERR']/10.0,df['L36_MAXI_ERR']/10.0],
	label='MAXI/GSC (2-10 keV)',fmt='o',
	color="#0041ff",markersize=2.8,
	markerfacecolor="#0041ff",markeredgecolor="#0041ff")
axes.errorbar(df_nicer['MJDOBS'],
	Lx,yerr=[Lx_err_min,Lx_err_max],
	label='NICER (1-10 keV)',fmt='s',
	color="#FF0000",markersize=4.2,
	markerfacecolor="#FF0000",markeredgecolor="#FF0000")

for Ncycle in range(915,930):
	T_MJD = 43076.32 + Ncycle * (16.57794 - 4.01e-5 * Ncycle)
	axes.vlines(T_MJD,plot_ymin,plot_ymax,colors='k',linestyles='dashed')
	print(T_MJD)

axes.legend(loc='upper left',shadow=False)
axes.set_xlim(plot_mjd_start,plot_mjd_stop)
# MJD 50133 --> 1996-02-20 00:00:00.000 UTC
# MJD 58390 --> 2018-09-29 00:00:00.000 UTC
#axes.set_ylim(0.0,4.0)	
#axes.set_xlabel('MJD (day)')			
#axes.set_ylabel('X-ray intensity (Crab unit)')
#axes2 = axes.twiny()
#date1 = datetime.datetime(1996,2,20)
#date2 = datetime.datetime(2018,9,29)
#axes2.set_xlim(date1,date2)
#axes2.set_xticks([0.0,0.5])
axes.set_ylim(plot_ymin,plot_ymax)
axes.set_xlabel('MJD (day)')			
axes.set_ylabel('X-ray Luminosity (10$^{37}$ erg s$^{-1}$)')

axes2 = axes.twiny()
axes2.set_xlim(datetime_start,datetime_stop)
axes2.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))

plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('out/cirx1_asm_maxi_short_luminosity.pdf',dpi=300)


