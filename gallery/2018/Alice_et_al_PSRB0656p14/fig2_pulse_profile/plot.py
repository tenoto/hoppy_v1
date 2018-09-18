#!/usr/bin/env python

import pandas as pd 
import astropy.io.fits as pyfits 

df_radio = pd.read_csv('data/PSRJ0659+1414_2PC_radioprofile.csv')
#print(df_radio)
# Phase_Min,Phase_Max,Norm_Intensity

df_gamma = pd.read_csv('data/PSRJ0659+1414_2PC_gammaprofile.csv')
#print(df_gamma)

hdu_band1 = pyfits.open('data/psrb0656_merge_0mpu7_cl_o1u20cut_0p25to0p5keV.fht')
hdu_band1_xwid = 1/float(hdu_band1['PROFILE'].header['NBINPHAS'])

hdu_band2 = pyfits.open('data/psrb0656_merge_0mpu7_cl_o1u20cut_0p5to0p7keV.fht')
hdu_band2_xwid = 1/float(hdu_band2['PROFILE'].header['NBINPHAS'])

hdu_band3 = pyfits.open('data/psrb0656_merge_0mpu7_cl_o1u20cut_0p7to2p0keV.fht')
hdu_band3_xwid = 1/float(hdu_band2['PROFILE'].header['NBINPHAS'])

f = open('psrb0656_profile.qdp','w')
f.write('READ SERR 1 2\n')
f.write('@psrb0656_profile.pco\n')
for i in range(2):
	for index, row in df_radio.iterrows():
		dump = '%.8f 0.00390 %.8f 0.0\n' % (
			row['Phase_Min']+float(i),
			row['Norm_Intensity'])
		print(dump); f.write(dump)

f.write('NO NO NO NO\n')
for i in range(2):
	for index, row in df_gamma.iterrows():
		dump = '%.8f %.8f %.8f %.8f\n' % (
			0.5*(row['Phase_Min']+row['Phase_Max'])+float(i),
			0.5*(-row['Phase_Min']+row['Phase_Max']),
			row['GT100_WtCounts'],row['Unc_GT100_WtCounts'])
		print(dump); f.write(dump)

data = hdu_band1['PROFILE'].data
f.write('NO NO NO NO\n')
for i in range(len(data)):
	dump = '%.8f %.8f %.8f %.8f\n' % (
		data[i]['PULSE_PHASE'],0.5*hdu_band1_xwid,
		data[i]['RATE'],data[i]['RATE_ERR'])
	print(dump); f.write(dump)

data = hdu_band2['PROFILE'].data
f.write('NO NO NO NO\n')
for i in range(len(data)):
	dump = '%.8f %.8f %.8f %.8f\n' % (
		data[i]['PULSE_PHASE'],0.5*hdu_band2_xwid,
		data[i]['RATE'],data[i]['RATE_ERR'])
	print(dump); f.write(dump)

data = hdu_band3['PROFILE'].data
f.write('NO NO NO NO\n')
for i in range(len(data)):
	dump = '%.8f %.8f %.8f %.8f\n' % (
		data[i]['PULSE_PHASE'],0.5*hdu_band3_xwid,
		data[i]['RATE'],data[i]['RATE_ERR'])
	print(dump); f.write(dump)
f.close()

f = open('psrb0656_profile.pco','w')
dump = """
skip on 
time off
lab rotate 
lab csize 1.1
lab f 
yplo off 
lab pos y 2.8
r x 0.0 2.0

win 1 
loc 0.2 0.8 0.6 1.0
yplo on 1 
err off 1 
line on 1 
lab nx off
r y -0.1 1.2
lab y Radio intensity 

win 2 
loc 0.2 0.64 0.6 0.84
yplo on 2 
err on 2 
line on 2 
lab nx off
r y 10 75
lab y Counts

win 3
loc 0.2 0.48 0.6 0.68
yplo on 3
err on 3
line on 3 
lab nx off
r y 5.5 9.0 
lab y Cnts sec\\u-1\\d 

win 4
loc 0.2 0.32 0.6 0.52
yplo on 4
err on 4
line on 4 
lab nx off
r y 1.8 2.3
lab y Cnts sec\\u-1\\d 

win 5
loc 0.2 0.16 0.6 0.36
yplo on 5
err on 5
line on 5 
lab nx off
line step on 2,3,4,5
lab nx on
r y 0.8 1.4
lab y Cnts sec\\u-1\\d 

lab x Pulse Phase
"""
f.write(dump)
f.close()
