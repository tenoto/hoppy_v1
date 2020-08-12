#!/usr/bin/env python

import hoppy.xspec.xspec as xspec

import glob 
import os
import pandas as pd

indir = 'out/bkgd_rxte/'
rmffile = 'resp/nicer-rmf6s-teamonly-array50_wo14_34.rmf'
arffile = 'resp/nicer-consim135p-teamonly-array50_wo14_34.arf'
energy_bands = {"etot":[0.2,10.0],"e1":[0.2,0.4],"e2":[0.4,0.8],"e3":[0.8,1.6],"e4":[1.6,3.0],"e5":[3.0,6.0],"e6":[6.0,12.0],"e7":[12.0,15.0]}
outcsv = 'out/bkgd_rxte/bkgd_rxte_rates.csv'
flog = 'out/bkgd_rxte/bkgd_rxte_rates.log'
record_keywords = ["OBSID","DATEOBS","DATEEND","EXPOSURE","OBJECT","MJDOBS","MJD_DATEOBS","MJD_DATEEND"]

# Init
dict_out = {}
for keyword in record_keywords:
	dict_out[keyword] = []
for name in energy_bands:
	emin = energy_bands[name][0]
	emax = energy_bands[name][1]	
	keyword = 'rate_%s_%.1fto%.1fkeV' % (name,emin,emax)
	dict_out[keyword] = []
	keyword = 'rate_%s_%.1fto%.1fkeV_err' % (name,emin,emax)
	dict_out[keyword] = []

cmd = 'rm -f %s' % outcsv
print(cmd);os.system(cmd)

# Loop
f = open(flog,'w')
for obsid_dir in glob.glob('%s/*' % indir):
	obsid = os.path.basename(obsid_dir)
	src_pha = '%s/spec/ni%s_3c50_tot.pi' % (obsid_dir,obsid)
	bkg_pha = '%s/spec/ni%s_3c50_bkg.pi' % (obsid_dir,obsid)	

	if not os.path.exists(src_pha) or not os.path.exists(bkg_pha):
		f.write('%s error (file can not be found)' % obsid)
		continue

	cmd = 'rm -rf tmp_xspec_rate'
	os.system(cmd)

	xspec_pha = xspec.XspecPha(src_pha,outdir='tmp_xspec_rate',
		backgrnd=bkg_pha,rmffile=rmffile,arffile=arffile)
	xspec_pha.get_phafile_property()
	xspec_pha.show_property()	

	for keyword in record_keywords:
		dict_out[keyword].append(getattr(xspec_pha,keyword))

	for name in energy_bands:
		emin = energy_bands[name][0]
		emax = energy_bands[name][1]
		rate, rate_error = xspec_pha.get_rate_and_error(emin,emax)
		keyword = 'rate_%s_%.1fto%.1fkeV' % (name,emin,emax)
		dict_out[keyword].append(rate)
		keyword = 'rate_%s_%.1fto%.1fkeV_err' % (name,emin,emax)
		dict_out[keyword].append(rate_error)

	cmd = 'rm -rf tmp_xspec_rate'
	os.system(cmd)	

# end 
for name in energy_bands:
	emin = energy_bands[name][0]
	emax = energy_bands[name][1]	
	keyword = 'rate_%s_%.1fto%.1fkeV' % (name,emin,emax)
	record_keywords.append(keyword)
	keyword = 'rate_%s_%.1fto%.1fkeV_err' % (name,emin,emax)
	record_keywords.append(keyword)	

print(dict_out)
df = pd.DataFrame(dict_out,columns=record_keywords)
df.to_csv(outcsv)

f.close()



