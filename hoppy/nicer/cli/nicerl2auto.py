#!/usr/bin/env python

import os 
import glob 
import astropy.io.fits as pyfits

def run_nicerl2(obsid,indir):
	fversion = os.popen('fversion').read().rstrip()
	mkffile = '%s/auxil/ni%s.mkf.gz' % (indir,obsid)
	mkfhdu = pyfits.open(mkffile)
	object_name = mkfhdu[1].header['OBJECT']

	print(object_name,obsid)

	outdir = 'reprocess/%s/%s' % (fversion,object_name)

	if os.path.exists('%s/%s' % (outdir,obsid)):
		print('output directory (%s) has already existed.' % obsid)
		return -1

	if not os.path.exists(outdir):
		cmd = 'mkdir -p %s;\n' % (outdir)
	else:
		cmd = ''

	cmd += 'cp -r %s %s/;\n' % (indir,outdir)
	print(cmd);os.system(cmd)

	current_dir = os.getcwd()
	os.chdir(outdir)
	cmd = 'nicerl2 indir=%s clobber=YES' % (obsid)
	print(cmd);os.system(cmd)
	os.chdir(current_dir)	

for line in open('/Users/enoto/work/drbv1/soft/git/hoppy/data/nicer/nicerl2_input.lst'):
	obsid_text = line.split()[0]
	if obsid_text in ['#','%']:
		print('skip')
		continue 
	for indir in glob.glob('data/nicer/data/obs/*/%s' % obsid_text):
		obsid = os.path.basename(indir)
		print(obsid,indir)
		run_nicerl2(obsid,indir)





