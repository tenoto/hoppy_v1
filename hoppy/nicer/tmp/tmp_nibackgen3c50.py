#!/usr/bin/env python

import os 
import glob


indir = "/Users/enoto/work/drbv1/reporitory/heasarc/reprocess/nicer/24Jun2019_V6.26.1/"
outdir = "/Users/enoto/work/nicer/project/bgd/190806_nibackgen3C50_chk"

os.chdir(outdir)

num = 0
for obsid_path in glob.glob('%s/BKGD_RXTE_?/*' % indir):
	parent_path = os.path.dirname(obsid_path)
	sname = os.path.basename(parent_path)
	obsid = os.path.basename(obsid_path)

	outdir = 'out/%s/%s' % (sname,obsid)

	if os.path.exists(outdir):
		print("obsid %s has already existed. skip." % outdir)
	else:
		cmd  = 'mkdir -p %s' % outdir
		print(cmd);os.system(cmd)

		totspec = 'ni%s_nibackgen3C50_tot.pi' % obsid
		bkgspec = 'ni%s_nibackgen3C50_bgd.pi' % obsid	

		cmd  = 'nibackgen3C50 '
		cmd += 'rootdir=%s ' % parent_path
		cmd += 'obsid=%s ' % obsid
		cmd += 'bkgdir=%s ' % os.getenv('NICER_BKGDIR')
		cmd += 'totspec=%s ' % totspec
		cmd += 'bkgspec=%s ' % bkgspec
		print(cmd);os.system(cmd)

		cmd  = 'mv %s %s %s' % (totspec,bkgspec,outdir)
		print(cmd);os.system(cmd)
	num += 1 


