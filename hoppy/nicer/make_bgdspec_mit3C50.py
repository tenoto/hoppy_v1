#!/usr/bin/env python

#__name__    = 'make_bgdspec_mit3C50.py'
__author__  = 'Teru Enoto (for Rons model with Mikes application)'
__version__ = '1.01'
__date__    = '2018 July 11'
__version_note__ = 'exposure map normalized by grand_total_exposure (ML comment)'
__header_keyword__ =  '3C50_kyoto_v%s'  % __version__

"""
see the document 'NICER XTI Background Model 3C50'
Description and User Instructions
(BGMod_3C50.tar)
v.3 2018 0531

with Mikes application 
"""
import os
import sys
import subprocess
import numpy as np 
import astropy.io.fits as pyfits
import pandas as pd 
import argparse
from datetime import datetime 

BGD_MODEL_VERSION = "BGMod_3C50"
BGD_LIBRARY_PATH = "/Users/enoto/Dropbox/enoto/research/nicer/data/mitbgd/%s" % BGD_MODEL_VERSION
BG_GROUP_3C50_TABLE  = "%s/bg_group_3C50.table" % BGD_LIBRARY_PATH
DAY_NOISE_3C50_TABLE = "%s/day_noise_3C50.table" % BGD_LIBRARY_PATH
NUM_of_MPU = 7
NUM_of_FPM = 8
NUM_of_IBG = 5
NUM_of_HREJ = 6
NUM_of_NZ = 15
RMFFILE = '%s/nicer_v1.02.rmf' % os.getenv('NICER_RESP_PATH')
#ARFFILE = '%s/ni_xrcall_onaxis_v1.02.arf' % os.getenv('NICER_RESP_PATH')
ARFFILE = '%s/ni_xrcall_onaxis_v1.02_wo_detid14_34.arf' % os.getenv('NICER_RESP_PATH')

IBG52_HREJ52_pha_pairs = [[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2],[3,1],[3,2],[3,3],[4,1],[4,2],[4,3],[5,1],[5,2],[5,3],[5,4]]

def get_lc_rate(lcfits):
	"""
	Return the average of the 'Rate' column
	"""
	hdu = pyfits.open(lcfits)
	return np.mean(hdu['RATE'].data['RATE'])

class ExposureMap_IBGvsHREJ():
	"""
	Exposure map of IBG52 (vertical) and HREJ52 (horizontal)
	"""
	def __init__(self):
		self.matrix_ibg52_hrej52 = np.zeros((NUM_of_HREJ,NUM_of_IBG))
		self.matrix_ibg52_hrej52_rawexp = np.zeros((NUM_of_HREJ,NUM_of_IBG))
		self.ground_total_exposure = 0.0

	def identiy_bin(self,IBG52,HREJ52):
		if IBG52 < 0.10:
			if HREJ52<0.15: 
				i = 0; j = 0;
			elif 0.15<=HREJ52<0.40:
				i = 0; j = 1;				
			elif 0.40<=HREJ52<1.00:
				i = 0; j = 2;								
			elif 1.00<=HREJ52<2.00:
				i = 0; j = 3;												
			elif 2.00<=HREJ52<5.00:
				i = 0; j = 4;																
		elif 0.10 <= IBG52 < 0.20:
			if HREJ52<0.18: 
				i = 1; j = 0;
			elif 0.18<=HREJ52<0.40:
				i = 1; j = 1;
			elif 0.40<=HREJ52<1.00:
				i = 1; j = 2;
			elif 1.00<=HREJ52<2.00:
				i = 1; j = 3; 
			elif 2.00<=HREJ52<5.00:
				i = 1; j = 4;
		elif 0.20 <= IBG52 < 0.40:
			if HREJ52<0.20: 
				i = 2; j = 0;
			elif 0.20<=HREJ52<0.40:
				i = 2; j = 1; 
			elif 0.40<=HREJ52<1.00:
				i = 2; j = 2 ;
			elif 1.00<=HREJ52<2.00:
				i = 2; j = 3;
			elif 2.00<=HREJ52<5.00:
				i = 2; j = 4;
		elif 0.40 <= IBG52 < 1.00:
			if HREJ52<0.20: 
				i = 3; j = 0;
			elif 0.20<=HREJ52<0.40:
				i = 3; j = 1;
			elif 0.40<=HREJ52<0.90:
				i = 3; j = 2;
			elif 0.90<=HREJ52<2.00:
				i = 3; j = 3;
			elif 2.00<=HREJ52<5.00:
				i = 3; j = 4;
		elif 1.00 <= IBG52 < 3.00:
			if HREJ52<0.20: 
				i = 4; j = 0;
			elif 0.20<=HREJ52<0.40:
				i = 4; j = 1;
			elif 0.40<=HREJ52<1.00:
				i = 4; j = 2; 
			elif 1.00<=HREJ52<2.30:
				i = 4; j = 3;
			elif 2.30<=HREJ52<5.00:
				i = 4; j = 4;
		elif 3.00 <= IBG52 < 10.00:
			if HREJ52<0.20: 
				i = 5; j = 0;
			elif 0.20<=HREJ52<0.55:
				i = 5; j = 1;
			elif 0.55<=HREJ52<1.00:
				i = 5; j = 2;
			elif 1.00<=HREJ52<2.00:
				i = 5; j = 3;
			elif 2.00<=HREJ52<5.00:				
				i = 5; j = 4;

		return i, j 

	def add(self,IBG52,HREJ52,exposure,NFPM):
		i, j = self.identiy_bin(IBG52,HREJ52)

		target_phaname = 'bg_group_3C50_ngt_%d%d.pha' % (i+1,j+1)
		for line2 in open(BG_GROUP_3C50_TABLE):
			cols2 = line2.split()
			phaname = cols2[0]
			if phaname == target_phaname:
				IBG_REF = float(cols2[1])
				break 

		print(i,j,exposure,IBG52,NFPM)
		weight = exposure * IBG52 / IBG_REF * NFPM / 52.0
		self.matrix_ibg52_hrej52[i][j] += weight		
		self.matrix_ibg52_hrej52_rawexp[i][j] += exposure
		self.ground_total_exposure += float(exposure)
		return i,j

	def show(self):
		print(self.matrix_ibg52_hrej52)

	def show_rawexp(self):
		print(self.matrix_ibg52_hrej52_rawexp)

	def get_grand_exposure(self,norm_type='ground_total_exposure_v1.01'):
		if norm_type == 'ground_total_exposure_v1.01':
			out = self.ground_total_exposure
		elif norm_type == 'original_v1.00':
			out = sum(sum(self.matrix_ibg52_hrej52))
		else:
			sys.stderr.write('error: no corresponding norm_type...\n')
			exit()
		return out

	def get_normalized_exposuremap(self):
		self.normalized_matrix_ibg52_hrej52 =  self.matrix_ibg52_hrej52 / self.get_grand_exposure()
		return self.normalized_matrix_ibg52_hrej52	

class ExposureMap_NZ():
	"""
	Exposure map of NZ52 
	"""
	def __init__(self):
		self.hist_nz52 = np.zeros(NUM_of_NZ)
		self.hist_nz52_rawexp = np.zeros(NUM_of_NZ)
		self.ground_total_exposure = 0.0

	def add(self,NZ52,exposure,NFPM):
		if NZ52 < 190.0: 
			k = 0
		elif 190.0 <= NZ52 < 210.0:
			k = 1
		elif 210.0 <= NZ52 < 250.0:			
			k = 2
		elif 250.0 <= NZ52 < 300.0:
			k = 3
		elif 300.0 <= NZ52 < 400.0:
			k = 4
		elif 400.0 <= NZ52 < 500.0:
			k = 5
		elif 500.0 <= NZ52 < 600.0:
			k = 6
		elif 600.0 <= NZ52 < 700.0:
			k = 7
		elif 700.0 <= NZ52 < 800.0:
			k = 8
		elif 800.0 <= NZ52 < 950.0:
			k = 9
		elif 950.0 <= NZ52 < 1100.0:
			k = 10
		elif 1100.0 <= NZ52 < 1250.0:
			k = 11
		elif 1250.0 <= NZ52 < 1400.0:
			k = 12
		elif 1400.0 <= NZ52 < 1600.0:
			k = 13
		elif 1600.0 <= NZ52:
			k = 14

		target_phaname = 'day_noise_nz%02d.pha' % k
		if k == 0 or k == 14:
			weight = 0.0
		else:
			for line2 in open(DAY_NOISE_3C50_TABLE):
				cols2 = line2.split()
				phaname = cols2[0]
				print(cols2)
				if phaname == target_phaname:
					NZ52_REF = float(cols2[1])
					break 
			weight = exposure * NZ52 / NZ52_REF * NFPM / 52.0
		self.hist_nz52_rawexp[k] += exposure
		self.hist_nz52[k] += weight 
		self.ground_total_exposure += float(exposure)
		return k

	def show(self):
		print(self.hist_nz52)

	def show_rawexp(self):
		print(self.hist_nz52_rawexp)

	def get_grand_exposure(self,norm_type='ground_total_exposure_v1.01'):
		if norm_type == 'ground_total_exposure_v1.01':
			out = self.ground_total_exposure
		elif norm_type == 'original_v1.00':
			out = sum(self.hist_nz52)
		else:
			sys.stderr.write('error: no corresponding norm_type...\n')
			exit()
		return out

	def get_normalized_exposuremap(self):
		if self.get_grand_exposure() == 0.0:
			self.normalized_hist_nz52 = self.hist_nz52
		else:
			self.normalized_hist_nz52 = self.hist_nz52 / self.get_grand_exposure()
		return self.normalized_hist_nz52

def add_gti_divide(phafile,gtilist):
	tmp_gtilist = 'tmp_%s' % os.path.basename(gtilist)
	cmd = 'rm -f %s' % tmp_gtilist
	print(cmd);os.system(cmd)
	f = open(tmp_gtilist,'w')
	for line in open(gtilist):
		f.write(line.replace('# ',''))
	f.close()

	df = pd.read_csv(tmp_gtilist, sep=' ')
	print(df)

	hdu_append_list = [
		pyfits.Column(name='START', format='1D', array=df['START'].values.tolist()),
		pyfits.Column(name='STOP', format='1D', array=df['STOP'].values.tolist()),
		pyfits.Column(name='ObsID', format='10A', array=df['ObsID'].values.tolist()),		
		pyfits.Column(name='GTI_NUM', format='I', array=df['ID'].values.tolist()),			
		pyfits.Column(name='Exposure', format='1D', array=df['EXP(s)'].values.tolist()),		
		pyfits.Column(name='NUM_FPM', format='I', array=df['NFPM'].values.tolist()),		
		pyfits.Column(name='IBG', format='1D', array=df['IBG(cps)'].values.tolist()),
		pyfits.Column(name='HREJ', format='1D', array=df['HREJ(cps)'].values.tolist()),						
		pyfits.Column(name='NZ', format='1D', array=df['NZ(cps)'].values.tolist()),						
		pyfits.Column(name='IBG52', format='1D', array=df['IBG52'].values.tolist()),								
		pyfits.Column(name='HREJ52', format='1D', array=df['HREJ52'].values.tolist()),
		pyfits.Column(name='NZ52', format='1D', array=df['NZ52'].values.tolist()),	
		pyfits.Column(name='i_IBG', format='I', array=df['i_IBG'].values.tolist()),	
		pyfits.Column(name='j_HREJ', format='I', array=df['j_HREJ'].values.tolist()),
		pyfits.Column(name='k_NZ', format='I', array=df['k_NZ'].values.tolist())
		]

	hdul = pyfits.open(phafile)
	hdu = pyfits.BinTableHDU.from_columns(hdu_append_list)
	hdul.append(hdu)
	hdul.writeto(phafile,overwrite=True)
	cmd = 'fparkey GTI_DIVIDE %s+%d EXTNAME add=yes' % (phafile,len(hdul)-1)
	print(cmd);os.system(cmd)
	cmd = 'rm -f %s' % tmp_gtilist
	print(cmd);os.system(cmd)

if __name__ == '__main__':
	# ==============================
	# Check HEASoft environment
	# ==============================
	cmd = 'which fhelp'
	resp = subprocess.Popen(cmd.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	if len(resp.stdout.readlines()) == 0:
		sys.stderr.write('error: set HEASOFT (e.g., heainit)\n')
		exit()

	# ==============================
	# Get input parameters 
	# ==============================
	help_message = """
	(example) %s ufafile clfile 
	""" % sys.argv[0]

	parser = argparse.ArgumentParser(description=help_message)
	#parser.add_argument('obsid_path', metavar='obsid_path',type=str,        
	#	help='input obsid_path or @input_obsid_path.lst for multiple obsids.')
	parser.add_argument('ufafile', metavar='ufafile',type=str,        
		help='input ufafile or @ufafile.lst for multiple obsids.')
	parser.add_argument('clfile', metavar='clfile',type=str,        
		help='input clfile or @clfile.lst for multiple obsids.')
	parser.add_argument('--outdir',dest='outdir',action='store',
		help='output directory', default='mitbgd') 
	parser.add_argument('--prefix',dest='prefix',action='store',
		help='prefix of the output files', default='mitbgd') 
	parser.add_argument('--recreate', action='store_true', 
		help='flag to recreate of the output directory.') 
	parser.add_argument('--exclude',dest='exclude',action='store',
		help='excluding detid (e.g., --exclude 14,34', default='14,34') 
	parser.add_argument('--tbin',dest='tbin',action='store',
		help='GTI time scale (default 100 sec)', default=100.0) 
	parser.add_argument('--tlimit',dest='tlimit',action='store',
		help='GTI time scale lower threshold (default 60 sec)', default=60.0) 
	parser.add_argument('--title',dest='title',action='store',
		help='Title for plotting', default=None) 
	parser.add_argument('--lctbin',dest='lctbin',action='store',
		help='light curve time bin', default=20.0) 
	parser.add_argument('--lcemin',dest='lcemin',action='store',
		help='light curve energy low', default=0.2) 
	parser.add_argument('--lcemax',dest='lcemax',action='store',
		help='light curve energy high', default=15.0) 
	args = parser.parse_args()
	print(args)

	# ==============================
	# Prepare output directory 
	# ==============================
	outdir = args.outdir
	outdir_log = '%s/log' % outdir 
	if args.recreate:
		cmd = 'rm -rf %s; mkdir -p %s' % (outdir,outdir_log)
	elif not os.path.exists(outdir):
		cmd = 'mkdir -p %s' % outdir_log
	else:
		sys.stderr.write('output directory has already exited: %s\n' % outdir)
		quit()
	print(cmd);os.system(cmd)

	prefix = args.prefix

	# ==============================
	# Set up observation path as inputs 
	# ==============================
	ufafile_list = []
	if args.ufafile[0] is '@':
		for line in open(args.ufafile[1:]):
			ufafile_list.append(line.strip().strip('/'))
	else:
		ufafile_list.append(args.ufafile)
	print("ufafile_list: %s" % ufafile_list)

	clfile_list = []
	if args.clfile[0] is '@':
		for line in open(args.clfile[1:]):
			clfile_list.append(line.strip().strip('/'))
	else:
		clfile_list.append(args.clfile)
	print("clfile_list: %s" % clfile_list)

	# ==============================
	# Prepare input ufa / cl event list file 
	# ==============================
	fname_ufaevt_list = '%s/%s_ufaevt.lst' % (outdir_log,prefix)
	f = open(fname_ufaevt_list,'w')
	for ufaevt in ufafile_list:
		f.write(ufaevt+'\n')
	f.close()

	fname_clevt_list = '%s/%s_clevt.lst' % (outdir_log,prefix)
	f = open(fname_clevt_list,'w')
	for clevt in clfile_list:
		f.write(clevt+'\n')
	f.close()

	# ==============================
	# Set screening of the exlucded DET_ID 
	# ==============================
	exclude_detid_list = []
	for i in args.exclude.split(','):
		exclude_detid_list.append('%02d' % int(i))
	print('exclude_detid_list: %s\n' % exclude_detid_list)

	num_of_fpm = 56 - 4 - len(exclude_detid_list) 

	# ==============================
	# Make GTI information file (divided with tbin)
	# ==============================
	#outdir_gti = '%s/gti' % outdir
	#cmd = 'rm -rf %s; mkdir -p %s' % (outdir_gti,outdir_gti)
	#print(cmd);os.system(cmd)

	fname_org_gti_list  = '%s/%s_original_gti.lst' % (outdir_log,prefix)
	f = open(fname_org_gti_list,'w')
	for clevt in clfile_list:
		obsid = os.path.basename(clevt).split('_')[0].split('ni')[-1]
		hdu = pyfits.open(clevt)
		for row in hdu['GTI'].data:
			gti_exp = row['STOP'] - row['START']		
			message = '%.1f %.1f %.1f %s' % (row['START'],row['STOP'],gti_exp, obsid)
			f.write(message+'\n')
	f.close()

	gti_list = []
	fname_gti_list1 = '%s/%s_divided_gti1.lst' % (outdir_log,prefix)
	f = open(fname_gti_list1,'w')
	for clevt in clfile_list:
		obsid = os.path.basename(clevt).split('_')[0].split('ni')[-1]
		hdu = pyfits.open(clevt)
		print(clevt,obsid)
		for row in hdu['GTI'].data:
			gti_exp = row['STOP'] - row['START']
			if gti_exp < args.tlimit:
				message = '%.1f %.1f %.1f (skip)' % (row['START'],row['STOP'],gti_exp)
				print(message)
				continue
			ngti_subdiv = int(gti_exp/float(args.tbin))
			message = '%.1f %.1f %.1f %d' % (row['START'],row['STOP'],gti_exp,ngti_subdiv)
			print(message)		
			for i in range(ngti_subdiv+1):
				print(row['START'],row['STOP'],i,args.tbin)
				tstart = row['START'] + float(i) * float(args.tbin)
				tstop  = row['START'] + float(i+1) * float(args.tbin)
				if tstop > row['STOP']:
					tstop = row['STOP']
				tgti = tstop - tstart
				if tgti > args.tlimit:
					nid = len(gti_list)
					gti_list.append([nid,tstart,tstop,tgti,obsid])
					message = '  %d %.1f %.1f %.1f %s' % (nid,tstart,tstop,tgti,obsid)
					f.write('%.1f %.1f %d %s %.1f\n' % (tstart,tstop,nid,obsid,tgti))
				else:
					message = '  -- %.1f %.1f %.1f (skip)' % (tstart,tstop,tgti)
				print(message)
	f.close()			
	print('GTI: %s' % gti_list)

	if len(gti_list) == 0:
		sys.stderr.write('no GTI...\n')
		exit()

	# ==============================
	# Slow-only ufa file 
	# ==============================
	outdir_evtfile = '%s/evtsel' % outdir
	cmd = 'mkdir -p %s' % outdir_evtfile
	print(cmd);os.system(cmd)

	# 4.1 Make trimmed, slow-chain ufa files to support down-stream queries 
	fname_ufaslow = '%s/%s_ufaslow.evt' % (outdir_evtfile,prefix)
	cmd  = 'niextract-events '
	cmd += 'filename="@%s[PI=-10:1800,EVENT_FLAGS=bxxx1x000]" ' % fname_ufaevt_list
	cmd += 'eventsout="%s" gti=GTI ' % fname_ufaslow
	print(cmd);os.system(cmd)

	cmd  = 'gzip %s' % fname_ufaslow
	print(cmd);os.system(cmd)
	fname_ufaslow += '.gz'
	print(fname_ufaslow)

	# ==============================
	# Trumpet-selected ufa file 
	# ==============================
	# 4.2 Make trumpet-selected ufa files to screen FPMs & evaluate IBG, NZ

	fname_ufatsel = '%s/%s_ufatsel.evt' % (outdir_evtfile,prefix)
	cmd  = 'fselect %s %s ' % (fname_ufaslow,fname_ufatsel)
	cmd += '"ISNULL(PI_RATIO) || PI_RATIO < (1.1+120/PI)" '
	print(cmd);os.system(cmd)

	cmd  = 'gzip %s' % fname_ufatsel
	print(cmd);os.system(cmd)
	fname_ufatsel += '.gz'
	print(fname_ufatsel)

	# ==============================
	# Hatched-rejected ufa file 
	# ==============================
	# 4.3 Make hatchet-rejected ufa files to evaluate HREJ
	fname_ufahrej = '%s/%s_ufahrej.evt' % (outdir_evtfile,prefix)
	cmd  = 'fselect %s %s ' % (fname_ufaslow,fname_ufahrej)
	cmd += '"PI_RATIO >= 1.54" '
	print(cmd);os.system(cmd)

	cmd  = 'gzip %s' % fname_ufahrej
	print(cmd);os.system(cmd)
	fname_ufahrej += '.gz'
	print(fname_ufahrej)

	# ==============================
	# Merge clean cl events 
	# ==============================
	# 4.5 FPM Screening the cl event lists prior to extracting spectra per GTI
	fname_clevt_merge = '%s/%s_clmerge.evt' % (outdir_evtfile,prefix)
	cmd  = 'ftmerge @%s %s' % (fname_clevt_list,fname_clevt_merge)
	print(cmd);os.system(cmd)

	cmd  = 'gzip %s' % fname_clevt_merge
	print(cmd);os.system(cmd)
	fname_clevt_merge += '.gz'

	# ==============================
	# Screen DET_ID (exclude noisy modules) for cleam events 
	# ==============================
	# screen of detid 
	fname_clevt_screen = '%s/%s_clmerge_screen.evt' % (outdir_evtfile,prefix)
	if len(exclude_detid_list) == 0:
		cmd = 'cp %s %s ' % (fname_clevt_merge, fname_clevt_screen)
	else:
		cmd  = 'fselect %s %s ' % (fname_clevt_merge,fname_clevt_screen)
		cmd += '"'
		for detid in exclude_detid_list:
			cmd += '(DET_ID != %s)' % (detid)
			if exclude_detid_list.index(detid) < len(exclude_detid_list)-1:
				cmd += ' && '
		cmd += '"'
	print(cmd);os.system(cmd)

	cmd  = 'gzip %s' % fname_clevt_screen
	print(cmd);os.system(cmd)
	fname_clevt_screen += '.gz'

	# ==============================
	# Screen DET_ID (exclude noisy modules) for ufa events [Trumpet-selected]
	# ==============================
	# 4.6 FPM Screening ufa event lists prior to extracting IBG, HREJ, NZ
	# do 4.2 fselect command with addition "(DET_ID != 14) && (DET_ID != 34)"
	# do 4.3 fselect command with addition "(DET_ID != 14) && (DET_ID != 34)"

	fname_ufatsel_screen = '%s/%s_ufatsel_screen.evt' % (outdir_evtfile,prefix)
	if len(exclude_detid_list) == 0:
		cmd = 'cp %s %s ' % (fname_ufatsel, fname_ufatsel_screen)
	else:
		cmd  = 'fselect %s %s ' % (fname_ufatsel,fname_ufatsel_screen)
		cmd += '"'
		for detid in exclude_detid_list:
			cmd += '(DET_ID != %s)' % (detid)
			if exclude_detid_list.index(detid) < len(exclude_detid_list)-1:
				cmd += ' && '
		cmd += '"'
	print(cmd);os.system(cmd)

	cmd  = 'gzip %s' % fname_ufatsel_screen
	print(cmd);os.system(cmd)
	fname_ufatsel_screen += '.gz'

	# ==============================
	# Screen DET_ID (exclude noisy modules) for ufa events [Hatched-rejected]
	# ==============================
	fname_ufahrej_screen = '%s/%s_ufahrej_screen.evt' % (outdir_evtfile,prefix)
	if len(exclude_detid_list) == 0:
		cmd = 'cp %s %s ' % (fname_ufahrej,fname_ufahrej_screen)
	else:
		cmd  = 'fselect %s %s ' % (fname_ufahrej,fname_ufahrej_screen)
		cmd += '"'
		for detid in exclude_detid_list:
			cmd += '(DET_ID != %s)' % (detid)
			if exclude_detid_list.index(detid) < len(exclude_detid_list)-1:
				cmd += ' && '
		cmd += '"'
	print(cmd);os.system(cmd)

	cmd  = 'gzip %s' % fname_ufahrej_screen
	print(cmd);os.system(cmd)
	fname_ufahrej_screen += '.gz'

	# ==============================
	# Make cleaned spectrum / IBG, HREJ, NZ light curves for each GTI
	# ==============================
	# Make spectra and evaluate parameters for BG Mdeling
	# 6.1 xselect: extract spectra, per GTI, from event lists
	# for selected FPMs to obtain target spectra
	# (from *cl.evt.gz or *cl50.evt.gz or selection-specialized files)
	outdir_gtichk = '%s/gtichk' % outdir
	cmd = 'rm -rf %s; mkdir -p %s' % (outdir_gtichk,outdir_gtichk)
	print(cmd);os.system(cmd)

	fname_input_to_script = '%s/%s_gti_input_to_mit_script.lst' % (outdir_gtichk,prefix)
	fin = open(fname_input_to_script,'w')
	fname_gti_list2 = '%s/%s_divided_gti1_add_ratechk.lst' % (outdir_gtichk,prefix)
	flog = open(fname_gti_list2,'w')
	dump = '# START STOP ID ObsID EXP(s) NFPM IBG(cps) HREJ(cps) NZ(cps) IBG52 HREJ52 NZ52'
	flog.write(dump+'\n')
	for line in open(fname_gti_list1):
		cols = line.split()
		tstart = float(cols[0])
		tstop  = float(cols[1])
		gtinum = int(cols[2])
		obsid = cols[3]
		texp   = float(cols[4])
		subdir = '%s/gti%05d' % (outdir_gtichk,gtinum)
		cmd = 'rm -rf %s; mkdir -p %s' % (subdir,subdir)
		print(cmd);os.system(cmd)

		fname_clevt_screen_gti_evt = '%s/%s_clmerge_screen_gti%d.evt' % (subdir,prefix,gtinum)
		fname_clevt_screen_gti_pha = fname_clevt_screen_gti_evt.replace('.evt','.pha')
		#cmd  = 'niextract-events '
		#cmd += 'filename="%s[TIME=%.1f:%.1f]" ' % (fname_clevt_screen,tstart,tstop)
		#cmd += 'eventsout="%s" gti=GTI ' % fname_clevt_screen_gti_evt
		#print(cmd);os.system(cmd)

		cmd = """
xselect <<EOF
xsel
read event %s .
yes
filter time scc 
%.4f,%.4f
x 
set phaname PI
extract spectrum
save spectrum
%s
exit
no
EOF
mv xsel_timefile.asc xselect.log %s
""" % (fname_clevt_screen,tstart,tstop,fname_clevt_screen_gti_pha,subdir)
		print(cmd);os.system(cmd)

		fname_ufatsel_screen_gti_ibg_lc = '%s/%s_ufatsel_screen_gti%d_ibg.lc' % (subdir,prefix,gtinum)
		cmd = """
xselect <<EOF
xsel
read event %s .
yes
filter time scc 
%.4f,%.4f
x 
filter pha_cutoff 1500 1700
set binsize 1
show filter 
extract curve
save curve
%s
exit
no
EOF
mv xsel_timefile.asc xselect.log %s
""" %  (fname_ufatsel_screen,tstart,tstop,fname_ufatsel_screen_gti_ibg_lc,subdir)
		print(cmd);os.system(cmd)
		gti_ibg_average = get_lc_rate(fname_ufatsel_screen_gti_ibg_lc)

		fname_ufahrej_screen_gti_hrej_lc = '%s/%s_ufahrej_screen_gti%d_hrej.lc' % (subdir,prefix,gtinum)
		cmd = """
xselect <<EOF
xsel
read event %s .
yes
filter time scc 
%.4f,%.4f
x 
filter pha_cutoff 300 1800
set binsize 1
show filter 
extract curve
save curve
%s
exit
no
EOF
mv xsel_timefile.asc xselect.log %s
""" %  (fname_ufahrej_screen,tstart,tstop,fname_ufahrej_screen_gti_hrej_lc,subdir)
		print(cmd);os.system(cmd)
		gti_hrej_average = get_lc_rate(fname_ufahrej_screen_gti_hrej_lc)

		fname_ufatsel_screen_gti_nz_lc = '%s/%s_ufatsel_screen_gti%d_nz.lc' % (subdir,prefix,gtinum)
		cmd = """
xselect <<EOF
xsel
read event %s .
yes
filter time scc 
%.4f,%.4f
x 
filter pha_cutoff 0 20
set binsize 1
show filter 
extract curve
save curve
%s
exit
no
EOF
rm -f xsel_timefile.asc xselect.log
mv xsel_timefile.asc xselect.log %s
""" %  (fname_ufatsel_screen,tstart,tstop,fname_ufatsel_screen_gti_nz_lc,subdir)
		print(cmd);os.system(cmd)
		gti_nz_average = get_lc_rate(fname_ufatsel_screen_gti_nz_lc)

		dump = '%s %d %.4f %.4f %.4f' % (fname_clevt_screen_gti_pha,num_of_fpm,
			gti_ibg_average,gti_hrej_average,gti_nz_average)
		fin.write(dump+'\n')

		dump  = '%.1f %.1f %d %s %.1f ' % (tstart,tstop,gtinum,obsid,texp)
		dump += '%d %.4f %.4f %.4f ' % (num_of_fpm,gti_ibg_average,gti_hrej_average,gti_nz_average)
		dump += '%.4f %.4f %.4f ' % (
			52.0 / float(num_of_fpm) * gti_ibg_average,
			52.0 / float(num_of_fpm) * gti_hrej_average,
			52.0 / float(num_of_fpm) * gti_nz_average
			)
		flog.write(dump+'\n')
	fin.close()
	flog.close()


	# ==============================
	# Add GTI information into exposure map 
	# ==============================
	expmap_IBGvsHREJ = ExposureMap_IBGvsHREJ()
	expmap_NZ = ExposureMap_NZ()

	fname_gti_list3 = '%s/%s_divided_gti2_expmap_selected.lst' % (outdir_gtichk,prefix)
	flog = open(fname_gti_list3,'w')
	dump = '# START STOP ID ObsID EXP(s) NFPM IBG(cps) HREJ(cps) NZ(cps) IBG52 HREJ52 NZ52 i_IBG j_HREJ k_NZ'
	flog.write(dump+'\n')
	for line in open(fname_gti_list2):
		cols = line.split()
		if cols[0] == '#':
			continue 
		texp   = float(cols[4])
		NFPM = float(cols[5])	
		IBG52  = float(cols[9])
		HREJ52 = float(cols[10])
		NZ52   = float(cols[11])
		print(IBG52,HREJ52,NZ52,texp,NFPM)
		if not (IBG52 <= 10.0 and HREJ52 <= 5.0 and NZ52 <= 1600.0):
			print("skip...")
			continue
		i, j = expmap_IBGvsHREJ.identiy_bin(IBG52,HREJ52)
		if not [i,j] in IBG52_HREJ52_pha_pairs:
			print("skip...")
			continue
		expmap_IBGvsHREJ.add(IBG52,HREJ52,texp,NFPM)
		k = expmap_NZ.add(NZ52,texp,NFPM)
		print(i,j,k)
		flog.write(line.strip() + ' %d %d %d\n' % (i, j, k))
	expmap_IBGvsHREJ.show()
	expmap_IBGvsHREJ.show_rawexp()
	expmap_NZ.show()
	expmap_NZ.show_rawexp()
	flog.close()

	# ==============================
	# Normalize the exposure map 
	# ==============================
	norm_expmap_IBGvsHREJ = expmap_IBGvsHREJ.get_normalized_exposuremap() 
	norm_expmap_NZ = expmap_NZ.get_normalized_exposuremap()
	print(norm_expmap_IBGvsHREJ)
	print(norm_expmap_NZ)

	# ==============================
	# Extract source spectrum 
	# ==============================
	outdir_prod = '%s/product' % outdir
	cmd = 'mkdir -p %s' % outdir_prod
	print(cmd);os.system(cmd)

	fname_gti_list3_fits = '%s/%s_divided_gti2_expmap_selected.gti' % (outdir_gtichk,prefix)
	cmd  = 'fconv_txt2gti.py '
	cmd += '-i %s ' % fname_gti_list3
	cmd += '-o %s ' % fname_gti_list3_fits
	print(cmd);os.system(cmd)

	fname_clevt_screen_gti = '%s/%s_clscr_gtisel.evt' % (outdir_prod,prefix)
	cmd = 'nicerclean '
	cmd += 'infile=%s ' % fname_clevt_screen
	cmd += 'outfile=%s ' % fname_clevt_screen_gti
	cmd += 'gtifile=%s ' % fname_gti_list3_fits
	print(cmd);os.system(cmd)

	fname_clpha_screen_gti = '%s/%s_clscr_gtisel.pha' % (outdir_prod,prefix)
	cmd = 'fxselect_extract_spec.py '
	cmd += '-i %s ' % fname_clevt_screen_gti
	cmd += '-o %s ' % fname_clpha_screen_gti
	cmd += '-r %s ' % RMFFILE
	cmd += '-a %s ' % ARFFILE
	print(cmd);os.system(cmd)

	fname_cllc_screen_gti = '%s/%s_clscr_gtisel.lc' % (outdir_prod,prefix)
	cmd = 'fxselect_extract_curve.py '
	cmd += '-i %s ' % fname_clevt_screen_gti
	cmd += '-o %s ' % fname_cllc_screen_gti
	cmd += '-t %.1f ' % float(args.lctbin)
	cmd += '-d %.1f -u %.1f ' % (float(args.lcemin),float(args.lcemax))
	print(cmd);os.system(cmd)

	fname_cllc_screen_gti_ps = '%s/%s_clscr_gtisel_lc.ps' % (outdir_prod,prefix)
	cmd = 'lcurve 1 %s ' % fname_cllc_screen_gti
	cmd += 'window="-" dtnb=20 nbint=INDEF <<EOF\n'
	cmd += 'null.flc\n' 
	cmd += 'yes\n'
	cmd += '/xw\n'
	cmd += 'lwid 5 \n'
	cmd += 'la f %.1f-%.1f keV %.1f s bin\n' % (
		float(args.lcemin),float(args.lcemax),float(args.lctbin))
	cmd += 'hard %s/cps\n' % fname_cllc_screen_gti_ps
	cmd += 'exit\n'
	cmd += 'EOF\n'
	print(cmd);os.system(cmd)

	cmd = 'ps2pdf %s' % fname_cllc_screen_gti_ps
	print(cmd);os.system(cmd)
	cmd = 'rm -f %s' % fname_cllc_screen_gti_ps
	print(cmd);os.system(cmd)
	cmd = 'mv %s %s' % (os.path.basename(fname_cllc_screen_gti_ps).replace('.ps','.pdf'),outdir_prod)
	print(cmd);os.system(cmd)
	cmd = 'rm -f null.flc'
	print(cmd);os.system(cmd)

	# ==============================
	# mathpha 
	# ==============================
	cmd = 'ln -s %s .' % BGD_LIBRARY_PATH
	print(cmd);os.system(cmd)

	expr_list = []
	bgdpha_exposure = 0.0 # v1.01 
	for i in range(NUM_of_HREJ):
		for j in range(NUM_of_IBG):
			bg_pha = '%s/bg_group_3C50_ngt_%d%d.pha' % (BGD_MODEL_VERSION,i+1,j+1)
			print(i,j,norm_expmap_IBGvsHREJ[i][j],bg_pha)
			if norm_expmap_IBGvsHREJ[i][j] > 0.0:
				expr_list.append("'%s' * %.6f " % (bg_pha,norm_expmap_IBGvsHREJ[i][j]))
				hdu = pyfits.open(bg_pha)
				bgdpha_exposure += norm_expmap_IBGvsHREJ[i][j] * float(hdu['SPECTRUM'].header['EXPOSURE'])

	for k in range(1,NUM_of_NZ-1):
		noise_pha = '%s/day_noise_nz%02d.pha' % (BGD_MODEL_VERSION,k)
		print(k,norm_expmap_NZ[k],noise_pha)
		if norm_expmap_NZ[k] > 0.0:
			expr_list.append("'%s' * %.6f" % (noise_pha,norm_expmap_NZ[k]))
			#hdu = pyfits.open(noise_pha)
			#bgdpha_exposure += norm_expmap_NZ[k] * float(hdu['SPECTRUM'].header['EXPOSURE'])
	print(expr_list)		

	fname_mathpha_expr = '%s/%s_mathpha_expr.txt' % (outdir_gtichk,prefix)
	f = open(fname_mathpha_expr,'w')
	for i in expr_list:
		f.write('%s' % i)
		if expr_list.index(i) < len(expr_list)-1:
			f.write(' + ')
	f.close()

	cmd = 'ln -s %s .' % fname_mathpha_expr
	print(cmd);os.system(cmd)

	fname_bgd_pha = '%s/%s_%s.pha' % (outdir_prod,prefix,BGD_MODEL_VERSION)
	cmd  = 'mathpha '
	cmd += 'expr=@%s ' % os.path.basename(fname_mathpha_expr)
	cmd += 'units=R outfil="%s" ' % os.path.basename(fname_bgd_pha)
	cmd += 'exposure=%.3f ' % bgdpha_exposure
	cmd += 'errmeth=gaussian properr=yes ncomments=0 areascal=NULL clobber=yes'
	print(cmd);os.system(cmd)

	cmd = 'mv %s %s' % (os.path.basename(fname_bgd_pha),os.path.dirname(fname_bgd_pha))
	print(cmd);os.system(cmd)

	cmd = 'rm -f %s;' % os.path.basename(fname_mathpha_expr)
	print(cmd);os.system(cmd)

	#hdu = pyfits.open(fname_clpha_screen_gti)
	#src_exposure = float(hdu['SPECTRUM'].header['EXPOSURE'])
	for extnum in [1]:
		#cmd  = 'fparkey %.7f %s+%d EXPOSURE;' % (src_exposure,fname_bgd_pha,extnum)
		cmd  = 'fparkey %s %s+%d BGDVERSN add=yes;' % (__header_keyword__,fname_bgd_pha,extnum)
		cmd += 'fparkey %s %s+%d RESPFILE add=yes;' % (RMFFILE,fname_bgd_pha,extnum)		
		cmd += 'fparkey %s %s+%d ANCRFILE add=yes;' % (ARFFILE,fname_bgd_pha,extnum)				
		print(cmd);os.system(cmd)

	cmd = 'rm -f %s' % BGD_MODEL_VERSION
	print(cmd);os.system(cmd)

	# ==============================
	# add GTI 
	# ==============================
	add_gti_divide(fname_clpha_screen_gti,fname_gti_list3)
	add_gti_divide(fname_bgd_pha,fname_gti_list3)

	# ==============================
	# xspec show 
	# ==============================
	fname_read_xcm = '%s/%s_read.xcm' % (outdir_prod,prefix)
	f = open(fname_read_xcm,'w')
	dump = 'data 1 %s\n' % fname_clpha_screen_gti 
	dump += 'back 1 %s\n' % fname_bgd_pha
	dump += 'data 2 %s\n' % fname_clpha_screen_gti
	dump += 'data 3 %s\n' % fname_bgd_pha
	dump += 'resp 3 %s\n' % RMFFILE
	dump += 'arf 3 %s\n' % ARFFILE
	dump += 'ignore 1-3:**-0.2 15.0-**\n'
	dump += 'setplot energy\n'
	f.write(dump)
	f.close()

	hdu = pyfits.open(fname_clpha_screen_gti)
	target  = hdu[0].header['OBJECT']
	dateobs = hdu[0].header['DATE-OBS']
	exposure = hdu[1].header['EXPOSURE']
	title = '%s' % target
	subtitle = '%s %s (%.1f s)' % (obsid,dateobs,exposure)

	fname_pha_ps = fname_clpha_screen_gti.replace('.pha','_pha.ps')
	cmd = 'xspec<<EOF\n'
	cmd += 'setplot device tmp/null\n'			
	cmd += '@%s\n' % fname_read_xcm
	cmd += 'cpd /xw\n'
	#cmd += 'setp rebin 3 40 1\n'
	cmd += 'ipl ld\n'
	#cmd += 'ipl d\n'
	cmd += 'r x 0.2 15.0\n'
	#cmd += 'r y -8.0 8.0\n'
	cmd += 'lwid 5 on 1..100\n'
	cmd += 'lwid 5 \n'
	if args.title == None:
		cmd += 'la t %s \n' % title 
		cmd += 'la f %s \n' % subtitle
	else:
		cmd += 'la t %s \n' % args.title 
	cmd += 'time off\n'
	cmd += 'hard %s/cps\n' % fname_pha_ps 
	cmd += 'quit\n'
	cmd += 'quit\n'
	print(cmd);os.system(cmd)

	cmd = 'ps2pdf %s' % fname_pha_ps
	print(cmd);os.system(cmd)

	fname_pha_pdf = os.path.basename(fname_pha_ps.replace('.ps','.pdf'))
	cmd = 'mv %s %s' % (fname_pha_pdf, outdir_prod)
	print(cmd);os.system(cmd)

	cmd = 'rm -f %s ' % fname_pha_ps
	print(cmd);os.system(cmd)




