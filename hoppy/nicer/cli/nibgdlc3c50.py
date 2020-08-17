#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import numpy as np 

import hoppy.xspec.xspec as xspec

__author__ = 'Teruaki Enoto'
__version__ = '0.02'
# v0.02 : 2020-08-17 
# v0.01 : 2020-08-08 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('nibgdlc3c50.py',
#		usage='%(prog)s ',
		description="""
"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--clfile', type=str, default=None,
		help='clfile (see fhelp nibackgen3C50)')	
	parser.add_argument('--ufafile', type=str, default=None,
		help='ufafile (see fhelp nibackgen3C50)')	
	parser.add_argument('--tbinfile', type=str, default=None,
		help='time bin file (csv format)')		
	parser.add_argument('--bkgidxdir', type=str, default=None,
		help='bkgidxdir (see fhelp nibackgen3C50)')	
	parser.add_argument('--bkglibdir', type=str, default=None,
		help='bkglibdir (see fhelp nibackgen3C50)')	
	parser.add_argument('--gainepoch', type=str, default=None,
		help='gainepoch (see fhelp nibackgen3C50)')		
	parser.add_argument('--dtmin', type=float, default=None,
		help='dtmin (see fhelp nibackgen3C50)')	
	parser.add_argument('--dtmax', type=float, default=None,
		help='dtmax (see fhelp nibackgen3C50)')	
	parser.add_argument('--hbgcut', type=float, default=None,
		help='hbgcut (see fhelp nibackgen3C50)')	
	parser.add_argument('--fpmofflist', type=str, default=None,
		help='fpmofflist (see fhelp nibackgen3C50)')			
	parser.add_argument('--rmffile', type=str, default=None, 
		help='rmffile ')	
	parser.add_argument('--arffile', type=str, default=None, 
		help='arffile ')	
	parser.add_argument('--ratebands',metavar='ratebands',type=str,
		default=[[0.4,6.0],[1.0,10.0]],
		help='rate energy bands (list) example:0.8-6.0,2.0-10.0 .') 		
	parser.add_argument('--outcsvfile', type=str, default=None,
		help='output file (csv format)')														
	return parser

def nibgdlc3c50(args):
	print("\n %s" % (sys._getframe().f_code.co_name))

	print(args)

	if type(args.ratebands) == list:
		ratebands = args.ratebands
	elif type(args.ratebands) == str:
		ratebands = xspec.string_to_list(args.ratebands)
	print(ratebands)

	outdir = 'tmp_nibgdlc3c50'
	cmd = 'rm -rf %s;mkdir -p %s' % (outdir,outdir)
	print(cmd);os.system(cmd)

	clfile_basename = os.path.splitext(os.path.basename(args.clfile))[0]
	ufafile_basename = os.path.splitext(os.path.basename(args.ufafile))[0]

	df = pd.read_csv(args.tbinfile)
	print(df)

	add_column_names_basename = [
		"rate_sub_cps","rate_sub_err",
		"rate_src_cps","rate_src_err",
		"rate_bkg_cps","rate_bkg_err"]

	add_column_names = []
	for band in ratebands:
		emin = band[0]
		emax = band[1]
		for name in add_column_names_basename:
			add_column_names.append("%s_%.1f_%.1fkeV" % (name,emin,emax))

	add_dictionary = {}
	for keyword in add_column_names:
		add_dictionary[keyword] = []

	for index, row in df.iterrows():
		print(row['START_MET'],row['STOP_MET'],row['Exposure_sec'])

		subdir = '%s/gti_%05d' % (outdir,index)
		cmd = 'rm -rf %s;mkdir -p %s' % (subdir,subdir)
		print(cmd);os.system(cmd)

		fgti_txt = '%s/tmp_gti%05d.txt' % (subdir,index)
		ftmp = open(fgti_txt,'w')
		ftmp.write('%.10f	%.10f\n' % (row['START_MET'],row['STOP_MET']))
		ftmp.close()		

		fgti_fits = '%s/tmp_gti%05d.gti' % (subdir,index)
		cmd = 'fconv_txt2gti.py -i %s -o %s ' % (fgti_txt,fgti_fits)
		print(cmd);os.system(cmd)
		
		fgti_clevt = '%s/tmp_%s_gti%05d.evt' % (subdir,clfile_basename,index)
		cmd = 'niextract-events %s %s timefile=%s' % (args.clfile,fgti_clevt,fgti_fits)
		print(cmd);os.system(cmd)

		fgti_ufaevt = '%s/tmp_%s_gti%05d.evt' % (subdir,ufafile_basename,index)
		cmd = 'niextract-events %s %s timefile=%s' % (args.ufafile,fgti_ufaevt,fgti_fits)
		print(cmd);os.system(cmd)		

		totspec_name = 'tmp_%s_gti%05d_3c50_tot' % (clfile_basename,index)
		bkgspec_name = 'tmp_%s_gti%05d_3c50_bkg' % (clfile_basename,index)

		# prepare a script for nicerl2 for each self.obsid
		fcmd = '%s/tmp_%s_3c50.sh' % (subdir,clfile_basename)
		flog = '%s/tmp_%s_3c50.log' % (subdir,clfile_basename)
		f = open(fcmd,'w')
		dump  = '#!/bin/sh -f\n'
		dump += 'nibackgen3C50 ' 
		dump += 'rootdir=\'NONE\' ' 
		dump += 'obsid=\'NONE\' '
		dump += 'bkgidxdir=\'%s\' ' % args.bkgidxdir
		dump += 'bkglibdir=\'%s\' ' % args.bkglibdir
		dump += 'gainepoch=\'%s\' ' % args.gainepoch
		dump += 'calevtdir=\'NONE\' '
		dump += 'clfile=\'%s\' ' % fgti_clevt
		dump += 'ufafile=\'%s\' ' % fgti_ufaevt
		dump += 'totspec=\'%s\' ' % totspec_name
		dump += 'bkgspec=\'%s\' ' % bkgspec_name
		dump += 'dtmin=%.1f dtmax=%.1f hbgcut=%.1f ' % (args.dtmin,args.dtmax,args.hbgcut)
		dump += 'clobber=yes '
		dump += 'chatter=5 '
		dump += '> %s 2>&1 ' % flog
		dump += '\n'
		f.write(dump)
		f.close()	

		# run the script.
		cmd  = 'chmod +x %s\n' % fcmd
		cmd += './%s' % fcmd
		print(cmd);os.system(cmd)

		cmd = 'mv %s*.pi %s*.pi %s' % (totspec_name,bkgspec_name,subdir)
		print(cmd);os.system(cmd)			

		totspec = '%s/%s.pi' % (subdir,totspec_name)
		bkgspec = '%s/%s.pi' % (subdir,bkgspec_name)

		if not os.path.exists(totspec):
			for keyword in add_dictionary:
				add_dictionary[keyword].append(np.nan)
			continue 

		for band in ratebands:
			emin = band[0]
			emax = band[1]

			cmd = 'rm -rf tmp_xspec_rate;\n'
			os.system(cmd)

			xspec_pha = xspec.XspecPha(totspec,
				outdir='tmp_xspec_rate',backgrnd=bkgspec,
				rmffile=args.rmffile,arffile=args.arffile)
			rate_sub, rate_sub_error = xspec_pha.get_rate_and_error(emin,emax)
			print("Rate: %.5f +/- %.5f (%.2f-%.2f keV)" % (rate_sub,rate_sub_error,emin,emax))
			cmd = 'rm -rf tmp_xspec_rate;\n'
			os.system(cmd)
			keyword = 'rate_sub_cps_%.1f_%.1fkeV' % (emin,emax)
			add_dictionary[keyword].append(rate_sub)
			keyword = 'rate_sub_err_%.1f_%.1fkeV' % (emin,emax)			
			add_dictionary[keyword].append(rate_sub_error)	

			xspec_pha = xspec.XspecPha(totspec,
				outdir='tmp_xspec_rate',backgrnd=None,
				rmffile=args.rmffile,arffile=args.arffile)
			rate_src, rate_src_error = xspec_pha.get_rate_and_error(emin,emax)
			print("Rate: %.5f +/- %.5f (%.2f-%.2f keV)" % (rate_src,rate_src_error,emin,emax))
			cmd = 'rm -rf tmp_xspec_rate;\n'
			os.system(cmd)
			keyword = 'rate_src_cps_%.1f_%.1fkeV' % (emin,emax)
			add_dictionary[keyword].append(rate_src)
			keyword = 'rate_src_err_%.1f_%.1fkeV' % (emin,emax)			
			add_dictionary[keyword].append(rate_src_error)		

			xspec_pha = xspec.XspecPha(bkgspec,
				outdir='tmp_xspec_rate',backgrnd=None,
				rmffile=args.rmffile,arffile=args.arffile)
			rate_bkg, rate_bkg_error = xspec_pha.get_rate_and_error(emin,emax)
			print("Rate: %.5f +/- %.5f (%.2f-%.2f keV)" % (rate_bkg,rate_bkg_error,emin,emax))
			cmd = 'rm -rf tmp_xspec_rate;\n'
			os.system(cmd)
			keyword = 'rate_bkg_cps_%.1f_%.1fkeV' % (emin,emax)
			add_dictionary[keyword].append(rate_bkg)
			keyword = 'rate_bkg_err_%.1f_%.1fkeV' % (emin,emax)			
			add_dictionary[keyword].append(rate_bkg_error)		

	cmd = 'rm -rf %s' % (outdir)
	print(cmd);os.system(cmd)	

	df_add = pd.DataFrame(add_dictionary,columns=(add_column_names))
	df_new = df.join([df_add])
	df_new.to_csv(args.outcsvfile)


def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	nibgdlc3c50(args)

if __name__=="__main__":
	main()