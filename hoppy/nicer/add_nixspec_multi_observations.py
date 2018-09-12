#!/usr/bin/env python

import os 
import sys 
import argparse
import pandas as pd 

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='add_nixspec_multi_observations.py',
		usage='add_nixspec_multi_observations.py xspec_multifit_fit.csv',
		description='add spectral files from a csv file.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'csvfile',metavar='csvfile',type=str,        
		help='Input csv file')
	parser.add_argument(
		'--outdir',metavar='--outdir',type=str,default="addspec",
		help='Output directory (default=addspec).')
	parser.add_argument(
		'--outname',metavar='--outname',type=str,default="nicer",
		help='Output file basename (default=nicer).')	
	parser.add_argument(
		'--rmffile',metavar='--rmffile',type=str,
		default="%s/nicer_v1.02.rmf" % os.getenv('NICER_RESP_PATH'),
		help='rmf file.')
	parser.add_argument(
		'--arffile',metavar='--arffile',type=str,
		default="%s/ni_xrcall_onaxis_v1.02_wo_detid14_34.arf" % os.getenv('NICER_RESP_PATH'),
		help='arf file.')	
	args = parser.parse_args()	

	print(args)

	cmd = 'rm -rf %s; mkdir -p %s' % (args.outdir,args.outdir)
	print(cmd);os.system(cmd)

	df = pd.read_csv(args.csvfile)

	fname_srcpha_list = '%s/%s_src.lst' % (args.outdir,args.outname)
	f = open(fname_srcpha_list,'w')
	for index, row in df.iterrows():
		f.write(row['phafile']+'\n')	
	f.close()
	cmd = 'faddpha.py %s --addtype C' % fname_srcpha_list
	print(cmd);os.system(cmd)
	merged_srcpha = fname_srcpha_list.replace('.lst','_merge.pha')
	for extnum in [1]:
		cmd  = 'fparkey %s %s+%d RESPFILE add=yes;' % (args.rmffile,merged_srcpha,extnum)		
		cmd += 'fparkey %s %s+%d ANCRFILE add=yes;' % (args.arffile,merged_srcpha,extnum)				
		print(cmd);os.system(cmd)

	fname_bgdpha_list = '%s/%s_bgd.lst' % (args.outdir,args.outname)
	f = open(fname_bgdpha_list,'w')
	for index, row in df.iterrows():
		f.write(row['backgrnd']+'\n')	
	f.close()
	cmd = 'faddpha.py %s --addtype R' % fname_bgdpha_list
	print(cmd);os.system(cmd)
	merged_bgdpha = fname_bgdpha_list.replace('.lst','_merge.pha')	

	for extnum in [1]:
		cmd  = 'fparkey %s %s+%d RESPFILE add=yes;' % (args.rmffile,merged_bgdpha,extnum)		
		cmd += 'fparkey %s %s+%d ANCRFILE add=yes;' % (args.arffile,merged_bgdpha,extnum)				
		print(cmd);os.system(cmd)

