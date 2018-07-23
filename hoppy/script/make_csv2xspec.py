#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import sys 
import argparse
import yaml
import numpy as np 
import pandas as pd 

class CSVtoXSPEC():
	def __init__(self,filelist,yamlfile):
		self.filelist = filelist
		self.yamlfile = yamlfile

		if not os.path.exists(self.filelist):
			sys.stderr.write('file %s does not exist.' % self.filelist)
			exit()
		if not os.path.exists(self.yamlfile):
			sys.stderr.write('file %s does not exist.' % self.yamlfile)
			exit()			

		self.param = yaml.load(open(self.yamlfile))

	def make_csv2xspec(self,outcsvfile):
		df = pd.read_csv(self.filelist,names=('data_id','phafile','backgrnd','rmffile','arffile'))
		df["modelxcm"] = [self.param["modelxcm"] for i in range(len(df))]
		df["binminsig"] = [self.param["binminsig"] for i in range(len(df))]
		df["binmaxbin"] = [self.param["binmaxbin"] for i in range(len(df))]		
		df["fitemin"] = [self.param["fitemin"] for i in range(len(df))]				
		df["fitemax"] = [self.param["fitemax"] for i in range(len(df))]						
		df["plotxmin"] = [self.param["plotxmin"] for i in range(len(df))]								
		df["plotxmax"] = [self.param["plotxmax"] for i in range(len(df))]										
		df["plotymin"] = [self.param["plotymin"] for i in range(len(df))]												
		print(df)

		"""
		csvfile = '%s/%s.csv' % (param['outdir'],param['outdir'].replace('/','_'))
COLUMNS = ['data_id','src_pha','bgd_pha','rmffile','arffile',
	'model_xcm','bin_minsig','bin_maxbin','fit_emin','fit_emax','param_fit','param_fix_for_flux',
	'bgdmodel_err_qdp'
	]
df_out = pd.DataFrame(csv_data_list,columns=COLUMNS)
df_out.to_csv(csvfile)
		"""

if __name__=="__main__":

	sys.stdout.write('\n... prepare csvfile for xspec fittings ...\n')

	parser = argparse.ArgumentParser(
		prog=__file__,
		usage='python %s filelist yamlfile outcsvfile',
		description='prepare csvfile for xspec fittings',
		epilog='',
		add_help=True,
		)
	parser.add_argument(
		'filelist',metavar='filelist',type=str,
		help='observation file list.') 
	parser.add_argument(
		'yamlfile',metavar='yamlfile',type=str,
		help='parameter file in the yaml format.') 	
	parser.add_argument(
		'outcsvfile',metavar='outcsvfile',type=str,
		help='output csvfile') 		
	args = parser.parse_args()	
	print(args)

	csv = CSVtoXSPEC(args.filelist,args.yamlfile)
	csv.make_csv2xspec(args.outcsvfile)