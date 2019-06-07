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
		df["plotymax"] = [self.param["plotymax"] for i in range(len(df))]
		df["plotymineeuf"] = [self.param["plotymineeuf"] for i in range(len(df))]												
		df["plotymaxeeuf"] = [self.param["plotymaxeeuf"] for i in range(len(df))]																														
		df["ploty2min"] = [self.param["ploty2min"] for i in range(len(df))]																
		df["ploty2max"] = [self.param["ploty2max"] for i in range(len(df))]																		
		#df["flagrun"] = [1 for i in range(len(df))]																				
		print(df)

		outdir = os.path.dirname(outcsvfile)
		if not os.path.exists(outdir):
			cmd = 'mkdir -p %s' % outdir
			print(cmd);os.system(cmd)

		df.to_csv(outcsvfile)

if __name__=="__main__":

	sys.stdout.write('\n... prepare csvfile for xspec fittings ...\n')

	parser = argparse.ArgumentParser(
		prog=__file__,
		usage='python %s filelist yamlfile outcsvfile' % __file__,
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