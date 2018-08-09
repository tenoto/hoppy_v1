#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import sys 
import argparse
import yaml 
import pandas as pd 

from hoppy.xspec.xspec import XspecPha, XspecFitLog

class MonitoringManager():
	def __init__(self,csvfile,yamlfile):
		self.csvfile = csvfile 
		self.yamlfile = yamlfile

		if not os.path.exists(self.csvfile):
			sys.stderr.write('error: file %s does not exist.\n' % self.csvfile)
			exit()			
		if not os.path.exists(self.yamlfile):
			sys.stderr.write('error: file %s does not exist.\n' % self.yamlfile)
			exit()						
		print("csvfile: %s" % self.csvfile)
		print("yamlfile: %s" % self.yamlfile)

		self.df    = pd.read_csv(self.csvfile)
		self.df = self.df.drop(columns='Unnamed: 0')
		self.param = yaml.load(open(self.yamlfile))

		self.outcsvfile = '%s/%s' % (
			self.param['outdir'],
			os.path.basename(self.csvfile).replace('.csv','_fit.csv'))
		print("out csvfile: %s" % self.outcsvfile)

	def run(self):		
		add_column_names = ["grp_pha","OBSID","DATEOBS","DATEEND","EXPOSURE",
			"MJD_DATEOBS","MJD_DATEEND",
			"OBJECT","TELESCOP","INSTRUME","MJDOBS","TSTART","TSTOP",
			"fxcm_fit","title",
			"chisquare","reduced_chisquare","dof","probability"]
		add_values_matrix = [[] for i in range(len(add_column_names))]

		add_column_names_rate = []
		for ebands in self.param['ratebands']:
			add_column_names_rate.append("rate_%sto%skeV" % (str(ebands[0]).replace('.','p'),str(ebands[1]).replace('.','p')))
			add_column_names_rate.append("rate_%sto%skeV_err" % (str(ebands[0]).replace('.','p'),str(ebands[1]).replace('.','p')))			
		add_values_matrix_rate = [[] for i in range(len(add_column_names_rate))]

		add_column_names_flux = []
		for ebands in self.param['fluxbands']:
			add_column_names_flux.append("flux_%sto%skeV" % (str(ebands[0]).replace('.','p'),str(ebands[1]).replace('.','p')))
			add_column_names_flux.append("flux_%sto%skeV_err_min" % (str(ebands[0]).replace('.','p'),str(ebands[1]).replace('.','p')))			
			add_column_names_flux.append("flux_%sto%skeV_err_max" % (str(ebands[0]).replace('.','p'),str(ebands[1]).replace('.','p')))						
		add_values_matrix_flux = [[] for i in range(len(add_column_names_flux))]

		add_column_names_par = []
		for parnum in self.param['parerrnum']:
			add_column_names_par.append("par_%d" % parnum)
			add_column_names_par.append("par_%d_err_min" % parnum)
			add_column_names_par.append("par_%d_err_max" % parnum)
		add_values_matrix_par = [[] for i in range(len(add_column_names_par))]

		exposure_list = []
		for index, dataset in self.df.iterrows():
			outdir = '%s/%s' % (self.param['outdir'],dataset['data_id'])
			xspec_pha = XspecPha(
				phafile=dataset['phafile'],
				outdir=outdir,
				backgrnd=dataset['backgrnd'],
				rmffile=dataset['rmffile'],
				arffile=dataset['arffile'],
				modelxcm=dataset['modelxcm'],
				binminsig=dataset['binminsig'],binmaxbin=dataset['binmaxbin'],
				fitemin=dataset['fitemin'],fitemax=dataset['fitemax'],
				plotxmin=dataset['plotxmin'],plotxmax=dataset['plotxmax'],
				plotymin=dataset['plotymin'],plotymax=dataset['plotymax'],
				ploty2min=dataset['ploty2min'],ploty2max=dataset['ploty2max'],
				ratebands=self.param['ratebands'],
				fluxbands=self.param['fluxbands'],
				parerrnum=self.param['parerrnum'])
			xspec_pha.run()
			for key, value in vars(xspec_pha).iteritems(): # python 2 
				if key in add_column_names:
					i = add_column_names.index(key)
					add_values_matrix[i].append(value)

			for i in range(len(self.param['ratebands'])):
				add_values_matrix_rate[i*2].append(xspec_pha.ratelist[i][0])
				add_values_matrix_rate[i*2+1].append(xspec_pha.ratelist[i][1])				

			for i in range(len(self.param['fluxbands'])):
				add_values_matrix_flux[i*3].append(xspec_pha.fluxlist[i][0])
				add_values_matrix_flux[i*3+1].append(xspec_pha.fluxlist[i][1])				
				add_values_matrix_flux[i*3+2].append(xspec_pha.fluxlist[i][2])	

			for i in range(len(self.param['parerrnum'])):
				add_values_matrix_par[i*3].append(xspec_pha.parerrorlist[i][0])
				add_values_matrix_par[i*3+1].append(xspec_pha.parerrorlist[i][1])				
				add_values_matrix_par[i*3+2].append(xspec_pha.parerrorlist[i][2])	

			print(add_values_matrix)

		add_dictionary = {}
		for key in add_column_names:
			i = add_column_names.index(key)
			add_dictionary[key] = add_values_matrix[i]
		for i in range(len(add_column_names_rate)):
			add_dictionary[add_column_names_rate[i]] = add_values_matrix_rate[i] 
		for i in range(len(add_column_names_flux)):
			add_dictionary[add_column_names_flux[i]] = add_values_matrix_flux[i] 			
		for i in range(len(add_column_names_par)):
			add_dictionary[add_column_names_par[i]] = add_values_matrix_par[i] 						
		self.df_add = pd.DataFrame(add_dictionary,
			columns=(add_column_names+add_column_names_rate+add_column_names_flux+add_column_names_par))
		self.df_new = self.df.join([self.df_add])
		self.df_new.to_csv(self.outcsvfile)

if __name__=="__main__":

	sys.stdout.write('\n... run multiple xspec fitting ...\n')

	parser = argparse.ArgumentParser(
		prog=__file__,
		usage='python %s csvfile yamlfile' % __file__,
		description='xspec fitting of multiple observations.',
		epilog='',
		add_help=True,
		)
	parser.add_argument(
		'csvfile',metavar='csvfile',type=str,
		help='observation list in the csv format.') 
	parser.add_argument(
		'yamlfile',metavar='yamlfile',type=str,
		help='parameter file in the yaml format.') 	
	args = parser.parse_args()	
	print(args)

	manager = MonitoringManager(args.csvfile,args.yamlfile)
	manager.run()