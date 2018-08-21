#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 August 21'
__version__ = '0.01'

"""
http://xte.mit.edu/asmlc/ASM.html
"""

import os 
import sys 
import yaml
import pandas as pd 
import argparse 

import matplotlib as mpl
import matplotlib.pyplot as plt 

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = '12'
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['xtick.top'] = 'True'
mpl.rcParams['ytick.right'] = 'True'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['axes.grid'] = 'True'
mpl.rcParams['axes.xmargin'] = '.05' #'.05'
mpl.rcParams['axes.ymargin'] = '.05'
mpl.rcParams['savefig.facecolor'] = 'None'
mpl.rcParams['savefig.edgecolor'] = 'None'
mpl.rcParams['savefig.bbox'] = 'tight'

ASM_SUMBAND_INTENSITY_CRAB = 75

class ASMCurve():

	def __init__(self,fyaml_param):
		self.fyaml_param = fyaml_param
		if not os.path.exists(self.fyaml_param):
			sys.stderr.write('file %s does not exist.' % self.fyaml_param)

		self.param = yaml.load(open(self.fyaml_param))

	def read(self):
		for line in open(self.param['input_asm_file']):
			cols = line.split()
			if cols[1] == 'COLUMNS':
				words = line.replace('% COLUMNS : ','').replace('\n','').split(', ')

		words.insert(0, "MJD")
		self.df = pd.read_csv(self.param['input_asm_file'],sep=' ',skiprows=5,names=words)
		print(self.df)

		self.df['crab_intensity'] = self.df['Sum Band Intensity']/ASM_SUMBAND_INTENSITY_CRAB
		self.df['crab_intensity_uncertainty'] = self.df['Sum Band Uncertainties']/ASM_SUMBAND_INTENSITY_CRAB		

	def dump_crab_intensity(self):
		outdir = os.path.dirname(self.param['output_asm_crabintensity'])
		cmd = 'rm -rf %s; mkdir -p %s' % (outdir,outdir)
		print(cmd);os.system(cmd)

		self.df.to_csv(self.param['output_asm_crabintensity'],
			columns=['MJD','crab_intensity','crab_intensity_uncertainty'],
			index=False)

	def plot(self):
		fig, axes = plt.subplots(1,1,#sharex=True,sharey=False,
			figsize=(self.param['panel_size'][0],self.param['panel_size'][1]))		
		axes.errorbar(
			self.df['MJD'],self.df['crab_intensity'],
			yerr=[self.df['crab_intensity_uncertainty'],self.df['crab_intensity_uncertainty']],
			fmt='o',
			color='r',
			markersize=1.0,
#			markerfacecolor=datagroup['markerfacecolor'],
#			markeredgecolor=datagroup['markeredgecolor'])
			)

		#axes.legend(loc=self.param['legend_location'],shadow=False)
		axes.set_title(self.param['title'])
		axes.set_xlim(self.param['xlim'])
		axes.set_ylim(self.param['ylim'])	
		axes.set_xlabel(self.param['xlabel'])			
		axes.set_ylabel(self.param['ylabel'])
		plt.subplots_adjust(wspace=0, hspace=0)
		plt.savefig(self.param['outpdf'])

	def run(self):
		self.read()
		self.dump_crab_intensity()
		self.plot()

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='plot_asm_light_curve.py',
		usage='plot_asm_light_curve.py fyaml_param',
		description='Plotting the ASM light curve.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'fyaml_param',metavar='fyaml_param',type=str,        
		help='yaml file for input parameters.')
	args = parser.parse_args()	
	
	asmlc = ASMCurve(args.fyaml_param)
	asmlc.run()
