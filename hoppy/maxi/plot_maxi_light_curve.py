#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 August 21'
__version__ = '0.01'

"""
http://134.160.243.88/top/index.html
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

MAXI_CRAB_RATE = [3.3,1.7,1.2,0.35]
# picked up by eye, http://134.160.243.88/star_data/J0534+220/J0534+220_00055058g_lc_all.gif

class MAXICurve():
	def __init__(self,fyaml_param):
		"""
		http://134.160.243.88/top/lc_readme.txt
		"""
		self.fyaml_param = fyaml_param
		if not os.path.exists(self.fyaml_param):
			sys.stderr.write('file %s does not exist.' % self.fyaml_param)

		self.param = yaml.load(open(self.fyaml_param))

	def read(self):
		words = [
			'MJDcenter',
			'2-20keV','2-20keV_err',
			'2-4keV','2-4keV_err',
			'4-10keV','4-10keV_err',
			'10-20keV','10-20keV_err']
		self.df = pd.read_csv(self.param['input_maxi_file'],sep=' ',names=words)
		print(self.df)

		self.df['crab_intensity_2-20keV'] = self.df['2-20keV']
		self.df['crab_intensity_2-20keV_err'] = self.df['2-20keV_err']		
		self.df['crab_intensity_2-4keV'] = self.df['2-4keV']
		self.df['crab_intensity_2-4keV_err'] = self.df['2-4keV_err']				
		self.df['crab_intensity_4-10keV'] = self.df['4-10keV']
		self.df['crab_intensity_4-10keV_err'] = self.df['4-10keV_err']	
		self.df['crab_intensity_10-20keV'] = self.df['10-20keV']
		self.df['crab_intensity_10-20keV_err'] = self.df['10-20keV_err']	

	def dump_crab_intensity(self):
		outdir = os.path.dirname(self.param['output_maxi_crabintensity'])
		cmd = 'rm -rf %s; mkdir -p %s' % (outdir,outdir)
		print(cmd);os.system(cmd)

		self.df.to_csv(self.param['output_maxi_crabintensity'],
			columns=['MJDcenter','crab_intensity_2-20keV','crab_intensity_2-20keV_err'],
			index=False)

	def plot(self):
		fig, axes = plt.subplots(1,1,#sharex=True,sharey=False,
			figsize=(self.param['panel_size'][0],self.param['panel_size'][1]))		
		axes.errorbar(
			self.df['MJDcenter'],self.df['crab_intensity_2-20keV'],
			yerr=[self.df['crab_intensity_2-20keV_err'],self.df['crab_intensity_2-20keV_err']],
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
		prog='plot_maxi_light_curve.py',
		usage='plot_maxi_light_curve.py fyaml_param',
		description='Plotting the MAXI light curve.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'fyaml_param',metavar='fyaml_param',type=str,        
		help='yaml file for input parameters.')
	args = parser.parse_args()	
	
	maxilc = MAXICurve(args.fyaml_param)
	maxilc.run()
