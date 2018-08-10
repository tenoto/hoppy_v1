#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import sys
import yaml 
import numpy as np 
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

class Plotter():
	def __init__(self,fyaml_setup):
		self.fyaml_setup = fyaml_setup

		if not os.path.exists(self.fyaml_setup):
			sys.stderr.write('error: yaml file %s does not exist.' % self.fyaml_setup)
			quit()
		self.param = yaml.load(open(self.fyaml_setup))

	def draw_x1y1_multigroup(self):
		fig, axes = plt.subplots(1,1,#sharex=True,sharey=False,
			figsize=(self.param['panel_size'][0],self.param['panel_size'][1]))

		for datagroup in self.param['datagroup']:
			did  = datagroup['id']
			dfile = datagroup['file']
			dtype = datagroup['type']
			xcolum = datagroup['xcolumn']
			xcolumn_error = datagroup['xcolumn_error']
			ycolum = datagroup['ycolumn']
			ycolumn_error = datagroup['ycolumn_error']			

			if dtype == 'dat':
				df = pd.read_csv(dfile,sep=' ')
			elif dtype == 'csv':
				df = pd.read_csv(dfile)
			x = df[xcolum].values			
			if xcolumn_error != "None":
				xe = df[xcolumn_error].values
			else:
				xe = np.zeros(len(x))
			y = df[ycolum].values
			if ycolumn_error != "None":
				ye = df[ycolumn_error].values
			else:
				ye = np.xeros(len(y))
			print(x,y,xe,ye)
			axes.errorbar(
				x,y,xerr=[xe,xe],yerr=[ye,ye],
				label=datagroup['label'],
				fmt=datagroup['fmt'],
				color=datagroup['color'],
				markersize=datagroup['markersize'],
				markerfacecolor=datagroup['markerfacecolor'],
				markeredgecolor=datagroup['markeredgecolor'])

		axes.legend(loc=self.param['legend_location'],shadow=False)
		axes.set_title(self.param['title'])
		axes.set_xlim(self.param['xlim'])
		axes.set_ylim(self.param['ylim'])	
		axes.set_xlabel(self.param['xlabel'])			
		axes.set_ylabel(self.param['ylabel'])
		plt.subplots_adjust(wspace=0, hspace=0)
		outdir = os.path.dirname(self.param['outpdf'])
		if not os.path.exists(outdir):
			cmd = 'mkdir -p %s' % outdir
			print(cmd);os.system(cmd)
		plt.savefig(self.param['outpdf'])

if __name__=="__main__":

	sys.stdout.write('... plot x1y1 plot for multi group ...\n')

	parser = argparse.ArgumentParser(
		prog='plot_x1y1_multigroup.py',
		usage='plot_x1y1_multigroup.py plot_x1y1_multigroup_setup.yaml',
		description='Plot data points for multiple groups.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'fyaml_setup',metavar='fyaml_setup',type=str,
		help='input yamlfile for setups.') 		
	args = parser.parse_args()	
	
	plotter = Plotter(args.fyaml_setup)
	plotter.draw_x1y1_multigroup()

