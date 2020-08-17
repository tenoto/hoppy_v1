#!/usr/bin/env python

import os
import re
import sys
import argparse
import pandas as pd 

import matplotlib.pylab as plt 

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-17 : original version


def get_parser():
	parser = argparse.ArgumentParser('nibgdlc3c50_plot.py',
		description="""
Plot the output csvfile from nibgdlc3c50.py
"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--infile', type=str, default=None, required=True,
		help='input csv file')
	parser.add_argument('--outfile', type=str, default=None, 
		help='output pdf file')	
	return parser

def nibgdlc3c50_plot(args):
	print("\n %s" % (sys._getframe().f_code.co_name))

	if args.outfile == None:
		args.outfile = os.path.splitext(os.path.basename(args.infile))[0] + '.pdf'

	df = pd.read_csv(args.infile)
	print(df)

	keyword_main = []
	for keyword in df:
		print(keyword)
		if re.match('rate_sub_cps', keyword) != None:
			keyword_main.append(keyword)
	print(len(keyword_main))

	tcenter = 0.5 * (df['START_MET'] + df['STOP_MET'])-df['START_MET'][0]

	fig, axs = plt.subplots(len(keyword_main),1,
		figsize=(11.69,8.27), sharex=True, gridspec_kw={'hspace': 0})	
	for i in range(len(keyword_main)):
		print(i,keyword_main[i],keyword_main[i].replace('cps','err'))
		ylabel = keyword_main[i].split('rate_sub_cps_')[-1].replace('_','-').replace('keV',' keV')

		axs[i].errorbar(tcenter,df[keyword_main[i]],
		yerr=df[keyword_main[i].replace('cps','err')],
		marker='o',ls='')
		axs[i].set_ylabel(ylabel)
		axs[i].axhline(0.0,ls='--',color='k')
	axs[-1].set_xlabel('Time (sec)')
	fig.align_ylabels(axs)
	plt.tight_layout(pad=2)
	plt.rcParams["font.family"] = "serif"
	plt.rcParams["mathtext.fontset"] = "dejavuserif"	
	plt.savefig(args.outfile)


def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	nibgdlc3c50_plot(args)

if __name__=="__main__":
	main()