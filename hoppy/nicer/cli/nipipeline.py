#!/usr/bin/env python

import argparse

from hoppy.nicer.nicer import NicerElf

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-03 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('nipipeline.py',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="""
Run the nicerl2 command to a sequence of ObsIDs. 
		"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command.')
	parser.add_argument('--setup_yamlfile', '-y', type=str, required=True, 
		help='setup yamlfile for the analysis paraemters.')
	parser.add_argument('--obsid_lstfile', '-o', type=str, required=True, 
		help='obsid list for the analyses')

	return parser

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	nielf = NicerElf(setup_yamlfile=args.setup_yamlfile,obsid_lstfile=args.obsid_lstfile)
	nielf.make_directory()
	nielf.run_nicerl2()	
	nielf.run_nibackgen3C50()
	nielf.plot_lightcurve()
	nielf.run_barycorr()
	nielf.devide_to_gti()

if __name__=="__main__":
	main()