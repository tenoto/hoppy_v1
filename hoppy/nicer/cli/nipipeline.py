#!/usr/bin/env python

import argparse

from hoppy.nicer.nicer import NicerManager

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

	niman = NicerManager(setup_yamlfile=args.setup_yamlfile,obsid_lstfile=args.obsid_lstfile)
	niman.make_directory()
	niman.run_nicerl2()	
	niman.run_nibackgen3C50()
	niman.plot_lightcurve()
	niman.run_barycorr()
	niman.fit_avespec()	
	niman.devide_to_segment()
	niman.fit_of_segment()
	niman.devide_to_block()
	niman.fit_of_block()
	niman.extract_segment_table()

if __name__=="__main__":
	main()