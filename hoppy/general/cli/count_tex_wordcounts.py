#!/usr/bin/env python

import os
import argparse

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2021-07-01 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('count_tex_wordcounts.py',
		usage='%(prog)s texfile ',
		description="""
Show word counts in the tex file between '%%detex-abstract-begin' and '%%detex-abstract-end'
"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument(dest='texfile',type=str,help='tex file.')	
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')	
	parser.add_argument('--abstract_begin',type=str,default='detex-abstract-begin',
		help='keyword for the start of the abstract (default=detex-abstract-begin).')	
	parser.add_argument('--abstract_end',type=str,default='detex-abstract-end',
		help='keyword for the end of the abstract (default=detex-abstract-end).')	
	parser.add_argument('--main_begin',type=str,default='detex-main-begin',
		help='keyword for the start of the main body (default=detex-main-begin).')	
	parser.add_argument('--main_end',type=str,default='detex-main-end',
		help='keyword for the end of the main body (default=detex-main-end).')					
	return parser


def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	print("----- abstract -----")
	cmd = "sed -n '/%s/,/%s/p' %s " % (args.abstract_begin,args.abstract_end,args.texfile)
	cmd += "| detex "
	cmd += "| wc"
	print(cmd);os.system(cmd)

	print("----- main body -----")
	cmd = "sed -n '/%s/,/%s/p' %s " % (args.main_begin,args.main_end,args.texfile)
	cmd += "| detex "
	cmd += "| wc"
	print(cmd);os.system(cmd)

if __name__=="__main__":
	main()