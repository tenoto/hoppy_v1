#!/usr/bin/env python

import argparse

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-01 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('CommandName',
		usage='%(prog)s -a arg_file -b arg_value [-l arg_option]',
		description='Description of this command.'
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--inputfile', '-a', metavar='filename', type=str, required=True, 
		help='input filename (required)')	
	parser.add_argument('--inputvalue', '-b', metavar='value', type=float, required=True, 
		help='input float value (required)')		
	parser.add_argument('--loption', '-l', type=optional_command,
		help='run optional command')
	return parser

def optional_command(arg_option):
	print("Optional command is executed with a input '%s'" % arg_option)

def main(args=None):
	"""
	Main entry point for your project.
	Args:
		args : list
		A of arguments as if they were input in the command line. 
		Leave it None to use sys.argv.
    """
	parser = get_parser()
	args = parser.parse_args(args)

	print("argument filename: %s" % args.inputfile)
	print("argument value: %s" % args.inputvalue)

if __name__=="__main__":
	main()