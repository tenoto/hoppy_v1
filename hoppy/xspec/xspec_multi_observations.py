# -*- coding: utf-8 -*-


if __name__=="__main__":

	sys.stdout.write('\n... run multiple xspec fitting ...\n')

	parser = argparse.ArgumentParser(
		prog=__file__,
		usage='python %s csvfile yamlfile',
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
