#!/usr/bin/env python

import os 
import argparse
import astropy.io.fits as fits

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='extract_xspec_rate.py',
		usage='extract_xspec_rate.py pha rmf arf emin emax',
		description='Extract count rate from a xspec file.',
		epilog='',
		add_help=True,
		)

	parser.add_argument('pha',metavar='pha',type=str, help='xspec pha file.')
	parser.add_argument('rmf',metavar='rmf',type=str, help='xspec rmf file.')	
	parser.add_argument('arf',metavar='arf',type=str, help='xspec arf file.')		
	parser.add_argument('emin',metavar='emin',type=float, help='energy min (keV).')			
	parser.add_argument('emax',metavar='emax',type=float, help='energy max (keV).')				
	parser.add_argument('--keyword',metavar='keyword',type=str,default='RT01',help='energy max (keV).')	

	args = parser.parse_args()	
	#print(args.emin)			
	#print(args.emax)				

	cmd  = 'rm -f tmp_xspec.log tmp_rate.txt;\n'
	cmd += "xspec <<EOF\n"
	cmd += "data 1 %s\n" % args.pha
	cmd += "resp 1 %s\n" % args.rmf
	cmd += "arf 1 %s\n" % args.arf
	cmd += "setplot energy\n"
	cmd += "ignore 1:**-%.1f,%.1f-**\n" % (args.emin,args.emax)
	cmd += 'pl ld\n'
	cmd += 'log tmp_xspec.log\n'
	cmd += 'show rate\n'
	cmd += 'log none\n'
	cmd += 'exit\n'
	cmd += 'EOF\n'
	#print(cmd);
	os.system(cmd)

	f = open('tmp_rate.txt','w')
	cmd  = 'grep \'Net\' tmp_xspec.log | awk \'{print $7}\' > tmp_rate.txt;'
	cmd += 'grep \'Net\' tmp_xspec.log | awk \'{print $9}\' >> tmp_rate.txt;'	
	#print(cmd);
	os.system(cmd)

	with open('tmp_rate.txt', mode='rt', encoding='utf-8') as f:
		read_data = list(f)
	rate  = float(read_data[0])
	error = float(read_data[1])
	print("Count rate: %.6f" % rate)
	print("Rate error: %.6f" % error)
	cmd  = 'rm -f tmp_xspec.log tmp_rate.txt;\n'
	#print(cmd);
	os.system(cmd)

	cmd  = 'fparkey %s %s %s add=yes;\n' % (rate, args.pha, args.keyword)
	cmd += 'fparkey %s %s %s add=yes;\n' % (error, args.pha, ('%s_ERR' % args.keyword))
	cmd += 'fparkey %s %s %s add=yes;\n' % (args.emin, args.pha, ('%sEMIN' % args.keyword))	
	cmd += 'fparkey %s %s %s add=yes;\n' % (args.emax, args.pha, ('%sEMAX' % args.keyword))		
	print(cmd);os.system(cmd)

