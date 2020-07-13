#!/usr/bin/env python

import os 
import argparse
import astropy.io.fits as pyfits

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='extract_gti_fits.py',
		usage='extract_gti_fits.py input.fits gtinum',
		description='Extracct GTI as a fits file (start from 1)',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'input_fits',metavar='input_fits',type=str,        
		help='Input fits file.')
	parser.add_argument(
		'gtinum',metavar='gtinum',type=int,        
		help='GTI number (int, start from 1).')	
	args = parser.parse_args()	

	outfits = '%s_gti%003d.evt' % (os.path.splitext(os.path.basename(args.input_fits))[0],args.gtinum)
	gtifits = '%s_gti%003d.gti' % (os.path.splitext(os.path.basename(args.input_fits))[0],args.gtinum)
	if os.path.exists(outfits) or os.path.exists(gtifits):
		print("Error: output files %s or %s have already existed." % (outfits,gtifits))
		quit()

	hdu = pyfits.open(args.input_fits)
	extname_list = []

	for i in range(len(hdu)):
		try:
			extname_list.append(hdu[i].header['EXTNAME'])
		except:
			pass
			#print('skip the extension...')
	if "GTI" in extname_list:
		target_gti = "GTI"
	elif "STDGTI" in extname_list:
		target_gti = "STDGTI"

	maximum_gtinum = len(hdu[target_gti].data)
	if args.gtinum > maximum_gtinum:
		print("Input GTI number %d should be smaller than %d." % (args.gtinum,maximum_gtinum))
		quit()

	gti_start = hdu[target_gti].data[args.gtinum-1]['START']
	gti_stop  = hdu[target_gti].data[args.gtinum-1]['STOP']
	print(gti_start,gti_stop)

	cmd = 'rm -f tmp_gti.txt xselect.log'
	print(cmd);os.system(cmd)

	fname = "tmp_gti.txt"
	f = open(fname,'w')
	dump = "%.8f %.8f\n" % (gti_start,gti_stop)
	f.write(dump)
	f.close()

	cmd = 'fconv_txt2gti.py -i tmp_gti.txt -o %s \n' % gtifits
	print(cmd);os.system(cmd)

	cmd  = 'xselect<<EOF\n'
	cmd += 'xsel\n'
	cmd += 'read event %s ./\n' % args.input_fits
	cmd += 'yes\n'
	cmd += 'filter time file\n'
	cmd += '%s\n' % gtifits
	cmd += 'extract event\n'
	cmd += 'save event %s\n' % outfits
	cmd += 'yes\n'
	cmd += 'exit\n'
	cmd += 'no\n'
	cmd += 'exit\n'	
	cmd += 'EOF\n'
	print(cmd);os.system(cmd)	

	cmd = 'rm -f tmp_gti.txt xselect.log xsel_timefile.asc'
	print(cmd);os.system(cmd)
