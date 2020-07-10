#!/usr/bin/env python

import sys 
import argparse
import astropy.io.fits as pyfits

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='group_gti.py',
		usage='group_gti.py input.fits outgti.txt',
		description='Grouping GTI to output gti text file',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'input_fits',metavar='input_fits',type=str,        
		help='Input fits file.')
	parser.add_argument(
		'output_gtitxt',metavar='output_gtitxt',type=str,        
		help='Output GTI text file.')	
	parser.add_argument(
		'--gti_separation_threshold',metavar='--gti_separation_threshold',type=float,        
		default=20.0,
		help='GTI separation threshold (sec) [default=20.0].')		
	args = parser.parse_args()	

	hdu = pyfits.open(args.input_fits)
	extname_list = []

	for i in range(len(hdu)):
		try:
			extname_list.append(hdu[i].header['EXTNAME'])
		except:
			pass
			#print('skip the extension...')

	if 'GTI' in extname_list:
		gtiname = 'GTI'
	elif 'STDGTI' in extname_list:
		gtiname = 'STDGTI'
	else:
		sys.stderr.write('input file does not have GTI extension.\n')

	modified_gti_list = []
	num_of_gti = len(hdu[gtiname].data)
	print("# GTI-num: START STOP Exposure(s) TDIFF New-GTI-num")
	for i in range(num_of_gti):
		gti_start = float(hdu[gtiname].data[i]['START'])
		gti_stop  = float(hdu[gtiname].data[i]['STOP'])
		gti_exp   = gti_stop - gti_start

		if i == 0:
			grpgti_tstart = gti_start
			grpgti_tstop  = gti_stop
			print("%03d:  %.3f %.3f %.3f -- 0" % (i,gti_start,gti_stop,gti_exp))
			continue

		gti_tdiff = gti_start - float(hdu['GTI'].data[i-1]['STOP'])

		if i+1 == num_of_gti:
			if gti_tdiff <= args.gti_separation_threshold:
				grpgti_tstop = gti_stop 
				modified_gti_list.append([grpgti_tstart,grpgti_tstop])
				print("%03d:  %.3f %.3f %.3f %.3f %d" % (i,gti_start,gti_stop,gti_exp,gti_tdiff,len(modified_gti_list)))			
			else:
				modified_gti_list.append([grpgti_tstart,grpgti_tstop])
				print("%03d:  %.3f %.3f %.3f %.3f %d" % (i,gti_start,gti_stop,gti_exp,gti_tdiff,len(modified_gti_list)))			
				grpgti_tstart = gti_start
				grpgti_tstop  = gti_stop
				modified_gti_list.append([grpgti_tstart,grpgti_tstop])		
			break 

		if gti_tdiff <= args.gti_separation_threshold:
			grpgti_tstop = gti_stop
			print("%03d:  %.3f %.3f %.3f %.3f %d" % (i,gti_start,gti_stop,gti_exp,gti_tdiff,len(modified_gti_list)))
			continue
			
		modified_gti_list.append([grpgti_tstart,grpgti_tstop])
		print("%03d:  %.3f %.3f %.3f %.3f %d" % (i,gti_start,gti_stop,gti_exp,gti_tdiff,len(modified_gti_list)))
		grpgti_tstart = gti_start		
		grpgti_tstop = gti_stop

	f = open(args.output_gtitxt,'w')
	f.write('# input evt: %s\n' % args.input_fits)
	f.write('# gti_separation_threshold: %.3f (s)\n' % args.gti_separation_threshold)
	for gti in modified_gti_list:
		print(gti)
		start = gti[0]
		stop = gti[1]
		f.write('%.6f %.6f\n' % (start,stop))
	f.close()		



