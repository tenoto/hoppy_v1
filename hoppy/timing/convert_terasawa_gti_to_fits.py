#!/usr/bin/env python

__author__  = 'Teruaki Enoto'
__version__ = '0.01'
__date__    = '2018 September 28'

import os 
import sys 
import pandas as pd
import astropy.io.fits as fits
import argparse

class RadioGTIbyTerasawa():
	def __init__(self,textfile):
		self.textfile = textfile
		sys.stdout.write('input gti file: %s\n' % self.textfile)
		if not os.path.exists(self.textfile):
			sys.stderr.write('file %s does not exist.' % self.textfile)
			quit()

	def convert2fitsfile(self,output_gti=None):
		"""
		see http://docs.astropy.org/en/stable/io/fits/usage/table.html 
		original columns 
		nGTI     
		UTstart      
		UTend   
		sodStart     
		sodEnd         
		MJDobsStart           
		MJDobsEND        
		TDBstart          
		TDBend         
		MJDtdbStart           
		MJDtdbEND        
		duration           
		accum 
		"""
		if output_gti == None:
			output_gti = '%s.fits' % os.path.splitext(os.path.basename(self.textfile))[0]

		self.df = pd.read_csv(self.textfile,header=1,delim_whitespace=True)
		cols = []
		cols.append(fits.Column(name='nGTI',format='I',array=self.df['nGTI']))		
		cols.append(fits.Column(name='UTstart',format='D',array=self.df['UTstart'],unit='hhmmss.xxx'))
		cols.append(fits.Column(name='UTend',format='D',array=self.df['UTend'],unit='hhmmss.xxx'))
		cols.append(fits.Column(name='sodStart',format='D',array=self.df['sodStart'],unit='sec'))	
		cols.append(fits.Column(name='sodEnd',format='D',array=self.df['sodEnd'],unit='sec'))
		cols.append(fits.Column(name='MJDobsStart',format='D',array=self.df['MJDobsStart'],unit='day'))
		cols.append(fits.Column(name='MJDobsEND',format='D',array=self.df['MJDobsEND'],unit='day'))	
		cols.append(fits.Column(name='TDBstart',format='D',array=self.df['TDBstart'],unit='sec'))
		cols.append(fits.Column(name='TDBend',format='D',array=self.df['TDBend'],unit='sec'))
		cols.append(fits.Column(name='MJDtdbStart',format='D',array=self.df['MJDtdbStart'],unit='day'))				
		cols.append(fits.Column(name='MJDtdbEND',format='D',array=self.df['MJDtdbEND'],unit='day'))
		cols.append(fits.Column(name='duration',format='D',array=self.df['duration'],unit='sec'))	
		cols.append(fits.Column(name='accum',format='D',array=self.df['accum'],unit='sec'))
		cols.append(fits.Column(name='BARY_START',format='D',array=self.df['MJDtdbStart'],unit='day'))				
		cols.append(fits.Column(name='BARY_STOP',format='D',array=self.df['MJDtdbEND'],unit='day'))

		hdu_primary = fits.PrimaryHDU()
		hdu_table = fits.BinTableHDU.from_columns(cols,name='GTI')
		hdulist = fits.HDUList([hdu_primary,hdu_table])
		hdulist.writeto(output_gti,overwrite=True)


if __name__=="__main__":

	usage = """
This script convert a text-format GTI file provided by Terasawa into fits format.
"""

	parser = argparse.ArgumentParser(
		prog='convert_terasawa_gti_to_fits.py',
		usage='convert_terasawa_gti_to_fits.py terasawa_gti.txt',
		description=usage,epilog='',add_help=True)
	parser.add_argument('input_gti',metavar='input_gti',type=str,help='Input GTI text file.')
	parser.add_argument('--output_gti',metavar='output_gti',type=str,help='Output GTI fitsfile.',default=None)
	args = parser.parse_args()	

	gti = RadioGTIbyTerasawa(args.input_gti)
	gti.convert2fitsfile(args.output_gti)

