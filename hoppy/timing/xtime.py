#!/usr/bin/env python

__author__  = 'Teruaki Enoto'
__version__ = '0.01'
__date__    = '2018 September 14'

import sys 
import argparse
from astropy.time import Time 

FERMI_MET_ORIGIN   = Time('2001-01-01T00:00:00.000',format='isot',scale='utc')
NUSTAR_MET_ORIGIN  = Time('2010-01-01T00:00:00.000',format='isot',scale='utc')
RXTE_MET_ORIGIN    = Time('1994-01-01T00:00:00.000',format='isot',scale='utc')
SUZAKU_MET_ORIGIN  = Time('2000-01-01T00:00:00.000',format='isot',scale='utc')
SWIFT_MET_ORIGIN   = Time('2001-01-01T00:00:00.000',format='isot',scale='utc')
NICER_MET_ORIGIN   = Time('2014-01-01T00:00:00.000',format='isot',scale='utc')

class XTime():
	def __init__(self,time_value,format,scale,flag_display=False):
		self.time_value  = time_value
		self.format = format
		self.scale  = scale
		self.flag_display = flag_display

		if self.time_value.isdigit():
			self.time = Time(float(self.time_value),format=self.format,scale=self.scale)
		else:
			self.time = Time(str(self.time_value),format=self.format,scale=self.scale)

		if self.flag_display:
			dump = '%s (format=%s scale=%s)\n' % (
				self.time_value,self.format,self.scale)		
			sys.stdout.write(dump)

	def get_tt(self,format='isot'):
		if format == 'isot':
			return self.time.tt.isot
		elif format == 'jd':
			return self.time.tt.jd	
		elif format == 'mjd':
			return self.time.tt.mjd		
		elif format == 'yday':
			return self.time.tt.yday
		else:
			sys.stderr.write('error: no corresponding format [isot, jd, mjd, yday].')

	def get_utc(self,format='isot'):
		if format == 'isot':
			return self.time.utc.isot
		elif format == 'jd':
			return self.time.utc.jd	
		elif format == 'mjd':
			return self.time.utc.mjd		
		elif format == 'yday':
			return self.time.utc.yday
		else:
			sys.stderr.write('error: no corresponding format [isot, jd, mjd, yday].')	

	def get_met_chandra(self):
		return self.time.cxcsec

	def get_met_fermi(self):
		return self.time.gps - FERMI_MET_ORIGIN.gps

	def get_met_nustar(self):
		return self.time.gps - NUSTAR_MET_ORIGIN.gps

	def get_met_rxte(self):
		return self.time.gps - RXTE_MET_ORIGIN.gps

	def get_met_suzaku(self):
		return self.time.gps - SUZAKU_MET_ORIGIN.gps

	def get_met_swift(self):
		return self.time.gps - SWIFT_MET_ORIGIN.gps

	def get_met_nicer(self):
		return self.time.gps - NICER_MET_ORIGIN.gps		

if __name__=="__main__":

	usage = """
script 
"""

	parser = argparse.ArgumentParser(
		prog='xtime.py',
		usage='xtime.py input_time input_format input_scale',
		description=usage,epilog='',add_help=True)
	parser.add_argument('input_time',metavar='input_time',type=str,help='Input time value')
	parser.add_argument('input_format',metavar='input_format',type=str,
		help='Input time format : isot, jd, mjd, yday, met')	
	parser.add_argument('input_scale',metavar='input_scale',type=str,
		help='Input time scale : utc, tt, met_[nicer]')		
	args = parser.parse_args()	

	print(args)
	#xtime = XTime("2014-01-01T00:00:00","isot","utc")
	xtime = XTime(args.input_time,args.input_format,args.input_scale)

	dump  = "----- Calendar Time Formats -----\n"
	dump += "ISO8601_TT : %s (TT)\n" % xtime.get_tt(format='isot')
	dump += " JD_TT     : %.8f (TT) \n" % xtime.get_tt(format='jd')
	dump += "MJD_TT     :   %.8f (TT)\n" % xtime.get_tt(format='mjd')
	dump += "DOY_TT     :   %s (TT)\n" % xtime.get_tt(format='yday')
	dump += "\n"
	dump += "ISO8601_UTC: %s (UTC)\n" % xtime.get_utc(format='isot')
	dump += " JD_UTC    : %.8f (UTC) \n" % xtime.get_utc(format='jd')
	dump += "MJD_UTC    :   %.8f (UTC) \n" % xtime.get_utc(format='mjd')
	dump += "DOY_UTC   :   %s (UTC)\n" % xtime.get_utc(format='yday')
	dump += "\n"
	dump += "----- Mission-Specific Time Formats (Misson Elapsed Time, NET) -----\n"
	dump += "Fermi seconds sicne 2001.0 UTC (decimal)     : %.6f\n" % xtime.get_met_fermi()
	dump += "NuSTAR seconds since 2010.0 UTC (decimal)    : %.6f\n" % xtime.get_met_nustar()
	dump += "RXTE seconds since 1994.0 UTC (decimal)      : %.8f\n" % xtime.get_met_rxte()
	dump += "Suzaku seconds since 2000.0 UTC (decimal)    : %.6f\n" % xtime.get_met_suzaku()
	dump += "Swift seconds since 2001.0 UTC (decimal): %.8f\n" % xtime.get_met_swift()
	dump += "XMM/Chandra seconds since 1998.0 TT (decimal) : %.6f\n" % xtime.get_met_chandra()
	dump += "NICER seconds since 2014.0 UTC (decimal)     : %.6f\n" % xtime.get_met_nicer() 
	#if options.flag_time_correction:
	#	dump += "Caution : TIMEZERO correction is included.\n"
	print(dump)
