#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import sys
import argparse
import astropy.io.fits as pyfits 

__author__  = 'Teruaki Enoto'
__version__ = '1.01'
__date__    = '2018 September 12'

def fplot_pulseprofile(inputfits,outfits=None,
		expr=None,nbin=20,colname="PULSE_PHASE",
		bgdrate=0.0,title=''):
	dump  = "inputfits: %s\n" % inputfits
	if outfits == None:
		if os.path.splitext(inputfits)[-1] == '.evt':
			outfits = os.path.basename(inputfits).replace('.evt','_pls.fits')
		elif os.path.splitext(inputfits)[-1] == '.fits':
			outfits = os.path.basename(inputfits).replace('.fits','_pls.fits')
	dump += "outfits: %s\n" % outfits
	dump += "expr: %s\n" % expr 
	dump += "nbin: %d\n" % nbin
	dump += "colname: %s\n" % colname
	dump += "bgdrate: %.4e\n" % bgdrate
	dump += 'title: %s\n' % title
	print(dump)

	if os.path.exists(outfits):
		sys.stderr.write('error: file %s has already existed.\n' % outfits)
		exit()
	outdir = os.path.dirname(outfits)
	if not os.path.exists(outdir):
		cmd = 'mkdir -p %s' % outdir
		print(cmd);os.system(cmd)

	cmd = 'rm -f tmp00_xselectmerge.evt;\n'
	if inputfits[0] == '@':
		sys.stdout.write('input : file list mode.\n')
		inputfits_list_file = inputfits[1:]
		if not os.path.exists(inputfits_list_file):
			sys.stderr.write('error: file %s does not exist.\n' % inputfits_list_file)
			exit()
		#cmd += 'ftmerge @%s tmp00_ftmerge.evt' % inputfits_list_file
		#print(cmd);os.system(cmd)
		cmd  = 'xselect <<EOF\n'
		cmd += 'xsel\n'
		cmd += 'read event %s .\n' % inputfits
		cmd += 'yes\n'
		cmd += 'no\n'
		cmd += 'extract event\n'
		cmd += 'save event\n'
		cmd += 'tmp00_xselectmerge.evt\n'
		cmd += 'no\n'
		cmd += 'exit\n'
		cmd += 'no\n'
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)
		cmd = 'rm -f xselect.log\n'
		print(cmd);os.system(cmd)
	else:
		if not os.path.exists(inputfits):
			sys.stderr.write('error: file %s does not exist.\n' % inputfits)
			exit()
		cmd += 'ln -s %s tmp00_xselectmerge.evt' % inputfits
		print(cmd);os.system(cmd)

	cmd = 'rm -f tmp01_fselect.evt;\n'
	if expr != None:
		cmd += 'fselect tmp00_xselectmerge.evt tmp01_fselect.evt "%s"' % (expr)
	else:
		cmd += 'ln -s tmp00_xselectmerge.evt tmp01_fselect.evt'
	print(cmd);os.system(cmd)


	binwidth = 1.0/float(nbin)
	cmd  = 'rm -f tmp02_fhisto.evt;\n'
	cmd += 'fhisto tmp01_fselect.evt tmp02_fhisto.evt '
	cmd += '%s %.5f ' % (colname, binwidth)
	cmd += 'lowval=0.0 highval=1.0 outcolx=%s outcoly=COUNTS outcolz=ERROR extname=PROFILE' % colname
	print(cmd);os.system(cmd)

	cmd  = 'rm -f tmp03_fhisto_copy.evt;\n'
	cmd += 'fcalc clobber=yes tmp02_fhisto.evt tmp03_fhisto_copy.evt '
	cmd += '%s %s+1 ' % (colname,colname)
	print(cmd);os.system(cmd)

	cmd  = 'rm -f tmp04_ftmerge.evt;\n'
	cmd += 'ftmerge tmp02_fhisto.evt,tmp03_fhisto_copy.evt tmp04_ftmerge.evt'
	print(cmd);os.system(cmd)

	hdu = pyfits.open('tmp00_xselectmerge.evt')
	exposure = float(hdu[1].header['EXPOSURE'])
	cmd  = 'fparkey %.5f tmp04_ftmerge.evt[1] EXPOSURE add=yes;\n' % exposure
	cmd += 'fcalc tmp04_ftmerge.evt %s RATE "COUNTS/(#EXPOSURE/((float) %d ))"\n' % (outfits,nbin)
	cmd += 'fcalc %s %s RATE_ERR "sqrt(COUNTS)/(#EXPOSURE/((float) %d ))" clobber=YES;' % (outfits,outfits,nbin)
	print(cmd);os.system(cmd)

	hdu = pyfits.open(outfits)
	rate_list = hdu['PROFILE'].data['RATE'][0:nbin]
	max_rate = max(rate_list)
	min_rate = min(rate_list)
	bgdsub_max_rate = max_rate - bgdrate
	bgdsub_min_rate = min_rate - bgdrate
	pulsed_fraction = (bgdsub_max_rate-bgdsub_min_rate)/(bgdsub_max_rate+bgdsub_min_rate)
	cmd  = 'fparkey %.5f %s[1] MAX_WBGD add=yes comm="Maximum rate (cps) before background subtraction.";\n' % (max_rate,outfits)
	cmd += 'fparkey %.5f %s[1] MIN_WBGD add=yes comm="Minimum rate (cps) before background subtraction.";\n' % (min_rate,outfits)
	cmd += 'fparkey %.5f %s[1] MAX_RATE add=yes comm="Background-subtracted maximum rate (cps).";\n' % (bgdsub_max_rate,outfits)
	cmd += 'fparkey %.5f %s[1] MIN_RATE add=yes comm="Background-subtracted minimum rate (cps).";\n' % (bgdsub_min_rate,outfits)
	cmd += 'fparkey %.5f %s[1] PULFRACT add=yes comm="Pulsed fraction after the background subtraction.";\n' % (pulsed_fraction,outfits)
	cmd += 'fparkey %d %s[1] NBINPHAS add=yes comm="Number of phase bins";\n' % (nbin,outfits)
	cmd += 'fparkey %.5e %s[1] BACKRATE add=yes comm="Subtracted background rate";\n' % (bgdrate,outfits)
	print(cmd);os.system(cmd)
	
	subtitle = 'BGD=%.2e (cps) Max=%.2e Min=%.2e Pulsed_Fraction=%.1f%%' % (bgdrate,
		bgdsub_max_rate,bgdsub_min_rate,pulsed_fraction*100.0)
	psfile = '%s.ps' % os.path.splitext(outfits)[0]
	cmd  = 'fplot %s %s "RATE[RATE_ERR]" - /xw @ ' % (outfits,colname)
	cmd += '<<EOF\n'
	cmd += 'line step on\n'
	cmd += 'time off\n'
	cmd += 'lwid 5 on\n'
	cmd += 'lwid 5\n'
	cmd += 'la t %s\n' % title
	cmd += 'la f %s\n' % subtitle 
	cmd += 'we %s\n' % os.path.splitext(os.path.basename(psfile))[0]
	cmd += 'hard %s/cps\n' % psfile
	cmd += 'quit\n'
	cmd += 'EOF\n'
	print(cmd);os.system(cmd)

	cmd = 'mv %s.{qdp,pco} %s' % (os.path.splitext(os.path.basename(psfile))[0],outdir)
	print(cmd);os.system(cmd)

	pdffile = '%s.pdf' % os.path.splitext(os.path.basename(psfile))[0]
	cmd  = 'ps2pdf %s ;' % (psfile)
	if outdir == '':
		cmd += 'mv %s . ' % pdffile
	else:	
		cmd += 'mv %s %s ' % (pdffile,outdir)
	print(cmd);os.system(cmd)

	for tmpfile in ['tmp00_xselectmerge.evt','tmp01_fselect.evt','tmp02_fhisto.evt','tmp03_fhisto_copy.evt','tmp04_ftmerge.evt',psfile]:
		cmd = 'rm -f %s' % tmpfile
		print(cmd);os.system(cmd	)

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='fplot_pulseprofile.py',
		usage='fplot_pulseprofile.py input.fits',
		description='Plot a pulse profile from fits file.',
		epilog='',
		add_help=True)

	parser.add_argument(
		'inputfits',metavar='inputfits',type=str,
		help='Input fitsfile') 	
	parser.add_argument(
		'--outfits',metavar='outfits',type=str, default=None,
		help='Output fitsfile name') 				
	parser.add_argument(
		'--expr',metavar='expr',type=str, default=None,
		help='Filter expression for fselect.') 		
	parser.add_argument(
		'--nbin',metavar='nbin',type=int, default=20,
		help='Folding number of bins.') 
	parser.add_argument(
		'--colname',metavar='colname',type=str,default='PULSE_PHASE',
		help='Phase column name of fits file.') 
	parser.add_argument(
		'--bgdrate',metavar='bgdrate',type=float,default=0.0,
		help='Background rate (cps).') 	
	parser.add_argument(
		'--title',metavar='title',type=str,default='Pulse Profile',
		help='Title for output figure.') 								
	args = parser.parse_args()	
	
	fplot_pulseprofile(args.inputfits,outfits=args.outfits,
		expr=args.expr,nbin=args.nbin,colname=args.colname,
		bgdrate=args.bgdrate,title=args.title)

