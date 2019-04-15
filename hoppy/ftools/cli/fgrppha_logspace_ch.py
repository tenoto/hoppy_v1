#!/usr/bin/env python

__author__ = 'Teru Enoto'
__version__ = '1.00'
# v1.00 : update to click (2019-04-15)

import os 
import click
import json 
import numpy
import astropy.io.fits as pyfits
 
@click.command(help='binning pha file in log space. ')
@click.argument('fname_inpha',type=click.Path(exists=True))
@click.argument('fname_outpha')
@click.option('--binning_ch', default=[], help='binning channel, e.g., --binning_ch="[[20,130,12],[131,300,3],[301,1000,10]]"')
@click.option('--chatter',type=int,default=1,help='mount of verbosity of the task. For chatter=0 no output is printed. For chatter=2, debugging output is printed.')

def fgrppha_logscpace_ch(fname_inpha,fname_outpha,binning_ch,chatter):
	if binning_ch == []:
		click.echo('--binning_ch option is needed. (e.g., --binning_ch="[[20,130,12],[131,300,3],[301,1000,10]]")')
		quit()
	try:
		binning_ch = json.loads(binning_ch)   
	except ValueError:
		click.echo('--binning_ch option is needed. (e.g., --binning_ch="[[20,130,12],[131,300,3],[301,1000,10]]")')
		quit()

	basename = os.path.splitext(os.path.basename(fname_outpha))[0]
	shfile   = '%s.sh' % basename
	f = open(shfile,'w')
	dump  = '#!/bin/sh -f\n\n'
	dump += 'grppha << EOF\n'
	dump += '%s\n' % fname_inpha
	dump += '%s\n' % fname_outpha
	for binset in binning_ch:
		pi_min, pi_max, nbin = binset
		tmp_float_pi = numpy.logspace(
			numpy.log10(pi_min-1), numpy.log10(pi_max), 
			nbin, base=10)
		for i in range(len(tmp_float_pi)-1):
			p0 = int(round(tmp_float_pi[i]))+1
			p1 = int(round(tmp_float_pi[i+1]))
			binsize = p1 - p0 + 1 
			dump += 'group %d %d %d\n' % (p0, p1, binsize)
	dump += 'exit\n'
	click.echo(dump)
	f.write(dump)
	f.close()

	cmd = 'chmod +x %s' % shfile
	os.system(cmd)
	os.system('./%s' % shfile)

if __name__ == '__main__':
    fgrppha_logscpace_ch()

