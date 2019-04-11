#!/usr/bin/env python

__author__ = 'Teru Enoto'
__version__ = '1.00'
# v1.00 : original version (2019-04-11)

import os 
import click
import astropy.io.fits as pyfits
 
@click.command(help='filter event fitsfile by the energy column (emin_keV <= energy < emax_keV)') # (1)
@click.argument('inevt',type=click.Path(exists=True))
@click.argument('outevt')
@click.argument('emin_kev',type=float)
@click.argument('emax_kev',type=float)
@click.option('--chatter',type=int,default=1,help='mount of verbosity of the task. For chatter=0 no output is printed. For chatter=2, debugging output is printed.')
def fselect_energy(inevt, outevt, emin_kev, emax_kev, chatter):
	msg  = 'input event file (inevt) : {inevt}\n'.format(inevt=inevt)
	msg += 'output event file (outevt): {outevt}\n'.format(outevt=outevt)
	msg += 'minimum energy in keV (emin_kev): {emin_kev}\n'.format(emin_kev=emin_kev)	
	msg += 'maximum energy in keV (emax_kev): {emax_kev}'.format(emax_kev=emax_kev)		
	if chatter >= 2:
		click.echo(msg)

	if os.path.exists(outevt):
		click.echo('Error: output event file has already existed ({outevt})'.format(outevt=outevt))
		quit()	

	hdu = pyfits.open(inevt)
	num_of_evt_in = len(hdu['EVENTS'].data)
	if num_of_evt_in == 0:
		click.echo('Skip: No events in fits file {}'.format(inevt))
		quit()

	kev_to_pi = None		
	telescope  = hdu[0].header['TELESCOP']
	instrument = hdu[0].header['INSTRUME']
	if telescope == 'NICER' and instrument == 'XTI':
		kev_to_pi = 1000.0/10.0 # Conversion from PI to keV (PI is in units of 10 eV)

	if chatter >= 2:		
		click.echo('telescope={telescope},instrument={instrument},kev_to_pi={kev_to_pi}'.format(telescope=telescope,instrument=instrument,kev_to_pi=kev_to_pi))
	pi_min = int(kev_to_pi * emin_kev)
	pi_max = int(kev_to_pi * emax_kev)
	cmd  = 'fselect %s %s ' % (inevt,outevt)
	cmd += '"(PI >= %d).and.(PI < %d)" ' % (pi_min,pi_max)
	if chatter >=1:
		print(cmd);
	os.system(cmd)

	if chatter >=1:
		hdu_out = pyfits.open(outevt)
		num_of_evt_out = len(hdu_out['EVENTS'].data)
		filtered_rate = float(num_of_evt_out)/float(num_of_evt_in)
		click.echo('number of filered events: {} ({:.3f}%% out of {})'.format(num_of_evt_out,100.0*filtered_rate,num_of_evt_in))

if __name__ == '__main__':
    fselect_energy()