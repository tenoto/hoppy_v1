#!/usr/bin/env python

__author__ = 'Teru Enoto'
__version__ = '1.00'
# v1.00 : original version (2019-04-12)

import os 
import click
import pandas

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = '14'
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['xtick.top'] = 'True'
mpl.rcParams['ytick.right'] = 'True'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['axes.grid'] = 'True'
mpl.rcParams['axes.xmargin'] = '.05' #'.05'
mpl.rcParams['axes.ymargin'] = '.05'
mpl.rcParams['savefig.facecolor'] = 'None'
mpl.rcParams['savefig.edgecolor'] = 'None'
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['xtick.labelsize']=8
mpl.rcParams['ytick.labelsize']=8

"""
###---format of files named:  $SRC.dtype_etype_ref..phot

     These are photometry files, extractions from spectra processed in various ways

     dtype can be :  cms(clean,merged,slow-chain(raw))
     	       	     bg(background model file)
		     bgsub(BG-subtracted == cms - BG)
		     ufa(specialized)
     etype = extraction type can be:
     	     	     cl      cleaned events, 52 operating FPMs
     	     	     cl50    cleaned events, 50 FPMs, exclusing 14,34
     ref can be:  cms  (clean,merged,slow-chain measurements for operating FPMs)
     	     	  norm (normalized to 52 FPMs, not the number used

##--column   quantity			unit
1	filename (basis of extractions)	string(contains GTI seq. #)
2	exposure 	   		s
3	MET_day				days (midpoint of exposure

4	soft_1 (0.2-0.3 keV)		c/s 
5	soft_2 (0.3-0.4 keV)		c/s 
6	A_band (0.4-1.0 keV)		c/s 
7	B_band (1-2 keV)		c/s 
8	C_band (2-4 keV)		c/s 
9	D_band (4-12 keV)		c/s 
10	"in-band" (A+B+C+D = 0.4-12)	c/s
11	HBG (13-15 keV)			c/s
12	MJD (=MET_day + 56658.0)	day  ##--missing on older files

##########################

###---note that for bgsub type photometry files, HBG can be used to diagnose cases 
###---where there is a shape match between the true and predicted backgrund spectrum.
###---
###---accordingly, a file "$file.fphot" is a HBG-filtered version of "$file.phot"
###---    the filter condition is: awk '(($11 > -0.5) && ($11 < 0.5))' 
###---
###---Files named $file.ffphot are a more stringently filtered version of $file.phot
###---    the filter condition is: awk '(($11 > -0.1) && ($11 < 0.1))' 
###---This is a more appropriate filter for sources with in-band flux < 1000 c/s.
###---
###---HBG is beyond the effective area curve of the XTI optics, and so HBG should be 
###---near zero for BG-subtracted pha spectra, except for the "straight-through" 
###---hard photons...which should have very low count rates except for rare
###---bright-hard sources.....like the outurst maximum time range of swiftj0243+61.
###---
###---In summary, use $file.ffphot for all faint sources and be prepared to use 
###---$file.fphot, if col_10 > 1000 c/s.  Simply look at the correlation between 
###---col.10 and col.11 in the phot file to determine what is safe.....
###---fitering is intended to eliminate outliers, not to mask a portion of an outburst
###---
"""

MITBGD3C50_PHOTFILE_COLNAMES = ["filename","exposure","MET_day","soft_1","soft_2","A_band","B_band","C_band","D_band","in-band","HBG","MJD"]
 
@click.command(help='filter event fitsfile by the energy column (emin_keV <= energy < emax_keV)') # (1)
@click.argument('photfile',type=click.Path(exists=True))
@click.option('--chatter',type=int,default=1,help='mount of verbosity of the task. For chatter=0 no output is printed. For chatter=2, debugging output is printed.')
def plot_mitbgd3c50_photfile(photfile, chatter):
	click.echo('plotting {photfile}'.format(photfile=photfile))
	outpdf = '%s.pdf' % os.path.basename(os.path.splitext(photfile)[0])

	df = pandas.read_csv(photfile,delim_whitespace=True,names=MITBGD3C50_PHOTFILE_COLNAMES)
	print(df["MJD"])

	plt.clf()
	fig, axes = plt.subplots(8,1,figsize=(6,15),sharex=True)
	fig.subplots_adjust(hspace=0)
	axes[0].plot(df["MJD"],df["soft_1"],'o',
		ms=4,color='#000000',markerfacecolor='#ffffff',lw=1.5)
	xmin,xmax = axes[0].get_xlim()	
	axes[0].hlines(0.0,xmin,xmax,linestyle='dashed',linewidth=1.5)
	axes[0].set_xlim(xmin,xmax)
	axes[0].set_ylabel("soft_1 (0.2-0.3 keV) c/s",fontsize=10)

	axes[1].plot(df["MJD"],df["soft_2"],'o',
		ms=4,color='#000000',markerfacecolor='#ffffff',lw=1.5)
	xmin,xmax = axes[0].get_xlim()	
	axes[1].hlines(0.0,xmin,xmax,linestyle='dashed',linewidth=1.5)
	axes[1].set_xlim(xmin,xmax)
	axes[1].set_ylabel("soft_2 (0.3-0.4 keV) c/s",fontsize=10)

	axes[2].plot(df["MJD"],df["A_band"],'o',
		ms=4,color='#000000',markerfacecolor='#ffffff',lw=1.5)
	xmin,xmax = axes[0].get_xlim()	
	axes[2].hlines(0.0,xmin,xmax,linestyle='dashed',linewidth=1.5)
	axes[2].set_xlim(xmin,xmax)
	axes[2].set_ylabel("A_band (0.4-1.0 keV) c/s",fontsize=10)

	axes[3].plot(df["MJD"],df["B_band"],'o',
		ms=4,color='#000000',markerfacecolor='#ffffff',lw=1.5)
	xmin,xmax = axes[0].get_xlim()	
	axes[3].hlines(0.0,xmin,xmax,linestyle='dashed',linewidth=1.5)
	axes[3].set_xlim(xmin,xmax)
	axes[3].set_ylabel("B_band (1-2 keV) c/s",fontsize=10)

	axes[4].plot(df["MJD"],df["C_band"],'o',
		ms=4,color='#000000',markerfacecolor='#ffffff',lw=1.5)
	xmin,xmax = axes[0].get_xlim()	
	axes[4].hlines(0.0,xmin,xmax,linestyle='dashed',linewidth=1.5)
	axes[4].set_xlim(xmin,xmax)
	axes[4].set_ylabel("C_band (2-4 keV) c/s",fontsize=10)

	axes[5].plot(df["MJD"],df["D_band"],'o',
		ms=4,color='#000000',markerfacecolor='#ffffff',lw=1.5)
	xmin,xmax = axes[0].get_xlim()	
	axes[5].hlines(0.0,xmin,xmax,linestyle='dashed',linewidth=1.5)
	axes[5].set_xlim(xmin,xmax)
	axes[5].set_ylabel("D_band (4-12 keV) c/s",fontsize=10)

	axes[6].plot(df["MJD"],df["in-band"],'o',
		ms=4,color='#000000',markerfacecolor='#ffffff',lw=1.5)
	xmin,xmax = axes[0].get_xlim()	
	axes[6].hlines(0.0,xmin,xmax,linestyle='dashed',linewidth=1.5)
	axes[6].set_xlim(xmin,xmax)
	axes[6].set_ylabel("in-band (A+B+C+D) c/s",fontsize=10)

	axes[7].plot(df["MJD"],df["HBG"],'o',
		ms=4,color='#000000',markerfacecolor='#ffffff',lw=1.5)
	xmin,xmax = axes[0].get_xlim()	
	axes[7].hlines(0.0,xmin,xmax,linestyle='dashed',linewidth=1.5)
	axes[7].set_xlim(xmin,xmax)
	axes[7].set_ylabel("HBG (13-15 keV) c/s",fontsize=10)
	axes[7].set_xlabel("MJD",fontsize=10)

	fig.align_ylabels()
	plt.savefig(outpdf)


if __name__ == '__main__':
    plot_mitbgd3c50_photfile()


