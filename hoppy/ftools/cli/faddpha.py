#!/usr/bin/env python

__author__ = 'Teru Enoto'
__version__ = '2.00'
# v2.00 : update to click (2019-04-11)
# v1.00 : original version 

import os 
import click
import astropy.io.fits as pyfits
 
@click.command(help='Add pha files based on the list.')
@click.argument('fname_inlst',type=click.Path(exists=True))
@click.argument('fname_outpha')
@click.option('--addtype',type=str,default='C',help='addition type C or R.')
@click.option('--chatter',type=int,default=1,help='mount of verbosity of the task. For chatter=0 no output is printed. For chatter=2, debugging output is printed.')
def faddpha(fname_inlst,fname_outpha,addtype,chatter):
	inpha_lst = []
	for line in open(fname_inlst):
		inpha_lst.append(line.split()[0])

	fname_expr = 'tmp_expr.txt'
	cmd = 'rm -f %s' % fname_expr
	print(cmd);os.system(cmd)

	exposure_list = []
	for inpha in inpha_lst:
		hdu = pyfits.open(inpha)
		exposure_list.append(float(hdu['SPECTRUM'].header['EXPOSURE']))
	sum_of_exposure = sum(exposure_list)
	ave_of_exposure = 0.0
	for i in range(len(exposure_list)):
		ave_of_exposure += exposure_list[i]/sum_of_exposure * exposure_list[i]

	f = open(fname_expr,'w')
	for i in range(len(inpha_lst)):
		cmd = 'ln -s %s tmp_%s' % (inpha_lst[i],os.path.basename(inpha_lst[i]))
		print(cmd);os.system(cmd)
		if i < len(inpha_lst)-1:
			if addtype == 'C':
				f.write("'tmp_%s' + " % os.path.basename(inpha_lst[i]))
			elif addtype == 'R':
				f.write("'tmp_%s' * %.6f + " % (os.path.basename(inpha_lst[i]), exposure_list[i]/sum_of_exposure))
		else:
			if addtype == 'C':		
				f.write("'tmp_%s'" % os.path.basename(inpha_lst[i]))			
			elif addtype == 'R':
				f.write("'tmp_%s' * %.6f" % (os.path.basename(inpha_lst[i]), exposure_list[i]/sum_of_exposure))
	f.close()	

	print(fname_expr)
	if addtype == 'C':
		cmd  = 'mathpha <<EOF\n'
		cmd += '@%s\n' % fname_expr
		cmd += 'C\n'
		cmd += '%s\n' % fname_outpha
		cmd += 'CALC\n'
		cmd += 'NULL\n'
		cmd += '0\n'
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)
	elif addtype == 'R':
		cmd  = 'mathpha '
		cmd += 'expr=@%s ' % fname_expr
		cmd += 'units=R outfil="%s" ' % fname_outpha
		cmd += 'exposure=NULL '
		#cmd += 'exposure=%.3f ' % ave_of_exposure
		cmd += 'errmeth=gaussian properr=yes ncomments=0 areascal=NULL clobber=yes'
		print(cmd);os.system(cmd)

	cmd = 'rm -f %s' % fname_expr
	print(cmd);os.system(cmd)
	for i in range(len(inpha_lst)):
		cmd = 'rm -f tmp_%s' % (os.path.basename(inpha_lst[i]))
		print(cmd);os.system(cmd)

if __name__ == '__main__':
    faddpha()

