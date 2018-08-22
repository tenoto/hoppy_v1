#!/usr/bin/env python
# add pha files 

import sys
import os 
import math
import argparse
import astropy.io.fits as pyfits 

__author__  = "Teru Enoto"
__version__ = "1.0.0"
__email__   = "teruaki.enoto@gmail.com"
__status__  = "test"

#if len(sys.argv) != 2:
#	sys.stderr.write("usage: %s input_pha.lst\n" % sys.argv[0])
#	quit()

help_message = """
(example) %s input_pha.lst 
""" % sys.argv[0]

parser = argparse.ArgumentParser(description=help_message)
parser.add_argument('fname_inlst', metavar='fname_inlst',type=str, help='Input pha file list.')
parser.add_argument('--addtype', metavar='addtype',type=str, default='C',help='addition type C or R.')
args = parser.parse_args()
print(args)

#fname_inlst = sys.argv[1]
fname_inlst = args.fname_inlst
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
		if args.addtype == 'C':
			f.write("'tmp_%s' + " % os.path.basename(inpha_lst[i]))
		elif args.addtype == 'R':
			f.write("'tmp_%s' * %.6f + " % (os.path.basename(inpha_lst[i]), exposure_list[i]/sum_of_exposure))
	else:
		if args.addtype == 'C':		
			f.write("'tmp_%s'" % os.path.basename(inpha_lst[i]))			
		elif args.addtype == 'R':
			f.write("'tmp_%s' * %.6f" % (os.path.basename(inpha_lst[i]), exposure_list[i]/sum_of_exposure))
f.close()	

fname_outpha = '%s_merge.pha' % os.path.splitext(sys.argv[1])[0]

if args.addtype == 'C':
	cmd  = 'mathpha <<EOF\n'
	cmd += '@%s\n' % fname_expr
	cmd += 'C\n'
	cmd += '%s\n' % fname_outpha
	cmd += 'CALC\n'
	cmd += 'NULL\n'
	cmd += '0\n'
	cmd += 'EOF\n'
	print(cmd);os.system(cmd)
elif args.addtype == 'R':
	cmd  = 'mathpha '
	cmd += 'expr=@%s ' % fname_expr
	cmd += 'units=R outfil="%s" ' % fname_outpha
	cmd += 'exposure=%.3f ' % ave_of_exposure
	cmd += 'errmeth=gaussian properr=yes ncomments=0 areascal=NULL clobber=yes'
	print(cmd);os.system(cmd)

cmd = 'rm -f %s' % fname_expr
print(cmd);os.system(cmd)
for i in range(len(inpha_lst)):
	cmd = 'rm -f tmp_%s' % (os.path.basename(inpha_lst[i]))
	print(cmd);os.system(cmd)

