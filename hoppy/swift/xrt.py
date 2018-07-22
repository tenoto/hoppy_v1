# -*- coding: utf-8 -*-

import os
import sys  
import yaml 
import glob 
import argparse

class XRTObsID():
	"""
	Swift XRT ObsID path 
	"""
	def __init__(self,obsid_path,fyaml_param):
		self.obsid_path = obsid_path
		self.fyaml_param = fyaml_param
		self.obsid = os.path.basename(self.obsid_path)

	def check_input_files(self):
		if not os.path.exists(self.obsid_path):
			sys.stderr.write('[error] obsid_path does not exist: %s\n' % self.obsid_path)			
			quit()				
		if not os.path.exists(self.fyaml_param)	:
			sys.stderr.write('[error] fyaml_param does not exist: %s\n' % self.fyaml_param)			
			quit()				
		try:
			self.param = yaml.load(open(self.fyaml_param))
		except:
			sys.stderr.write('[error] fyaml_param can not be loaded.: %s\n' % self.fyaml_param)			
			quit()				
		self.outdir = self.param['outdir']
		self.outdir_obsid = '%s/%s' % (self.outdir,self.obsid)

		sys.stdout.write('===============================\n')
		sys.stdout.write('XRTObsID\n')		
		sys.stdout.write('obsid_path: %s\n' % self.obsid_path)				
		sys.stdout.write('outdir: %s\n' % self.outdir)										
		sys.stdout.write('fyaml_param: %s\n' % self.fyaml_param)						
		print(self.param)		
		sys.stdout.write('===============================\n')		

	def run_xrtpipeline(self):
		"""
		# ref. http://www.swift.ac.uk/XRT.shtml
		"""
		self.flog_xrtpipeline = "sw%s_xrtpipeline.log" % (self.obsid)
		self.fscript_xrtpipeline = "sw%s_xrtpipeline.sh" % (self.obsid)		
		cmd  = 'xrtpipeline '
		cmd += 'indir=\"%s\" ' % self.obsid_path
		cmd += 'outdir=%s ' % self.outdir_obsid
		cmd += 'steminputs=sw%s ' % self.obsid
		cmd += 'stemoutputs=DEFAULT '
		cmd += 'srcra=%.5f srcdec=%.5f ' % (float(self.param['SOURCE_RA']), float(self.param['SOURCE_DEC'])) 
		cmd += '>& %s ' % self.flog_xrtpipeline
		print(cmd)

		if not os.path.exists("%s/%s" % (self.outdir_obsid,self.flog_xrtpipeline)):
			f = open(self.fscript_xrtpipeline,'w')
			f.write('#!/bin/sh -f\n')
			f.write(cmd)
			f.close()

			cmd = 'chmod +x %s; ./%s' % (self.fscript_xrtpipeline,self.fscript_xrtpipeline)
			print(cmd);os.system(cmd)
			cmd = 'mv %s %s %s' % (self.flog_xrtpipeline,self.fscript_xrtpipeline,self.outdir_obsid)
			print(cmd);os.system(cmd)
		self.flog_xrtpipeline = "%s/%s" % (self.outdir_obsid,self.flog_xrtpipeline)
		self.fscript_xrtpipeline = "%s/%s" % (self.outdir_obsid,self.fscript_xrtpipeline)

		self.flag_xrtpipeline_status = False
		if not os.path.exists(self.flog_xrtpipeline):
			self.flag_xrtpipeline_status = False
		else:
			for line in open(self.flog_xrtpipeline).readlines():
				if line.find("Exit with no errors") == 20:
					self.flag_xrtpipeline_status = True
		print("flag_xrtpipeline_status:%s" % self.flag_xrtpipeline_status)

		if not self.flag_xrtpipeline_status: 
			sys.stderr.write('[error] xrtpipeline error.\n')
			exit()

	def get_rmffile(self,mode):
		rmffile = None
		for line in open(self.flog_xrtpipeline).readlines():
			if line.find('response') != -1:
				if line.find(mode) != -1:
					rmffile = line.split()[1]
		return rmffile

	def make_symboliclink(self,file):
		pwd = os.getcwd()
		os.chdir(self.outdir_obsid)
		os.system('ln -s %s .' % file)
		os.chdir(pwd)

	def check_XRTmode(self):
		self.xpc_mode = False 
		self.xwt_mode = False
		print(glob.glob('%s/sw%sxpcw?posr.pha' % (self.outdir_obsid,self.obsid)))
		print(glob.glob('%s/sw%sxpcw?po_cl.evt' % (self.outdir_obsid,self.obsid)))
		if len(glob.glob('%s/sw%sxpcw?posr.pha' % (self.outdir_obsid,self.obsid))) == 1:
			self.xpc_mode = True

			self.xpc_default_clevt = glob.glob('%s/sw%sxpcw?po_cl.evt' % (self.outdir_obsid,self.obsid))[0]
			self.xpc_default_pha   = glob.glob('%s/sw%sxpcw?posr.pha' % (self.outdir_obsid,self.obsid))[0]
			self.xpc_default_arf   = glob.glob('%s/sw%sxpcw?posr.arf' % (self.outdir_obsid,self.obsid))[0]	
			self.xpc_default_rmf   = self.get_rmffile(mode='xpc')	
			self.xpc_default_img   = glob.glob('%s/sw%sxpcw?po_sk.img' % (self.outdir_obsid,self.obsid))[0]
			#self.xpcXspecPha = XspecPha(self.xpc_default_pha)
			#self.xpcXspecPha.setProperties()
			self.make_symboliclink(self.xpc_default_rmf)

		if len(glob.glob('%s/sw%sxwtw?stsr.pha' % (self.outdir_obsid,self.obsid))) == 1:
			self.xwt_mode = True

			self.xwt_default_clevt = glob.glob('%s/sw%sxwtw?st_cl.evt' % (self.outdir_obsid,self.obsid))[0]
			self.xwt_default_pha   = glob.glob('%s/sw%sxwtw?stsr.pha' % (self.outdir_obsid,self.obsid))[0]		
			self.xwt_default_arf   = glob.glob('%s/sw%sxwtw?stsr.arf' % (self.outdir_obsid,self.obsid))[0]		
			self.xwt_default_rmf   = self.get_rmffile(mode='xwt')	
			#self.xwtXspecPha = XspecPha(self.xwt_default_pha)
			#self.xwtXspecPha.setProperties()
			self.make_symboliclink(self.xwt_default_rmf)

		print(vars(self))

	def make_ring_region(self, outerPix=120.0, innerPix=70.0, ringtype='src'):
		PIX2ARCSEC = 2.38 #   (2.38 arcsec / 1 pix)
		outerArcsec = PIX2ARCSEC * outerPix
		innerArcsec = PIX2ARCSEC * innerPix

		regfname = "%s/sw%s_%sring_%sto%spx.reg" % (self.outdir_obsid, self.obsid, 
			ringtype, str(innerPix).replace('.','p'), str(outerPix).replace('.','p'))

		reg =  '# Region file format: DS9 version 4.1 \n'
		reg += '# Filename: %s \n' % regfname
		reg += 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n'
		reg += 'fk5 \n'
		reg += 'circle(%.5f, %.5f, %.5f")\n' % (float(self.param['SOURCE_RA']), float(self.param['SOURCE_DEC']), float(outerArcsec))
		reg += '-circle(%.5f, %.5f, %.5f")' % (float(self.param['SOURCE_RA']), float(self.param['SOURCE_DEC']), float(innerArcsec))

		print reg
		fw = open(regfname, 'w')
		fw.write(reg)
		fw.close()
		return regfname

	def run_xrtproduct_region(self, inevtfile, regfile, outbase):
		logfname = regfile.split(".reg")[0] + "_xrtproduct.log"
		cmd  = 'xrtproducts '
		cmd += 'infile=%s ' % inevtfile		
		cmd += 'regionfile=%s ' % regfile 
		cmd += 'stemout=%s ' % outbase
		cmd += 'outdir=%s ' % self.outdir_obsid
		cmd += 'expofile=NONE correctlc=no '
		cmd += 'clobber=yes binsize=-99 > %s' % logfname
		print cmd 
		os.system(cmd)
		return logfname

	def make_bgdreg_ring(self, outerPix=120.0, innerPix=70.0):
		self.xpc_bgdreg = self.make_ring_region(outerPix, innerPix, ringtype='bgd')
		self.xpc_bgdreg_basename = self.xpc_bgdreg.split(".reg")[0] 
		self.run_xrtproduct_region(self.xpc_default_clevt, self.xpc_bgdreg, os.path.basename(self.xpc_bgdreg_basename))
		self.xpc_bgdring_pha = glob.glob('%s*.pha' % self.xpc_bgdreg_basename)[0]

	def check_pileup(self):
		img = glob.glob("%s/sw%sxpcw?po_sk.img" % (self.outdir_obsid, self.obsid))[0]
		evt = glob.glob("%s/sw%sxpcw?po_cl.evt" % (self.outdir_obsid, self.obsid))[0]
		center = commands.getoutput("fimgstat %s INDEF INDEF | grep \"The location of maximum\" | awk \'{print $10}\'" % img)
		x_center = float(center.split(",")[0].split("(")[1])
		y_center = float(center.split(",")[1].split(")")[0])
		print x_center, y_center
		offset = 25.0
		x_out = x_center - offset
		y_out = y_center - offset

		cmd  = 'xselect <<EOF \n' 
		cmd += 'xsel \n'
		cmd += 'read event %s ./ \n' % evt
		cmd += 'yes \n'
		cmd += 'extract image \n'
		cmd += 'save image %s \n' % "psfchk.img"
		cmd += '\$XIMAGE \n'
		cmd += 'read %s \n' % "psfchk.img"
		cmd += 'cpd /xtk \n'
		cmd += 'disp \n'
		cmd += 'back \n'
		cmd += 'psf \n'
		cmd += '%.1f \n' % x_center
		cmd += '%.1f \n' % y_center
		cmd += '%.1f \n' % x_out
		cmd += '%.1f \n' % y_out
		cmd += 'col off 1 2 3 4 6 \n'
		cmd += 'rescale x 15 \n'
		cmd += 'model king \n'
		cmd += '5.8, -1 \n'
		cmd += '1.55, -1 \n'
		cmd += '1 \n'
		cmd += 'fit \n'
		cmd += 'rescale x \n'
		cmd += 'rescale y 0.02 20 \n'
		cmd += 'fit plot 100 \n'
		cmd += 'col 2 on 5 \n'
		cmd += 'lwid 5 \n'
		cmd += 'hard psfchk.ps/cps \n'
		cmd += 'exit \n'
		cmd += 'exit \n'
		cmd += 'exit \n'
		cmd += 'no \n'
		cmd += 'EOF \n'
		print cmd
		os.system(cmd)
		if os.path.isfile("psfchk.ps"):
			os.system("ps2pdf psfchk.ps")

	def run(self):
		self.check_input_files()
		self.run_xrtpipeline()
		self.check_XRTmode()
		#if self.xpc_mode:
		#	self.make_bgdreg_ring()

if __name__=="__main__":

	sys.stdout.write('... run a single ObsID process (xrtpipeline) ...\n')

	parser = argparse.ArgumentParser(
		prog='xrt.py',
		usage='python xrt.py obsid_path parameter.yaml',
		description='A series of interfaces for Swift XRT pipelines.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'obsid_path',metavar='obsid_path',type=str,
		help='input ObsID path') 
	parser.add_argument(
		'fyaml_param',metavar='fyaml_param',type=str,
		help='input yamlfile for process parameters') 		
	args = parser.parse_args()	
	
	xrtobsid = XRTObsID(args.obsid_path,args.fyaml_param)
	xrtobsid.run()

