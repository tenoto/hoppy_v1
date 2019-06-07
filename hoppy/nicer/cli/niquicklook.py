#!/usr/bin/env python

import os 
import sys 
import numpy as np
import astropy.io.fits as pyfits
from astropy.stats import LombScargle
import matplotlib.pyplot as plt 

plt.style.use('https://raw.githubusercontent.com/tenoto/macsettings/master/supermongo.mplstyle')
plt.ioff()

def wget_nicer_cleanevt(obsid,yyyy_mm):
	clevt_path = 'https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/%s/%s/xti/event_cl/ni%s_0mpu7_cl.evt.gz' % (yyyy_mm,obsid,obsid)
	cmd = 'wget -q -nH --no-check-certificate --cut-dirs=8 -r -l0 -c -N -np -R \'index*\' -erobots=off --retr-symlinks %s' % clevt_path
	filename = 'ni%s_0mpu7_cl.evt.gz' % obsid
	print(cmd);os.system(cmd)
	return filename

def wget_nicer_orbfile(obsid,yyyy_mm):
	orbfile_path = 'https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/%s/%s/auxil/ni%s.orb.gz' % (yyyy_mm,obsid,obsid)
	cmd = 'wget -q -nH --no-check-certificate --cut-dirs=8 -r -l0 -c -N -np -R \'index*\' -erobots=off --retr-symlinks %s' % orbfile_path
	filename = 'ni%s.orb.gz' % obsid
	print(cmd);os.system(cmd)
	return filename

def keV2pi(keV,instrume):
	if instrume == 'XTI' or instrume == 'xti':
		a = 100.0
		b = 0.0
	pi = int(a * keV + b)
	return pi
	
def pi2keV(pi,instrume):
	if instrume == 'XTI' or instrume == 'xti':
		a = 1/100.0
		b = 0.0
	keV = a * pi + b
	return keV

def get_lightcurve(times,tbin=1.):
	nbin = int(round(times[-1] - times[0])/float(tbin))
	ncounts, bins_edge = np.histogram(times-times[0],bins=nbin)
	times = 0.5*(bins_edge[1:] + bins_edge[:-1])	
	terr = 0.5*(bins_edge[1:] - bins_edge[:-1])	
	nerr = np.sqrt(ncounts)
	return times, ncounts, terr, nerr

def plot(x,y,xerr=None,yerr=None,xlabel='X',ylabel='Y',
	xlim=None,ylim=None,xlog=False,ylog=False,
	legend='',title='',pdfname='hoge.pdf'):
	plt.clf()
	plt.figure(figsize=(9,6))
	plt.errorbar(x,y,xerr=xerr,yerr=yerr,drawstyle='steps-mid',marker='',color='k',ecolor='k',rasterized=True)
	plt.xlabel(xlabel);plt.ylabel(ylabel);plt.title(title)
	if xlim != None: plt.xlim(xlim)
	if ylim != None: plt.ylim(ylim)	
	if xlog: plt.xscale('log')
	if ylog: plt.yscale('log')	
	plt.legend([legend])	
	plt.savefig(pdfname)
	return pdfname

class EVTFILE():
	def __init__(self,filename):
		try:
			self.hdu = pyfits.open(filename)
			print('-- %s is successfully loaded.' % filename)
			self.filename = filename
			self.set_property()
		except IndexError:
			print('-- Usage: EVTFILE(filename)')
		except IOError:
			print('-- Error: %s cannot be opened.' % filename)	 
			 
	def set_property(self):
		self.telescop = self.hdu[0].header['TELESCOP'] 
		self.obs_id = self.hdu[0].header['OBS_ID'] 
		self.date_obs = self.hdu[0].header['DATE-OBS'] 
		self.date_end = self.hdu[0].header['DATE-END'] 
		self.instrume = self.hdu[0].header['INSTRUME']
		self.object = self.hdu[0].header['OBJECT']
		self.ra_obj = self.hdu[0].header['RA_OBJ']
		self.dec_obj = self.hdu[0].header['DEC_OBJ']
		self.exposure = self.hdu['EVENTS'].header['EXPOSURE']

		self.nevents = len(self.hdu['EVENTS'].data)	 
		self.ngtis	 = len(self.hdu['GTI'].data)
		if self.exposure > 0.0:
			self.rate = float(self.nevents)/float(self.exposure)
		else:
			self.rate = 0.0

		self.title = '%s Obs#%s %s' % (self.object,self.obs_id,self.date_obs)
		self.pdf_list = []

	def show_property(self):
		print("FILENAME: %s" % self.filename)
		print("TELESCOP: %s" % self.telescop)
		print("INSTRUME: %s" % self.instrume)
		print("OBJECT:	 %s" % self.object)		
		print("OBD_ID:	 %s" % self.obs_id)		
		print("DATE-OBS: %s" % self.date_obs)
		print("DATE-END: %s" % self.date_end)
		print("EXPOSURE: %.1f (s)" % self.exposure)
		print("NUM_OF_EVENTS: %d" % self.nevents)
		print("NUM_OF_GTI: %d" % self.ngtis)
		print("RATE: %.2f (cps)" % self.rate)
		
	def get_curve_counts(self,lower_keV,upper_keV,time_bin):
		lower_pi = keV2pi(lower_keV,instrume=self.instrume)
		upper_pi = keV2pi(upper_keV,instrume=self.instrume)		
		filter_energy = np.logical_and(self.hdu['EVENTS'].data['PI'] >= lower_pi, self.hdu['EVENTS'].data['PI'] <= upper_pi)
		time_series = self.hdu['EVENTS'].data['TIME'][filter_energy]
		duration = time_series[-1]-time_series[0]
		nbin = int(float(duration) / float(time_bin))
		binned_counts, binned_edges, patches = plt.hist(time_series-time_series[0],bins=nbin,histtype='step')
		binned_centers = 0.5*(binned_edges[1:] + binned_edges[:-1])
		#plt.errorbar(binned_centers,binned_counts,yerr=np.sqrt(binned_counts),marker='.',drawstyle='steps-mid',color='k',linecolor='k')
		return binned_centers,binned_counts

	def plot_curve(self,lower_keV,upper_keV,time_bin,color='k'):
		binned_centers,binned_counts = self.get_curve_counts(lower_keV,upper_keV,time_bin)
		fig, ax = plt.subplots(figsize=(12,8))		
		plt.errorbar(binned_centers,binned_counts/time_bin,yerr=np.sqrt(binned_counts)/time_bin,
								 label='%.1f-%.1f keV' % (lower_keV,upper_keV),
								 marker='',drawstyle='steps-mid',color=color[i],ecolor=color[i])
		ax.set_xlim(binned_centers[0]-0.5*time_bin,binned_centers[-1]+0.5*time_bin)
		ax.set_xlabel('Time (sec)')
		ax.set_ylabel('Count rate (cps)')
		ax.legend(loc='upper left',scatterpoints=1,fontsize=18)
		pdfname = 'ni%s_%.1f_%.1fkeV_%dbin.pdf' % (self.obs_id,lower_keV,upper_keV,time_bin)
		plt.savefig(pdfname)		
		self.pdf_list.append(pdfname)
		
	def plot_curves(self,ebands,time_bin,color=['#FF0000','#E8A317','#488AC7','k','k']):
		curves = []
		for eband in ebands:
			curves.append(self.get_curve_counts(eband[0],eband[1],time_bin))
		print(curves)
		
		fig, ax = plt.subplots(figsize=(12,8))		
		for i in range(len(curves)):			
			plt.errorbar(curves[i][0],curves[i][1]/time_bin,
				yerr=np.sqrt(curves[i][1])/time_bin,
				label='%.1f-%.1f keV' % (ebands[i][0],ebands[i][1]),
				marker='',drawstyle='steps-mid',color=color[i],ecolor=color[i])			
		title = '%s (T=%.1fs, bin=%.1fs)' % (self.title,self.exposure,time_bin)
		ax.set_xlim(curves[0][0][0]-0.5*time_bin,curves[0][0][-1]+0.5*time_bin)
		ax.set_xlabel('Time (sec)')
		ax.set_ylabel('Count rate (cps)')
		ax.set_title(title,fontsize=18)
		ax.legend(loc='upper left',scatterpoints=1,fontsize=18)
		pdfname = 'ni%s_%dbands_%dbin.pdf' % (self.obs_id,len(ebands),time_bin)
		plt.savefig(pdfname)
		self.pdf_list.append(pdfname)
		
	def plot_curves_default(self):
		self.plot_curves([[0.3,3.0],[8.0,15.0],[0.2,0.3]],10.0)
		
	def get_spectrum(self,energy_keV,bins=2**6,lower_keV=0.2,upper_keV=15.0,
		exposure=1.0,flag_exposure=False,pdfname='spec.pdf'):
		logbins = np.logspace(np.log10(lower_keV),np.log10(upper_keV),bins)
		binned_counts, binned_edges, patches = plt.hist(energy_keV,bins=logbins,log=True)
		binned_centers = 0.5*(binned_edges[1:] + binned_edges[:-1])
		binned_widths = 0.5*(+binned_edges[1:] - binned_edges[:-1])

		if flag_exposure:
			sepctrum = binned_counts/(2.0*binned_widths)/exposure
			spectrum_error = np.sqrt(binned_counts)/(2.0*binned_widths)/exposure
		else:
			sepctrum = binned_counts/(2.0*binned_widths)
			spectrum_error = np.sqrt(binned_counts)/(2.0*binned_widths)			
		return binned_centers,sepctrum,binned_widths,spectrum_error

	def set_rebinned_gti(self,separation_threshold=100.0):
		bingti_start = []
		bingti_stop = []	
		bingti_exposure = []	
		exposure = 0.0 
		for i in range(self.ngtis):
			if i == 0: 
				new_start = self.hdu['GTI'].data[0]['START']
				exposure += self.hdu['GTI'].data[i]['STOP']-self.hdu['GTI'].data[i]['START']
				continue
			separation = float(self.hdu['GTI'].data[i]['START']-self.hdu['GTI'].data[i-1]['STOP'])
			if separation >= 100.0:
				bingti_start.append(new_start)
				bingti_stop.append(self.hdu['GTI'].data[i-1]['STOP'])
				bingti_exposure.append(exposure)
				new_start = self.hdu['GTI'].data[i]['START']
				exposure	= self.hdu['GTI'].data[i]['STOP']-self.hdu['GTI'].data[i]['START']				
			else:
				exposure += self.hdu['GTI'].data[i]['STOP']-self.hdu['GTI'].data[i]['START']				

			if i == self.ngtis-1:
				bingti_start.append(new_start)
				bingti_stop.append(self.hdu['GTI'].data[i]['STOP'])		
				bingti_exposure.append(exposure)				
		self.bingti_start = np.array(bingti_start)		
		self.bingti_stop = np.array(bingti_stop)				
		self.bingti_exposure = np.array(bingti_exposure)

	def plot_spectrum(self):
		x,y,xerr,yerr = self.get_spectrum(
			pi2keV(self.hdu['EVENTS'].data['PI'],instrume='XTI'),
			exposure=self.exposure,flag_exposure=True)
		pdfname = 'ni%s_avespec.pdf' % (self.obs_id)

		title = '%s (T=%.1fs)' % (self.title,self.exposure)		
		plot(x,y,xerr=xerr,yerr=yerr,xlabel='Energy (keV)',ylabel='Count/sec/keV',
			xlim=(0.2,15.0),title=title,
			xlog=True,ylog=True,legend='NICER XTI',pdfname=pdfname)
		self.pdf_list.append(pdfname)

	def plot_spectrum_gti(self):
		self.set_rebinned_gti(self)
		tzero = self.hdu['EVENTS'].data['TIME'][0] 

		pdfname = 'ni%s_gtispec.pdf' % (self.obs_id)
		self.pdf_list.append(pdfname)
		cmd	= 'gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite '
		cmd += '-sOutputFile=%s ' % pdfname
		for i in range(len(self.bingti_start)):
			tstart = self.bingti_start[i]
			tstop = self.bingti_stop[i]	
			texp = self.bingti_exposure[i]		

			mask_time = np.logical_and(self.hdu['EVENTS'].data['TIME'] >= tstart,self.hdu['EVENTS'].data['TIME'] < tstop)

			x,y,xerr,yerr = self.get_spectrum(
				pi2keV(self.hdu['EVENTS'].data['PI'][mask_time],instrume='XTI'),
				exposure=texp,flag_exposure=True)
			pdfname = 'tmp_gti%02d_spec.pdf' % i

			title = '%s GTI-%02d (%.1f...%.1f Texp=%.1f s)' % (self.title,i,tstart-tzero,tstop-tzero,texp)
			plot(x,y,xerr=xerr,yerr=yerr,xlabel='Energy (keV)',ylabel='Count/sec/keV',
				xlim=(0.2,15.0),ylim=(0.02,500.0),title=title,
				xlog=True,ylog=True,legend='NICER XTI',pdfname=pdfname)
			cmd += '%s ' % pdfname
		print(cmd);os.system(cmd)
		os.system('rm -f tmp_gti*_spec.pdf')

	def pdf_merge(self):
		pdfname = 'ni%s_ql.pdf' % (self.obs_id)		
		cmd	= 'gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite '
		cmd += '-sOutputFile=%s ' % pdfname
		for pdf in self.pdf_list:
			cmd += '%s ' % pdf 
		print(cmd);os.system(cmd)
		for pdf in self.pdf_list:
			os.system('rm -f %s' % pdf)

	def run_barycorr(self,orbfile,ra=None,dec=None,ephem=None):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		if ra == None: ra = self.ra_obj
		if dec == None: dec = self.dec_obj
		self.baryevt = "%s_bary.evt" % os.path.splitext(self.filename.replace('.gz',''))[0]
		self.barylog = "%s.log" % os.path.splitext(self.baryevt)[0]
		cmd	= 'barycorr ra=%.6f dec=%.6f ' % (ra,dec)
		cmd += 'infile=%s ' % self.filename
		cmd += 'outfile=%s ' % self.baryevt
		cmd += 'orbitfiles=%s ' % orbfile
		if ephem == None:
			cmd += 'refframe=ICRS	ephem=JPLEPH.430 '
		else:
			cmd += '%s ' % ephem 
		cmd += ">& %s " % self.barylog
		print(cmd); os.system(cmd)

	def plot_LombScargle(self,lower_keV=0.8,upper_keV=5.0,tbin=0.5):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		lower_pi = keV2pi(lower_keV,instrume=self.instrume)
		upper_pi = keV2pi(upper_keV,instrume=self.instrume)		

		self.hdu_bary =	pyfits.open(self.baryevt)
		evttimes = self.hdu_bary['EVENTS'].data['TIME']		
		filter_energy = np.logical_and(self.hdu_bary['EVENTS'].data['PI'] >= lower_pi, self.hdu_bary['EVENTS'].data['PI'] <= upper_pi)
		time_series = self.hdu_bary['EVENTS'].data['TIME'][filter_energy]

		lc_x,lc_y,lc_xerr,lc_yerr=get_lightcurve(time_series,tbin=tbin)
		title = '%s Tbin=%1.fs' % (self.title,tbin)
		pdfname = 'ni%s_lc_for_lombscargle.pdf' % obsid
		plot(lc_x,lc_y,lc_xerr,lc_yerr,xlabel='Time (sec)',ylabel='Count/bin',
			title=title,pdfname=pdfname)
		self.pdf_list.append(pdfname)

		frequency, power = LombScargle(lc_x,lc_y).autopower(method='fast')
		fmax=frequency[np.argmax(power)]
		title = '%s Tbin=%.1fs' % (self.title,tbin)
		pdfname = 'ni%s_lombscargle.pdf' % obsid
		plot(frequency,power,xlabel='Frequency (Hz)',ylabel='Power',
			xlog=True,ylog=True,title=title,pdfname=pdfname)		
		self.pdf_list.append(pdfname)

def run(obsid,yyyy_mm,flag_download=True,flag_pulsesearch=True):
	outdir = 'ql%s' % obsid
	cmd = 'rm -rf %s; mkdir -p %s' % (outdir,outdir)
	print(cmd);os.system(cmd)
	os.chdir(outdir)

	if flag_download:
		clevt_filename = wget_nicer_cleanevt(obsid,yyyy_mm)
	else:
		clevt_filename = 'ni%s_0mpu7_cl.evt.gz' % obsid
	clevt = EVTFILE(clevt_filename)
	clevt.show_property()
	clevt.plot_curves_default()
	clevt.plot_spectrum()
	#clevt.plot_spectrum_gti()
	if flag_pulsesearch:
		orbfile = wget_nicer_orbfile(obsid,yyyy_mm)
		clevt.run_barycorr(orbfile)
		clevt.plot_LombScargle()
	clevt.pdf_merge()			
	os.chdir('../')

if __name__ == '__main__':
	obsid = sys.argv[1]
	yyyy_mm = sys.argv[2]
	run(obsid,yyyy_mm)
	


