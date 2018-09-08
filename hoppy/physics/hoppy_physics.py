def get_luminosity(flux,d_kpc):
	luminosity = 1.2e+32 * (flux / 1e-12) * d_kpc**2 # erg/s
	return luminosity

def get_spindown_luminosity(period, pdot, moment_of_inertia=1.0):
	"""
	moment of inertia normalized to I = 1e+45 g cm2
	"""
	spindown_luminosity = 3.94e+35 * moment_of_inertia / period**3 * (pdot/1e-11)
	return spindown_luminosity

def get_norm_sf(mean,std,value,flag_plot=True, plot_nstd=5):
	import scipy.stats
	import numpy as np
	import matplotlib.pyplot as plt 

	if flag_plot:
		xmin = min(mean - plot_nstd * std, 0.9*value)
		xmax = max(mean + plot_nstd * std, 1.1*value)

		fig, axes = plt.subplots(1,1,figsize=(7,5))
		x = np.linspace(xmin,xmax,100)
		ymax = max(scipy.stats.norm.pdf(x,mean,std))
		title = "Survival probability of  %.3e at %.3e" % (scipy.stats.norm.sf(value,mean,std),value)
		label = "Gaussian (mean=%.3e,std=%.3e)" % (mean,std)
		axes.plot(x,scipy.stats.norm.pdf(x,mean,std),color='k',label=label)
		axes.arrow(value, 0.6*ymax, 0.0, -0.1, head_width=0.05, head_length=0.1, fc='r', ec='r')
		axes.legend(loc='upper left')
		axes.set_xlim(xmin,xmax)
		axes.set_ylim(-0.1*ymax,1.3*ymax)
		axes.set_xlabel('x')
		axes.set_ylabel('Probability')	
		axes.set_title(title)
		outpdf = 'gauss_m%s_s%s_v%s.pdf' % (str('%.3e' % mean),str('%.3e' % std),str('%.3e' % value))
		plt.savefig(outpdf)

	return scipy.stats.norm.sf(value,mean,std)