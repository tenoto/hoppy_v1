def get_luminosity(flux,d_kpc):
	luminosity = 1.2e+32 * (flux / 1e-12) * d_kpc**2 # erg/s
	return luminosity

def get_spindown_luminosity(period, pdot, moment_of_inertia=1.0):
	"""
	moment of inertia normalized to I = 1e+45 g cm2
	"""
	spindown_luminosity = 3.94e+35 * moment_of_inertia / period**3 * (pdot/1e-11)
	return spindown_luminosity