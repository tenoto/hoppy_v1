def get_luminosity(flux,d_kpc):
	luminosity = 1.2e+32 * (flux / 1e-12) * d_kpc**2 
	return luminosity
