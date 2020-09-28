import numpy as np




H0=70. #km/s/Mpc
Omega_m=0.27
Omega_l=0.73


Omega_m = 0.307
Omega_l = 0.693
H0 = 67.8

c=299792

##############
# Commoving radial distance X of a source at redshift z



from scipy.integrate import quad



def kpc_2_arcsec(z):

	def X(a):
		t1=(c/H0)/np.sqrt(a*Omega_m+a**2*(1-Omega_m-Omega_l)+a**4*Omega_l)
		return t1



	a0=1./(1+z)
	a1=1
	X = quad(X,a0,a1)
	dA=X[0]/(1+z)
	dL=X[0]*(1+z)
	scale=1*dA*1e3/206265.

	#scale[kpc/arcsec]
	return scale
