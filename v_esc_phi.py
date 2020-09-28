import numpy as np
import matplotlib.pylab as plt

from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker



from NFW_potential import M_vir
from concentration import concetration
from virial_radius import r_vir
from scipy.interpolate import interp1d


Om_mat_cero = 0.307
Om_lambda_cero = 0.693


from angular_distance import kpc_2_arcsec


def g(c):
	g=1/(np.log(1+c)-c/(1+c))
	return g


def v_escape_halo(r,Mstar,z):

	Mvir=M_vir(10**Mstar)
	c= concetration(Mvir)
	R_vir= r_vir(10**Mvir,0.016,Om_mat_cero,Om_lambda_cero)
	#R_vir is in Kpc

	#
	# r must be in arcsec
	#
	G=(4.30091e-3)*(1e-3)#kpc*Msun-1 (km/s)**2
	Mv = 10**Mvir
	Rv = R_vir


	scale = kpc_2_arcsec(z)# [kpc/arcsec]
	Rvir_2_arcsec = Rv*(1./scale)

	#r = r
	if r == 0:
		phi = -(c*g(c)*G*Mv/Rv)
	else:
		s= np.divide(r,Rvir_2_arcsec)
		v_virial_square = (G*Mv/Rv)
		phi = -g(c)*(G*Mv/Rv)*np.log(1+c*s)/s
	v_esc=2*abs(phi)
	return R_vir,v_esc


