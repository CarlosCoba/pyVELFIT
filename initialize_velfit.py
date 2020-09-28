import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits

from axis import AXIS 
from fit_rotcurve import Rturn_Vmax


from diskfit import VELFIT 
from diskfit import RC_sigma_rnd2 
from diskfit import RC_emcee


from kinematic_centre_vsys import KC
from CBAR import colorbar as cb


c = 3e5

def SN(flux,eflux,sn):
	A = np.divide(flux,eflux)
	A[A<sn] = np.nan
	return np.divide(A,A)


def rotcur(galaxy, vel_map,evel_map,SN,z_star,PA,INC,X0,Y0,pixel_scale,Mstar,Reff,vary_PA,vary_INC,vary_XC,vary_YC,vary_VSYS,delta, rstart, rfinal, ring_space, frac_pixel, r_back, r_bar_max,vmode, save_plots = 1 ):
	"""
	vary = [Vrot,Vrad,PA,INC,XC,YC,VSYS,theta,Vtan]
	"""
	[ny,nx] = vel_map.shape
	ext = [nx/2., -nx/2,-ny/2.,ny/2.]


	#if evel_map == None:
	#	
	#	evel_map = np.ones((ny,nx))
	#else:
	#
	#	evel_map[evel_map>SN]=np.nan



	mask_vel=np.divide(evel_map,evel_map)
	vel_ha = vel_map*mask_vel
	vel_ha[vel_ha==0]=np.nan

	XK,YK,VSYS,e_vsys = KC(vel_map,X0,Y0,pixel_scale)
	VSYS = z_star*c
	guess = [50,0,PA,INC,X0,Y0,VSYS]



	def first():

		delta_arcsec = int(2/1)

		vary = [True,True,vary_PA,vary_INC,vary_XC,vary_YC,vary_VSYS,True,True]
		sigma = []
		PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL = VELFIT(vel_ha,evel_map,guess,vary, sigma, delta, rstart, rfinal, ring_space, frac_pixel,r_back, r_bar_max,pixel_scale,vmode)
		



		return PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL




	PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL = first()

	#
	# I havent calculated the errors! 20.08.2020
	#


	e_PA,e_VSYS,e_INC = 1,1,1


	def emcee():
		vary = [ True, True, True, True,False, False, True]
		guess = [60,0,pa_med,inc_med,XK,YK,vsys_med]
		sigma = [0,0,5*e_pa,5*e_inc,0,0,5*e_vsys]
		delta_arcsec = int(2/1)

		print("sigma:", sigma)

		a= RC_emcee(vel_ha,e_vel_ha,nrings,guess,vary,sigma,mode = "radial", delta=delta_arcsec,ring = "pixel",iter = 5,pixel_scale = pixel_size)
		return a


	#c = emcee()

	#
	## Write output
	#
	from write_table import write
	if vmode == "circular" or vmode == "radial":
		kin_params_table = "ana_kin_model_MaNGA_%s.csv"%vmode
		kin_params = [galaxy,XC,YC,PA,INC,VSYS]
		write(kin_params,kin_params_table,column = False)

	if vmode == "bisymmetric":
		kin_params_table = "ana_kin_model_MaNGA_%s.csv"%vmode
		kin_params = [galaxy,XC,YC,PA,INC,VSYS,THETA]
		write(kin_params,kin_params_table,column = False)



	e_vr,e_Vrad,e_pa,e_inc,e_vsys,e_Vtan,e_theta = 0,0,0,0,0,0,0


	save_file = save_plots
	from plot_models import plot_kin_models
	plot_kin_models(galaxy,vmode,vel_ha,R,Vrot,Vrad,Vtan, VSYS, MODEL, ext,plot = 0, save = save_file)

	from kin_model_fits import save_model
	if vmode == "circular": VMODELS = [Vrot]
	if vmode == "radial": VMODELS = [Vrot,Vrad]
	if vmode == "bisymmetric": VMODELS = [Vrot,Vrad,Vtan]


	s = save_model(galaxy,vmode,R,VMODELS,PA,INC,XC,YC,VSYS,save=save_file)

	from fit_rotation_curve_fits import fit_rotcur
	Vmax,Rturn = fit_rotcur(galaxy,vmode)

	from resolved_escape_vel import resolved_vescape
	resolved_vescape(galaxy,PA,INC,X0,Y0,Mstar,Reff,Vmax,Rturn,nx,ny,pixel_scale,z_star)


	#return Vmax,Rturn



 
