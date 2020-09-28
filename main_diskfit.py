import numpy as np
import sys
from initialize_velfit import rotcur


"""
INPUT:
galaxy = 	object name.
vel_map = 	velocity map in km/s.
evel_map = 	error map in km/s.
SN = 		Maximum value in the error map allowed for computing the rotation curve.
		Better results are obtained when SN has low values 15-20 km/s.
		if SN == None then the whole velocity map will be used for estimating the rotation curve
z_star	=	Redshift of the object. 
		This value will be passed as z*c for the guess in Vsys.
PA	= 	Position Angle of the major axis			
INC	=	Inclination.
X0,Y0	=	Pixel coordenates of the photometric center. 
pixel_scale	Size of the pixel in arcsecs

This two are for computin the escape velocity
Mstar	=	Stellar mass
Reff	=	Effective radius in arcsec

ADDITIONAL
The following directories need to be created in the same path
./plots
./fits
./vmax_rturn
./v_esc_norm_Re
"""


#
# RUNNING THE CODE:
#
#rotcur(galaxy, vel_map,evel_map,SN,z_star,PA,INC,X0,Y0,pixel_scale,vary_PA= True,vary_INC=True,vary_XC = True,vary_YC = True,vary_VSYS = True,delta = 2, rstart=2, rfinal = 55, ring_space = 3, frac_pixel = 1/3., r_back = 20, r_bar_max = 25,vmode = "circular")


#
# EXAMPLE
#

fits_file = "vels.fits"
fits_file = "NGC_3198_NA_MOM1_THINGS.FITS"

from astropy.io import fits
data = fits.getdata(fits_file)*1e-3
data = data[0,0,:,:]

[ny,nx] = data.shape
pixel_scale = 1
pixel_scale = 4.166666768E-04*3600


vel_ha = data
e_vel_ha= np.ones((ny,nx))

galaxy = "simulation"
galaxy = "THINGS"
XC,YC = 130.50,126.50
XC,YC = 513,510
vel_map,evel_map,SN,z_star,PA,INC,X0,Y0,pixel_scale,Mstar,Reff = vel_ha,e_vel_ha,1e3,110/3e5,45,60,130.50,126.50,pixel_scale,10,10
vel_map,evel_map,SN,z_star,PA,INC,X0,Y0,pixel_scale,Mstar,Reff = vel_ha,e_vel_ha,1e3,576/3e5,22,67,XC,YC,pixel_scale,10,10

try: 
	#Maybe change the ring widths?
	#rotcur(galaxy, vel_map,evel_map,SN,z_star,PA,INC,X0,Y0,pixel_scale,Mstar,Reff,vary_PA= True,vary_INC=True,vary_XC = True,vary_YC = True,vary_VSYS = True,delta = 3, ring_space = 6,rstart=3, rfinal = 140, frac_pixel = 1/3., r_back = 0, r_bar_max = 110,vmode = "bisymmetric")

	rotcur(galaxy, vel_map,evel_map,SN,z_star,PA,INC,X0,Y0,pixel_scale,Mstar,Reff,vary_PA= True,vary_INC=True,vary_XC = True,vary_YC = True,vary_VSYS = True,delta = 8, ring_space = 16,rstart=3, rfinal = 1e3, frac_pixel = 1/3., r_back = 0, r_bar_max = 110,vmode = "bisymmetric")

#except(ValueError,AttributeError,TypeError,ZeroDivisionError,OSError,AssertionError, IndexError,UnboundLocalError):
except(1):
	pass

