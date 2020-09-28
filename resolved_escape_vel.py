import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp1d
from astropy.io import fits
from v_escape_Eq import vel_escape


def Rings(xy_mesh,pa,inc,x0,y0,pixel_scale):
	PA,inc=(pa)*np.pi/180,inc*np.pi/180
	(x,y) = xy_mesh

	X = (- (x-x0)*np.sin(PA) + (y-y0)*np.cos(PA))
	Y = (- (x-x0)*np.cos(PA) - (y-y0)*np.sin(PA))

	R= np.sqrt(X**2+(Y/np.cos(inc))**2)
	return R*pixel_scale





def resolved_vescape(galaxy,PA,INC,X0,Y0,Mstar,Reff,v_max,R_turn,nx,ny,pixel_scale,redshift):
	# Reff in arcsec
	# Mstar in log
	# v_max in km/s
	# Rturn in arcsec
	try:
		#print(PA,INC,X0,Y0,Mstar,Reff,v_max,R_turn,pixel_scale)

		x = np.arange(0, nx, 1)
		y = np.arange(0, ny, 1)

		XY_mesh = np.meshgrid(x,y,sparse=True)
		R = Rings(XY_mesh, PA,INC,X0,Y0,pixel_scale)
		

		# Prepare for estimate resolved V_escape
		if np.isfinite(v_max) == False: 1/0



		escape = np.zeros((ny,nx))
		for j in range(ny):
			for ii in range(nx):
				r = R[j][ii]
				escape_at_Re = vel_escape(Reff,v_max,R_turn,Mstar,redshift)
				escape[j][ii] = vel_escape(r,v_max,R_turn,Mstar,redshift)/escape_at_Re

		

		escape[escape == 0] = np.nan

		from write_table import write
		table = [galaxy,escape_at_Re]
		write(table,"v_esc_at_Re.MaNGA.csv",column = False)


		hdu = fits.ImageHDU()
		hdu.data = escape
		hdu.header['FILE0'] = 'v_escape/v_escape_at_Re'
		hdu.header['UNITS0'] = "dimensionless"
		hdu.header['FILE1'] = "v_esc at Re"
		hdu.header['VESC_RE'] = escape_at_Re
		hdu.header['UNITS1'] = "km/s"
		hdu.writeto("./v_esc_norm_Re/%s.escape_vel.fits"%galaxy,clobber=True)

	except(ZeroDivisionError,ValueError, IndexError):
		escape = np.zeros((ny,nx))
		hdu = fits.PrimaryHDU()
		hdu.data = escape
		hdu.writeto("./v_esc_norm_Re/%s.escape_vel.fits"%galaxy,clobber=True)

