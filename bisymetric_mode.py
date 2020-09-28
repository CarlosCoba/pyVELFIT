import numpy as np
import matplotlib.pylab as plt
import lmfit
import sys
import matplotlib
from lmfit import Model
from lmfit import Parameters, fit_report, minimize
from astropy.stats import sigma_clip
from scipy.interpolate import interp1d
from matplotlib.colors import Normalize
#from numpy.linalg import LinAlgError

from vlos_bisym import Vlos

import cmap_califa
import cmap_vfield
califa=cmap_vfield.CALIFA()



from model_params import M_radial 
from model_params import Rings 
from pixel_params import pixels
import fit_params
from fit_params import fit



def bisym_mod(vel, evel, guess0, vary, rstart, rfinal, ring_space, frac_pixel, r_back, delta, pixel_scale, r_bar_max):
		plt.imshow(vel)
		plt.show()

		vrot0,vr20,pa0,inc0,x0,y0,vsys0,vtan,theta_b = guess0
		vmode = "bisymmetric"
		[ny,nx] = vel.shape


		"""

		 					BYSIMETRIC MODEL


		"""


		chisq_global = 1e10
		PA, INC, XC,YC,VSYS = 0,0,0,0,0
		Vrot, Vrad, Vsys,Vtan = [],[],[],[]
		R = 0

		if r_back < ring_space:
			r_back = ring_space


		#for jj in range(0,r_back,ring_space):
		for jj in range(1):
			for it in range(5):
				guess = [vrot0,vr20,pa0,inc0,x0,y0,vsys0,0,theta_b]

				#theta_b = np.arctan(np.tan(theta_b*np.pi/180-pa0*np.pi/180)/np.cos(inc0*np.pi/180))*180/np.pi

				xi_sq_array = np.asarray([])
				N = np.asarray([])


				vrot_model,vrad_model,vtan_model = np.asarray([]),np.asarray([]),np.asarray([])
				los_vel = np.array([])
				e_los_vel = np.array([])
				x_pos = np.array([])
				y_pos = np.array([])

				los_vel = []
				e_los_vel = []
				x_pos = []
				y_pos = []
				xy_pos = []
				r = []


				for j in np.arange(rstart,rfinal-jj,ring_space):

					XY_mesh, vel_val, e_vel, f_pixel = pixels(vel,evel,guess,ringpos = j, delta=delta,ring = "ARCSEC",pixel_scale=pixel_scale)
					npixels = len(XY_mesh[0])


					if j < 10: 
						f_pixel = 1


					if f_pixel > frac_pixel:


						if j < r_bar_max:
							# Create bisymetric model
							w_sys,w_rot,w_rad,w_tan,vrot,vrad,vsys,vtan = M_radial(XY_mesh,0,0,pa0,inc0,x0,y0,vsys0,0,theta_b,vel_val-vsys0,e_vel,vmode)
							vrot_model,vrad_model,vtan_model = np.append(vrot_model,vrot),np.append(vrad_model,vrad),np.append(vtan_model,vtan)

						else:
							# Create ciruclar model
							w_rot,w_rad,vrot,vrad = M_radial(XY_mesh,0,0,pa0,inc0,x0,y0,vsys0,0,theta_b,vel_val-vsys0,e_vel,vmode= "circular")
							vrot_model,vrad_model,vtan_model = np.append(vrot_model,vrot),np.append(vrad_model,0),np.append(vtan_model,0)

						r.append(j)
						los_vel.append(vel_val)
						e_los_vel.append(e_vel)
						x_pos.append(XY_mesh[0])
						y_pos.append(XY_mesh[1])
						xy_pos.append(XY_mesh)



				N = len(los_vel)

				guess = [vrot_model,vrad_model,pa0,inc0,x0,y0,vsys0, vtan_model,theta_b]

				vrot , vrad, vsys0,  pa0, inc0 , x0, y0, vtan, theta_b, xi_sq, n_data = fit("vsys",los_vel,e_los_vel,xy_pos,guess,vary,vmode,fit_method = "Powell", sigma = [], r_bar_max = r_bar_max)
				print(pa0, inc0 , x0, y0, theta_b,xi_sq)

				if 1.1*xi_sq < chisq_global:

					PA, INC, XC,YC,VSYS,THETA = pa0, inc0, x0, y0,vsys0,theta_b
					Vrot = vrot
					Vrad = vrad
					Vtan = vtan
					chisq_global = xi_sq
					R = r

					MODEL = np.zeros((ny,nx)) 

					for k in range(N):
						for mm,nn in zip(xy_pos[k][0],xy_pos[k][1]): 
							MODEL[nn][mm] = Vlos(mm,nn,Vrot[k],Vrad[k],PA,INC,XC,YC,VSYS,Vtan[k],THETA,vmode) - VSYS




		MODEL[MODEL == 0] = np.nan
		#plt.imshow(MODEL, vmin = -260, vmax = 260,  origin = "l", cmap = califa)
		#plt.show()



		"""

		#los_vel,e_los_vel,x_pos,y_pos,xy_pos = [],[],[],[],[]
		vrad_,vrot_,vtan_ = [],[],[]
		R2 = []
		ring_space = 1
		for j in np.arange(1,35,ring_space):

				xy_pos = []
				los_vel = []
				e_los_vel = []
				los_vel = np.array([])
				XY_mesh, vel_val, e_vel, f_pixel = pixels(vel,evel,guess,ringpos = j, delta=0.5,ring = "ARCSEC",pixel_scale=pixel_scale)
				npixels = len(XY_mesh[0])

				xy_pos.append(XY_mesh)
				#los_vel.append(vel_val)
				los_vel = [vel_val]
				#los_vel = np.append(los_vel,vel_val)
				e_los_vel.append(e_vel)

				if f_pixel > 0.5:

					# Create model
					w_sys,w_rot,w_rad,w_tan,vrot,vrad,vsys,vtan = M_radial(XY_mesh,0,0,PA,INC,XC,YC,VSYS,0,0,vel_val-VSYS,e_vel,vmode)
					vrot_model,vrad_model,vtan_model = np.append(vrot_model,vrot),np.append(vrad_model,vrad),np.append(vtan_model,vtan)
					R2.append(j)


					vary = [True,True,False,False,False,False,False,True,False]
					guess = [[vrot],[vrad],PA,INC,XC,YC,VSYS,[vtan],THETA]


					vel_val, e_vel, XY_mesh = np.array([vel_val]), np.array([e_vel]),[XY_mesh]



					vrot0 , vrad0, vsys,  pa, inc , x0, y0, vtan0, theta, xi_sq, n_data = fit("vsys",vel_val,e_vel,xy_pos,guess,vary,vmode,fit_method = "Powell", sigma = [], r_bar_max = r_bar_max)
					vrot_.append(vrot0)
					vtan_.append(vtan0)
					vrad_.append(vrad0)
	

		print(PA, INC, XC,YC,VSYS,theta, chisq_global)
		"""

		#plt.plot(R,Vrot,"b-", label = "min Xi")
		#plt.plot(R,Vrad,"g-")
		#plt.plot(R,Vtan,"r-")

		#plt.plot(R2,vrot_,"b--", label = "fix params")
		#plt.plot(R2,vrad_,"g--")
		#plt.plot(R2,vtan_,"r--")
		#plt.legend()
		#plt.show()

		R = np.asarray(R)
		Vtan = np.asarray(Vtan)
		Vrad = np.asarray(Vrad)
		Vrot = np.asarray(Vrot)

		return PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL
