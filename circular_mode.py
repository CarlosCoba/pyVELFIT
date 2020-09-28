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
from pixel_params import pixels
import fit_params
from fit_params import fit




def circ_mod(vel, evel, guess0, vary, rstart = 2, rfinal = 55, ring_space = 2, frac_pixel = 0.7, r_back = 20, delta = 1, pixel_scale = 0.2, r_bar_max = 20):

		vrot0,vr20,pa0,inc0,x0,y0,vsys0,vtan,theta_b = guess0
		vmode = "circular"
		[ny,nx] = vel.shape

		"""

		 					CIRCULAR MODEL


		"""

		chisq_global = 1e4
		PA, INC, XC,YC,VSYS = 0,0,0,0,0
		Vrot, Vrad, Vsys,Vtan = [],[],[],[]
		R = 0

		if r_back < ring_space:
			r_back = ring_space


		for jj in range(0,r_back,ring_space):
			for it in range(5):
				guess = [vrot0,vr20,pa0,inc0,x0,y0,vsys0,0,theta_b]


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

						# Create model
						try:
							w_rot,w_rad,vrot,vrad = M_radial(XY_mesh,0,0,pa0,inc0,x0,y0,vsys0,0,0,vel_val-vsys0,e_vel,vmode)
							vrot_model,vrad_model = np.append(vrot_model,vrot),np.append(vrad_model,vrad)



						except(np.linalg.LinAlgError):
							vrot_model,vrad_model = np.append(vrot_model,100),np.append(vrad_model,0)



						r.append(j)

						los_vel.append(vel_val)
						e_los_vel.append(e_vel)
						x_pos.append(XY_mesh[0])
						y_pos.append(XY_mesh[1])
						xy_pos.append(XY_mesh)




				#vary = [True,True,True,True,True,True,True]
				N = len(los_vel)
				guess = [vrot_model,vrad_model,pa0,inc0,x0,y0,vsys0, vtan_model,theta_b]
				vrot , vsys0,  pa0, inc0, x0, y0, xi_sq, n_data = fit("vsys",los_vel,e_los_vel,xy_pos,guess,vary,vmode, sigma = [],fit_method = "Powell")


				if inc0>85: 
					pa0 = guess0[2]
					x0,y0 = guess0[4],guess0[5]
					inc0 = 55
					frac_pixel = 0.6
					xi_sq = 1e5


				if 1.1*xi_sq < chisq_global:

					PA, INC, XC,YC,VSYS,THETA = pa0, inc0, x0, y0,vsys0,theta_b
					R = r
					Vrot = vrot
					chisq_global = xi_sq


					MODEL = np.zeros((ny,nx)) 

					for k in range(N):
						for mm,nn in zip(xy_pos[k][0],xy_pos[k][1]): 
							MODEL[nn][mm] = Vlos(mm,nn,Vrot[k],0,PA,INC,XC,YC,VSYS,0,0,vmode) - VSYS



		MODEL[MODEL == 0] = np.nan



		#plt.imshow(MODEL, origin = "l", cmap = califa)
		#plt.show()

		Vrot = np.array(Vrot)
		R = np.array(R)


		#plt.plot(R,Vrot,"k-")
		#plt.plot(R,Vrot,"ko")
		#plt.show()

		return PA,INC,XC,YC,VSYS,0,R,Vrot,0*Vrot,0*Vrot,MODEL


