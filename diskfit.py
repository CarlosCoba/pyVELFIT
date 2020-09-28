import numpy as np
import matplotlib.pylab as plt
import random


from circular_mode import circ_mod
from radial_mode import rad_mod
from bisymetric_mode import bisym_mod

#from emcee_fit import EMCEE
def VELFIT(vel,evel,guess0,vary,sigma, delta,  rstart, rfinal, ring_space, frac_pixel, r_back, r_bar_max, pixel_scale, vmode):

	vrot0,vr20,pa0,inc0,x0,y0,vsys0 = guess0
	guess0 = [vrot0,vr20,pa0,inc0,x0,y0,vsys0,0,50]
	n_sigma = len(sigma)


	if n_sigma != 0:
		e_vrot,e_vr2,e_pa,e_inc,e_x0,e_y0,e_vsys = sigma




	if vmode == "circular":

		PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL  = circ_mod(vel, evel, guess0, vary, rstart =  rstart, rfinal = rfinal, ring_space = ring_space, frac_pixel = frac_pixel, r_back = r_back, delta = delta, pixel_scale = pixel_scale, r_bar_max = r_bar_max)

		return PA,INC,XC,YC,VSYS,0,R,Vrot,0*Vrot,0*Vrot,MODEL



	if vmode == "radial":
		PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL = rad_mod(vel, evel, guess0, vary, rstart =  rstart, rfinal = rfinal, ring_space = ring_space, frac_pixel = frac_pixel, r_back = r_back, delta = delta, pixel_scale = pixel_scale, r_bar_max = r_bar_max)
		return PA,INC,XC,YC,VSYS,0,R,Vrot,Vrad,0*Vrot,MODEL



	
	if vmode == "bisymmetric":
		PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL = bisym_mod(vel, evel, guess0, vary, rstart =  rstart, rfinal = rfinal, ring_space = ring_space, frac_pixel = frac_pixel, r_back = r_back, delta = delta, pixel_scale = pixel_scale, r_bar_max = r_bar_max)

		return PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL




import random

def RC_sigma_rnd2(vel,evel,nrings,guess0,vary,sigma = [],vmode = "rotation", delta=1,ring = "pixel",rstart = 4,iter = 5, ring_position = [],pixel_scale = 1 ):


	vrot0,vr20,pa0,inc0,X0,Y0,vsys0,vt20,theta0 = guess0
	guess = [vrot0,vr20,pa0,inc0,X0,Y0,vsys0,vt20,theta0]
	guess_copy = np.copy(guess)

	n_sigma = len(sigma)

	if n_sigma != 0:
		e_vrot,e_vr2,e_pa,e_inc,e_x0,e_y0,e_vsys,e_vt20,e_theta0 = sigma


	[ny,nx] = vel.shape

	x = np.arange(0, nx, 1)
	y = np.arange(0, ny, 1)
	mesh = np.meshgrid(x,y,sparse=True)



	from pixel_params import pixels
	import fit_params
	from fit_params import fit
	from fit_params import fit_polynomial
	from fit_params import fit_linear


	vary = [True,True,True,True,False,False,True,True,True]

	vr1,vr21,vsys1,pa1,inc1 = np.asarray([]),np.asarray([]),np.asarray([]),np.asarray([]),np.asarray([])
	vr0,vr20,vsys0,pa0,inc0,r0,vtan,theta0 = [],[],[],[],[],[],[],[]
	r = np.asarray([])

	mode = "bisymetric"
	for n_iter in range(iter):
		for i in ring_position:
			sigma = []
			sol = np.asarray(guess)
			sol[0] =sol[0]+10*random.uniform(-1,1)
			sol[1] =sol[1]+10*random.uniform(-1,1)
			sol[2] =sol[2]+e_pa*random.uniform(-1,1)
			sol[3] =sol[3]+e_inc*random.uniform(-1,1)
			sol[4] =sol[4]#+(2/pixel_scale)*random.uniform(-1,1)
			sol[5] =sol[5]#+(2/pixel_scale)*random.uniform(-1,1)
			sol[6] =sol[6]+e_vsys*random.uniform(-1,1)
			sol[7] =sol[7]+10*random.uniform(-1,1)
			sol[8] =sol[8]+e_theta0*random.uniform(-1,1)

			try:


				XY_mesh, vel_val, e_vel, f_pixel = pixels(vel,evel,guess,ringpos = i, delta=1,ring = "ARCSEC",pixel_scale=pixel_scale)
				vary = [True,True,True,True,True,False,False,True,True]

				if f_pixel > 0.30:

					Vrot , Vrad, Vsys,  pa, inc , Vtan, theta, xi_sq, n_data = fit("vsys",vel_val,e_vel,XY_mesh,guess,vary,vmode,fit_method = "nelder-mead", sigma = [])
					if 1>0:

						vr0.append(Vrot)
						vr20.append(Vrad)
						vtan.append(Vtan)
						theta0.append(theta)
						vsys0.append(Vsys)
						pa0.append(pa)
						inc0.append(inc)
						r0.append(i)
					else:

						vr0.append(np.nan)
						vr20.append(np.nan)
						vtan.append(Vtan)
						theta0.append(theta)
						vsys0.append(np.nan)
						pa0.append(np.nan)
						inc0.append(np.nan)
						r0.append(i)


			except(TypeError,ZeroDivisionError,ValueError):
			#except(1):
				pass

	r_sort, pa_sort = zip(*sorted(zip(r0, pa0)))
	r_sort, inc_sort = zip(*sorted(zip(r0, inc0)))
	r_sort, vr_sort = zip(*sorted(zip(r0, vr0))) 
	r_sort, vsys_sort = zip(*sorted(zip(r0, vsys0))) 
	r_sort, vr20_sort = zip(*sorted(zip(r0, vr20))) 

	r_sort, vtan_sort = zip(*sorted(zip(r0, vtan))) 
	r_sort, theta_sort = zip(*sorted(zip(r0, theta0))) 


	max_r = np.nanmax(r_sort)
	max_r = int(max_r)

	inc_sort = np.asarray(inc_sort)
	pa_sort = np.asarray(pa_sort)
	r_final = np.asarray(r_sort)
	vr_sort = np.asarray(vr_sort)
	vsys_sort = np.asarray(vsys_sort)
	vr20_sort = np.asarray(vr20_sort)
	vtan_sort = np.asarray(vtan_sort)
	theta_sort = np.asarray(theta_sort)
	r_sort = np.asarray(r_sort)

	R = np.unique(r_sort)


	pa_mean_ring = np.asarray([])
	inc_mean_ring = np.asarray([])
	vr_mean_ring = np.asarray([])
	vsys_mean_ring = np.asarray([])
	vr20_mean_ring = np.asarray([])
	r_mean_ring = np.asarray([])

	vtan_mean_ring = np.asarray([])
	theta_mean_ring = np.asarray([])


	pa_e = np.asarray([])
	inc_e = np.asarray([])
	vr_e = np.asarray([])
	vsys_e = np.asarray([])
	vr20_e = np.asarray([])

	vtan_e = np.asarray([])
	theta_e = np.asarray([])

	#print(len(r_sort),len(ring_position))
	#for kk in range(max_r+1):
	for kk in R:
		mask = r_sort == kk
		pa_i = pa_sort[mask]
		inc_i = inc_sort[mask]
		vr_i = vr_sort[mask]
		vsys_i = vsys_sort[mask]
		vr20_i = vr20_sort[mask]

		vtan_i = vtan_sort[mask]
		theta_i = theta_sort[mask]

		if len(pa_i)>0:
			vr_mean_ring = np.append(vr_mean_ring,np.nanmean(vr_i))
			vsys_mean_ring = np.append(vsys_mean_ring,np.nanmean(vsys_i))
			vr20_mean_ring = np.append(vr20_mean_ring,np.nanmean(vr20_i))
			pa_mean_ring = np.append(pa_mean_ring,np.nanmean(pa_i))
			inc_mean_ring = np.append(inc_mean_ring,np.nanmean(inc_i))
			r_mean_ring = np.append(r_mean_ring,kk)
			vtan_mean_ring = np.append(vtan_mean_ring,np.nanmean(vtan_i))
			theta_mean_ring = np.append(vr20_mean_ring,np.nanmean(theta_i))

			vr_e = np.append(vr_e,np.nanstd(vr_i))
			vr20_e = np.append(vr20_e,np.nanstd(vr20_i))
			pa_e = np.append(pa_e,np.nanstd(pa_i))
			inc_e = np.append(inc_e,np.nanstd(inc_i))
			vtan_e = np.append(vtan_e,np.nanstd(vtan_i))
			theta_e = np.append(theta_e,np.nanstd(theta_i))

		#else:
		#	vr_mean_ring = np.append(vr_mean_ring,20)
		#	pa_mean_ring = np.append(pa_mean_ring,median_pa)
		#	inc_mean_ring = np.append(inc_mean_ring,median_inc)
		#	r_mean_ring = np.append(r_mean_ring,kk)





	return vr_e,vr20_e,pa_e,inc_e,vsys_e,vtan_e,theta_e
	#return vr0,vr20,vsys0,pa0,inc0,r




def RC_sigma_rnd(vel,evel,nrings,guess0,vary,sigma = [],mode = "rotation", delta=1,ring = "pixel",rstart = 4,iter = 5, ring_position = [] ,pixel_scale = 1):

	rstart = int(rstart/pixel_scale)
	vrot0,vr20,pa0,inc0,X0,Y0,vsys0 = guess0
	guess = [vrot0,vr20,pa0,inc0,X0,Y0,vsys0]
	guess_copy = np.copy(guess)

	n_sigma = len(sigma)

	if n_sigma != 0:
		e_vrot,e_vr2,e_pa,e_inc,e_x0,e_y0,e_vsys = sigma


	[ny,nx] = vel.shape

	x = np.arange(0, nx, 1)
	y = np.arange(0, ny, 1)
	mesh = np.meshgrid(x,y,sparse=True)



	from pixel_params import pixels
	import fit_params
	from fit_params import fit
	from fit_params import fit_polynomial
	from fit_params import fit_linear


	vary = [True,True,True,True,False,False,True]
	#ring_position = np.arange(rstart,nrings,2) 
	#ring_position = np.arange(rstart,nrings,4)

	vr1,vr21,vsys1,pa1,inc1 = np.asarray([]),np.asarray([]),np.asarray([]),np.asarray([]),np.asarray([])
	vr1,vr21,vsys1,pa1,inc1,r1 = [],[],[],[],[],[]
	r = np.asarray([])

	for n_iter in range(iter):
	#for i in ring_position:
		vr0,vr20,vsys0,pa0,inc0 =  np.asarray([]),np.asarray([]),np.asarray([]),np.asarray([]),np.asarray([])
		#for n_iter in range(iter):
		for i in ring_position:
			sigma = []
			sol = np.asarray(guess)
			sol[0] =sol[0]+10*random.uniform(-1,1)
			sol[1] =sol[1]+10*random.uniform(-1,1)
			sol[2] =sol[2]+e_pa*random.uniform(-1,1)
			sol[3] =sol[3]+e_inc*random.uniform(-1,1)
			sol[4] =sol[4]#+(2/pixel_scale)*random.uniform(-1,1)
			sol[5] =sol[5]#+(2/pixel_scale)*random.uniform(-1,1)
			sol[6] =sol[6]+e_vsys*random.uniform(-1,1)

			try:
				XY_mesh, vel_val, e_vel, f_pixel = pixels(vel,evel,guess,ringpos = i, delta=4,ring = "pixel",pixel_scale=pixel_scale)
				if f_pixel > 0.10:

					Vrot , Vr2, Vsys,  pa, inc ,xi = fit("vsys",vel_val,e_vel,XY_mesh,sol,vary,mode,sigma)
					if Vrot>10 and Vrot<vmaxrot:

						vr0 = np.append(vr0,Vrot)
						vr20 = np.append(vr20,Vr2)
						vsys0 = np.append(vsys0,Vsys)
						pa0 = np.append(pa0,pa)
						inc0 = np.append(inc0,inc)
						if n_iter == 0:
							r1.append(i)

					else:

						vr0 = np.append(vr0,Vrot)
						vr20 = np.append(vr20,Vr2)
						vsys0 = np.append(vsys0,Vsys)
						pa0 = np.append(pa0,pa)
						inc0 = np.append(inc0,inc)
						if n_iter == 0:
							r1.append(i)
				#else: print(i,"aqtuiiii")

			except(TypeError,ZeroDivisionError,ValueError):
				pass

		#if len(vr0) > 0 and  len(vr20)>0 and len(vsys0) >0 and len(pa0) >0 and len(inc0) >0:
		#if 1>0:


		#vr1=np.append(vr1,vr0)
		#vr21 = np.append(vr21,vr20)
		#vsys1 = np.append(vsys1,vsys0)
		#inc1 = np.append(inc1,inc0)
		#pa1 = np.append(pa1,pa0)

			#if n_iter == 0:
			#	r1.append(i)

		if len(vr0) !=0:
			vr1.append(vr0)
			vr21.append(vr20)
			vsys1.append(vsys0)
			inc1.append(inc0)
			pa1.append(pa0)

	return vr1,vr21,vsys1,pa1,inc1,r1
	#return vr0,vr20,vsys0,pa0,inc0,r






def RC_emcee(vel,evel,nrings,guess0,vary,sigma = [],mode = "rotation", delta=1,ring = "pixel",rstart = 4,iter = 5, pos = 2,pixel_scale = 1):

	rstart = int(rstart/pixel_scale)
	vrot0,vr20,pa0,inc0,X0,Y0,vsys0 = guess0
	guess = [vrot0,vr20,pa0,inc0,X0,Y0,vsys0]
	guess_copy = np.copy(guess)

	n_sigma = len(sigma)

	if n_sigma != 0:
		e_vrot,e_vr2,e_pa,e_inc,e_x0,e_y0,e_vsys = sigma


	[ny,nx] = vel.shape

	x = np.arange(0, nx, 1)
	y = np.arange(0, ny, 1)
	mesh = np.meshgrid(x,y,sparse=True)



	from pixel_params import pixels
	from emcee_fit import EMCEE
	import fit_params
	from fit_params import fit
	from fit_params import fit_polynomial
	from fit_params import fit_linear


	vary = [True,True,True,True,False,False,True]
	ring_position = np.arange(rstart,nrings,2) 



	vr1,vr21,vsys1,pa1,inc1 = [],[],[],[],[]
	r = []
	for i in ring_position:
		vr0,vr20,vsys0,pa0,inc0 = [],[],[],[],[]



		for n_iter in range(iter):
			#sigma = []
			sol = np.asarray(guess)
			sol[0] =sol[0]+10*random.uniform(-1,1)
			sol[1] =sol[1]+10*random.uniform(-1,1)
			sol[2] =sol[2]+e_pa*random.uniform(-1,1)
			sol[3] =sol[3]+e_inc*random.uniform(-1,1)
			sol[4] =sol[4]+(2/pixel_scale)*random.uniform(-1,1)
			sol[5] =sol[5]+(2/pixel_scale)*random.uniform(-1,1)
			sol[6] =sol[6]+e_vsys*random.uniform(-1,1)

			try:
				XY_mesh, vel_val, e_vel, f_pixel = pixels(vel,evel,sol,ringpos = i, delta=4,ring = "pixel",pixel_scale=pixel_scale)
				if f_pixel > 0.50:
					res = EMCEE(XY_mesh, vel_val, e_vel,guess,sigma)
			except(TypeError,ZeroDivisionError,ValueError):pass
			#except(1):pass

	#return vr1,vr21,vsys1,pa1,inc1,r
	return 0






