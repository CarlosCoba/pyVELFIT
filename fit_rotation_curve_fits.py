import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits
from scipy.optimize import curve_fit


from write_table import write
import errno
import os


from figure import  fig_ambient
from axis import AXIS


"""
Rotation Curve Models
"""



# From Bouche + 2015
def exp(r_sky,v_max,R_turn):
	return v_max*(1-np.exp(-r_sky/R_turn))

def arctan(r_sky,v_max,R_turn):
	return (2/np.pi)*v_max*np.arctan(r_sky/R_turn)

def tanh(r_sky,v_max,R_turn):
	return v_max*np.tanh(r_sky/R_turn)

#From Haeun Chung, to model the rising RC
def tanh_linear(r_sky,v_max,R1,R2):
		R_turn = R1
		# R2 can be negatica in case of descending RCs
		v=v_max*(np.tanh(r_sky/R_turn) + r_sky*R2 )
		return v
#From Barrera-Ballesteros
def multi_param(r_sky,v_max,R_turn,gamma):
		beta = 0
		x=R_turn/r_sky
		A = (1+x)**beta
		B = (1+x**gamma)**(1./gamma)
		v=v_max*A/B
		return v

# Bertola 1991
def bertola(r_sky,v_max,k,gamma):
	v = v_max*r_sky/(r_sky**2 + k**2)**(gamma/2.)
	return v



def fit_RC(r_array,v_array,method):
	if method == "exp":
		try:
			popt, pcov = curve_fit(exp, r_array, v_array, bounds=([0,0], [350.,60]))
			vmax,r_turn = popt
			best_fit = exp(r_array,vmax,r_turn)
		except(RuntimeError):
			vmax,r_turn,best_fit = 0,0,r_array*0
		return vmax,r_turn,best_fit

	if method == "arctan":
		try:
			popt, pcov = curve_fit(arctan, r_array, v_array, bounds=([0,0], [350.,60]))
			vmax,r_turn = popt
			best_fit = arctan(r_array,vmax,r_turn)
		except(RuntimeError):
			vmax,r_turn,best_fit = 0,0,r_array*0
		return vmax,r_turn,best_fit

	if method == "tanh":
		try:
			popt, pcov = curve_fit(tanh, r_array, v_array, bounds=([0,0], [350.,60]))
			vmax,r_turn = popt
			best_fit = tanh(r_array,vmax,r_turn)
		except(RuntimeError):
			vmax,r_turn,best_fit = 0,0,r_array*0
		return vmax,r_turn,best_fit

	if method == "tanh-linear":
		try:
			# R2 can be negative, common values 0.010 < R2 < 0.070
			popt, pcov = curve_fit(tanh_linear, r_array, v_array, bounds=([0,0,0], [350.,60,0.1]))
			vmax,r_turn,R2 = popt
			best_fit = tanh_linear(r_array,vmax,r_turn,R2)
		except(RuntimeError):
			vmax,r_turn,R2,best_fit = 0,0,0,r_array*0
		return vmax,r_turn,R2,best_fit

	if method == "multi-param":
		try:
			popt, pcov = curve_fit(multi_param, r_array, v_array, bounds=([0,0,-1], [350.,60,12]))
			vmax, r_turn,gamma = popt
			best_fit = multi_param(r_array,vmax,r_turn,gamma)
		except(RuntimeError):
			vmax,r_turn,gamma,best_fit = 0,0,0,r_array*0
		return vmax,r_turn,gamma,best_fit


	if method == "bertola":
		try:
			popt, pcov = curve_fit(bertola, r_array, v_array, bounds=([0,0,1.], [350.,60,3/2.]))
			vmax, k,gamma = popt
			best_fit = bertola(r_array,vmax,k,gamma)
		except(RuntimeError):
			vmax, k,gamma,best_fit = 0,0,0,r_array*0
		return vmax, k,gamma,best_fit





model = ["exp", "arctan", "tanh","tanh-linear","multi-param","bertola"]


axes = fig_ambient(fig_size = (4,4), ncol = 1, nrow = 1, left = 0.15, right = 0.95, top = 0.95, hspace = 0, bottom = 0.15, wspace = 0 )
ax = axes[0]

def fit_rotcur(galaxy,vmode):
	try:

		data = fits.getdata("./fits/%s.%s_model.fits"%(galaxy,vmode))

		R = data[0]
		V = data[1]
		mask = V>5
		R=R[mask]
		V=V[mask]
		min_vel,max_vel = int(np.nanmin(V)),int(np.nanmax(V))
		temp = []
		v_max_temp = []
		R_turn_temp = []
		for i in model:
			best_params = fit_RC(R,V,i)
			data = best_params[:-1]
			Vmax,Rturn = data[0],data[1]
			if Vmax > 340: Vmax = 0 
			if Vmax == 0: Vmax = np.nan 
			temp.append(Vmax)
			temp.append(Rturn)
			v_max_temp.append(Vmax)
			R_turn_temp.append(Rturn)

			ax.plot(R,V,"k-")
			ax.plot(R,V,"ks")
			ax.plot(R,best_params[-1],label = i)


		mean_V,median_V,std_V = np.nanmedian(v_max_temp),np.nanmean(v_max_temp),np.nanstd(v_max_temp)
		mean_R,median_R,std_R = np.nanmedian(R_turn_temp),np.nanmean(R_turn_temp),np.nanstd(R_turn_temp)
		table = [galaxy, temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],temp[7],temp[8],temp[9],temp[10],temp[11],mean_V,median_V,std_V,mean_R,median_R,std_R]
		write(table,"manga_vmax_velfit_1.5.%s.csv"%vmode,column = False)

		plt.legend( fontsize = "xx-small")
		ax.set_xlabel("r (arcsec)", fontsize = 12,labelpad = 0)
		ax.set_ylabel('V$_\mathrm{ROT}$ (km/s)',fontsize=12,labelpad = 0)
		ax.set_ylim(-50,max_vel+50)
		AXIS(ax)
		plt.savefig("./vmax_rturn/%s.fit_rotcur.%s.png"%(galaxy,vmode),dpi = 300)
		#plt.show()
		plt.cla()
		return median_V,median_R

	except(IOError, OSError):
		table = [galaxy, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
		write(table,"manga_vmax_velfit_1.5.%s.csv"%vmode,column = False)
		return 1,1	
		pass
