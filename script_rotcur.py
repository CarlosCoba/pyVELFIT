import numpy as np
from califa_bisymetric import rotcur


califa_id = np.genfromtxt("../tables/get_proc_elines_CALIFA.clean.csv",usecols =  0, dtype = str, delimiter = ",")
z_Star = np.genfromtxt("../tables/get_proc_elines_CALIFA.clean.csv",usecols =  14, dtype = float, delimiter = ",")
Center = np.genfromtxt("../tables/get_proc_elines_CALIFA.clean.csv",usecols =  [167,168], dtype = float, delimiter = ",")

name = np.genfromtxt("../tables/isophotal_gsdss_califa.txt",usecols =  1, dtype = str)
P_A = np.genfromtxt("../tables/isophotal_gsdss_califa.txt",usecols =  4, dtype = float)
INCL =  np.genfromtxt("../tables/isophotal_gsdss_califa.txt",usecols =  5, dtype = float)




n = len(califa_id)
#for i in range(597,n):
for i in range(1):

	try:
		galaxy = califa_id[i]
		galaxy = "NGC3600"


		index = np.where(califa_id == galaxy)[0][0]
		index2 = np.where(name == galaxy)[0][0]

		z_star = z_Star[index]
		X0,Y0 = Center[index]
		PA = P_A[index2]
		INC = INCL[index2]
		print(PA,INC)


		#i = index
		galaxy,XK,YK,median_pa,sigma_pa,median_inc,sigma_inc,median_vsys,sigma_vsys,vmax, r_turn,beta,gamma = rotcur(galaxy,z_star,PA,INC,X0,Y0)
		print(i,galaxy,XK,YK,median_pa,sigma_pa,median_inc,sigma_inc,median_vsys,sigma_vsys,vmax, r_turn,beta,gamma)


	except(ValueError,AttributeError,TypeError,ZeroDivisionError,OSError,AssertionError, IndexError,UnboundLocalError):
	#except(1):


		try:
			z_star = z_Star[index]
			X0,Y0 = Center[index]
			PA = -P_A[index2]
			INC = INCL[index2]

			galaxy,XK,YK,median_pa,sigma_pa,median_inc,sigma_inc,median_vsys,sigma_vsys,vmax, r_turn,beta,gamma = rotcur(galaxy,z_star,PA,INC,X0,Y0)
			print(i,galaxy,XK,YK,median_pa,sigma_pa,median_inc,sigma_inc,median_vsys,sigma_vsys,vmax, r_turn,beta,gamma)


		except(AttributeError,ValueError,TypeError,OSError,AssertionError,UnboundLocalError):
		#except(1):
			print(i,galaxy,0,0,0,0,0,0,0,0,0,0,0,0)

