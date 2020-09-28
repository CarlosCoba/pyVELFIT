import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp1d


###
# From Aldo Rodriguez-Puebla +2016,MNRAS 462, 893-916 (2016), Fig. 19a
###

log10_Mvir = np.genfromtxt("SMD_Planck_cvir_0.00000.txt",usecols=0,dtype=float,delimiter=",") 

log10_cvir= np.genfromtxt("SMD_Planck_cvir_0.00000.txt",usecols=1,dtype=float,delimiter=",") 


#Hay algo mas. En esta tabla las masas de los halos estan en terminos de h 
#asi que para usar una constante de Hubble de H0 = 70, h=0.7, tienes que hacer 
#lo siguiente: corrige las masas de halos usando que  log10 Mvir  - log10(h). 

h=0.7
log10_Mvir_cor= log10_Mvir-np.log10(h)






def concetration(log_Mh):
	x,y=log10_Mvir,10**log10_cvir
	f = interp1d(x,y)
	return f(log_Mh)
	


#plt.plot(log10_Mvir,10**log10_cvir)
#plt.show()



