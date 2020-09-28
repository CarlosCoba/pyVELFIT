import numpy as np
import matplotlib.pylab as plt

#
# Obtain the M_halo from the stellar mass.
#


#
# Constants valid to z<0.1
#

log_M1=12.58
M0=10**(10.90)
beta=0.48
delta=0.29
gamma=1.52

def M_vir(M_star):
	return log_M1+beta*(np.log10(M_star/M0))+((M_star/M0)**delta/(1+(M_star/M0)**(-gamma)))-0.5




x= np.linspace(7,12,100)
x_log=10**x

#plt.plot(x,M_vir(x_log))
#plt.show()

