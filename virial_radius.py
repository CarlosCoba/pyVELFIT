import numpy as np
############################Constants##########################

G = 4.299E-9    #Gravitational constant Mpc Msol**-1 (km/s)**2



H0 = 100        #Today's Hubble constant km/s/Mpc



Myear = 1E6



a_year_in_seconds = np.pi * 1E7 #units of seconds



speed_of_light = 3E10 #cm/s



Msolar = 1.99E33 #units of grams





# In[83]:





#############Bolshoi-Planck Cosmological paramters##############

Om_mat_cero = 0.307



Om_lambda_cero = 0.693



Om_baryons = 0.048



sigma8  =  0.829



h  =  0.678





# In[84]:





######################################Cosmology#######################################

def Om_m(Om_mat_cero,Om_lambda_cero,z):

    #Omega Matter as a function of redshift

    return Om_mat_cero * ( 1. + z )**3/( Om_lambda_cero + Om_mat_cero * ( 1. + z )**3 )



def Om_l(Om_mat_cero,Om_lambda_cero,z):        #Omega lambda

    #Omega Lambda as a function of redshift

    return Om_lambda_cero / ( Om_lambda_cero + Om_mat_cero * ( 1. + z )**3 )



def H(Om_mat_cero,Om_lambda_cero,z):

    #Hubble parameters as a function of redshift

    #output km/s/Mpc

    return H0 * np.sqrt( Om_lambda_cero + Om_mat_cero * ( 1. + z )**3 )



def rho_crit(Om_mat_cero,Om_lambda_cero,z):

    #Critical density as a function of redshift

    #output M_sun / Mpc**3

    return 3 * H(Om_mat_cero,Om_lambda_cero,z)**2 / 8 / np.pi /G



def rho_m(Om_mat_cero,Om_lambda_cero,z):

    return Om_m(Om_mat_cero,Om_lambda_cero,z)*rho_crit(Om_mat_cero,Om_lambda_cero,z)



def del_h(z):

    

    #def DELTA(z):

        

    x= Om_mat_cero * (1+z)**3 / (Om_mat_cero * (1+z)** 3+ Om_lambda_cero)-1

        

    return (18 * np.pi * np.pi + 82 * x - 39 * x * x) / (1+x)

    

def rho_crit(Om_mat_cero,Om_lambda_cero,z):

    #Critical density as a function of redshift to denote a flat universe

    #if this density is bigger, the universe has a positve curvature and recolapses

    #if this density is negatiev then the univerese is opened and it expands for ever

    #output M_sun / Mpc**3

    

    return 3 * H(Om_mat_cero,Om_lambda_cero,z)**2 / 8 / np.pi /G



def r_vir(Mvir, z, Om_mat_cero, Om_lambda_cero):

    #limiting radius/ virial radius: where the mean matter density is in the halo

    #input M_sun, and takes Mpc**3

    #output Kpc

    

    rho = rho_crit(Om_mat_cero,Om_lambda_cero,z) * Om_m(Om_mat_cero,Om_lambda_cero,z)

    

    d_h = del_h(z)

    

    rh = Mvir / ( 4/3 * np.pi * rho * del_h(z) )



    rh = rh **0.3333

    

    rh = rh * 1000



    return rh #* 1000 / h
