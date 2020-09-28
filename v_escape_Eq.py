import numpy as np
from v_esc_phi import v_escape_halo

def vel_escape(r,v_max,R_turn,Mstar,z):#,ang_dist):

	# r must be in arcsec
	# Rturn must be in arcsec

	#r = r*ang_dist
	#R_turn = R_turn*ang_dist

	v_esc_in_square=((v_max/R_turn)**2)*(R_turn-r)**2
	v_esc_out_square= v_escape_halo(r,Mstar,z)[1]
	Rvir, v_halo= v_escape_halo(r,Mstar,z)

	# Assiming an isothermal sphere the escape velocity is
	#v_esc_out_square = 2*v_max**2*(np.log(Rvir/r)+1)


	v_escape=[]

	if r<R_turn:
		v_esc_square=v_esc_in_square + v_esc_out_square
		if np.isfinite(v_esc_square)==True:
			v_escape.append(v_esc_square)
		#else:
		#	v_escape.append(np.nan)

	else:
		v_esc_square=v_esc_out_square

		if np.isfinite(v_esc_square)==True:
			v_escape.append(v_esc_square)
		#else:
		#	v_escape.append(np.nan)

	v_escape=np.asarray(v_escape)

	return np.sqrt(v_esc_square)



