import numpy as np
import matplotlib.pylab as plt
from matplotlib.gridspec import GridSpec
from axis import AXIS 
from CBAR import colorbar as cb
#import sys
#sys.path.append("/usr/local/bin/cmaps")
import cmap_califa
import cmap_vfield
califa=cmap_vfield.CALIFA()


def plot_kin_models(galaxy,vmode,vel_ha,R,Vrot,Vrad,Vtan, VSYS, MODEL, ext,plot = 0, save = 1):

	mask_MODEL = np.divide(MODEL,MODEL)
	fig=plt.figure(figsize=(6,2))
	gs2 = GridSpec(1, 3)
	gs2.update(left=0.06, right=0.62,top=0.83,hspace=0.01,bottom=0.15,wspace=0.)


	ax=plt.subplot(gs2[0,0])
	ax1=plt.subplot(gs2[0,1])
	ax2=plt.subplot(gs2[0,2])


	im0 = ax.imshow(vel_ha - VSYS,cmap = califa, origin = "lower",vmin = -250,vmax = 250, aspect = "auto", extent = ext, interpolation = "nearest")
	im2 = ax1.imshow(MODEL,cmap = califa, origin = "lower", aspect = "auto", vmin = -250,vmax = 250, extent = ext, interpolation = "nearest")

	residual = (vel_ha*mask_MODEL - VSYS)- MODEL
	im2 = ax2.imshow(residual,cmap = califa, origin = "lower", aspect = "auto",vmin = -50,vmax = 50, extent = ext, interpolation = "nearest")


	AXIS(ax,tickscolor = "k")
	AXIS(ax1,tickscolor = "k",remove_yticks= True)
	AXIS(ax2,tickscolor = "k",remove_yticks= True)



	ax.set_ylabel(r'$\Delta$ Dec (pix)',fontsize=8,labelpad=0)
	ax.set_xlabel(r'$\Delta$ RA (pix)',fontsize=8,labelpad=0)
	ax1.set_xlabel(r'$\Delta$ RA (pix)',fontsize=8,labelpad=0)
	ax2.set_xlabel(r'$\Delta$ RA (pix)',fontsize=8,labelpad=0)

	ax.text(0.05,0.9, "VLOS", fontsize = 7, transform = ax.transAxes)
	ax1.text(0.05,0.9,"MODEL", fontsize = 7, transform = ax1.transAxes)
	ax2.text(0.05,0.9,"RESIDUAL",fontsize = 7, transform = ax2.transAxes)


	ax.set_facecolor('#e8ebf2')
	ax1.set_facecolor('#e8ebf2')
	ax2.set_facecolor('#e8ebf2')


	gs2 = GridSpec(1, 1)
	gs2.update(left=0.68, right=0.995,top=0.83,bottom=0.15)
	ax3=plt.subplot(gs2[0,0])

	ax3.plot(R,Vrot, color = "k",linestyle='-', alpha = 0.6, label = "V$_\mathrm{circ}$")
	ax3.scatter(R,Vrot, color = "k",s = 5, marker = "s")
	#ax3.errorbar(R,Vrot, yerr=e_vr, fmt='s', color = "k",markersize = 3, label = "V_\mathrm{circ}")


	ax3.plot(R,Vrad, color = "orange",linestyle='-', alpha = 0.6, label = "V$_\mathrm{rad}$")
	ax3.scatter(R,Vrad, color = "orange",s = 5, marker = "s")
	#ax3.errorbar(R,Vrad, yerr=e_Vrad, fmt='s', color = "orange",markersize = 3)


	ax3.plot(R,Vtan, color = "skyblue",linestyle='-', alpha = 0.6, label = "V$_\mathrm{tan}$")
	ax3.scatter(R,Vtan, color = "skyblue",s = 5, marker = "s")
	#ax3.errorbar(R,Vtan, yerr=e_Vtan, fmt='s', color = "skyblue",markersize = 3)

	ax3.legend(loc = "center", fontsize = 6.5, bbox_to_anchor = (0, 1, 1, 0.1), ncol = 3, frameon = False)

	vels = [Vrot, Vrad, Vtan]
	max_vel,min_vel = int(np.nanmax(vels)),int(np.nanmin(vels)) 

	ax3.set_ylim(min_vel-50,max_vel+50)
	ax3.plot([0,np.nanmax(R)],[0,0],color = "k",linestyle='-', alpha = 0.6)
	ax3.set_xlabel('r (arcsec)',fontsize=8,labelpad = 0)
	ax3.set_ylabel('V$_\mathrm{ROT}$ (km/s)',fontsize=8,labelpad = 0)
	ax3.set_facecolor('#e8ebf2')

	AXIS(ax3,tickscolor = "k")


	cb(im0,ax,orientation = "horizontal", colormap = califa, bbox= (0,1.12,1,1),width = "100%", height = "5%",label_pad = -23, label = "(km/s)",font_size=7)
	cb(im2,ax2,orientation = "horizontal", colormap = califa, bbox= (0,1.12,1,1),width = "100%", height = "5%",label_pad = -23, label = "(km/s)",font_size=7)


	if save == 1 and plot == 1:
		plt.savefig("./plots/kin_%s_model_%s.png"%(vmode,galaxy),dpi = 300)
		plt.show()
		plt.clf()
	else:

		if plot == 1:
			plt.show()
			plt.clf()
		if save == 1:
			plt.savefig("./plots/kin_%s_model_%s.png"%(vmode,galaxy),dpi = 300)
			plt.clf()




