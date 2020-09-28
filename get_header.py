from astropy.io import fits
import sys
nargs=len(sys.argv)

if (nargs==3):
	fits_file = sys.argv[1]
	ext = int(sys.argv[2])
else: 
	print("USE: get_header.py fits_file extension")
	exit()


hdu = fits.open(fits_file)
#data = hdu[3].data
hdr = hdu[ext].header

print(hdr)


