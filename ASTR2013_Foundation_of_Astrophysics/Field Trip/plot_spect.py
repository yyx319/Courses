import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
plt.ion()
infile = 'astro_course/20190402/without_11_12/b/T2m3wb-20190902.110619-0021.p11.fits'
infile = 'astro_course/20190402/without_11_12/r/T2m3wr-20190902.114636-0024.p11.fits'

#Read in data, and find wavelength scale.
cube = pyfits.getdata(infile)
header = pyfits.getheader(infile)
wave = header['CRVAL3'] + np.arange(cube.shape[0])*header['CDELT3']

#Create an image...
cube_skysub = cube.copy()
for wave_ix in range(cube.shape[0]):
    #Subtract off the median of the cube, excluding the 
    cube_skysub[wave_ix] -= np.median(cube_skysub[wave_ix,1:,1:])

#Now display the image.
im_skysub = np.median(cube_skysub, axis=0)
plt.figure(1)
plt.clf()
plt.imshow(im_skysub, vmin=0, vmax=np.percentile(im_skysub,96))

#Display a spectrum of a single spaxel
plt.figure(2)
plt.clf()
plt.plot(wave, cube_skysub[:,14,13])
