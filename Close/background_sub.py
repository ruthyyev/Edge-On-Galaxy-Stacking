import sep
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np

#%matplotlib inline

#------------------------------------------------------------------------------
#IMPORT REQUIRED CSV DATA FILE
#------------------------------------------------------------------------------

sources = np.loadtxt('sample.csv', dtype=str, delimiter=',', skiprows = 1)

source_names = sources[:,0]
source_ra = sources[:,1]
source_dec = sources[:,2]
source_t = sources[:,3]
source_d25 = sources[:,6]
source_inclination = sources[:,8]
source_dist = sources[:,16]
source_pos_angle = sources[:,23]

#------------------------------------------------------------------------------
#IMPORT REQUIRED MAPS + BACKGROUND SUBTRACT
#------------------------------------------------------------------------------

for num in range(len(source_names)):
    name = source_names[num]
    
    for wavelength in [250, 350, 500]:
        #READ IN THE REQUIRED DATA:
        data = fits.getdata('Data/SPIRE/'+name+'_SPIRE_'+str(wavelength)+'.fits')
        
        #SET THE MEAN AND STANDARD DEVIATION IN THE DATA
        m, s = np.mean(data), np.std(data)
        
        #MAKE SURE THE DATA IS READABLE TO SEP MODULE, GET BACKGROUND:
        data = data.byteswap().newbyteorder()
        bkg = sep.Background(data)
        
        #GET A 'GLOBAL' MEAN AND NOISE OF THE IMAGE BACKGROUND
        #print(bkg.globalback)
        #print(bkg.globalrms)
        bkg_rms = bkg.rms()
        print(bkg_rms)
        
        data_sub = data - bkg.globalrms
        
        fits.writeto('Data/Background_sub/'+name+'_SPIRE_'+str(wavelength)+'_subtracted.fits', data_sub, clobber = True)
        
        
        
        

"""

rcParams['figure.figsize'] = [10., 8.]

data = fits.getdata('Data/SPIRE/NGC0891_SPIRE_250.fits')

m, s = np.mean(data), np.std(data)
#plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-s, vmax=m+s, origin='lower')
#plt.colorbar();

data = data.byteswap().newbyteorder()
bkg = sep.Background(data)

# get a "global" mean and noise of the image background:
print(bkg.globalback)
print(bkg.globalrms)

# evaluate background as 2-d array, same size as original image
bkg_image = bkg.back()
# bkg_image = np.array(bkg) # equivalent to above


# show the background
#plt.imshow(bkg_image, interpolation='nearest', cmap='gray', origin='lower')
#plt.colorbar();


# evaluate the background noise as 2-d array, same size as original image
bkg_rms = bkg.rms()

# show the background noise
#plt.imshow(bkg_rms, interpolation='nearest', cmap='gray', origin='lower')
#plt.colorbar();


data_sub = data - bkg.globalrms
plt.imshow(data_sub, cmap = 'gray')
#plt.imshow(data, cmap = 'gray')#, interpolation='nearest', cmap='gray', vmin=m-s, vmax=m+s, origin='lower')
#plt.colorbar();
plt.show()

fits.writeto('background_test.fits', data_sub, clobber = True)

"""

