import sep
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import progressbar
import montage_wrapper
import reproject
from reproject import reproject_exact
import os
from astropy.nddata.utils import Cutout2D

os.chdir('/Users/c1541417/Documents/Edge-On2')
#np.set_printoptions(threshold=np.inf)

bar = progressbar.ProgressBar()

#------------------------------------------------------------------------------
#IMPORT REQUIRED CSV DATA FILE
#------------------------------------------------------------------------------
#Here we want to import the csv file containing information on all sources, 
#specifying individual data columns


sources = np.loadtxt('sample.csv', dtype=str, delimiter=',', skiprows = 1)

source_names = sources[:,0]
source_ra = sources[:,1]
source_dec = sources[:,2]
source_dist = sources[:,16]
source_pos_angle = sources[:,22]

#convert imported strings to floats:
dist_list = source_dist.astype(np.float)
pos_angle1 = source_pos_angle.astype(np.float)

#correct for position angle:
pos_angle = [360. - x for x in pos_angle1]

#width is the size of the image that montage outputs (degrees):
width = 5.


x_position = []
y_position = []
cutouts = []

#------------------------------------------------------------------------------
#IMPORT REQUIRED MAPS + BACKGROUND SUBTRACT
#------------------------------------------------------------------------------
#Here we will use the 'sep' python package to subtract the background from 
#each maps, and output the new maps with their old headers


for num in bar(range(len(source_names))):
    name = source_names[num]
    ra = source_ra[num]
    dec = source_dec[num]
    angle = pos_angle[num]
    
    for wavelength in [250, 350, 500]:
    
          
      
        #read in the required data:
        data = fits.getdata('Data/SPIRE/'+name+'_SPIRE_'+str(wavelength)+'.fits')
        header = fits.getheader('Data/SPIRE/'+name+'_SPIRE_'+str(wavelength)+'.fits')
        
        #set the mean and standard deviation in the data:
        m, s = np.mean(data), np.std(data)
        
        #make sure the data is readable to the sep module, get background:
        data = data.byteswap().newbyteorder()
        bkg = sep.Background(data)
        
        #get a 'global' mean and noise of the data, remove from the data:
        bkg_rms = bkg.rms()
        data_sub = data - bkg.globalrms
        
        #write the background subracted maps out to new fits files (ACTIVATE):
        fits.writeto('Data/Background_sub/'+name+'_SPIRE_'+str(wavelength)+'_subtracted.fits', data_sub, header = header, clobber = True)
        
    
#------------------------------------------------------------------------------
#REGRID THE IMAGES SO THEY ARE THE SAME KPC PER PIXEL
#------------------------------------------------------------------------------
#Here we will convert each pixel to the same numer of kiloparsecs, so that 
#individual galaxies can be compared. We will use montage_wrapper to re-grid 
#and reproject the maps

        pix_size = header['CDELT2']
        pix_size_arcsec = pix_size * 3600. #converts pixel size from deg to "
        lin_dist = ((dist_list * pix_size_arcsec) / 206265. ) * 1000.
        #find the largest scale in sample, here it is lin_dist[31]
    
        
#3,19,31,60,63,65,73 are furthest away = worse res.
#at this point, should we remove furthest away objects? test the difference!!


        max_lin_dist = lin_dist[31] #the furthest away galaxy (with most kpc per pixel)
        #print max_lin_dist
        scale_factor = max_lin_dist / lin_dist #find ratio to change by
        new_pix_size = pix_size_arcsec * scale_factor #new pixel size for images!
        #print(new_pix_size)

  
        montage_wrapper.mHdr(str(ra)+', '+str(dec), width, 'Data/Rescaled/'+name+'_SPIRE_'+str(wavelength)+'_rescaled_kpc.txt', pix_size = new_pix_size[num], rotation = angle)

        montage_wrapper.wrappers.reproject('Data/Background_sub/'+name+'_SPIRE_'+str(wavelength)+'_subtracted.fits', 'Data/Rescaled/'+name+'_SPIRE_'+str(wavelength)+'_rescaled_kpc.fits', header='Data/Rescaled/'+name+'_SPIRE_'+str(wavelength)+'_rescaled_kpc.txt', silent_cleanup = True)
        
   
#-------------------------------------------------------------------------------
#MAKE CUT-OUTS OF THE DATA - THE GALAXY AT THE CENTRE
#-------------------------------------------------------------------------------



    for wavelength in [250]: #change according to wavelength! need to sort this out...
        data_cut = fits.open('Data/Rescaled/'+name+'_SPIRE_'+str(wavelength)+'_rescaled_kpc.fits')

    #get the galaxy central pixels from the header files
    x_coords = data_cut[0].header['CRPIX1']
    y_coords = data_cut[0].header['CRPIX2']

    #append information to empty lists
    x_position.append(x_coords)
    y_position.append(y_coords)

    x_pos = x_position[num]
    y_pos = y_position[num]

    fin_data = data_cut[0].data     



    #make the cutouts
    cutouts = Cutout2D(fin_data, (x_pos, y_pos), (60, 60), mode = 'partial', fill_value = 0.)

    #save the resulting cutouts as fits files
    fits.writeto('Data/Cutouts/'+name+'_'+str(wavelength)+'_cutout.fits', cutouts.data,  clobber = True)
      
#WORKS UP TO HERE 17:30 21/03/17

#-------------------------------------------------------------------------------
#STACK THE IMAGES
#-------------------------------------------------------------------------------

outputs = []

for wavelength in [250, 350, 500]:

    image_list = []

    for name in source_names:

        try:
            image_list.append(fits.getdata('Data/Cutouts/'+name+'_'+str(wavelength)+'_cutout.fits'))
        except:
            print 'could not find '+name+''

    image_list = [i for i in image_list if i.shape == (60, 60)]
    image_list = [np.ma.masked_invalid(i) for i in image_list]

    try:
        outputs.append(np.ma.sum(image_list, axis=0))
    except:
        pass

    plt.figure()
    plt.imshow(outputs[-1])

    outfile = 'Results/'+str(wavelength)+'_cutout_stacked.fits'
    fits.writeto(outfile, outputs[-1].data, clobber = True)


#all working, next step is to test the script with furthest galaxies removed
#to see how much better the resulting maps look
#also! - sort out header files, i.e. use mHdr/project? 