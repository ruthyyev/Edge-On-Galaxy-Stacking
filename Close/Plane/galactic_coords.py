#------------------------------------------------------------------------------
#IMPORT REQUIRED MODULES
#------------------------------------------------------------------------------
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
from numpy import genfromtxt
from astropy.table import Table
import numpy as np
from matplotlib import pyplot as plt
import pylab as pl
import math
import seaborn as sns


sns.set_style("whitegrid")

#os.chdir('/Users/c1541417/Documents/Edge-On2/')

#------------------------------------------------------------------------------
#READ IN THE DATA
#------------------------------------------------------------------------------
ra1 = np.loadtxt('/Users/c1541417/Documents/Edge-On2/Plane/sample_ra.csv')
dec1 = np.loadtxt('/Users/c1541417/Documents/Edge-On2/Plane/sample_dec.csv')


x = SkyCoord(ra1*u.degree, dec1*u.degree, frame='icrs')

#convert from RA in degrees to galactic coords:
y = x.transform_to(Galactic)

y.l.wrap_at(180*u.degree, inplace=True)

#------------------------------------------------------------------------------
#READ AND WRITE THE TABLE DATA
#------------------------------------------------------------------------------

t = Table([y.l, y.b], names=('L', 'B'), meta={'name': 'first table'})
t.write('/Users/c1541417/Documents/Edge-On2/Plane/galactic_coords_no_plane.csv', format='csv')

#------------------------------------------------------------------------------
#IMPORT TABLE DATA
#------------------------------------------------------------------------------
"""Where L is galactic longitude and B is galactic latitude"""

z = genfromtxt('/Users/c1541417/Documents/Edge-On2/Plane/galactic_coords_no_plane.csv', delimiter=',')
t = Table(z, names = ('L', 'B'))
z = t['L'][0:]
y = t['B'][0:]

#wrap-around galactic longitude vales:
x = z-180

#------------------------------------------------------------------------------
#PLOT THE CELESTIAL EQUATOR IN GALACTIC COORDINATES
#------------------------------------------------------------------------------
degtorad = math.pi/180.
alpha = np.arange(-180,180.,1.)
alpha *= degtorad

#------------------------------------------------------------------------------
#SCIENCE STUFFS
#------------------------------------------------------------------------------
#From Meeus, Astronomical algorithms (with delta = 0):
x1 = np.sin(192.25*degtorad - alpha)
x2 = np.cos(192.25*degtorad - alpha)*np.sin(27.4*degtorad)
yy = np.arctan2(x1, x2)
longitude = 303*degtorad - yy 
x3 = np.cos(27.4*degtorad) * np.cos(192.25*degtorad - alpha)
latitude  = np.arcsin(x3)

#------------------------------------------------------------------------------
#PUT THE ANGLES IN THE CORRECT DIRECTION
#------------------------------------------------------------------------------
for i in range(0,len(alpha)):
    if longitude[i] > 2.*math.pi:
        longitude[i] -= 2.*math.pi
    longitude[i] -= math.pi
    latitude[i] = -latitude[i]

#avoid a line in the middle of the plot (the curve must not loop):
for i in range(0,len(longitude)-1):
    if (longitude[i] * longitude[i+1] < 0 and longitude[i] > 170*degtorad and longitude[i+1] < -170.*degtorad):
        indice = i
        break

#array is put in increasing longitude:
longitude2 = np.zeros(len(longitude))
latitude2 = np.zeros(len(latitude))
longitude2[0:len(longitude)-1-indice] = longitude[indice+1:len(longitude)]
longitude2[len(longitude)-indice-1:len(longitude)] = longitude[0:indice+1]
latitude2[0:len(longitude)-1-indice] = latitude[indice+1:len(longitude)]
latitude2[len(longitude)-indice-1:len(longitude)] = latitude[0:indice+1]

xrad = x * degtorad
yrad = y * degtorad

#------------------------------------------------------------------------------
#PLOT THE DATA
#------------------------------------------------------------------------------

fig2 = pl.figure(2)
ax1 = fig2.add_subplot(111, projection="mollweide")

ax1.scatter(xrad,yrad, c = 'darkslategray')
ax1.plot([-math.pi, math.pi], [0,0],'cadetblue')
ax1.plot([0,0],[-math.pi, math.pi], 'cadetblue')

#plot celestial equator:
plt.plot(longitude2,latitude2, c = 'orangered')

plt.title("Edge-on Late-Type Sources")
plt.grid(True)


#plt.show()

pl.savefig('/Users/c1541417/Documents/Edge-On2/Plane/high_inc_sample_sky.pdf')
