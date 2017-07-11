import numpy as np
from scipy.special import sph_harm
import mayavi.mlab as mlab

pn = 5

phi, theta = np.mgrid[0:np.pi:101j, 0:2*np.pi:101j]

x = np.sin(phi)*np.cos(theta)
y = np.sin(phi)*np.sin(theta)
z = np.cos(phi)

c = np.sqrt(2)

mlab.figure(size=(1000, 1000))
for n in np.arange(pn+1):
    for m in np.arange(-n,n+1):

        if m==0:
           s = sph_harm(0,n,theta,phi).real
        if m > 0:
           s = c*sph_harm(m,n,theta,phi).real
        if m < 0:
           s = c*sph_harm(-m,n,theta,phi).imag

        mlab.mesh(abs(s)*x-1.5*n, abs(s)*y+1.5*m, abs(s)*z,scalars=s,colormap='Spectral')

mlab.view(azimuth=45.0, elevation=45.0, distance=35.0)
mlab.savefig('sph_harm.png')
mlab.show()
