import numpy as np
from scipy.special import sph_harm
import mayavi.mlab as mlab
import sh

pn = 5

theta, phi = np.mgrid[0:np.pi:101j, 0:2*np.pi:101j]

x = np.sin(theta)*np.cos(phi)
y = np.sin(theta)*np.sin(phi)
z = np.cos(theta)


mlab.figure(size=(1200, 900))
for n in np.arange(pn+1):
    for m in np.arange(-n,n+1):

        s = np.zeros([101,101],np.float64)
        for i in np.arange(101):
            t = theta[i][0]
            for j in np.arange(101):
                p = phi[0][j]
                s[i][j]= sh.sh_value.ylm(n,m,t,p)

        mlab.mesh(abs(s)*x-1.5*n, abs(s)*y+1.5*m, abs(s)*z,scalars=s, colormap='Spectral')


mlab.view(azimuth=45.0, elevation=45.0, distance=30.0)
mlab.savefig('../data/SphericalHarmonics.png')
mlab.show()
