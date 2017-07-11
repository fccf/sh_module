import matplotlib.pyplot as plt
import numpy as np
import sh

lmax = 5


theta = np.arange(0,np.pi,np.pi/100)
x = np.cos(theta)
nn = np.size(x);
sty = ['-r','-c','-g','-m','-b','-k']

plt.figure(figsize=(12, 9))

for l in np.arange(0,lmax+1):
    pl = np.zeros([nn],np.float64)
    for i in np.arange(nn):
        xi = x[i]
        pl[i] = sh.sh_value.plm(l,0,xi)
    plt.plot(x,pl,sty[l],label='p'+str(l))

plt.title('Legendre polynomial')
plt.xlabel('x')
plt.ylabel('pl(x)')
plt.legend(loc='best')
plt.savefig('../data/LegendrePolynomial.png')
plt.show()
