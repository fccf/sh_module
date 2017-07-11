import matplotlib.pyplot as plt
import numpy as np
import sh

lmax = 4


theta = np.arange(0,np.pi,np.pi/100)
x = np.cos(theta)
nn = np.size(x);

sty = ['-r','-c','--c','-g','--g','-.g','-m','--m','-.m','-..m','-b','--b','-.b','-..b','-*b']

plt.figure(figsize=(16, 12))

k = 0
for l in np.arange(0,lmax+1):
    for m in np.arange(0,l+1):
        plm = np.zeros([nn],np.float64)
        for i in np.arange(nn):
            xi = x[i]
            plm[i] = sh.sh_value.plm(l,m,xi)*sh.sh_value.alm(l,m)*np.sqrt(4.0*np.pi/(2.0*l+1.0))
        plt.plot(x,plm,sty[k],label='p'+str(l)+str(m))
        k = k+1


plt.title('Associate Legendre polynomial(normalized)')
plt.xlabel('x')
plt.ylabel('plm(x)')
plt.legend(bbox_to_anchor=(1.01,0.89),loc='best')
plt.savefig('../data/AssociateLegendrePolynomial.png')
plt.show()
