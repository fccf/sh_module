"""
Plot the sparsity pattern of mesh
"""
from __future__ import division
import matplotlib.pyplot as pl
import numpy as np

a1 = np.loadtxt(file('../app/a1.arr'))
ax = np.loadtxt(file('../app/ax.arr'))
ay = np.loadtxt(file('../app/ay.arr'))
az = np.loadtxt(file('../app/az.arr'))

tp = 2
pl.spy(a1)
pl.title('A1')
pl.savefig('../data/A1-tp'+str(tp)+'.png')
pl.show()


pl.spy(ax)
pl.title('Ax')
pl.savefig('../data/Ax-tp'+str(tp)+'.png')
pl.show()


pl.spy(ay)
pl.title('Ay')
pl.savefig('../data/Ay-tp'+str(tp)+'.png')
pl.show()


pl.spy(az)
pl.title('Az')
pl.savefig('../data/Az-tp'+str(tp)+'.png')
pl.show()
