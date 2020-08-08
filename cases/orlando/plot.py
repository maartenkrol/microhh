from netCDF4 import Dataset
from pylab import *
import numpy as np

slist = 'ch4,h2o2,n2o5,hald,co,hcho,isopooh,isop,xo2,mvkmacr,isopao2,no2,no3,ch3o2,isopbo2,no,ho2,o3,oh'
specs = slist.split(',')
profs = []
error = []
with Dataset("orlando.default.0000000.nc") as x:
    zax = x.variables['z'][:]
    z = x.groups['default']
    for spec in specs:
        profs.append(z.variables[spec][:])
        spece = spec+'_2'
        error.append(z.variables[spece][:])

profs = np.array(profs)
error = np.array(error)

for i,spec in enumerate(specs):
    f,ax = subplots()
    for j in range(8):
        ax.plot(profs[i,j,:], zax)
    ax.errorbar(profs[i,-1,:],zax, xerr=sqrt(error[i,-1,:]), fmt='.')
    ax.set_xlabel(spec + '(ppb)')
    ax.set_ylabel('z (m)')
    f.savefig(spec+'_32_4e2.png')
