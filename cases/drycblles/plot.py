from netCDF4 import Dataset
from pylab import *
with Dataset("drycblles.default.0000000.nc") as x:
    zax = x.variables['z'][:]
    z = x.groups['default']
    no = z.variables['no'][:]
    no2 = z.variables['no2'][:]
    o3 = z.variables['o3'][:]
f,ax = subplots(1,3,sharey=True,figsize=(12,5))
for i,prof in enumerate(o3):
    ax[0].plot(prof,zax)
    ax[1].plot(no[i,:],zax)
    ax[2].plot(no2[i,:],zax)
ax[0].set_xlabel('O3 (ppb)')
ax[0].set_ylabel('z (m)')
ax[1].set_xlabel('NO (ppb)')
ax[2].set_xlabel('NO2 (ppb)')
f.savefig('test.png')
