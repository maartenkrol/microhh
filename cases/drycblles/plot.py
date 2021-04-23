from netCDF4 import Dataset
from pylab import *
with Dataset("drycblles.default.0000000.nc") as x:
    zax = x.variables['z'][:]
    z = x.groups['default']
    no = z.variables['no'][:]
    no2 = z.variables['no2'][:]
    o3 = z.variables['o3'][:]
    ho2 = z.variables['ho2'][:]
    rh = z.variables['rh'][:]
    oh = z.variables['oh'][:]
    ohe = z.variables['oh_2'][:]
    noe = z.variables['no_2'][:]
    no2e = z.variables['no2_2'][:]
    o3e = z.variables['o3_2'][:]
    ho2e = z.variables['ho2_2'][:]
    rhe = z.variables['rh_2'][:]
f,ax = subplots(1,6,sharey=True,figsize=(16,5))
ax[0].errorbar(o3[-1,:],zax, xerr=sqrt(o3e[-1,:]), fmt='.')
ax[1].errorbar(no[-1,:],zax, xerr=sqrt(noe[-1,:]), fmt='.')
ax[2].errorbar(no2[-1,:],zax, xerr=sqrt(no2e[-1,:]), fmt='.')
ax[3].errorbar(rh[-1,:],zax, xerr=sqrt(rhe[-1,:]), fmt='.')
ax[4].errorbar(ho2[-1,:],zax, xerr=sqrt(ho2e[-1,:]), fmt='.')
ax[5].errorbar(oh[-1,:],zax, xerr=sqrt(ohe[-1,:]), fmt='.')
ax[0].set_xlabel('O3 (ppb)')
ax[0].set_ylabel('z (m)')
ax[1].set_xlabel('NO (ppb)')
ax[2].set_xlabel('NO2 (ppb)')
ax[3].set_xlabel('RH (ppb)')
ax[4].set_xlabel('HO2 (ppb)')
ax[5].set_xlabel('OH (ppb)')
f.savefig('test_final.png')
