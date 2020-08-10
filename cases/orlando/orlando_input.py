import numpy as np
import netCDF4 as nc

float_type = "f8"
# float_type = "f4"

case = 'gcss' # Original RICO

# Get number of vertical levels and size from .ini file
with open('orlando.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

dthetadz = 0.003

# set the height
z  = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl   = np.zeros(z.size)
qt    = np.zeros(z.size)
u     = np.zeros(z.size)
ug    = np.zeros(z.size)
v     = np.zeros(z.size)
vg    = np.zeros(z.size)
wls   = np.zeros(z.size)
thlls = np.zeros(z.size)
qtls  = np.zeros(z.size)
ch4 = np.zeros(np.size(z))
h2o2 = np.zeros(np.size(z))
n2o5 = np.zeros(np.size(z))
hald = np.zeros(np.size(z))
co = np.zeros(np.size(z))
hcho = np.zeros(np.size(z))
isopooh = np.zeros(np.size(z))
isop = np.zeros(np.size(z))
xo2 = np.zeros(np.size(z))
mvkmacr = np.zeros(np.size(z))
isopao2 = np.zeros(np.size(z))
no2 = np.zeros(np.size(z))
no3 = np.zeros(np.size(z))
ch3o2 = np.zeros(np.size(z))
isopbo2 = np.zeros(np.size(z))
no = np.zeros(np.size(z))
ho2 = np.zeros(np.size(z))
o3 = np.zeros(np.size(z))

for k in range(kmax):

    # Liquid water potential temperature: same in GCSS and SS08
    if(z[k] < 740.):
        thl[k] = 297.9
    else:
        thl[k] = 297.9 + (317.0 - 297.9)/(4000. - 740.) * (z[k] - 740.) 

    if(case == 'gcss'):
        if(z[k] < 740.):
            qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
        elif(z[k] < 3260.):
            qt[k] = 13.8 + (2.4 - 13.8) / (3260. - 740.) * (z[k] - 740.) 
        else:
            qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.) 

    elif(case == 'ss08'):
        if(z[k] < 740.):
            qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
        elif(z[k] < 3260.):
            qt[k] = 13.8 + (4.4 - 13.8) / (3260. - 740.) * (z[k] - 740.) 
        else:
            qt[k] = 4.4 + (3.6 - 4.4)/(4000. - 3260.) * (z[k] - 3260.) 

    elif(case == 'test'):
        q0 = 18.
        q1 = 15.8
        if(z[k] < 740.):
            qt[k] = q0 + (q1 - q0) / 740. * z[k]
        elif(z[k] < 3260.):
            qt[k] = q1 + (2.4 - q1) / (3260. - 740.) * (z[k] - 740.) 
        else:
            qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.) 

    # Subsidence
    if(z[k] < 2260):
        wls[k] = -0.005 * (z[k] / 2260.)
    else:
        wls[k] = -0.005

    # U and V component wind
    u[k]  = -9.9 + 2.0e-3 * z[k]
    ug[k] = u[k]
    v[k]  = -3.8
    vg[k] = v[k]

    # Advective and radiative tendency thl
    thlls[k] = -2.5 / 86400.

    # Advective tendency qt
    if(z[k] < 2980):
        qtls[k] = -1.0 / 86400. + (1.3456/ 86400.) * z[k] / 2980.
    else:
        qtls[k] = 4e-6

# normalize profiles to SI
qt  /= 1000.
qtls/= 1000.

ch4[:] = 1800
h2o2[:] = 2.0
n2o5[:] = 1.0
hald[:] = 1.0
co[:] = 150.0
hcho[:] = 2.0
isopooh[:] = 1.0
isop[:] = 5.0
xo2[:] = 0.1
mvkmacr[:] = 2.0
isopao2[:] = 0.1
no2[:] = 1.0
no3[:] = 0.1
ch3o2[:] = 0.1
isopbo2[:] = 0.1
no[:] = 0.1
ho2[:] = 0.030
o3[:] = 30.0

"""
# well mixed profile with jump
h    = 1000.
dth  = 10.
dthz = 100.

for k in range(kmax):
    if(z[k] <= h - 0.5*dthz):
        th[k] = 300.
    elif(z[k] <= h + 0.5*dthz):
        th[k] = 300. + dth/dthz * (z[k]-(h-0.5*dthz))
    else:
        th[k] = 300. + dth + dthetadz*(z[k]-(h+0.5*dthz))
    thls[k] = 2.*(z[k]/zsize - 0.5) / 3600.
    wls [k] = -0.01*(z[k]/zsize)
"""

nc_file = nc.Dataset("orlando_input.nc", mode="w", datamodel="NETCDF4", clobber=False)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_thl   = nc_group_init.createVariable("thl"   , float_type, ("z"))
nc_qt    = nc_group_init.createVariable("qt"    , float_type, ("z"))
nc_u     = nc_group_init.createVariable("u"     , float_type, ("z"))
nc_ug    = nc_group_init.createVariable("ug"    , float_type, ("z"))
nc_v     = nc_group_init.createVariable("v"     , float_type, ("z"))
nc_vg    = nc_group_init.createVariable("vg"    , float_type, ("z"))
nc_wls   = nc_group_init.createVariable("w_ls"  , float_type, ("z"))
nc_thlls = nc_group_init.createVariable("thl_ls", float_type, ("z"))
nc_qtls  = nc_group_init.createVariable("qt_ls" , float_type, ("z"))

nc_z    [:] = z    [:]
nc_thl  [:] = thl  [:]
nc_qt   [:] = qt   [:]
nc_u    [:] = u    [:]
nc_ug   [:] = ug   [:]
nc_v    [:] = v    [:]
nc_vg   [:] = vg   [:]
nc_wls  [:] = wls  [:]
nc_thlls[:] = thlls[:]
nc_qtls [:] = qtls [:]

nc_ch4 = nc_group_init.createVariable("ch4", float_type, ("z"))
nc_h2o2 = nc_group_init.createVariable("h2o2", float_type, ("z"))
nc_n2o5 = nc_group_init.createVariable("n2o5", float_type, ("z"))
nc_hald = nc_group_init.createVariable("hald", float_type, ("z"))
nc_co = nc_group_init.createVariable("co", float_type, ("z"))
nc_hcho = nc_group_init.createVariable("hcho", float_type, ("z"))
nc_isopooh = nc_group_init.createVariable("isopooh", float_type, ("z"))
nc_isop = nc_group_init.createVariable("isop", float_type, ("z"))
nc_xo2 = nc_group_init.createVariable("xo2", float_type, ("z"))
nc_mvkmacr = nc_group_init.createVariable("mvkmacr", float_type, ("z"))
nc_isopao2 = nc_group_init.createVariable("isopao2", float_type, ("z"))
nc_no2 = nc_group_init.createVariable("no2", float_type, ("z"))
nc_no3 = nc_group_init.createVariable("no3", float_type, ("z"))
nc_ch3o2 = nc_group_init.createVariable("ch3o2", float_type, ("z"))
nc_isopbo2 = nc_group_init.createVariable("isopbo2", float_type, ("z"))
nc_no = nc_group_init.createVariable("no", float_type, ("z"))
nc_ho2 = nc_group_init.createVariable("ho2", float_type, ("z"))
nc_o3 = nc_group_init.createVariable("o3", float_type, ("z"))

nc_ch4[:] = ch4[:]
nc_h2o2[:] = h2o2[:]
nc_n2o5[:] = n2o5[:]
nc_hald[:] = hald[:]
nc_co[:] = co[:]
nc_hcho[:] = hcho[:]
nc_isopooh[:] = isopooh[:]
nc_isop[:] = isop[:]
nc_xo2[:] = xo2[:]
nc_mvkmacr[:] = mvkmacr[:]
nc_isopao2[:] = isopao2[:]
nc_no2[:] = no2[:]
nc_no3[:] = no3[:]
nc_ch3o2[:] = ch3o2[:]
nc_isopbo2[:] = isopbo2[:]
nc_no[:] = no[:]
nc_ho2[:] = ho2[:]
nc_o3[:] = o3[:]


nc_file.close()

ep = 287.04 / 461.5 

# Surface settings
def esat(T):
    c0 = 0.6105851e+03; c1 = 0.4440316e+02; c2 = 0.1430341e+01; c3 = 0.2641412e-01 
    c4 = 0.2995057e-03; c5 = 0.2031998e-05; c6 = 0.6936113e-08; c7 = 0.2564861e-11 
    c8 = -.3704404e-13 
    x  = max(-80.,T-273.15)
    return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

def qsat(p, T):
    return ep*esat(T)/(p-(1.-ep)*esat(T))

ps  = 101540.
SST = 299.8 
ths = SST / (ps/1.e5)**(287.04/1005.)
qs  = qsat(ps, SST) 
