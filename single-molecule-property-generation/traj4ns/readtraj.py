import netCDF4 as netcdf
from simtk import unit
import glob

#lentraj = 10000
#perckeep = 0.8
ncfiles = glob.glob('AlkEthOH_*.nc')

for i in ncfiles:
    #indkeep = int(lentraj*perckeep)
    data = netcdf.Dataset(i)
    xyz = data.variables['coordinates']
    #xyzn = unit.Quantity(xyz[-indkeep:], unit.angstroms)
    xyzn = unit.Quantity(xyz[:],unit.angstroms)
    if len(xyzn) > 2000:
        print i,len(xyzn)



