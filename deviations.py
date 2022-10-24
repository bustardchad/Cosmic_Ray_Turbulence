# plots slices of density for various times

import yt
from yt.units import dimensions
from yt.units import pc
from yt import YTQuantity
import numpy as np
yt.enable_parallelism()

# conversion factors
edenstocgs = 6.54e-11
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.
#eddy = (3.0856e21*0.667/5e6)/3.155e13
eddy = 0.667/(0.5*0.11)

devdensArr = []
devvelArr = []
devEcArr = []
ts = yt.DatasetSeries('../cr.out1.0004*',parallel=False)
for ds in ts.piter():
 dd = ds.all_data()
 time = ds.current_time.v/eddy
 stddev_density,avg_density = dd.quantities.weighted_variance("density",weight="ones")
# stddev_density,avg_density = dd.quantities.weighted_standard_deviation("density",weight="ones")
 print(avg_density)
 print(stddev_density)
 stddev_vel,avg_vel = dd.quantities.weighted_variance("velocity_magnitude",weight="ones")
 stddev_Ec,avg_Ec = dd.quantities.weighted_variance("Ec",weight="ones")
 devdensArr.append(stddev_density/avg_density)
 devvelArr.append(stddev_vel/avg_vel)
 devEcArr.append(stddev_Ec/avg_Ec)

print("delta rho / rho")
print(devdensArr)
print("delta v / v")
print(devvelArr)
print("delta Ec / Ec")
print(devEcArr)
