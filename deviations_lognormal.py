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
 avg_density = dd.quantities.weighted_average_quantity("density",weight="ones")
 avg_vel = dd.quantities.weighted_average_quantity("velocity_magnitude",weight="ones")
 avg_Ec = dd.quantities.weighted_average_quantity("Ec",weight="ones")
 def density_log(field, data):
   return np.log(data['density']/avg_density)

 ds.add_field(('gas', u'density_log'), function = density_log, units="")

 def vel_log(field, data):
   return np.log(data['velocity_magnitude']/avg_vel)

 ds.add_field(('gas', u'vel_log'), function = vel_log, units="")

 def Ec_log(field, data):
   return np.log(data['Ec']/avg_Ec)

 ds.add_field(('gas', u'Ec_log'), function = Ec_log, units="")

 stddev_density, avg_density = dd.quantities.weighted_variance("density_log",weight="ones")
 stddev_vel,avg_vel = dd.quantities.weighted_variance("vel_log",weight="ones")
 stddev_Ec,avg_Ec = dd.quantities.weighted_variance("Ec_log",weight="ones")
 devdensArr.append(stddev_density)
 devvelArr.append(stddev_vel)
 devEcArr.append(stddev_Ec)

print("delta rho ")
print(devdensArr)
print("delta v ")
print(devvelArr)
print("delta Ec ")
print(devEcArr)
