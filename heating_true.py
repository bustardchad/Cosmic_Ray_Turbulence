import matplotlib.pyplot as plt
import yt
from yt.units import erg, pc
from yt.units.yt_array import YTQuantity
import numpy as np
# conversion factors
edenstocgs = 6.54e-11
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.
 
yt.enable_parallelism()
 
vmax = 40.0

#Streaming energy lossesa
def heating(field,data):
          vb = (data['velocity_x']*data['magnetic_field_x'] + data['velocity_y']*data['magnetic_field_y'] + data['velocity_z']*data['magnetic_field_z'])*YTQuantity(1,"s/cm")/data['magnetic_field_magnitude']
          sigma1 = (data['alfven_speed']*YTQuantity(1,"s/cm"))*(1./(1./data['Sigma_adv1'] + 1./data['Sigma_diff1']))
          s1 = data['Fc1'] - (4./(3.))*data['Ec']*vb/vmax
          return np.abs(sigma1*s1)*8.0/(3e-5)
         # return np.abs(data['Fc1']/data['Ec'])
 
 
yt.add_field(("gas","heating"), function=heating,display_name="Streaming Energy Loss Rate",units="") 



times = []
totalCollArr = []  # array holding the collisional energy loss
totalDampArr = []  # array holding the streaming loss
totalCRLossArr = [] # array holding the total (collisional + streaming) loss
ts = yt.DatasetSeries('../cr.out1.0000*',parallel=False)
i = 0

for ds in ts.piter():
  ad = ds.all_data()
  print(ds.field_list)
  time = ds.current_time.v/Myr                       # store the time (with units removed)
  times.append(time)
 
  totalDamp = ad.quantities.weighted_average_quantity("heating",weight="ones")
  print(totalDamp)
  totalDampArr.append(totalDamp)
 
print("Collisionless")
print(totalDampArr)

print(np.mean(totalDampArr))
