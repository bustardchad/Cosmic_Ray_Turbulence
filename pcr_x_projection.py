# plots slices of density for various times

import yt
from yt.units import dimensions
from yt.units import pc
from yt import YTQuantity
yt.enable_parallelism()

# conversion factors
edenstocgs = 6.54e-11
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.
eddy = (3.0856e21/4e6)/3.155e13


def _pcr(field, data):
  return 0.333*data['Ec']*edenstocgs*YTQuantity(1,"erg/cm**3")

yt.add_field(('gas', u'pcr'), function = _pcr, units="g/(cm*s**2)",display_name=r"P$_{cr}$", dimensions=dimensions.pressure)



plotvar = 'pcr'
#varmax = 1.e-24
varmax = 1.e-15
#varmax = 2.3e-27
varmin = 1.e-16

ts = yt.DatasetSeries('../cr.out1*',parallel=10)
for ds in ts.piter():
# time = ds.current_time.v/eddy
 time = ds.current_time.v/Myr
 slc = yt.ProjectionPlot(ds, 'x', plotvar,weight_field='ones',fontsize=20)
 slc.set_zlim(plotvar, varmin, varmax)
 slc.set_cmap(field=plotvar, cmap='dusk')
 slc.set_xlabel('y')
 slc.set_ylabel('z')
# slc.annotate_title(r"t = %3.1f $\tau_{eddy}$" % time)
 slc.annotate_title(r"t = %3.1f Myr" % time)
 slc.save()
