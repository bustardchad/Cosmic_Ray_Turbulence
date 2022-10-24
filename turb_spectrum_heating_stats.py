import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.units.yt_array import YTQuantity
from yt.fields.api import ValidateParameter
import scipy.stats as stats
pUnit = YTQuantity(1, 'cm**2/s**2')
PEUnit = YTQuantity(1, 'cm/s**2')

Myr = 1.

vmax = 5.0

#Streaming energy lossesa
def heating(field,data):
          vb = (data['velocity_x']*data['magnetic_field_x'] + data['velocity_y']*data['magnetic_field_y'] + data['velocity_z']*data['magnetic_field_z'])*YTQuantity(1,"s/cm")/data['magnetic_field_magnitude']
          sigma1 = (data['alfven_speed']*YTQuantity(1,"s/cm"))*(1./(1./data['Sigma_adv1'] + 1./data['Sigma_diff1']))
          s1 = data['Fc1'] - (4./(3.))*data['Ec']*vb/vmax
          return np.abs(sigma1*s1)*8.0/(3e-5)
         # return np.abs(data['Fc1']/data['Ec'])


yt.add_field(("gas","heating"), function=heating,display_name="Streaming Energy Loss Rate",units="")


E0vec = []
E1vec = []
E2vec = []
E3vec = []
E4vec = []
E5vec = []
E6vec = []
E7vec = []
E8vec = []
E9vec = []
E10vec = []
E11vec = []
E12vec = []
E13vec = []
E14vec = []
E15vec = []
E16vec = []
E17vec = []
E18vec = []
E19vec = []
E20vec = []
log10E0vec = []
log10E1vec = []
log10E2vec = []
log10E3vec = []
log10E4vec = []
log10E5vec = []
log10E6vec = []
log10E7vec = []
log10E8vec = []
log10E9vec = []
log10E10vec = []
log10E11vec = []
log10E12vec = []
log10E13vec = []
log10E14vec = []
log10E15vec = []
log10E16vec = []
log10E17vec = []
log10E18vec = []
log10E19vec = []
log10E20vec = []
lambdaVec = []

lambda0 = "8 kpc"
lambda1 = "4 kpc"
lambda2 = "8/3 kpc"
lambda3 = "2 kpc"
lambda4 = "8/5 kpc"
lambda5 = "4/3 kpc"
lambda6 = "8/7 kpc"
lambda7 = "1 kpc"
lambda8 = "8/9 kpc"
lambda9 = "4/5 kpc"
lambda10 = "8/11 kpc"
lambda11 = "2/3 kpc"
lambda12 = "8/13 kpc"
lambda13 = "4/7 kpc"
lambda14 = "8/15 kpc"
lambda15 = "1/2 kpc"
lambda16 = "8/17 kpc"
lambda17 = "4/9 kpc"
lambda18 = "8/19 kpc"
lambda19 = "2/5 kpc"
lambda20 = "8/21 kpc"

def doit_alt(image,i):
 
     # a FFT operates on uniformly gridded data.  We'll use the yt
     # covering grid for this.
     npix = image.shape[0]
     power = 1.0
     fourier_image = np.fft.fftn(image**power)
     fourier_amplitudes = np.abs(fourier_image)**2
 
     kfreq = np.fft.fftfreq(npix) * npix
     kfreq3D = np.meshgrid(kfreq, kfreq, kfreq)
     knrm = np.sqrt(kfreq3D[0]**2 + kfreq3D[1]**2 + kfreq3D[2]**2)
 
     knrm = knrm.flatten()
     fourier_amplitudes = fourier_amplitudes.flatten()
 
     kbins = np.arange(0.5, npix//2+1, 1.)
     k = 0.5 * (kbins[1:] + kbins[:-1])
     Abins, _, _ = stats.binned_statistic(knrm, fourier_amplitudes,
                                      statistic = "mean",
                                      bins = kbins)
     Abins *= 4.*np.pi/3. * (kbins[1:]**3 - kbins[:-1]**3)
     return k, Abins

def doit(ds,i):

    # a FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this.

    max_level = ds.index.max_level
    print(max_level)
 #   ref = int(np.product(ds.ref_factors[0:max_level]))

    low = ds.domain_left_edge
    print(low)
   # dims = ds.domain_dimensions*ref
    dims = ds.domain_dimensions*int(1.0)
    print(dims)
    nx, ny, nz = dims
    
    nindex_rho = 1./2.
   # nindex_rho = 0.0

    #Kk = np.zeros( (nx//2+1,ny//2+1) )
    Kk = np.zeros( (nx//2+1, ny//2+1, nz//2+1))
   # print(nindex_rho, max_level, low, dims)
   # for vel in [("gas", "magnetic_field_x"),("gas", "magnetic_field_y")]:
   # for vel in [("gas", "magnetic_field_y")]:
   
    Kk = fft_comp(ds, ("gas","heating"), max_level, low, dims)
    # wavenumbers
    L = (ds.domain_right_edge - ds.domain_left_edge).d
    print(L)
    print(np.fft.rfftfreq(nx))
    kx = np.fft.fftfreq(nx)*nx/L[0]
    ky = np.fft.fftfreq(ny)*ny/L[1]
    kz = np.fft.fftfreq(nz)*nz/L[2]

    # physical limits to the wavenumbers
    kmin = np.min(1.0/L)
    kmax = np.min(0.5*dims/L)

    kbins = np.arange(0.5, nx//2+1, 1.)
    N = len(kbins)

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, kx, kx)
    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)
    knrm = k.flatten()
    E_spectrum, _, _ = stats.binned_statistic(knrm, Kk.flatten(),
                                      statistic = "mean",
                                      bins = kbins)
    E_spectrum *= 4.*np.pi/3. * (kbins[1:]**3 - kbins[:-1]**3)
    
  #  print("Kk:")
  #  print(Kk)
   # k = 0.5*(kbins[0:N-1] + kbins[1:N])
    k = 0.5 * (kbins[1:] + kbins[:-1])
    #k = kbins[0:N]
    print(k.shape)
    print(1/(k*0.001))
   # E_spectrum = E_spectrum[1:N]
    print(E_spectrum.shape)

    
    index = np.argmax(E_spectrum)
    kmax = k[index]
    print("Wavelength with highest energy density (in kpc): ")
    print(1.0/(kmax))
    Emax = E_spectrum[index]
    print("Emax: ")
    print(Emax)

    print("Energy spectrum: ")
    print(E_spectrum)

    return k, E_spectrum

def fft_comp(ds, iu, level, low, delta ):

    cube = ds.covering_grid(level, left_edge=low,
                            dims=delta,
                            fields=[iu])

    u = cube[iu].d

    nx, ny, nz = u.shape
    # the first half of the axes -- that's what we keep.  Our
    # normalization has an '8' to account for this clipping to one
    # octant.
   # ru = np.fft.fftn(rho**nindex_rho * u * u)[0:nx//2+1,0:ny//2+1,0]
   # ru = np.fft.fftn(rho**nindex_rho * u)[0:nx//2+1,0:ny//2+1,0]
   # ru = np.fft.fftn(rho**nindex_rho * u)[0:nx//2+1,0,0]
    ru = np.fft.fftn(u) # take sqrt so returned value has units of dE/dt x k^-1
   # ru = np.fft.rfft(rho**nindex_rho * u)[0:nx//2+1,0:ny//2+1,0]
   # print(u[0:nx//2+1,int(ny/2),0])
   # ru = np.fft.fft(rho**nindex_rho * u)[0:nx//2+1]
   # ru = np.fft.fftn(rho**nindex_rho * u)[0:nx//2+1]
   # ru = 4.0*ru/(nx*ny)
   # ru = 2.0*ru/(nx)
   # ru = 8.0*ru/(nx*ny*nz)
   # print(ru)


   # return np.abs(ru)**2  #power spectrum -- added for x and y components will give B^2
    return np.abs(ru)**2  # rho v^2
   # return np.abs(ru)  # 




ts = yt.DatasetSeries("../cr.out1.0007*")
i = 0
EkinVector = []
EkinrhoVector = []
EthermVector = []
EthermrhoVector = []
magEVector = []
timeVector = []
ratioVector = []
ratiorhoVector = []
totalEVector = []
gravPEVector = []
log10EkinrhoVector = []
teddyvec = []
ekvec = np.zeros(128)

for ds in ts:
    dd = ds.all_data()
    teddyvec.append(ds.current_time.v*3.155e13/(3.0856e21*0.666667/5e6))
    timeVector.append(ds.current_time.v/Myr)
    kvec,ekvecout = doit(ds,i)
    ekvec = np.vstack([ekvec, ekvecout]) # does the eigenmode analysis and spits out KE vs k
    i = i+1

a = ekvec[1:9,:]  # last few rows (times) of spectrum array

avg = np.mean(a,axis = 0)
mina = np.min(a,axis = 0)
maxa = np.max(a,axis = 0)


k = kvec[0:len(avg)]


print("pcr/pg = 1, beta = 1, res = 256: ")
print(avg)
print(mina)
print(maxa)
print("k: ")
print(k)

plt.loglog((k), avg*(k**2), 'bo-',label=r"$P_{CR}/P_{g} \sim 1, \beta \sim 1$")
plt.fill_between((k), mina*(k**2), maxa*(k**2), facecolor='blue', alpha=0.5)
plt.ylim(1E-5,4E-3)
#plt.xlim(3E-2,2E0)
plt.xlim(0.5,50)
#plt.xlabel(r"$\lambda / L$",fontsize=18)
plt.xlabel(r"$k$",fontsize=18)
plt.ylabel(r"FFT of CR Heating Rate",fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig("heating_spectrum_streaming_res256_pcpg1_varybeta.pdf")
plt.close()
