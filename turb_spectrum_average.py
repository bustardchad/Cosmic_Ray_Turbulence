import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.units.yt_array import YTQuantity
from yt.fields.api import ValidateParameter
pUnit = YTQuantity(1, 'cm**2/s**2')
PEUnit = YTQuantity(1, 'cm/s**2')

Myr = 1.

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
    for vel in [("gas", "velocity_x"), ("gas", "velocity_y"),
                ("gas", "velocity_z")]:
   
        Kk += 0.5*fft_comp(ds, ("gas", "density"), vel,
                           nindex_rho, max_level, low, dims)

    # wavenumbers
    L = (ds.domain_right_edge - ds.domain_left_edge).d
    print(L)
    print(np.fft.rfftfreq(nx))
    kx = np.fft.rfftfreq(nx)*nx/L[0]
    ky = np.fft.rfftfreq(ny)*ny/L[1]
    kz = np.fft.rfftfreq(nz)*nz/L[2]

    # physical limits to the wavenumbers
    kmin = np.min(1.0/L)
    kmax = np.min(0.5*dims/L)

    kbins = np.arange(kmin, kmax, kmin)
    N = len(kbins)

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

   # k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)
   # k = np.sqrt(kx3d**2)
   # k = np.sqrt(kx**2)
   # print(kx)
   # print("k:")
   # print(k)
   # print("k.flat:")
   # print(k.flat)
    whichbin = np.digitize(k.flat, kbins)
   # print("whichbin:")
   # print(whichbin)
    ncount = np.bincount(whichbin)
   # print("ncount:")
   # print(ncount)
    E_spectrum = np.zeros(len(ncount)-1)

   # for n in range(1,len(ncount)):
   #     E_spectrum[n-1] = np.sum(Kk.flat[whichbin==n])
    
    for n in range(0,len(ncount)-1):
        E_spectrum[n] = np.sum(Kk.flat[whichbin==n])
    
  #  print("Kk:")
  #  print(Kk)
#    k = 0.5*(kbins[0:N-1] + kbins[1:N])
    k = kbins[0:N]
    print(k.shape)
    print(1/(k*0.001))
    E_spectrum = E_spectrum[0:N]
    print(E_spectrum.shape)

    
    index = np.argmax(E_spectrum)
    kmax = k[index]
    print("Wavelength with highest energy density (in kpc): ")
    print(1.0/(kmax))
    Emax = E_spectrum[index]
    print("Emax: ")
    print(Emax)

    E0vec.append(E_spectrum[0]) # first mode 
    E1vec.append(E_spectrum[1]) # second mode
    E2vec.append(E_spectrum[2]) # second mode
    E3vec.append(E_spectrum[3]) # second mode
    E4vec.append(E_spectrum[4]) # second mode
    E5vec.append(E_spectrum[5]) # second mode
    E6vec.append(E_spectrum[6]) # second mode
    E7vec.append(E_spectrum[7]) # second mode
    E8vec.append(E_spectrum[8]) # second mode
    E9vec.append(E_spectrum[9]) # second mode
    E10vec.append(E_spectrum[10]) # second mode
    E11vec.append(E_spectrum[11]) # second mode
    E12vec.append(E_spectrum[12]) # second mode
    E13vec.append(E_spectrum[13]) # second mode
    E14vec.append(E_spectrum[14]) # second mode
    E15vec.append(E_spectrum[15]) # second mode
    E16vec.append(E_spectrum[16]) # second mode
    E17vec.append(E_spectrum[17]) # second mode
    E18vec.append(E_spectrum[18]) # second mode
    E19vec.append(E_spectrum[19]) # second mode
    E20vec.append(E_spectrum[20]) # second mode
    log10E0vec.append(np.log10(E_spectrum[0])) # second mode
    log10E1vec.append(np.log10(E_spectrum[1])) # second mode
    log10E2vec.append(np.log10(E_spectrum[2])) # second mode
    log10E3vec.append(np.log10(E_spectrum[3])) # second mode
    log10E4vec.append(np.log10(E_spectrum[4])) # second mode
    log10E5vec.append(np.log10(E_spectrum[5])) # second mode
    log10E6vec.append(np.log10(E_spectrum[6])) # second mode
    log10E7vec.append(np.log10(E_spectrum[7])) # second mode
    log10E8vec.append(np.log10(E_spectrum[8])) # second mode
    log10E9vec.append(np.log10(E_spectrum[9])) # second mode
    log10E10vec.append(np.log10(E_spectrum[10])) # second mode
    log10E11vec.append(np.log10(E_spectrum[11])) # second mode
    log10E12vec.append(np.log10(E_spectrum[12])) # second mode
    log10E13vec.append(np.log10(E_spectrum[13])) # second mode
    log10E14vec.append(np.log10(E_spectrum[14])) # second mode
    log10E15vec.append(np.log10(E_spectrum[15])) # second mode
    log10E16vec.append(np.log10(E_spectrum[16])) # second mode
    log10E17vec.append(np.log10(E_spectrum[17])) # second mode
    log10E18vec.append(np.log10(E_spectrum[18])) # second mode
    log10E19vec.append(np.log10(E_spectrum[19])) # second mode
    log10E20vec.append(np.log10(E_spectrum[20])) # second mode
    
    lambda0 = "{} pc".format(str((k[0]*0.001)**(-1)))
    lambda1 = "{} pc".format(str((k[1]*0.001)**(-1)))
    lambda2 = "{} pc".format(str((k[2]*0.001)**(-1)))
    lambda3 = "{} pc".format(str((k[3]*0.001)**(-1)))
    lambda4 = "{} pc".format(str((k[4]*0.001)**(-1)))
    lambda5 = "{} pc".format(str((k[5]*0.001)**(-1)))
    lambda6 = "{} pc".format(str((k[6]*0.001)**(-1)))
    lambda7 = "{} pc".format(str((k[7]*0.001)**(-1)))
    lambda8 = "{} pc".format(str((k[8]*0.001)**(-1)))
    lambda9 = "{} pc".format(str((k[9]*0.001)**(-1)))
    lambda10 = "{} pc".format(str((k[10]*0.001)**(-1)))
    lambda11 = "{} pc".format(str((k[11]*0.001)**(-1)))
    lambda12 = "{} pc".format(str((k[12]*0.001)**(-1)))
    lambda13 = "{} pc".format(str((k[13]*0.001)**(-1)))
    lambda14 = "{} pc".format(str((k[14]*0.001)**(-1)))
    lambda15 = "{} pc".format(str((k[15]*0.001)**(-1)))

    time = str(ds.current_time.v*3.155e13/(3.0856e21*0.666667/5e6))
    time = (time[:4]) if len(time) > 4 else time
    #t = "{} $\tau_{eddy}$".format(str(time))
    
    stri = str(i).zfill(4)
    num = "spectrum_{}.png".format(stri)
    
    plt.loglog((k)**(-1), E_spectrum*(k**2.0), 'bo',label='Simulation')
   # plt.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5")
   # plt.loglog((k)**(-1), (E_spectrum[5]*k[5]**2.0 + 2.e-4)*(k/k[5])**(-5./3. + 2.0), ls=":", color="0.5",label = r"k$^{-5/3}$")
    plt.loglog((k)**(-1), (E_spectrum[5]*k[5]**2.0 + 2.e-4)*(k/k[5])**(-6./3. + 2.0), ls="-.", color="0.5",label = r"k$^{-2}$")
   # plt.loglog(k, E_spectrum, 'bo')
   # plt.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5")
   # plt.ylim(1E-2,1E1)
    plt.ylim(1E-6,1E-3)
   # plt.xlim(1E0,3E1)
    plt.xlim(3E-2,3E0)
    plt.xlabel(r"$\lambda / L$",fontsize=18)
    plt.ylabel(r"E(k)d$^{3}$k",fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc="lower right")
    plt.tight_layout()
    #plt.title(r"$t/ \tau_{\rm eddy}$ = )
    plt.savefig(num,format = 'png')
    plt.close()
    print("k in function:")
    print(k)
    return k, E_spectrum

def fft_comp(ds, irho, iu, nindex_rho, level, low, delta ):

    cube = ds.covering_grid(level, left_edge=low,
                            dims=delta,
                            fields=[irho, iu])

    rho = cube[irho].d
    u = cube[iu].d

    nx, ny, nz = rho.shape
    # the first half of the axes -- that's what we keep.  Our
    # normalization has an '8' to account for this clipping to one
    # octant.
   # ru = np.fft.fftn(rho**nindex_rho * u * u)[0:nx//2+1,0:ny//2+1,0]
   # ru = np.fft.fftn(rho**nindex_rho * u)[0:nx//2+1,0:ny//2+1,0]
   # ru = np.fft.fftn(rho**nindex_rho * u)[0:nx//2+1,0,0]
    ru = np.fft.fftn(rho**nindex_rho * u)[0:nx//2+1,0:ny//2+1,0:nz//2+1]
   # ru = np.fft.rfft(rho**nindex_rho * u)[0:nx//2+1,0:ny//2+1,0]
   # print(u[0:nx//2+1,int(ny/2),0])
   # ru = np.fft.fft(rho**nindex_rho * u)[0:nx//2+1]
   # ru = np.fft.fftn(rho**nindex_rho * u)[0:nx//2+1]
   # ru = 4.0*ru/(nx*ny)
   # ru = 2.0*ru/(nx)
    ru = 8.0*ru/(nx*ny*nz)
   # print(ru)


   # return np.abs(ru)**2  #power spectrum -- added for x and y components will give B^2
    return np.abs(ru)**2  # rho v^2




ts = yt.DatasetSeries("../cr.out1.00*")
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
ekvec = np.zeros(63)

for ds in ts:
    dd = ds.all_data()
    teddyvec.append(ds.current_time.v*3.155e13/(3.0856e21*0.666667/5e6))
    timeVector.append(ds.current_time.v/Myr)
    kvec,ekvecout = doit(ds,i)
    ekvec = np.vstack([ekvec, ekvecout]) # does the eigenmode analysis and spits out KE vs k
    i = i+1

a = ekvec[10:16,:]  # last few rows (times) of spectrum array
print("a: ")
print(a)

avg = np.mean(a,axis = 0)
mina = np.min(a,axis = 0)
maxa = np.max(a,axis = 0)

k = kvec[0:len(avg)]

ts2 = yt.DatasetSeries("../../res128_crbeta1/cr.out1.00*")
i = 0
teddyvec = []
ekvec2 = np.zeros(63)

for ds2 in ts2:
    dd = ds2.all_data()
    teddyvec.append(ds2.current_time.v*3.155e13/(3.0856e21*0.666667/5e6))
    timeVector.append(ds2.current_time.v/Myr)
    kvec,ekvecout = doit(ds2,i)
    ekvec2 = np.vstack([ekvec2, ekvecout]) # does the eigenmode analysis and spits out KE vs k
    i = i+1

a2 = ekvec2[10:16,:]  # last few rows (times) of spectrum array
print("a2: ")
print(a2)

avg2 = np.mean(a2,axis = 0)
mina2 = np.min(a2,axis = 0)
maxa2 = np.max(a2,axis = 0)
#print(avg)

print(avg)
print(avg2)

plt.loglog((k)**(-1), avg*(k**2), 'bo-',label=r"$P_{CR}/P_{g} \sim 10^{-4}$")
plt.fill_between((k)**(-1), mina*(k**2), maxa*(k**2), facecolor='blue', alpha=0.5)
plt.loglog((k)**(-1), avg2*(k**2), 'go-',label=r"$P_{CR}/P_{g} \sim 1$")
plt.fill_between((k)**(-1), mina2*(k**2), maxa2*(k**2), facecolor='green', alpha=0.5)
plt.loglog((k)**(-1), (avg[5]*k[5]**2.0 + 2.e-4)*(k/k[5])**(-6./3. + 2.0), ls="-.", color="0.5",label = r"k$^{-2}$")
plt.ylim(1E-6,1E-3)
plt.xlim(3E-2,2E0)
plt.xlabel(r"$\lambda / L$",fontsize=18)
plt.ylabel(r"E(k)d$^{3}$k",fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig("spectrum_composite_10_16.pdf")
plt.close()
