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

    time = str(ds.current_time.v/Myr)
    time = (time[:4]) if len(time) > 4 else time
    t = "{} Myrs".format(str(time))
    
    stri = str(i).zfill(4)
    num = "spectrum_{}.png".format(stri)
    
    plt.loglog((k)**(-1), E_spectrum*(k**2.0), 'bo',label='Simulation')
   # plt.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5")
    plt.loglog((k)**(-1), (E_spectrum[5]*k[5]**2.0 + 2.e-3)*(k/k[5])**(-5./3.), ls=":", color="0.5",label = r"k$^{-5/3}$")
    plt.loglog((k)**(-1), (E_spectrum[5]*k[5]**2.0 + 2.e-3)*(k/k[5])**(-6./3.), ls="-.", color="0.5",label = r"k$^{-2}$")
   # plt.loglog(k, E_spectrum, 'bo')
   # plt.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5")
   # plt.ylim(1E-2,1E1)
    plt.ylim(1E-6,1E-3)
   # plt.xlim(1E0,3E1)
    plt.xlim(3E-2,3E0)
    plt.xlabel(r"$\lambda$ (kpc)")
    plt.ylabel(r"E(k)d$^{3}$k")
    plt.legend()
    plt.title(t)
    plt.savefig(num,format = 'png')
    plt.close()


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

for ds in ts:
    dd = ds.all_data()
   # print("Total Mass Flux above z = 2 kpc: ")
   # print(totalMassFlux.in_units('Msun/yr'))

    timeVector.append(ds.current_time.v/Myr)
    #print("Total mass in simulation box is %0.3e Msun" % \
    #      ((totalMass).in_units('Msun')))

    doit(ds,i) # does the eigenmode analysis and spits out KE vs k
    i = i+1

    # Create a 1 Mpc radius sphere, centered on the max density.
    #sp = ds.sphere("max", (4.0, "kpc"))


    # Use the total_quantity derived quantity to sum up the
    # values of the cell_mass and particle_mass fields
    # within the sphere.
    #baryon_mass = cut_dd.quantities.total_quantity(["cell_mass"])
    #print("Total mass above 2 kpc is %0.3e Msun" % \
    #      ((baryon_mass).in_units('Msun')))


growthRate = []
growthRatehat = []
growthTime = []


# Change these to get slope between different start and end points
start = 1
end = 8
#start = 100
#end = 200


p0 = np.polyfit(timeVector[start:end],log10E0vec[start:end],1)
growthRate.append(p0[0]*np.log(10))
growthTime.append(1./(p0[0]*np.log(10)))
p1 = np.polyfit(timeVector[start:end],log10E1vec[start:end],1)
growthRate.append(p1[0]*np.log(10))
growthTime.append(1./(p1[0]*np.log(10)))
p2 = np.polyfit(timeVector[start:end],log10E2vec[start:end],1)
growthRate.append(p2[0]*np.log(10))
growthTime.append(1./(p2[0]*np.log(10)))
p3 = np.polyfit(timeVector[start:end],log10E3vec[start:end],1)
growthRate.append(p3[0]*np.log(10))
growthTime.append(1./(p3[0]*np.log(10)))
p4 = np.polyfit(timeVector[start:end],log10E4vec[start:end],1)
growthRate.append(p4[0]*np.log(10))
growthTime.append(1./(p4[0]*np.log(10)))
p5 = np.polyfit(timeVector[start:end],log10E5vec[start:end],1)
growthRate.append(p5[0]*np.log(10))
growthTime.append(1./(p5[0]*np.log(10)))
p6 = np.polyfit(timeVector[start:end],log10E6vec[start:end],1)
growthRate.append(p6[0]*np.log(10))
growthTime.append(1./(p6[0]*np.log(10)))
p7 = np.polyfit(timeVector[start:end],log10E7vec[start:end],1)
growthRate.append(p7[0]*np.log(10))
growthTime.append(1./(p7[0]*np.log(10)))
p8 = np.polyfit(timeVector[start:end],log10E8vec[start:end],1)
growthRate.append(p8[0]*np.log(10))
growthTime.append(1./(p8[0]*np.log(10)))
p9 = np.polyfit(timeVector[start:end],log10E9vec[start:end],1)
growthRate.append(p9[0]*np.log(10))
growthTime.append(1./(p9[0]*np.log(10)))
p10 = np.polyfit(timeVector[start:end],log10E10vec[start:end],1)
growthRate.append(p10[0]*np.log(10))
growthTime.append(1./(p10[0]*np.log(10)))
p11 = np.polyfit(timeVector[start:end],log10E11vec[start:end],1)
growthRate.append(p11[0]*np.log(10))
growthTime.append(1./(p11[0]*np.log(10)))
p12 = np.polyfit(timeVector[start:end],log10E12vec[start:end],1)
growthRate.append(p12[0]*np.log(10))
growthTime.append(1./(p12[0]*np.log(10)))
p13 = np.polyfit(timeVector[start:end],log10E13vec[start:end],1)
growthRate.append(p13[0]*np.log(10))
growthTime.append(1./(p13[0]*np.log(10)))
p14 = np.polyfit(timeVector[start:end],log10E14vec[start:end],1)
growthRate.append(p14[0]*np.log(10))
growthTime.append(1./(p14[0]*np.log(10)))
p15 = np.polyfit(timeVector[start:end],log10E15vec[start:end],1)
growthRate.append(p15[0]*np.log(10))
growthTime.append(1./(p15[0]*np.log(10)))
p16 = np.polyfit(timeVector[start:end],log10E16vec[start:end],1)
growthRate.append(p16[0]*np.log(10))
growthTime.append(1./(p16[0]*np.log(10)))
p17 = np.polyfit(timeVector[start:end],log10E17vec[start:end],1)
growthRate.append(p17[0]*np.log(10))
growthTime.append(1./(p17[0]*np.log(10)))
p18 = np.polyfit(timeVector[start:end],log10E18vec[start:end],1)
growthRate.append(p18[0]*np.log(10))
growthTime.append(1./(p18[0]*np.log(10)))
p19 = np.polyfit(timeVector[start:end],log10E19vec[start:end],1)
growthRate.append(p19[0]*np.log(10))
growthTime.append(1./(p19[0]*np.log(10)))
p20 = np.polyfit(timeVector[start:end],log10E20vec[start:end],1)
growthRate.append(p20[0]*np.log(10))
growthTime.append(1./(p20[0]*np.log(10)))
#print("growth rate of 5.4 kpc mode:")
#print(growthRate)
#print("growth time of 5.4 kpc mode (Myrs):")
#print(growthTime)

kvec = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
# This part calculates k_hat and w_hat
 # Would need to check and change these if changing c from 0 to 1, for instance
for i in np.arange(len(kvec)):
#    kvechat.append(2*np.pi*7.714e20*kvec[i]/(8*3.0856e21))
   # growthRatehat.append(growthRate[i]*9.549e20/(1.05e6 * 3.155e13))
   # growthRatehat.append(growthRate[i]*7.714e20/(1.05e6 * 3.155e13))
    lambdaVec.append(50.0/kvec[i])

"""
plt.plot(kvec[0],growthRate[0])
plt.plot(kvec[1],growthRate[1])
plt.plot(kvec[2],growthRate[2])
plt.plot(kvec[3],growthRate[3])
plt.plot(kvec[4],growthRate[4])
plt.plot(kvec[5],growthRate[5])
plt.plot(kvec[6],growthRate[6])
plt.plot(kvec[7],growthRate[7])
plt.plot(kvec[8],growthRate[8])
plt.plot(kvec[9],growthRate[9])
plt.plot(kvec[10],growthRate[10],'ko')
plt.plot(kvec[11],growthRate[11],'ro')
plt.plot(kvec[12],growthRate[12],'bo')
plt.plot(kvec[13],growthRate[13],'go')
plt.plot(kvec[14],growthRate[14],'co')
plt.plot(kvec[15],growthRate[15],'mo')
"""
plt.plot(kvec,growthRate,'ko')
plt.xlabel("k")
plt.ylabel("Growth Rate (1/Myrs)")
plt.savefig("growthRates.png")
plt.close()

"""
plt.plot(kvechat,growthRatehat,'ko')
plt.xlabel("k")
plt.ylabel("$wH/a_{g}$")
plt.ylim(0,0.8)
plt.xlim(0,1.0)
plt.savefig("growthRateHats.png")
plt.close()


plt.plot(kvec[0],growthTime[0],'kx')
plt.plot(kvec[1],growthTime[1],'rx')
plt.plot(kvec[2],growthTime[2],'bx')
plt.plot(kvec[3],growthTime[3],'gx')
plt.plot(kvec[4],growthTime[4],'cx')
plt.plot(kvec[5],growthTime[5],'mx')
plt.plot(kvec[6],growthTime[6],color='tan')
plt.plot(kvec[7],growthTime[7],color='gray')
plt.plot(kvec[8],growthTime[8],color='lightcoral')
plt.plot(kvec[9],growthTime[9],color='orange')
plt.plot(kvec[10],growthTime[10],'ko')
plt.plot(kvec[11],growthTime[11],'ro')
plt.plot(kvec[12],growthTime[12],'bo')
plt.plot(kvec[13],growthTime[13],'go')
plt.plot(kvec[14],growthTime[14],'co')
plt.plot(kvec[15],growthTime[15],'mo')
"""
plt.semilogy(lambdaVec,growthTime,'ko')
plt.xlabel(r"$\lambda $ (pc)")
#plt.xlabel(r"$\lambda$ (kpc)")
plt.ylabel("Growth Time (Myrs)")
plt.ylim(0.01,10)
plt.savefig("growthTimes.png")
plt.close()

#print("kvechat:")
#print(kvechat)

#print("growthratehat:")
#print(growthRatehat)

#print("growthtimes:")
#print(growthTime)

ax = plt.subplot(111)
ax.semilogy(timeVector,E0vec,'k-',label=lambda0)
ax.semilogy(timeVector,E1vec,'r-',label=lambda1)
ax.semilogy(timeVector,E2vec,'b-',label=lambda2)
ax.semilogy(timeVector,E3vec,'g-',label=lambda3)
ax.semilogy(timeVector,E4vec,'c-',label=lambda4)
ax.semilogy(timeVector,E5vec,'m-',label=lambda5)
ax.semilogy(timeVector,E6vec, color='tan',label=lambda6)
ax.semilogy(timeVector,E7vec,color='gray',label=lambda7)
ax.semilogy(timeVector,E8vec,color='lightcoral',label=lambda8)
ax.semilogy(timeVector,E9vec,color='orange',label=lambda9)
ax.semilogy(timeVector,E10vec,'k-o',label=lambda10)
ax.semilogy(timeVector,E11vec,'r-o',label=lambda11)
ax.semilogy(timeVector,E12vec,'b-o',label=lambda12)
ax.semilogy(timeVector,E13vec,'g-o',label=lambda13)
ax.semilogy(timeVector,E14vec,'c-o',label=lambda14)
ax.semilogy(timeVector,E15vec,'m-o',label=lambda15)
chartBox=ax.get_position()
ax.set_position([chartBox.x0,chartBox.y0,chartBox.width*0.8,chartBox.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('Time (Myrs)')
plt.ylabel(r'$KE$(k)')
plt.ylim(1E-4, 1E0)
plt.savefig("ModeGrowth.pdf")
plt.close()
