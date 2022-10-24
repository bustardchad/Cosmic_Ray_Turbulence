import numpy as np
import matplotlib.pyplot as plt
# Athena++ history data
# [1]=time     [2]=dt       [3]=mass     [4]=1-mom    [5]=2-mom    [6]=3-mom    [7]=1-KE     [8]=2-KE     [9]=3-KE     [10]=tot-E   [11]=1-ME    [12]=2-ME    [13]=3-ME  
#time, ke1, ke2, ke3, tote, me1, me2, me3, ec  = np.loadtxt('../../../drive_x_4/cr.hst',skiprows = 2, usecols = (0,6,7,8,9,10,11,12,13),unpack=True)

#time_noStream, ke1_noStream, ke2_noStream, ke3_noStream, tote_noStream, me1_noStream, me2_noStream, me3_noStream, ec_noStream  = np.loadtxt('noStream/cr.hst',skiprows = 2, usecols = (0,6,7,8,9,10,11,12,13),unpack=True)

time_diff3e24, ke1_diff3e24, ke2_diff3e24, ke3_diff3e24, tote_diff3e24, me1_diff3e24, me2_diff3e24, me3_diff3e24, ec_diff3e24  = np.loadtxt('diff3e24/cr.hst',skiprows = 2, usecols = (0,6,7,8,9,10,11,12,13),unpack=True)

time_diff3e25, ke1_diff3e25, ke2_diff3e25, ke3_diff3e25, tote_diff3e25, me1_diff3e25, me2_diff3e25, me3_diff3e25, ec_diff3e25  = np.loadtxt('diff3e25/cr.hst',skiprows = 2, usecols = (0,6,7,8,9,10,11,12,13),unpack=True)

time_diff3e26, ke1_diff3e26, ke2_diff3e26, ke3_diff3e26, tote_diff3e26, me1_diff3e26, me2_diff3e26, me3_diff3e26, ec_diff3e26  = np.loadtxt('diff3e26/cr.hst',skiprows = 2, usecols = (0,6,7,8,9,10,11,12,13),unpack=True)

time_diff3e27, ke1_diff3e27, ke2_diff3e27, ke3_diff3e27, tote_diff3e27, me1_diff3e27, me2_diff3e27, me3_diff3e27, ec_diff3e27  = np.loadtxt('diff3e27/cr.hst',skiprows = 2, usecols = (0,6,7,8,9,10,11,12,13),unpack=True)

#time_diff3e28, ke1_diff3e28, ke2_diff3e28, ke3_diff3e28, tote_diff3e28, me1_diff3e28, me2_diff3e28, me3_diff3e28, ec_diff3e28  = np.loadtxt('diff3e28/cr.hst',skiprows = 2, usecols = (0,6,7,8,9,10,11,12,13),unpack=True)

edenstocgs = 6.54e-11
#cellVol = (250*3.0856e18)**3.0
cellVol = 0.05**3.0
ketot = []
metot = []
thermaletot = []
ectot = []
ketot_noStream = []
metot_noStream = []
thermaletot_noStream = []
ectot_noStream = []
ketot_diff3e24 = []
metot_diff3e24 = []
thermaletot_diff3e24 = []
ectot_diff3e24 = []

ketot_diff3e25 = []
metot_diff3e25 = []
thermaletot_diff3e25 = []
ectot_diff3e25 = []

ketot_diff3e26 = []
metot_diff3e26 = []
thermaletot_diff3e26 = []
ectot_diff3e26 = []

ketot_diff3e27 = []
metot_diff3e27 = []
thermaletot_diff3e27 = []
ectot_diff3e27 = []

ketot_diff3e28 = []
metot_diff3e28 = []
thermaletot_diff3e28 = []
ectot_diff3e28 = []
"""
for j in range(len(time)):
  keval = np.sqrt(ke1[j]**2 + ke2[j]**2 + ke3[j]**2)*edenstocgs/cellVol
  meval = np.sqrt(me1[j]**2 + me2[j]**2 + me3[j]**2)*edenstocgs/cellVol
  thermalval = (tote[j]*edenstocgs/cellVol) - keval - meval
  ecval = ec[j]*edenstocgs/cellVol
  ketot.append(keval)
  metot.append(meval)
  thermaletot.append(thermalval)
  ectot.append(ecval)

for j in range(len(time_noStream)):
  keval_noStream = np.sqrt(ke1_noStream[j]**2 + ke2_noStream[j]**2 + ke3_noStream[j]**2)*edenstocgs/cellVol
  ketot_noStream.append(keval_noStream)
  meval_noStream = np.sqrt(me1_noStream[j]**2 + me2_noStream[j]**2 + me3_noStream[j]**2)*edenstocgs/cellVol
  thermalval_noStream = (tote_noStream[j]*edenstocgs/cellVol) - keval_noStream - meval_noStream
  ecval_noStream = ec_noStream[j]*edenstocgs/cellVol
  metot_noStream.append(meval_noStream)
  thermaletot_noStream.append(thermalval_noStream)
  ectot_noStream.append(ecval_noStream)
""" 
for j in range(len(time_diff3e24)):
  keval_diff3e24 = np.sqrt(ke1_diff3e24[j]**2 + ke2_diff3e24[j]**2 + ke3_diff3e24[j]**2)*edenstocgs/cellVol
  ketot_diff3e24.append(keval_diff3e24)
  meval_diff3e24 = np.sqrt(me1_diff3e24[j]**2 + me2_diff3e24[j]**2 + me3_diff3e24[j]**2)*edenstocgs/cellVol
  thermalval_diff3e24 = (tote_diff3e24[j]*edenstocgs/cellVol) - keval_diff3e24 - meval_diff3e24
  ecval_diff3e24 = ec_diff3e24[j]*edenstocgs/cellVol
  metot_diff3e24.append(meval_diff3e24)
  thermaletot_diff3e24.append(thermalval_diff3e24)
  ectot_diff3e24.append(ecval_diff3e24)

for j in range(len(time_diff3e25)):
  keval_diff3e25 = np.sqrt(ke1_diff3e25[j]**2 + ke2_diff3e25[j]**2 + ke3_diff3e25[j]**2)*edenstocgs/cellVol
  ketot_diff3e25.append(keval_diff3e25)
  meval_diff3e25 = np.sqrt(me1_diff3e25[j]**2 + me2_diff3e25[j]**2 + me3_diff3e25[j]**2)*edenstocgs/cellVol
  thermalval_diff3e25 = (tote_diff3e25[j]*edenstocgs/cellVol) - keval_diff3e25 - meval_diff3e25
  ecval_diff3e25 = ec_diff3e25[j]*edenstocgs/cellVol
  metot_diff3e25.append(meval_diff3e25)
  thermaletot_diff3e25.append(thermalval_diff3e25)
  ectot_diff3e25.append(ecval_diff3e25)
 
for j in range(len(time_diff3e26)):
  keval_diff3e26 = np.sqrt(ke1_diff3e26[j]**2 + ke2_diff3e26[j]**2 + ke3_diff3e26[j]**2)*edenstocgs/cellVol
  ketot_diff3e26.append(keval_diff3e26)
  meval_diff3e26 = np.sqrt(me1_diff3e26[j]**2 + me2_diff3e26[j]**2 + me3_diff3e26[j]**2)*edenstocgs/cellVol
  thermalval_diff3e26 = (tote_diff3e26[j]*edenstocgs/cellVol) - keval_diff3e26 - meval_diff3e26
  ecval_diff3e26 = ec_diff3e26[j]*edenstocgs/cellVol
  metot_diff3e26.append(meval_diff3e26)
  thermaletot_diff3e26.append(thermalval_diff3e26)
  ectot_diff3e26.append(ecval_diff3e26)

for j in range(len(time_diff3e27)):
  keval_diff3e27 = np.sqrt(ke1_diff3e27[j]**2 + ke2_diff3e27[j]**2 + ke3_diff3e27[j]**2)*edenstocgs/cellVol
  ketot_diff3e27.append(keval_diff3e27)
  meval_diff3e27 = np.sqrt(me1_diff3e27[j]**2 + me2_diff3e27[j]**2 + me3_diff3e27[j]**2)*edenstocgs/cellVol
  thermalval_diff3e27 = (tote_diff3e27[j]*edenstocgs/cellVol) - keval_diff3e27 - meval_diff3e27
  ecval_diff3e27 = ec_diff3e27[j]*edenstocgs/cellVol
  metot_diff3e27.append(meval_diff3e27)
  thermaletot_diff3e27.append(thermalval_diff3e27)
  ectot_diff3e27.append(ecval_diff3e27)
"""
for j in range(len(time_diff3e28)):
  keval_diff3e28 = np.sqrt(ke1_diff3e28[j]**2 + ke2_diff3e28[j]**2 + ke3_diff3e28[j]**2)*edenstocgs/cellVol
  ketot_diff3e28.append(keval_diff3e28)
  meval_diff3e28 = np.sqrt(me1_diff3e28[j]**2 + me2_diff3e28[j]**2 + me3_diff3e28[j]**2)*edenstocgs/cellVol
  thermalval_diff3e28 = (tote_diff3e28[j]*edenstocgs/cellVol) - keval_diff3e28 - meval_diff3e28
  ecval_diff3e28 = ec_diff3e28[j]*edenstocgs/cellVol
  metot_diff3e28.append(meval_diff3e28)
  thermaletot_diff3e28.append(thermalval_diff3e28)
  ectot_diff3e28.append(ecval_diff3e28)
"""

#plt.plot(time,thermaletot,'k-',label = r'E$_{th}$',linewidth=3)
#plt.plot(time,ketot,label = r'Streaming',linewidth=1)
#plt.plot(time_noStream,ketot_noStream,label = r'Advection',linewidth=1)
plt.plot(time_diff3e24,ketot_diff3e24,label = r'$\kappa_{||} = 3$ x $10^{24}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e25,ketot_diff3e25,label = r'$\kappa_{||} = 3$ x $10^{25}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e26,ketot_diff3e26,label = r'$\kappa_{||} = 3$ x $10^{26}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e27,ketot_diff3e27,label = r'$\kappa_{||} = 3$ x $10^{27}$ cm$^{2}$/s',linewidth=1)
#plt.plot(time_diff3e28,ketot_diff3e28,label = r'$\kappa_{||} = 3$ x $10^{28}$ cm$^{2}$/s',linewidth=1)
#plt.plot(time,metot,'g-',label = r'E$_{B}$',linewidth=3)
#plt.plot(time,ectot,'m-',label = r'E$_{CR}$',linewidth=3)
plt.xlabel('Time (Myrs)',fontsize=18)
plt.ylabel('KE Energy Density',fontsize=18)
plt.legend()
plt.savefig('KE_EnergyDensity_Compare.pdf')
plt.close()

#plt.plot(time,thermaletot,label = r'Streaming',linewidth=1)
#plt.plot(time_noStream,thermaletot_noStream,label = r'Advection',linewidth=1)
plt.plot(time_diff3e24,thermaletot_diff3e24,label = r'$\kappa_{||} = 3$ x $10^{24}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e25,thermaletot_diff3e25,label = r'$\kappa_{||} = 3$ x $10^{25}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e26,thermaletot_diff3e26,label = r'$\kappa_{||} = 3$ x $10^{26}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e27,thermaletot_diff3e27,label = r'$\kappa_{||} = 3$ x $10^{27}$ cm$^{2}$/s',linewidth=1)
#plt.plot(time_diff3e28,thermaletot_diff3e28,label = r'$\kappa_{||} = 3$ x $10^{28}$ cm$^{2}$/s',linewidth=1)
#plt.plot(time,metot,'`g-',label = r'E$_{B}$',linewidth=3)
#plt.plot(time,ectot,'m-',label = r'E$_{CR}$',linewidth=3)
plt.xlabel('Time (Myrs)',fontsize=18)
plt.ylabel('Thermal Energy Density',fontsize=18)
plt.legend()
plt.savefig('Thermal_EnergyDensity_Compare.pdf')
plt.close()

#plt.plot(time,ectot,label = r'Streaming',linewidth=1)
#plt.plot(time_noStream,ectot_noStream,label = r'Advection',linewidth=1)
plt.plot(time_diff3e24,ectot_diff3e24,label = r'$\kappa_{||} = 3$ x $10^{24}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e25,ectot_diff3e25,label = r'$\kappa_{||} = 3$ x $10^{25}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e26,ectot_diff3e26,label = r'$\kappa_{||} = 3$ x $10^{26}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e27,ectot_diff3e27,label = r'$\kappa_{||} = 3$ x $10^{27}$ cm$^{2}$/s',linewidth=1)
#plt.plot(time_diff3e28,ectot_diff3e28,label = r'$\kappa_{||} = 3$ x $10^{28}$ cm$^{2}$/s',linewidth=1)
#plt.plot(time,metot,'g-',label = r'E$_{B}$',linewidth=3)
#plt.plot(time,ectot,'m-',label = r'E$_{CR}$',linewidth=3)
plt.xlabel('Time (Myrs)',fontsize=18)
plt.ylabel('CR Energy Density',fontsize=18)
plt.legend()
plt.savefig('CR_EnergyDensity_Compare.pdf')
plt.close()

#plt.plot(time,metot,label = r'Streaming',linewidth=1)
#plt.plot(time_noStream,metot_noStream,label = r'Advection',linewidth=1)
plt.plot(time_diff3e24,metot_diff3e24,label = r'$\kappa_{||} = 3$ x $10^{24}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e25,metot_diff3e25,label = r'$\kappa_{||} = 3$ x $10^{25}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e26,metot_diff3e26,label = r'$\kappa_{||} = 3$ x $10^{26}$ cm$^{2}$/s',linewidth=1)
plt.plot(time_diff3e27,metot_diff3e27,label = r'$\kappa_{||} = 3$ x $10^{27}$ cm$^{2}$/s',linewidth=1)
#plt.plot(time_diff3e28,metot_diff3e28,label = r'$\kappa_{||} = 3$ x $10^{28}$ cm$^{2}$/s',linewidth=1)
#plt.plot(time,metot,'g-',label = r'E$_{B}$',linewidth=3)
#plt.plot(time,ectot,'m-',label = r'E$_{CR}$',linewidth=3)
plt.xlabel('Time (Myrs)',fontsize=18)
plt.ylabel('Magnetic Energy Density',fontsize=18)
plt.legend()
plt.savefig('Magnetic_EnergyDensity_Compare.pdf')
plt.close()
