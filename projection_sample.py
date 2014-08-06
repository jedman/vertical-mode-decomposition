import numpy as np 
import netCDF4 as nc
import matplotlib.pyplot as plt
from vstats import * 
import normalmodes as vmd 

data= nc.Dataset('sgpoldwpg.nc')
z = getvar('z',data)
theta = getvar('theta',data)
w = getvar('w', data)
time = getvar('time',data)

tindex = np.where((time>132) & (time<134))
theta_sub = theta[tindex,:][0]
w_sub = w[tindex,:][0]
N_sq = vmd.brunt_vaisala(profile(theta),z) 
time_sub = time[tindex]

zi_full = vmd.make_zi_full(z) #make z interface vector including boundaries

# find eigenvalues and eigenvectors for lid at mean cold point tropopause (calculation note shown here) 

eigs, zp = vmd.get_vertical_modes(N_sq, zi_full, rigidlid = 16368)  

# get the phase speeds for each mode
speeds = vmd.get_phase_speeds(eigs)

#project each vertical slice (5min averages in this case) onto vertical eigenveectors
for i in range(0,w_sub.shape[0]):
     coeffs.append(vmd.get_projection_coeffs(w_sub[i,:], eigs, N_sq, zp))

# get the time series of projection coefficients for each vertical mode 
nmodes = 9 # look at the first 9 modes 
tcomp = np.zeros((nmodes,len(coeffs)))
for i, val in enumerate(coeffs):
     for j in range(0,nmodes):
         tcomp[j,i] = val[j]



# make periodogram 
from scipy import signal

f, pgram = signal.periodogram(tcomp[0,:], 1./(dt))
f2, pgram2 = signal.periodogram(tcomp[1,:], 1./(dt))
f3, pgram3 = signal.periodogram(tcomp[2,:], 1./(dt))
f4, pgram4 = signal.periodogram(tcomp[3,:], 1./(dt))

# get the expected WPG resonance frequency 
freq=[]
for speed in speeds: 
    freq.append(86400./(2*np.pi*128000./speed)) 

# plot the power spectra
fig, ax = plt.subplots(1)
ax.plot(f,pgram, marker='o',lw=3, label = 'Mode 1' )
ax.plot(f2,pgram2, marker='o',lw=3,label = 'Mode 2' )
ax.plot(f3,pgram3, marker='o', lw=3,label = 'Mode 3' )
#plot(f4,pgram4, marker='o', lw=3)

ax.legend(fontsize=14)
ax.set_xlim(0,10)
ax.set_ylim(0,0.8)
ax.vlines(freq[0:20], 0,0.8)
ax.tick_params(labelsize=14)
ax.set_xlabel('cycles per day', fontsize=14)
fig.savefig('power_spectra.pdf')




