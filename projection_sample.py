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

coeffs = []
#project each vertical slice (5min averages in this case) onto vertical eigenveectors
for i in range(0,w_sub.shape[0]):
     coeffs.append(vmd.get_projection_coeffs(w_sub[i,:], eigs, N_sq, zp))

# get the time series of projection coefficients for each vertical mode 
nmodes = 4 # look at the first 9 modes 
tcomp = np.zeros((nmodes,len(coeffs)))
for i, val in enumerate(coeffs):
     for j in range(0,nmodes):
         tcomp[j,i] = val[j]


### plot the time series of projection coefficients for each mode
fig, ax = plt.subplots(1) 
ax.plot(time_sub,tcomp[0,:], label= 'Mode 1', lw =2 )
ax.plot(time_sub,tcomp[1,:], label ='Mode 2', lw =2)
ax.plot(time_sub,tcomp[2,:], label = 'Mode 3', lw =2)
ax.plot(time_sub,tcomp[3,:], label = 'Mode 4', lw =2)
ax.legend() 
fig.savefig('mode_timeseries.pdf', bbox_inches=0) 

dt = (time[1]-time[0]) # 5 minutes (time has units of days)

from scipy import signal

## kwarg 'spectrum' produces a power spectrum with 
# units of m^2/s^2; it defaults to 'density,' which
# just scales the y axis by sampling frequency
f, pgram = signal.periodogram(tcomp[0,:], 1./(dt), scaling='spectrum')
f2, pgram2 = signal.periodogram(tcomp[1,:], 1./(dt), scaling='spectrum' )
f3, pgram3 = signal.periodogram(tcomp[2,:], 1./(dt), scaling='spectrum')
f4, pgram4 = signal.periodogram(tcomp[3,:], 1./(dt), scaling='spectrum')

freq=[]
for speed in speeds: 
    freq.append(86400./(2*np.pi*128000./speed)) 

fac = (128000.*2*np.pi)/86400. #scaling factor for frequency to wavespeeds in m/s

# Plot the power spectra
fig, ax = plt.subplots(1)
ax.plot(fac*f,pgram, marker='o',lw=3, label = 'Mode 1' )
ax.plot(fac*f2,pgram2, marker='v',lw=3,label = 'Mode 2' )
ax.plot(fac*f3,pgram3, marker='x', lw=3,label = 'Mode 3' )
ax.plot(fac*f4,pgram4, marker='^', lw=3, label = 'Mode 4')

ax.legend(fontsize=14)
ax.set_xlim(0,65)
#ax.set_ylim(0,3)
ax.vlines(speeds[0:20], 0,0.45, lw=2)
ax.tick_params(labelsize=14)
ax.set_xlabel('c (m/s)', fontsize=14)
fig.savefig('power_spectra.pdf', bbox_inches=0)


#Now compute the projectiosn of each mode for the whole time series

#the normalized eigenvectors were used to compute the projection coefficients
normedeigs=[] 
for index,eig in enumerate(eigs): 
    normedeigs.append(vmd.normalize(eigs,N_sq, zp, index))
levs = len(normedeigs[0])

mode1 = np.zeros( (levs,len(tindex[0])))
mode2 = np.zeros( (levs,len(tindex[0])))
mode3 = np.zeros( (levs,len(tindex[0])))
mode4 = np.zeros( (levs,len(tindex[0])))

modes = (mode1, mode2, mode3, mode4)
for i in range(len(time_sub)):
    for mindex, m in enumerate(modes): 
        m[:,i] = tcomp[mindex,i]*normedeigs[mindex]

fig = plt.figure(figsize=(10,4))
# subplot2grid(shape, loc, rowspan=1, colspan=1)
ax1 = plt.subplot2grid((2,3), (0,0), colspan=2)
ax2 = plt.subplot2grid((2,3), (1,0), colspan=2)
ax3 = plt.subplot2grid((2,3), (0,2), rowspan=2)
im1= ax1.pcolor(time_sub-120.23, zp/1000., modes[0], vmin=-1, vmax=1)
im2 = ax2.pcolor(time_sub-120.23, zp/1000., modes[1], vmin=-1, vmax=1)
cbar1 = plt.colorbar(im1, ax=ax1)
cbar2 = plt.colorbar(im2, ax=ax2)

cbar1.ax.tick_params(labelsize=11) 
cbar2.ax.tick_params(labelsize=11) 
ax1.set_ylim(0,16.2)
ax2.set_ylim(0,16.2)
ax2.set_xlim(11.8,13.7)
ax1.set_xlim(11.8,13.7)
ax2.set_xlabel('days', fontsize=12)
for ind,m in enumerate(normedeigs[0:2]):
    ax3.plot(m,zp/1000., lw =3, label='mode '+str(ind+1))
ax3.set_ylim(0,16.2)
ax3.set_xlim(-3,1.5)
ax3.legend(loc='lower left', fontsize=12)
ax = (ax1,ax2,ax3)
for a in ax: 
    a.tick_params(labelsize=12)
    a.set_ylabel('z(km)', fontsize=12)
ax3.locator_params(axis ='x',nbins=5)
fig.savefig('modedecomp12.png',bbox_inches=0)

f = plt.figure(figsize=(10,4))
# subplot2grid(shape, loc, rowspan=1, colspan=1)
ax1 = plt.subplot2grid((2,3), (0,0), colspan=2)
ax2 = plt.subplot2grid((2,3), (1,0), colspan=2)
ax3 = plt.subplot2grid((2,3), (0,2), rowspan=2)
im1= ax1.pcolor(time_sub-120.23, z/1000., w_sub.transpose(), vmin=-1, vmax=1)
im2 = ax2.pcolor(time_sub-120.23, zp/1000., modes[0]+modes[1], vmin=-1, vmax=1)
cbar1 = plt.colorbar(im1, ax=ax1)
cbar2 = plt.colorbar(im2, ax=ax2)

cbar1.ax.tick_params(labelsize=11) 
cbar2.ax.tick_params(labelsize=11) 
ax1.set_ylim(0,16.2)
ax2.set_ylim(0,16.2)
ax2.set_xlim(11.8,13.7)
ax1.set_xlim(11.8,13.7)
ax2.set_xlabel('days', fontsize=12)
for ind,m in enumerate(normedeigs[0:2]):
    ax3.plot(m,zp/1000., lw =3, label='mode '+str(ind+1))
ax3.set_ylim(0,16.2)
ax3.set_xlim(-3,1.5)
ax3.legend(loc='lower left', fontsize=12)
ax = (ax1,ax2,ax3)
for a in ax: 
    a.tick_params(labelsize=12)
    a.set_ylabel('z(km)', fontsize=12)
ax3.locator_params(axis ='x',nbins=5)

#ax1.set_title('$w$ from simulation (top) and first 2 modes (bottom)')
f.savefig('../../Research_notes/Waves/figures/modes12_and_w.png',dpi=100, bbox_inches=0)

