# this is a library for calculating the normal mode decomposition of vertical profiles of 
# atmospheric variables. 
#depends on numpy 
import numpy as np 

def brunt_vaisala(thetav_prof, z):
    '''construct N^2 = -g/thetav d/dz(thetav)  on the interface levels'''
    N_sq_prof = np.zeros(len(z)+1)
    ggr = 9.81
    thetav_int = scale2int(thetav_prof,z)
    for i in range(1,len(z)):
        N_sq_prof[i] = ggr/thetav_int[i] * (thetav_prof[i]-thetav_prof[i-1])/(z[i]-z[i-1])
    N_sq_prof[len(z)] = N_sq_prof[len(z)-1] 
    N_sq_prof[0] = N_sq_prof[1]
    return N_sq_prof 
def find_vertical_modes(N_sq, zi_full, rigidlid = []):
       '''solve the vertical structure equation d^2/dz^2 W = -\lambda N_sq W \
       and return set of eigenfunctions W_n and eigenvalues \lambda_n. 
        N_sq is the brunt-vaisala frequency, zi_full includes top and bottom interfaces;  
       kwarg *rigidlid* is the height in z of the rigid lid.''' 
    if(rigidlid):
        lid = rigidlid
    else: 
        lid = zi_full[-1]
        
    index = np.where((zi_full<lid) & (zi_full>0))
    zi = zi_full[index]
    Nsqgrd= N_sq[index]
    #set up matrix for interface levels 
    M = np.zeros((len(zi), len(zi))) 
    for i in range(1,len(zi)-1):
        M[i,i] =  -1./(zi[i+1]-zi[i]) - 1./(zi[i]-zi[i-1]) 
        M[i,i-1] = 1./(zi[i]-zi[i-1])
        M[i,i+1] = 1./(zi[i+1]-zi[i])
        M[i,:] = 1./(Nsqgrd[i])*2./(zi[i+1]-zi[i-1])*M[i,:]
        
    #boundary conditions
    M[0,0] = -1./(zi[1]-zi[0]) - 1./(zi[0]-0.)
    M[0,1] = 1./(zi[1]-zi[0]) 
    M[0,:] = 1./(Nsqgrd[0])*2./(zi[1]-0.)*M[0,:]
    M[-1,-1] = - 1./(lid - zi[-1]) - 1./(zi[-1]-zi[-2]) 
    M[-1,-2] = 1./(zi[-1]-zi[-2]) 
    M[-1,:] = 1./(Nsqgrd[-1])*2./(lid -zi[-2])*M[-1,:]
    
    c_w, Z_w = np.linalg.eig(-M) 
    return (c_w,Z_w, zi)

def innerproduct(EIGS,NSQ,ZP,val1,val2):
    ''' calculate the inner product of vertical normal modes val1 and val2 '''
    tot = EIGS[val1][1][0]*EIGS[val2][1][0]*NSQ[0]*(ZP[0])
    for i in range(1,len(ZP)):
        tot = tot + EIGS[val1][1][i]*EIGS[val2][1][i]*NSQ[i]*(ZP[i]-ZP[i-1])
    return tot

def normalize(EIGS,NSQ, ZP, val1):
    ''' ****THIS MIGHT BE THE WRONG THING TO DO****
    returns the normalized eigenfunction for val1''' 
    coeff = 1/np.sqrt(innerproduct(EIGS,NSQ,ZP, val1, val1))
    normalized_eig = coeff * EIGS[val1][1][:]
    return normalized_eig

def project_onto_vertical(var_profile, EIGS, NSQ, ZP,mode):
    '''project an observed profile onto eigenfunction *mode* '''
    eigvec=EIGS[mode][1][:]
    normed_eig = normalize(EIGS,NSQ,ZP, mode)
    eigvec = normed_eig
    coeff = eigvec[0]*var_profile[0]*NSQ[0]*(ZP[0])
    for i in range(1,len(ZP)):
        coeff = coeff + eigvec[i]*var_profile[i]*NSQ[i]*(ZP[i]-ZP[i-1])
    projection = coeff*eigvec
    return projection 

def make_zi_full(z):
    '''makes a vector of grid interfaces vector (including top and bottom)  from z on the scalar levels '''
    zi_full = np.zeros(len(z)+1)
    for i in range(1,len(z)):
        zi_full[i] = 0.5*(z[i]+z[i-1])
    zi_full[len(z)]= zi_full[len(z)-1] + (z[len(z)-1]-zi_full[len(z)-2])
    return zi_full 


