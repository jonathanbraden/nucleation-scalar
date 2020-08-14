import numpy as np
import matplotlib.pyplot as plt

fInd = { 'phi'        : 0,
         'phidot'     : 1,
         'G^2 (FD)'   : 2,
         'V(phi)'     : 3,
         'G^2 (SP)'   : 4,
         'm^2phi^2'   : 5
}

def load_fields(fName,nLat):
    return np.genfromtxt(fName).reshape((-1,nLat,6))

def compute_rho(dat):
    return 0.5*dat[:,:,fInd['phidot']]**2 + 0.5*dat[:,:,fInd['G^2 (SP)']] + dat[:,:,fInd['V(phi)']]

def compute_cos(dat):
    return np.cos(dat[:,:,fInd['phi']])

def compute_smooth_cos(dat,kcut):
    return np.cos(project_k_modes(dat[:,:,fInd['phi']],kcut))

def compute_rho_from_field(fld,dfld,lSize,pot,spec=True):
    """
    Compute energy density, including the derivatives from input field.
    
    Assumed shape of fld is (num_tsteps,num_lat,2)
    fld[:,:,0] : the field values
    fld[:,:,1] : the field derivative values

    lSize - Side length of simulation

    pot - Potential function

    spec - If True compute spectral (grad phi)^2
         - If False compute finite difference (grad phi)^2
    """
    if spec:
        return 0.5*dfld[:,:]**2 + 0.5*fourier_derivative(fld[:,:],lSize)**2 + pot(fld[:,:])
    else:
        return 0.5*dfld[:,:]**2 + 0.5*grad_squared_finite_diff(fld[:,:],lSize)**2 + pot(fld[:,:])

def fourier_laplacian(fld,lSize):
    nLat = fld.shape[1]; norm = 2.*np.pi*nLat/lSize
    fk = -np.fft.rfft(fld,axis=-1)*np.fft.rfftfreq(nLat)**2*norm**2
    return np.fft.irfft(fk,n=nLat)

# This has some small errors relative to my Fortran code.  Figure out if it's a problem or not.
def fourier_derivative(fld,lSize):
    """
    Compute gradient of fld on a grid of size lSize.

    fld is assumed to be of shape (num_tsteps,num_lat)
    """
    nLat = fld.shape[1]; norm = 2.*np.pi*nLat/lSize
    fk = 1j*np.fft.rfft(fld,axis=-1)*np.fft.rfftfreq(nLat)*norm
    return np.fft.irfft(fk,n=nLat)

def grad_squared_finite_diff(fld,lSize):
    """
    Compute (grad phi)^2 using a centered finite-difference approximation

    fld is assumed to be of shape (num_tsteps,num_lat)
    lSize - Side length of simulation volume
    """
    nLat = fld.shape[1]; dx = lSize/np.float64(nLat)
    ind_right = ( np.arange(nLat) + 1) % nLat
    ind_left = ( np.arange(nLat) - 1) % nLat
    df2 = (fld[:,ind_right]-fld[:,:])**2 + (fld[:,ind_left]-fld[:,:])**2
    return 0.5*df2 / dx**2
    
def project_k_modes(fld,kcut):
    """
    Smooth the input field by projecting modes at given kcut.
    Assumes fld is of shape (num_time,num_lat)
    """
    nLat = fld.shape[1]
    fk = np.fft.rfft(fld,axis=-1)
    fk[:,kcut:] = 0.
    return np.fft.irfft(fk,axis=-1,n=nLat)

def gaussian_smooth(fld,sig):
    """
    Do a Gaussian smoothing of the field.
    """
    return

if __name__=="__main__":
#    dat = load_fields('../fields.dat',1024)  # temporary for prototyping
    dat = load_fields('../fields-sg-l1.2-instanton.dat',1024)
#    dat = load_fields('../fields-sg-l1.3-instanton.dat',1024)
    rho = compute_rho(dat)
    
    f_1 = project_k_modes(dat[:,:,0],10)
    f_2 = project_k_modes(dat[:,:,0],15)
    f_3 = project_k_modes(dat[:,:,0],20)

    df_1 = project_k_modes(dat[:,:,1],10)
    df_2 = project_k_modes(dat[:,:,1],15)
    df_3 = project_k_modes(dat[:,:,1],20)

    pot = lambda x : np.cos(x)-1.+0.5*1.2**2*np.sin(x)**2
    r_1 = compute_rho_from_field(f_1,df_1,50.,pot)
    r_2 = compute_rho_from_field(f_2,df_2,50.,pot)
    r_3 = compute_rho_from_field(f_3,df_3,50.,pot)

    f,a = plt.subplots()
    a.contourf(r_1[:32])
    f1,a1 = plt.subplots()
    a1.contourf(0.5*df_1[:32]**2)
    f2,a2 = plt.subplots()
    a2.contourf(pot(f_1[:32]))
