import numpy as np   
from scipy.fftpack import fft, ifft, fftfreq

def ddy_period(f,il):

    df = np.real(ifft(il*fft(f,axis=1),axis=1))

    return df

def ddy_even(f,il):

    N = f.shape[1]
    fe = np.concatenate([f,f],axis=1)
    df = np.real(ifft(il*fft(fe,axis=1),axis=1))
    df = df[:,:N]
    
    return df

def ddy_odd(f,il):

    N = f.shape[1]
    fe = np.concatenate([f,-f],axis=1)
    df = np.real(ifft(il*fft(fe,axis=1),axis=1))
    df = df[:,:N]
    
    return df


def SPECTRAL_y(sim):

    # Set the wavenumber vectors  
    if sim.BCy =='periodic':
        ky = 2*np.pi*fftfreq(sim.Ny,d=sim.Ly/sim.Ny)
        sim.ky = ky.copy()
        sim.il = 1j*np.tile(ky.reshape((1,sim.Ny)),(sim.Nx,1))
    elif sim.BCy == 'walls':
        ky = 1j*np.pi*fftfreq(2*sim.Ny,d=sim.Ly/sim.Ny)
        sim.ky = ky.copy()
        sim.il = 1j*np.tile(ky.reshape((1,2*sim.Ny)),(sim.Nx,1))
    else:
        print "x boundary conditions must be from the list: periodic, walls"
        sys.exit()

    # Set the differentiation operators
    sim.ddy_period = ddy_period
    sim.ddy_even   = ddy_even
    sim.ddy_odd    = ddy_odd
