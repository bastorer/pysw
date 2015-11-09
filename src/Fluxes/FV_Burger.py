import numpy as np
import matplotlib.pyplot as plt
import sys

# Compute the flux terms
def fv_Burger_flux(sim):

    # Some storage fields
    the_flux = np.zeros(sim.sol.shape)

    # Loop through each layer and compute flux
    for ii in range(sim.Nz):
       
        u = sim.sol[sim.Iu,:,:,ii]

        # Compute the interpolations
        umx,upx = sim.ax(u)
        umy,upy = sim.ay(u)

        # Compute the winds
        EastWest = u.copy()

        EW = np.abs(EastWest + np.abs(EastWest)) > 0. # East Wind
        WW = np.abs(EastWest - np.abs(EastWest)) > 0. # West Wind

        # Compute fluxes

        ##
        ## u flux
        ## (u*u)_x
        ##

        FUp = 0.5*upx*upx
        FUm = np.roll(FUp,1,0)
        FDp = 0.5*umx*umx
        FDm = np.roll(FDp,-1,0)

        Fxp = EW*FUp - WW*FDp
        Fxm = EW*FUm - WW*FDm

        the_flux[sim.Iu,:,:,ii] =  -(Fxp - Fxm)/sim.dx[0]
        
    return the_flux

# Deal with the source terms
def fv_Burger_source(sim):
    
    # Some storage fields
    the_source = np.zeros(sim.sol.shape)

    return the_source

