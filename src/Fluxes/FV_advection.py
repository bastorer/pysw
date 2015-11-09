import numpy as np

def fv_advection_flux(sim):

    # Some storage fields
    the_flux = np.zeros(sim.sol.shape)

    # Loop through each layer and compute flux
    for ii in range(sim.Nz):

        u = sim.sol[sim.Iu,:,:,ii]
        h = sim.sol[sim.Ih,:,:,ii]
        
        Wind = u.copy()
        EW = np.abs(Wind + np.abs(Wind)) > 0.
        WW = np.abs(Wind - np.abs(Wind)) > 0.

        hxm, hxp = sim.ax(h)
        uxm, uxp = sim.ax(u)

        #  Flux Upwind/Down minus/plus (West/East)
        FUp = hxp*uxp
        FDp = hxm*uxm
        FUm = np.roll(FUp,1,0)
        FDm = np.roll(FDp,-1,0)
        
        Fxp = EW*FUp - WW*FDp
        Fxm = EW*FUm - WW*FDm

        the_flux[sim.Ih,:,:,ii] =  -(Fxp - Fxm)/sim.dx[0] 

    return the_flux

def fv_advection_source(sim):
    return np.zeros(sim.sol.shape)
