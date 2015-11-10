import numpy as np
import sys

def spectral_sw_flux(sim):

    # Some storage fields
    the_flux = np.zeros(sim.sol.shape)
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))

    # Difference in layer deformations gives layer depths
    hs[:,:,:sim.Nz] = sim.sol[sim.Ih,:,:,:sim.Nz] - sim.sol[sim.Ih,:,:,1:]

    # Loop through each layer and compute flux
    for ii in range(sim.Nz):

        ##
        ## For the moment, assume one layer
        ##
       
        h = hs[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.sol[sim.Iu,:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.sol[sim.Iv,:,:,ii].reshape((sim.Nx,sim.Ny))

        # Coriolis terms
        the_flux[sim.Iu,:,:,ii] += -sim.f0*v
        the_flux[sim.Iv,:,:,ii] +=  sim.f0*u

        # x fluxes
        if sim.Nx > 1:
            # Compute derivatives
            if sim.BCx == 'periodic':
                du = sim.ddx_period(u,sim.ik)
                dv = sim.ddx_period(v,sim.ik)
                dh = sim.ddx_period(h,sim.ik)
            elif sim.BCx == 'wall':
                du = sim.ddx_odd(u,sim.ik)
                dv = sim.ddx_even(v,sim.ik)
                dh = sim.ddx_even(h,sim.ik)

            # Compute fluxes
            the_flux[sim.Iu,:,:,ii] += -u*du - sim.gs[ii]*dh
            the_flux[sim.Iv,:,:,ii] += -u*dv    
            the_flux[sim.Ih,:,:,ii] += -h*du - u*dh

        # y fluxes
        if sim.Ny > 1:
            # Compute derivatives
            if sim.BCy == 'periodic':
                du = sim.ddy_period(u,sim.il)
                dv = sim.ddy_period(v,sim.il)
                dh = sim.ddy_period(h,sim.il)
            elif sim.BCy == 'wall':
                du = sim.ddy_even(u,sim.il)
                dv = sim.ddy_odd(v,sim.il)
                dh = sim.ddy_even(h,sim.il)

            # Compute fluxes
            the_flux[sim.Iu,:,:,ii] += -v*du
            the_flux[sim.Iv,:,:,ii] += -v*dv - sim.gs[ii]*dh
            the_flux[sim.Ih,:,:,ii] += -h*dv - v*dh

    return the_flux

def spectral_sw_source(sim):
    return np.zeros(sim.sol.shape)
