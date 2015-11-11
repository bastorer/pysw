import Differentiation as Diff
import numpy as np
import sys

def spectral_sw_flux(sim):

    # Some storage fields
    flux_u = np.zeros(sim.soln.u.shape)
    flux_v = np.zeros(sim.soln.v.shape)
    flux_h = np.zeros(sim.soln.h.shape)
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))

    # Difference in layer deformations gives layer depths
    hs[:,:,:sim.Nz] = sim.soln.h[:,:,:sim.Nz] - sim.soln.h[:,:,1:]

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = hs[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Coriolis terms
        flux_u[:,:,ii] += -sim.f0*v
        flux_v[:,:,ii] +=  sim.f0*u

        # x fluxes
        if sim.Nx > 1:
            # Compute derivatives
            if sim.geomx == 'periodic':
                du = sim.ddx_period(u,sim.ik)
                dv = sim.ddx_period(v,sim.ik)
                dh = sim.ddx_period(h,sim.ik)
            elif sim.geomx == 'wall':
                du = sim.ddx_odd(u,sim.ik)
                dv = sim.ddx_even(v,sim.ik)
                dh = sim.ddx_even(h,sim.ik)

            # Intra-layer dyanics
            flux_u[:,:,ii] += -u*du - sim.gs[ii]*dh
            flux_v[:,:,ii] += -u*dv    
            flux_h[:,:,ii] += -h*du - u*dh

        # y fluxes
        if sim.Ny > 1:
            # Compute derivatives
            if sim.geomy == 'periodic':
                du = sim.ddy_period(u,sim.il)
                dv = sim.ddy_period(v,sim.il)
                dh = sim.ddy_period(h,sim.il)
            elif sim.geomy == 'wall':
                du = sim.ddy_even(u,sim.il)
                dv = sim.ddy_odd(v,sim.il)
                dh = sim.ddy_even(h,sim.il)

            # Intra-layer dynamics
            flux_u[:,:,ii] += -v*du
            flux_v[:,:,ii] += -v*dv - sim.gs[ii]*dh
            flux_h[:,:,ii] += -h*dv - v*dh

    return flux_u, flux_v, flux_h

def spectral_sw_source(sim):
    return np.zeros(sim.soln.u.shape),np.zeros(sim.soln.v.shape),np.zeros(sim.soln.h.shape)

def spectral_sw(sim):
    sim.x_derivs = Diff.SPECTRAL_x
    sim.y_derivs = Diff.SPECTRAL_y
    sim.flux_function = spectral_sw_flux
    sim.source_function = spectral_sw_source
