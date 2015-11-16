import Differentiation as Diff
import numpy as np
import sys
from scipy.fftpack import fftn, ifftn, fftfreq

def spectral_sw_flux(sim):

    # Some storage fields
    sim.curr_flux.u *= 0.
    sim.curr_flux.v *= 0.
    sim.curr_flux.h *= 0.
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))

    # Difference in layer deformations gives layer depths
    hs[:,:,:sim.Nz] = sim.soln.h[:,:,:sim.Nz] - sim.soln.h[:,:,1:]

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = hs[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Coriolis terms
        sim.curr_flux.u[:,:,ii] +=  sim.f0*v
        sim.curr_flux.v[:,:,ii] += -sim.f0*u

        # x fluxes
        if sim.Nx > 1:
            # Compute derivatives
            du = sim.ddx_u(u,sim.ik)
            dv = sim.ddx_v(v,sim.ik)
            dh = sim.ddx_h(h,sim.ik)

            # Intra-layer dyanics
            sim.curr_flux.u[:,:,ii] += -u*du - sim.gs[ii]*dh
            sim.curr_flux.v[:,:,ii] += -u*dv    
            sim.curr_flux.h[:,:,ii] += -h*du - u*dh

        # y fluxes
        if sim.Ny > 1:
            # Compute derivatives
            du = sim.ddy_u(u,sim.il)
            dv = sim.ddy_v(v,sim.il)
            dh = sim.ddy_h(h,sim.il)

            # Intra-layer dynamics
            sim.curr_flux.u[:,:,ii] += -v*du
            sim.curr_flux.v[:,:,ii] += -v*dv - sim.gs[ii]*dh
            sim.curr_flux.h[:,:,ii] += -h*dv - v*dh

    return

def spectral_sw_source(sim):
    return 

def filter_periodic_periodic(sim):
    for ii in range(sim.Nz):
        sim.soln.u[:,:,ii] = ifftn(sim.sfilt*fftn(sim.soln.u[:,:,ii],axes=[0,1]),axes=[0,1]).real
        sim.soln.v[:,:,ii] = ifftn(sim.sfilt*fftn(sim.soln.v[:,:,ii],axes=[0,1]),axes=[0,1]).real
        sim.soln.h[:,:,ii] = ifftn(sim.sfilt*fftn(sim.soln.h[:,:,ii],axes=[0,1]),axes=[0,1]).real

def filter_periodic_wall(sim):
    for ii in range(sim.Nz):
        tmp = sim.soln.u[:,:,ii]
        tmp = np.concatenate([tmp, tmp[:,::-1]],axis=1)
        sim.soln.u[:,:,ii] = ifftn(sim.sfilt*fftn(tmp, axes=[0,1]),axes=[0,1]).real[:sim.Nx,:sim.Ny]

        tmp = sim.soln.v[:,:,ii]
        tmp = np.concatenate([tmp,-tmp[:,::-1]],axis=1)
        sim.soln.v[:,:,ii] = ifftn(sim.sfilt*fftn(tmp,axes=[0,1]),axes=[0,1]).real[:sim.Nx,:sim.Ny]
        
        tmp = sim.soln.h[:,:,ii]
        tmp = np.concatenate([tmp, tmp[:,::-1]],axis=1)
        sim.soln.h[:,:,ii] = ifftn(sim.sfilt*fftn(tmp,axes=[0,1]),axes=[0,1]).real[:sim.Nx,:sim.Ny]

def filter_wall_periodic(sim):
    for ii in range(sim.Nz):
        tmp = sim.soln.u[:,:,ii]
        tmp = np.concatenate([tmp,-tmp[::-1,:]],axis=0)
        sim.soln.u[:,:,ii] = ifftn(sim.sfilt*fftn(tmp, axes=[0,1]),axes=[0,1]).real[:sim.Nx,:sim.Ny]

        tmp = sim.soln.v[:,:,ii]
        tmp = np.concatenate([tmp, tmp[::-1,:]],axis=0)
        sim.soln.v[:,:,ii] = ifftn(sim.sfilt*fftn(tmp,axes=[0,1]),axes=[0,1]).real[:sim.Nx,:sim.Ny]
        
        tmp = sim.soln.h[:,:,ii]
        tmp = np.concatenate([tmp, tmp[::-1,:]],axis=0)
        sim.soln.h[:,:,ii] = ifftn(sim.sfilt*fftn(tmp,axes=[0,1]),axes=[0,1]).real[:sim.Nx,:sim.Ny]

def filter_wall_wall(sim):
    for ii in range(sim.Nz):
        tmp = sim.soln.u[:,:,ii]
        tmp = np.concatenate([np.concatenate([tmp,         -tmp[::-1,:]],axis=0), \
                              np.concatenate([tmp[:,::-1], -tmp[::-1,::-1]],axis=0)], \
                    axis=1)
        sim.soln.u[:,:,ii] = ifftn(sim.sfilt*fftn(tmp, axes=[0,1]),axes=[0,1]).real[:sim.Nx,:sim.Ny]

        tmp = sim.soln.v[:,:,ii]
        tmp = np.concatenate([np.concatenate([ tmp,          tmp[::-1,:]],axis=0), \
                              np.concatenate([-tmp[:,::-1], -tmp[::-1,::-1]],axis=0)], \
                    axis=1)
        sim.soln.v[:,:,ii] = ifftn(sim.sfilt*fftn(tmp,axes=[0,1]),axes=[0,1]).real[:sim.Nx,:sim.Ny]
        
        tmp = sim.soln.h[:,:,ii]
        tmp = np.concatenate([np.concatenate([tmp,         tmp[::-1,:]],axis=0), \
                              np.concatenate([tmp[:,::-1], tmp[::-1,::-1]],axis=0)], \
                    axis=1)
        sim.soln.h[:,:,ii] = ifftn(sim.sfilt*fftn(tmp,axes=[0,1]),axes=[0,1]).real[:sim.Nx,:sim.Ny]


def spectral_sw(sim):
    sim.x_derivs = Diff.SPECTRAL_x
    sim.y_derivs = Diff.SPECTRAL_y
    sim.flux_function = spectral_sw_flux
    sim.source_function = spectral_sw_source
    if sim.geomx == 'periodic':
        if sim.geomy == 'periodic':
            sim.apply_filter = filter_periodic_periodic
        elif sim.geomy == 'wall':
            sim.appy_filter = filter_periodic_wall
    elif sim.geomx == 'wall':
        if sim.geomy == 'periodic':
            sim.apply_filter = filter_wall_periodic
        elif sim.geomy == 'wall':
            sim.apply_filter = filter_wall_wall

