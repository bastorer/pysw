import numpy as np

def Euler(sim):
    flux_u, flux_v, flux_h = sim.flux()
    src_u,  src_v,  src_h  = sim.source()
    
    sim.soln.u += (flux_u+src_u)*sim.dt
    sim.soln.v += (flux_v+src_v)*sim.dt
    sim.soln.h += (flux_h+src_h)*sim.dt
    
    if sim.nfluxes > 0:
        sim.fluxes.u = [flux_u+src_u]
        sim.fluxes.v = [flux_v+src_v]
        sim.fluxes.h = [flux_h+src_h]
        sim.dts    = [sim.dt]
