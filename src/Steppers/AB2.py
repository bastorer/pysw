import numpy as np
from Euler import Euler
import matplotlib.pyplot as plt

def AB2(sim):
    if sim.nfluxes < 1:
        sim.nfluxes = 1

    if len(sim.fluxes.u) == 0:
        # Do euler

        Euler(sim)

    elif len(sim.fluxes.u) == 1:
        # Do AB2
        
        flux_u, flux_v, flux_h = sim.flux()
        src_u,  src_v,  src_h  = sim.source()

        w1 = sim.dt*(1. + 0.5*sim.dt/sim.dts[0])
        w2 = -0.5*sim.dt**2/sim.dts[0]
        
        sim.soln.u += w1*(flux_u+src_u) + w2*sim.fluxes.u[0]
        sim.soln.v += w1*(flux_v+src_v) + w2*sim.fluxes.v[0]
        sim.soln.h += w1*(flux_h+src_h) + w2*sim.fluxes.h[0]
        
        if sim.nfluxes == 1:
            sim.fluxes.u = [flux_u+src_u]
            sim.fluxes.v = [flux_v+src_v]
            sim.fluxes.h = [flux_h+src_h]
            sim.dts    = [sim.dt]
        else:
            sim.fluxes.u = [flux_u+src_u] + sim.fluxes.u
            sim.fluxes.v = [flux_v+src_v] + sim.fluxes.v
            sim.fluxes.h = [flux_h+src_h] + sim.fluxes.h
            sim.dts = [sim.dt] + sim.dts
