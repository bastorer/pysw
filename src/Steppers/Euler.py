import numpy as np

def Euler(sim):
    curr_flux = sim.flux() + sim.source()

    sim.sol += curr_flux*sim.dt
    
    if sim.nfluxes > 0:
        sim.fluxes = [curr_flux]
        sim.dts    = [sim.dt]
