import numpy as np
from Euler import Euler
import matplotlib.pyplot as plt

def AB2(sim):
    if sim.nfluxes < 1:
        sim.nfluxes = 1

    if len(sim.fluxes) == 0:
        # Do euler

        Euler(sim)

    elif len(sim.fluxes) == 1:
        # Do AB2
        
        curr_flux = sim.flux() + sim.source()

        """
        plt.figure()
        flx = sim.flux()
        plt.plot(flx[sim.Iu,:,0,0])
        src = sim.source()
        plt.plot(src[sim.Iu,:,0,0])
        plt.show()
        """

        w1 = sim.dt*(1. + 0.5*sim.dt/sim.dts[0])
        w2 = -0.5*sim.dt**2/sim.dts[0]
        
        sim.sol += w1*curr_flux + w2*sim.fluxes[0]
        
        if sim.nfluxes == 1:
            sim.fluxes = [curr_flux]
            sim.dts    = [sim.dt]
        else:
            sim.fluxes = [curr_flux] + sim.fluxes
            sim.fluxes = [sim.dt] + sim.dts
