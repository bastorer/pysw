import numpy as np
import sys
import matplotlib.pyplot as plt

def RK4(sim):

    # y^{n+1} = y^n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
    # k1 = flux(y^n)
    # k2 = flux(y^n + (dt/2)*k1)
    # k3 = flux(y^n + (dt/2)*k2)
    # k4 = flux(y^n + dt*k3)

    # Compute k1
    k1 = sim.flux() + sim.source()

    # Compute k2
    sim.sol += k1*sim.dt/2.
    k2 = sim.flux() + sim.source()

    # Compute k3
    sim.sol += -k1*sim.dt/2. + k2*sim.dt/2.
    k3 = sim.flux() + sim.source()

    # Compute k4
    sim.sol += -k2*sim.dt/2. + k3*sim.dt
    k4 = sim.flux() + sim.source()

    sim.sol += -k3*sim.dt + (sim.dt/6.)*(k1 + 2.*k2 + 2.*k3 + k4)
