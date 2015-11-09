import numpy as np
import sys
import matplotlib.pyplot as plt

def TVD_RK3(sim):

    dt = sim.dt
    un = sim.sol.copy()

    u1 = un + dt*(sim.flux() + sim.source())

    sim.sol = u1.copy()
    u2 = 0.75*un + 0.25*u1 + 0.25*dt*(sim.flux() + sim.source())

    sim.sol =  u2.copy()
    k3 = sim.flux() + sim.source()

    sim.sol = (1./3)*un + (2./3)*u2 + (2./3)*dt*k3
