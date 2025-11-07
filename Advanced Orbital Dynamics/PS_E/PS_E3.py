ps = "E3"
# Author: Lillian Shido
# Date: 11/6/2025

import sys
import pdb
import numpy as np
import pandas as pd
from math import pi, sqrt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from constants import r_Earth, mu_Earth
from methods import spatial_2bp_ode, calc_spatial_monodromy_half

# Circular orbit about Earth
e = 0
a = r_Earth + 2000
r = a
mu_2BP = mu_Earth
v0 = (mu_2BP/r)**(1/2) 

period = 2*pi/(mu_2BP/a**3)**(1/2)


IC = [
    r,0,0,0,v0,0,# IC states
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]
full_period_prop = solve_ivp(spatial_2bp_ode, [0, period], IC, args=(mu_2BP,r), rtol=1e-13,atol=1e-14)
half_period_prop = solve_ivp(spatial_2bp_ode, [0, period/2], IC, args=(mu_2BP,r), rtol=1e-13,atol=1e-14)
stm_half = half_period_prop.y[6:42,-1]
monodromy_full = full_period_prop.y[6:42,-1].reshape(6,6)
np.set_printoptions(threshold=sys.maxsize)
print(monodromy_full)
pdb.set_trace()
monodromy = calc_spatial_monodromy_half(stm_half)

fig, ax = plt.subplots()
# ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
# ax.yaxis.set_major_locator(ticker.MultipleLocator(50))
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
# ax.yaxis.set_minor_locator(ticker.MultipleLocator(25))
ax.plot(full_period_prop.y[0], full_period_prop.y[1], label='orbit')
ax.axis('equal')
ax.set(xlim=(-10000, 10000), ylim=(-10000,10000))
ax.legend(loc="upper center", fontsize=6)
plt.grid()
plt.xlabel("x [km]")
plt.ylabel("y [km]")
ax.tick_params(axis='both', which='major', labelsize=6)
ax.set_title(f'Circular orbit\n({ps}, Lillian Shido)')
plt.savefig(f'circular_orbit_{ps}.png', dpi=300, bbox_inches='tight')