ps = "E4 part c"
# Author: Lillian Shido
# Date: 11/9/2025

import pdb
import copy
import numpy as np
import pandas as pd
from math import pi, sqrt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
from scipy import stats
from great_tables import GT, md, html, style, loc, system_fonts
from pypalettes import load_cmap

from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, build_A_matrix_collinear, spatial_ode

import warnings
warnings.filterwarnings("ignore",category=PendingDeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)

# Determine the stable and unstable manifolds 
d_dim = 30 #[km]
d_val = d_dim/l_char #[non-dim]
tf = 1*pi
pos_tspan = [0,tf]
neg_tspan = [0,-tf]

A = build_A_matrix_collinear(mu, x_L1, y_L1, 0)
eigenvalues, eigenvectors = np.linalg.eig(A)
# Plot eigenvectors
fig1 = plt.figure()
ax1 = fig1.add_subplot()
# ax1.scatter(x_L1, y_L1, s=20, label="L1", color="red")
ax1.scatter(x_Moon, 0, s=20, label="Moon", color="gray")
for i in range(0,6):
    if i==0 or i==1:
        print(f"vector {i+1} i-comp: {eigenvectors[:,i][0].real[0,0]}")
        print(f"vector {i+1} y-comp: {eigenvectors[:,i][1].real[0,0]}")
        if i==0:
            label=f"$E^S$"
            linestyle="solid"
            label_point="$x_S^+$"
            label_manifold="$W_{{loc}}^S+$"
            color="blue"
        elif i==1:
            label=f"$E^U$"
            linestyle="dashed"
            label_point="$x_U^+$"
            label_manifold="$W_{{loc}}^U+$"
            color="red"
        x_eig = x_L1 + eigenvectors[:,i][0].real[0,0]
        y_eig = y_L1 + eigenvectors[:,i][1].real[0,0]
        nu_W = eigenvectors[:,i]/np.linalg.norm(eigenvectors[:,i][0:3])
        scaled_nu_W = d_val*nu_W
        x0 = x_L1 + scaled_nu_W[0].real[0,0]
        y0 = y_L1 + scaled_nu_W[1].real[0,0]
        vx0 = scaled_nu_W[3].real[0,0]
        vy0 = scaled_nu_W[4].real[0,0]
        # check the distance from eq
        check_distance = np.linalg.norm(np.array([[x0-x_L1],[y0-y_L1]]))*l_char
        print(f"Check distance!: {check_distance}")
        ax1.axline((x_L1,y_L1), (x_eig,y_eig), label=label,linestyle=linestyle)
        ax1.scatter(x0, y0, label=label_point, s=10,zorder=3)
        IC = [
            x0,y0,0,vx0,vy0,0,# IC states
            1,0,0,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,0,0,
            0,0,1,0,0,0,
            0,0,0,1,0,0,
            0,0,0,0,1,0,
            0,0,0,0,0,1
        ]
        pos_prop = solve_ivp(spatial_ode, pos_tspan, IC, args=(mu,), rtol=1e-12,atol=1e-14)
        ax1.plot(pos_prop.y[0], pos_prop.y[1],label=label_manifold+" +", color=color,zorder=2.5)
        neg_prop = solve_ivp(spatial_ode, neg_tspan, IC, args=(mu,), rtol=1e-12,atol=1e-14)
        ax1.plot(neg_prop.y[0], neg_prop.y[1],label=label_manifold+ "-", color=color,zorder=2.5)

plt.legend()
lim=1e-4
ax1.axis('equal')
# ax1.set(xlim=(-lim,lim),ylim=(-lim,lim))
ax1.set_xlabel("x [non-dim]")
ax1.set_ylabel("y [non-dim]")
plt.title(f"Stable and Unstable Manifolds ({ps}, Lillian Shido)")
plt.grid()
plt.savefig(f'manifolds_{ps}.png', dpi=300, bbox_inches='tight')
plt.show()






