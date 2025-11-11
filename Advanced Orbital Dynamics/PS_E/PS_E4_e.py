ps = "E4 part e"
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
from methods import system_properties, calc_L1, build_A_matrix_collinear, spatial_ode, calc_L2

import warnings
warnings.filterwarnings("ignore",category=PendingDeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)
x_L2, y_L2 = calc_L2(mu, a_Moon)

# Determine the stable and unstable manifolds 
d_dim = 0 #[km] from L1 point
d_val = d_dim/l_char #[non-dim]
print(d_val)
tf = .5*pi
pos_tspan = [0,tf]
neg_tspan = [0,-tf]

A = build_A_matrix_collinear(mu, x_L1, y_L1, 0)
eigenvalues, eigenvectors = np.linalg.eig(A)
# Plot eigenvectors
fig1 = plt.figure()
ax1 = fig1.add_subplot()
ax1.scatter(x_L1, y_L1, s=25, label="L1", color="green")
ax1.scatter(x_L2, y_L2, s=25, label="L2", color="green")
ax1.scatter(x_Moon, 0, s=30, label="Moon", color="gray")
ax1.scatter(x_Earth, 0, s=30, label="Earth", color="blue")
for i in range(0,6):
    if i==0 or i==1:
        print(f"vector {i+1} i-comp: {eigenvectors[:,i][0].real[0,0]}")
        print(f"vector {i+1} y-comp: {eigenvectors[:,i][1].real[0,0]}")
        if i==0:
            label=f"$E^S$"
            linestyle="solid"
            label_point="$x_S^+$"
            label_manifold="$W_{{loc}}^S$"
            color="blue"
        elif i==1:
            label=f"$E^U$"
            linestyle="dashed"
            label_point="$x_U^+$"
            label_manifold="$W_{{loc}}^U$"
            color="red"
        x_eig = x_L1 + eigenvectors[:,i][0].real[0,0]
        y_eig = y_L1 + eigenvectors[:,i][1].real[0,0]
        x0_pos = x_L1
        y0_pos = y_L1
        vx0_pos = eigenvectors[:,i][3].real[0,0]
        vy0_pos = eigenvectors[:,i][4].real[0,0]
        print(f"{x0_pos:.4f}")
        print(f"{y0_pos:.4f}")
        print(f"{vx0_pos:.4f}")
        print(f"{vy0_pos:.4f}")
        # check the distance from eq
        check_distance = np.linalg.norm(np.array([[x0_pos-x_L1],[y0_pos-y_L1]]))*l_char
        print(f"Check distance!: {check_distance}")
        ax1.axline((x_L1,y_L1), (x_eig,y_eig), label=label,linestyle=linestyle)
        # ax1.scatter(x0_pos, y0_pos, label=label_point, s=10,zorder=3)
        IC_pos = [
            x0_pos,y0_pos,0,vx0_pos,vy0_pos,0,# IC states
            1,0,0,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,0,0,
            0,0,1,0,0,0,
            0,0,0,1,0,0,
            0,0,0,0,1,0,
            0,0,0,0,0,1
        ]
        pos_prop = solve_ivp(spatial_ode, pos_tspan, IC_pos, args=(mu,), rtol=1e-12,atol=1e-14)
        ax1.plot(pos_prop.y[0], pos_prop.y[1],label=label_manifold+" +", color=color, linestyle=linestyle, zorder=2.5)
        neg_prop = solve_ivp(spatial_ode, neg_tspan, IC_pos, args=(mu,), rtol=1e-12,atol=1e-14)
        ax1.plot(neg_prop.y[0], neg_prop.y[1], color=color, linestyle=linestyle, zorder=2.5)
        # Now do the other sides of the manifolds
        x0_neg = x_L1 - eigenvectors[:,i][0].real[0,0]
        y0_neg = y_L1 - eigenvectors[:,i][1].real[0,0]
        vx0_neg = -eigenvectors[:,i][3].real[0,0]
        vy0_neg = -eigenvectors[:,i][4].real[0,0]
        IC_neg = [
            x0_neg,y0_neg,0,vx0_neg,vy0_neg,0,# IC states
            1,0,0,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,0,0,
            0,0,1,0,0,0,
            0,0,0,1,0,0,
            0,0,0,0,1,0,
            0,0,0,0,0,1
        ]
        pos_prop = solve_ivp(spatial_ode, pos_tspan, IC_neg, args=(mu,), rtol=1e-12,atol=1e-14)
        ax1.plot(pos_prop.y[0], pos_prop.y[1],label=label_manifold+" -", color=color, linestyle=linestyle, zorder=2.5)
        neg_prop = solve_ivp(spatial_ode, neg_tspan, IC_neg, args=(mu,), rtol=1e-12,atol=1e-14)
        ax1.plot(neg_prop.y[0], neg_prop.y[1], color=color, linestyle=linestyle, zorder=2.5)


plt.legend(ncol=3,framealpha=1)
lim=1e-4
ax1.axis('equal')
# ax1.set(xlim=(-lim,lim),ylim=(-lim,lim))
ax1.set_xlabel("x [non-dim]")
ax1.set_ylabel("y [non-dim]")
plt.title(f"Stable and Unstable Manifolds ({ps}, Lillian Shido)")
plt.grid()
plt.savefig(f'manifolds_{ps}.png', dpi=300, bbox_inches='tight')
plt.show()






