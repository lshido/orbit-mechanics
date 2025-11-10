ps = "E4 part b"
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
from methods import system_properties, calc_L1, build_A_matrix_collinear

import warnings
warnings.filterwarnings("ignore",category=PendingDeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)

A = build_A_matrix_collinear(mu, x_L1, y_L1, 0)
eigenvalues, eigenvectors = np.linalg.eig(A)
# Plot eigenvectors
fig1 = plt.figure()
ax1 = fig1.add_subplot()
ax1.scatter(x_L1, y_L1, s=20, label="L1", color="red")
for x in range(0,6):
    if x==0 or x==1:
        print(f"vector {x+1} x-comp: {eigenvectors[:,x][0].real[0,0]}")
        print(f"vector {x+1} y-comp: {eigenvectors[:,x][1].real[0,0]}")
        if x==0:
            # label=f"Eigenvector {x+1} (Stable)"
            label=f"$E^S$"
            linestyle="solid"
        elif x==1:
            # label=f"Eigenvector {x+1} (Stable)"
            label=f"$E^U$"
            linestyle="dashed"
        x_eig = x_L1 + eigenvectors[:,x][0].real[0,0]
        y_eig = y_L1 + eigenvectors[:,x][1].real[0,0]
        # ax1.plot([x_L1,x_eig],[y_L1,y_eig], label=label, linewidth=1.5)
        ax1.axline((x_L1,y_L1), (x_eig,y_eig), label=label,linestyle=linestyle)
plt.legend()
ax1.axis('equal')
ax1.set_xlabel("x [non-dim]")
ax1.set_ylabel("y [non-dim]")
# plt.title(f"Stable and Unstable Eigenvectors ({ps}, Lillian Shido)")
plt.title(f"Stable and Unstable Eigenspaces ({ps}, Lillian Shido)")
plt.grid()
# plt.savefig(f'Eigenvectors_{ps}.png', dpi=300, bbox_inches='tight')
plt.savefig(f'Eigenspaces_{ps}.png', dpi=300, bbox_inches='tight')

# Calc angle between the stable and unstable eigenvectors
b_vec = eigenvectors[:,0][0:2]
c_vec = eigenvectors[:,1][0:2]
norm_b = np.linalg.norm(b_vec)
norm_c = np.linalg.norm(c_vec)
try:
    dot = np.dot(np.asarray(b_vec).flatten(),np.asarray(c_vec).flatten()).real
except:
    pdb.set_trace()
angle = np.rad2deg(np.arccos(dot/(norm_b*norm_c)))
print(f"angle b/t Eig_1 and eig_2: {angle}")


