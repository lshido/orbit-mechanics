ps = "F1 part a"
# Author: Lillian Shido
# Date: 11/13/2025

import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, build_A_matrix_collinear

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)

A = build_A_matrix_collinear(mu, x_L1, y_L1, 0)
eigenvalues, eigenvectors = np.linalg.eig(A)
for i in range(0,6):
    print(f"lambda_{i+1}: {eigenvalues[i]:.5f}")
    print(f"eigenvectors_{i+1}:")
    print("\\begin{bmatrix}")
    try:
        for x in range(0,6): print(f"{eigenvectors[:,i].flatten()[0,x]:15.5f}\\\\")
    except:
        pdb.set_trace()
    print("\\end{bmatrix}")
    # Check if eigs match their eigvectors
    if np.isclose(A*eigenvectors[:,i],eigenvalues[i]*eigenvectors[:,i],atol=1e-8).all():
        print("PASS!")
    else:
        print("FAIL.")

print("\n".join([" ".join(f"{item:15.5f}" for item in row) for row in np.asarray(A)]))

fig1 = plt.figure()
ax1 = fig1.add_subplot()
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
# ax1.scatter(x_L1, y_L1, s=20, label="L1", color="red")
for x in range(0,6):
    if x==0 or x==1:
        print(f"vector {x+1} x-comp: {eigenvectors[:,x][3].real[0,0]}")
        print(f"vector {x+1} y-comp: {eigenvectors[:,x][4].real[0,0]}")
        if x==0:
            # label=f"Eigenvector {x+1} (Stable)"
            label=f"$E^S$"
            linestyle="solid"
        elif x==1:
            # label=f"Eigenvector {x+1} (Unstable)"
            label=f"$E^U$"
            linestyle="dashed"
        x_eig = x_L1 + eigenvectors[:,x][0].real[0,0]
        y_eig = y_L1 + eigenvectors[:,x][1].real[0,0]
        vx_eig = eigenvectors[:,x][3].real[0,0]
        vy_eig = eigenvectors[:,x][4].real[0,0]
        # ax1.plot([x_L1,x_eig],[y_L1,y_eig], label=label, linewidth=1.5)
        ax1.plot([0,vx_eig],[0,vy_eig], label=label, linewidth=1.5)
        # ax1.axline((x_L1,y_L1), (x_eig,y_eig), label=label+" pos",linestyle=linestyle)
        # ax1.axline((0,0), (vx_eig,vy_eig), label=label+" vel",linestyle=linestyle, color="green")
plt.legend(framealpha=1)
ax1.set(xlim=(-0.8,0.8))
ax1.axis('equal')
ax1.set_xlabel(r"$\dot{x}$"+"\n[non-dim]")
ax1.set_ylabel(r"$\dot{y}$"+"\n[non-dim]")
plt.title(f"Stable and Unstable Eigenvectors Projected onto the Velocity Space\n({ps}, Lillian Shido)")
# plt.title(f"Stable and Unstable Eigenspaces Projected onto the Velocity Space\n({ps}, Lillian Shido)")
plt.grid()
plt.savefig(f'Eigenvectors_velocity {ps}.png', dpi=300, bbox_inches='tight')
# plt.savefig(f'Eigenspaces_velocity{ps}.png', dpi=300, bbox_inches='tight')

# Calc angle between the stable and unstable eigenvectors projected on the velocity space
b_vec = eigenvectors[:,0][3:5]
c_vec = eigenvectors[:,1][3:5]
norm_b = np.linalg.norm(b_vec)
norm_c = np.linalg.norm(c_vec)
try:
    dot = np.dot(np.asarray(b_vec).flatten(),np.asarray(c_vec).flatten()).real
except:
    pdb.set_trace()
angle = np.rad2deg(np.arccos(dot/(norm_b*norm_c)))
print(f"angle b/t Eig_1 and eig_2: {angle}")

pdb.set_trace()