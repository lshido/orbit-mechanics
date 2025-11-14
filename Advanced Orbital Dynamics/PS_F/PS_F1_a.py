ps = "F1 part a"
# Author: Lillian Shido
# Date: 11/8/2025

import pdb
import numpy as np
import matplotlib.pyplot as plt

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
ax1.scatter(x_L1, y_L1, s=20, label="L1", color="red")
for x in range(0,6):
    if x==0 or x==1:
        print(f"vector {x+1} x-comp: {eigenvectors[:,x][0].real[0,0]}")
        print(f"vector {x+1} y-comp: {eigenvectors[:,x][1].real[0,0]}")
        if x==0:
            label=f"Eigenvector {x+1} (Stable)"
            # label=f"$E^S$"
            # linestyle="solid"
        elif x==1:
            label=f"Eigenvector {x+1} (Unstable)"
            # label=f"$E^U$"
            # linestyle="dashed"
        x_eig = x_L1 + eigenvectors[:,x][0].real[0,0]
        y_eig = y_L1 + eigenvectors[:,x][1].real[0,0]
        ax1.plot([x_L1,x_eig],[y_L1,y_eig], label=label, linewidth=1.5)
        # ax1.axline((x_L1,y_L1), (x_eig,y_eig), label=label,linestyle=linestyle)
plt.legend()
ax1.axis('equal')
ax1.set_xlabel("x [non-dim]")
ax1.set_ylabel("y [non-dim]")
plt.title(f"Stable and Unstable Eigenvectors ({ps}, Lillian Shido)")
# plt.title(f"Stable and Unstable Eigenspaces ({ps}, Lillian Shido)")
plt.grid()
plt.savefig(f'Eigenvectors_{ps}.png', dpi=300, bbox_inches='tight')
# plt.savefig(f'Eigenspaces_{ps}.png', dpi=300, bbox_inches='tight')

pdb.set_trace()