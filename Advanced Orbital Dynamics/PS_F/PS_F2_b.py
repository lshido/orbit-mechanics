ps = "F2 part b"
# Author: Lillian Shido
# Date: 11/18/2025

import pdb
import numpy as np
import pandas as pd
from math import pi
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import solve_ivp
from great_tables import GT, md, system_fonts
import altair as alt

from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, calc_initial_velocities, find_halfperiod, calc_Jacobi, spatial_ode, calc_poincare_exponents, calc_spatial_monodromy_half

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)

# Initial Conditions
xi = 0.01
eta = 0

# Distance from L1
xi_from_L1_dim = xi*l_char
eta_from_L1_dim = eta*l_char

# Calc starting guess values
xi_dot_0, eta_dot_0 = calc_initial_velocities(xi, eta, x_L1, y_L1, mu)
starting_x = x_L1 + xi
starting_y = y_L1 + eta
starting_xdot = xi_dot_0
ydot_guess = eta_dot_0

iterations, tf, arrival_states, converged_initial_states = find_halfperiod(starting_x, ydot_guess, mu, tolerance=1e-12)
period = 2*tf
# Spatial ICs, but for z=zdot=0
IC = [
    converged_initial_states[0], converged_initial_states[1], 0, converged_initial_states[2], converged_initial_states[3], 0,
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]
# Monodromy matrix from half period
print("Propagating the half-period")
print(f"period: {period}")
half_period_prop = solve_ivp(spatial_ode, [0, tf], IC, args=(mu,), rtol=1e-12,atol=1e-14)
stm_half = half_period_prop.y[6:42,-1].reshape(6,6)
monodromy = calc_spatial_monodromy_half(stm_half)
print("\n".join([" ".join(f"{item:15.5f}" for item in row) for row in np.asarray(monodromy)]))
df_monodromy = pd.DataFrame({
    'Determinant':[np.linalg.det(monodromy)],
    'Accuracy':[abs(np.linalg.det(monodromy))-1]
})

# Get eigs and check they match
eigenvalues, eigenvectors = np.linalg.eig(monodromy)
for i in range(0,6):
    print(f"Check eigenset {i+1}")
    if np.isclose(monodromy*eigenvectors[:,i],eigenvalues[i]*eigenvectors[:,i],atol=1e-8).all():
        print("PASS!")
    else:
        print("FAIL.")

# Print the eigs
for i in range(0,6):
    print(f"lambda_{i+1}: {eigenvalues[i]:.5f}")
    print(f"eigenvectors_{i+1}:")
    print("\\begin{bmatrix}")
    try:
        for x in range(0,6): print(f"{eigenvectors[:,i].flatten()[0,x]:15.5f}\\\\")
    except:
        pdb.set_trace()
    print("\\end{bmatrix}")

# Propagate by t1
t1 = 0.25*period
t1_prop = solve_ivp(spatial_ode, [0, t1], IC, args=(mu,), rtol=1e-12,atol=1e-14)
stm_t1 = t1_prop.y[6:42,-1].reshape(6,6)
print("STM @ t1")
print("\n".join([" ".join(f"{item:15.5f}" for item in row) for row in np.asarray(stm_t1)]))

# Calculate the eigenvectors at t1
t1_eigenvectors_list = []
for i in range(0,6):
    t1_eigenvectors_list.append(stm_t1 * eigenvectors[:,i])
t1_eigenvectors = np.array(t1_eigenvectors_list)

# Print the eigs
for i in range(0,6):
    print(f"eigenvectors_{i+1}:")
    print("\\begin{bmatrix}")
    try:
        for x in range(0,6): print(f"{t1_eigenvectors[i].flatten()[x]:15.5f}\\\\")
    except:
        pdb.set_trace()
    print("\\end{bmatrix}")

# Check their eigenvalues Av=λv
# for i in range(0,6):
#     Av = stm_t1 * t1_eigenvectors[i]
#     if nonzero_index = np.where(t1_eigenvector != 0)

eigenspace = pd.DataFrame({})
manifold = pd.DataFrame({})
for i in range(0,6):
    if i==0 or i==1:
        print(f"vector {i+1} i-comp: {t1_eigenvectors[0].real[0,0]}")
        print(f"vector {i+1} y-comp: {t1_eigenvectors[1].real[0,0]}")
        if i==0:
            name="Stable Eigenspace"
            label=f"$E^S$"
            linestyle="solid"
            label_point="$x_S^+$"
            label_manifold="Stable Half-Manifold"
            color="blue"
        elif i==1:
            name="Unstable Eigenspace"
            label=f"$E^U$"
            linestyle="dashed"
            label_point="$x_U^+$"
            label_manifold="Unstable Half-Manifold"
            color="red"
        scale = 1 # Extend eigenvector line out
        x_eig_max = x_L1 + scale*t1_eigenvectors[0].real[0,0]
        x_eig_min = x_L1 + scale*-t1_eigenvectors[0].real[0,0]
        y_eig_max = y_L1 + scale*t1_eigenvectors[1].real[0,0]
        y_eig_min = y_L1 + scale*-t1_eigenvectors[1].real[0,0]
        # nu_W = t1_eigenvectors/np.linalg.norm(t1_eigenvectors[0:3])
        # print(nu_W)
        # scaled_nu_W = d_val*nu_W
        # x0_pos = x_L1 + scaled_nu_W[0].real[0,0]
        # y0_pos = y_L1 + scaled_nu_W[1].real[0,0]
        # vx0_pos = scaled_nu_W[3].real[0,0]
        # vy0_pos = scaled_nu_W[4].real[0,0]
        # print(f"{x0_pos:.4f}")
        # print(f"{y0_pos:.4f}")
        # print(f"{vx0_pos:.4f}")
        # print(f"{vy0_pos:.4f}")
        # # check the distance from eq
        # check_distance = np.linalg.norm(np.array([[x0_pos-x_L1],[y0_pos-y_L1]]))*l_char
        # print(f"Check distance!: {check_distance}")
        eigenspace_data = pd.DataFrame({
            'name': name,
            'x':[x_eig_min,x_L1,x_eig_max],
            'y':[y_eig_min,y_L1,y_eig_max]
        })
        eigenspace = pd.concat([eigenspace, eigenspace_data], ignore_index=True)
pdb.set_trace()
# Build plot
x_min = 0.65
x_max = 1.05
y_lim = (x_max-x_min)/2
base = alt.Chart(eigenspace).mark_line(clip=True).encode(
    # x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    # y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N').title(None),
).properties(
    width=400,
    height=400,
    title=f"Stable and Unstable Eigenspaces at the fixed point ({ps}, Lillian Shido)"
)
base.save(f'eigenspaces_fixed_point_{ps}.png', ppi=200)