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

from symbols import eta_symbol
from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, calc_initial_velocities, find_halfperiod, calc_Jacobi, spatial_ode, calc_poincare_exponents, calc_spatial_monodromy_half

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[y_L1]})

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
x_fixed = converged_initial_states[0]
y_fixed = converged_initial_states[1]
vx_fixed = converged_initial_states[2]
vy_fixed = converged_initial_states[3]
fixed_point = pd.DataFrame({'name': [f"Fixed Point @ {eta_symbol}=0"],'x':[x_fixed],'y':[y_fixed]})

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
# For plotting purposes:
full_period_prop = solve_ivp(spatial_ode, [0, period], IC, args=(mu,), rtol=1e-12,atol=1e-14)
orbit = pd.DataFrame({
    'name':'orbit',
    't':full_period_prop.t,
    'x':full_period_prop.y[0],
    'y':full_period_prop.y[1]
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

eigenspace = pd.DataFrame({})
velocity_eigenspace = pd.DataFrame({})
manifold = pd.DataFrame({})
for i in range(0,6):
    if i==0 or i==1:
    # if i==0:
        if i==0:
            name="Unstable Eigendirection"
        elif i==1:
            name="Stable Eigendirection"
        scale = 0.03 # Extend eigenvector line out
        vscale = 0.01
        x_eig_max = x_fixed + scale*eigenvectors[:,i][0,0].real
        x_eig_min = x_fixed + scale*-eigenvectors[:,i][0,0].real
        y_eig_max = y_fixed+ scale*eigenvectors[:,i][1,0].real
        y_eig_min = y_fixed+ scale*-eigenvectors[:,i][1,0].real
        vx_eig_max = x_fixed+ vscale*eigenvectors[:,i][3,0].real
        vx_eig_min = x_fixed+ vscale*-eigenvectors[:,i][3,0].real
        vy_eig_max = y_fixed+ vscale*eigenvectors[:,i][4,0].real
        vy_eig_min = y_fixed+ vscale*-eigenvectors[:,i][4,0].real
        eigenspace_pos_data = pd.DataFrame({
            'name': f'{name} (+)',
            'x':[x_fixed,x_eig_max],
            'y':[y_fixed,y_eig_max]
        })
        eigenspace = pd.concat([eigenspace, eigenspace_pos_data], ignore_index=True)
        eigenspace_neg_data = pd.DataFrame({
            'name': f'{name} (-)',
            'x':[x_eig_min,x_fixed],
            'y':[y_eig_min,y_fixed]
        })
        eigenspace = pd.concat([eigenspace, eigenspace_neg_data], ignore_index=True)
        velocity_eigenspace_pos_data = pd.DataFrame({
            'name': f'{name} (+)',
            'vx':[x_fixed,vx_eig_max],
            'vy':[y_fixed,vy_eig_max]
        })
        velocity_eigenspace = pd.concat([velocity_eigenspace, velocity_eigenspace_pos_data], ignore_index=True)
        velocity_eigenspace_neg_data = pd.DataFrame({
            'name': f'{name} (-)',
            'vx':[vx_eig_min,x_fixed],
            'vy':[vy_eig_min,y_fixed]
        })
        velocity_eigenspace = pd.concat([velocity_eigenspace, velocity_eigenspace_neg_data], ignore_index=True)

# Build plot
x_min = 0.78
x_max = 0.88
y_lim = (x_max-x_min)/2
base = alt.Chart(orbit).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N').title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Stable and Unstable Eigendirections on the Position Space",f"at the Fixed point @ {eta_symbol}=0 ({ps}, Lillian Shido)"]
)

L1_loc = alt.Chart(L1).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L1'], range=['darkblue'])).title(None)
)

fixed_point_loc = alt.Chart(fixed_point).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=[f'Fixed Point @ {eta_symbol}=0'], range=['green'])).title(None)
)

eigendirections = alt.Chart(eigenspace).mark_line(clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N').title("Position Space Projections")
)

velocity_eigendirections = alt.Chart(velocity_eigenspace).mark_line(clip=True).encode(
    x='vx:Q',
    y='vy:Q',
    color=alt.Color('name:N').title("Velocity Space Projections")
)

final_pos = alt.layer(base, L1_loc, eigendirections, fixed_point_loc).resolve_scale(color='independent')

final_pos.save(f'eigenspaces_fixed_point_{ps}.png', ppi=200)

final_vel = alt.layer(base, L1_loc, velocity_eigendirections, fixed_point_loc).resolve_scale(color='independent')

final_vel.save(f'eigenspaces_fixed_point_velocity_{ps}.png', ppi=200)