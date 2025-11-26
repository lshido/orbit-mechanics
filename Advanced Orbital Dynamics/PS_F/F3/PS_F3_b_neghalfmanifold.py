ps = "F3 part b"
# Author: Lillian Shido
# Date: 11/24/2025

import pdb
import numpy as np
import pandas as pd
from math import atan2, sin
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import solve_ivp
from great_tables import GT, md, system_fonts
import altair as alt

from symbols import xi_symbol, eta_symbol
from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, calc_initial_velocities, find_halfperiod, calc_Jacobi, spatial_ode, calc_poincare_exponents, calc_spatial_monodromy_half

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[y_L1]})
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0]})
earth = pd.DataFrame({'name':["Earth"],'x':[x_Earth],'y':[0]})

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
fixed_point = pd.DataFrame({'label': [f"Fixed Point @ {eta_symbol}=0"],'name': ['Fixed Point'],'x':[x_fixed],'y':[y_fixed]})

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
    'name':'Periodic Orbit',
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

# Choose the stable eigenvalue
d_val = eigenvalues[1].real #[nondim]
d_dim = d_val*l_char #[km]
print(f'step off: {d_dim}[km]')
print(f'step off: {d_val}[non-dim]')

nu_W_stable = eigenvectors[:,1]/np.linalg.norm(eigenvectors[:,1][0:3])
print(nu_W_stable)
scaled_nu_W_stable = d_val*nu_W_stable
x0_pos = x_fixed - scaled_nu_W_stable[0].real[0,0]
y0_pos = y_fixed - scaled_nu_W_stable[1].real[0,0]
vx0_pos = vx_fixed - scaled_nu_W_stable[3].real[0,0]
vy0_pos = vy_fixed - scaled_nu_W_stable[4].real[0,0]
print(f"x0: {x0_pos:.4f}")
print(f"y0: {y0_pos:.4f}")
print(f"vx0: {vx0_pos:.4f}")
print(f"vy0: {vy0_pos:.4f}")
x_step_off = pd.DataFrame({
    'name':['Step-Off Point'],
    'x':[x0_pos],
    'y':[y0_pos],
    'vx':[vx0_pos],
    'vy':[vy0_pos],
})

# Propagate backward in time
IC_stable = [
    x0_pos, y0_pos, 0, vx0_pos, vy0_pos, 0,
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]
stable_prop = pd.DataFrame({})
reverse_prop = solve_ivp(spatial_ode, [0, -1.5*period], IC_stable, args=(mu,), rtol=1e-12,atol=1e-14)
stable_prop_data = pd.DataFrame({
    'name':'Backward Time',
    't':reverse_prop.t,
    'x':reverse_prop.y[0],
    'y':reverse_prop.y[1]
})
stable_prop = pd.concat([stable_prop, stable_prop_data], ignore_index=True)

# Build eigenvector df
scale = 0.1 # Extend eigenvector line out
vscale = 0.01
x_eig_max = x_fixed + scale*eigenvectors[:,1][0,0].real
x_eig_min = x_fixed + scale*-eigenvectors[:,1][0,0].real
y_eig_max = y_fixed + scale*eigenvectors[:,1][1,0].real
y_eig_min = y_fixed + scale*-eigenvectors[:,1][1,0].real
vx_eig_max = x_fixed + vscale*eigenvectors[:,1][3,0].real
vx_eig_min = x_fixed + vscale*-eigenvectors[:,1][3,0].real
vy_eig_max = y_fixed + vscale*eigenvectors[:,1][4,0].real
vy_eig_min = y_fixed + vscale*-eigenvectors[:,1][4,0].real
angle_max = np.rad2deg(atan2((y_eig_max - y_fixed),(x_eig_max - x_fixed)))
angle_min = np.rad2deg(atan2((y_eig_min - y_fixed),(x_eig_min - x_fixed)))
eigenspace = pd.DataFrame({})
eigenspace_pos_data = pd.DataFrame({
    'name': 'Stable Eigenspace (+)',
    'label':'Stable Eigenspace',
    'x':[x_fixed],
    'x2':[x_eig_max],
    'y':[y_fixed],
    'y2':[y_eig_max],
    'angle':[angle_max],
    'angle_wedge': [90-angle_max]
})
eigenspace = pd.concat([eigenspace, eigenspace_pos_data], ignore_index=True)
eigenspace_neg_data = pd.DataFrame({
    'name': 'Stable Eigenspace (-)',
    'label': 'Stable Eigenspace',
    'x':[x_fixed],
    'x2':[x_eig_min],
    'y':[y_fixed],
    'y2':[y_eig_min],
    'angle':[angle_min],
    'angle_wedge': [90-angle_min]
})
eigenspace = pd.concat([eigenspace, eigenspace_neg_data], ignore_index=True)

# Build Plots
orbit_plot = alt.Chart(orbit).mark_line(clip=True,strokeWidth=2).encode(
    x=alt.X('x:Q'),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q'),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N').title(None),
    order='t'
)

fixed_point_loc = alt.Chart(fixed_point).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N').title(None)
)

eigendirections = alt.Chart(eigenspace).mark_rule(strokeWidth=0.5,clip=True).encode(
    x='x:Q',
    y='y:Q',
    x2='x2:Q',
    y2='y2:Q',
    color=alt.Color('label:N', scale=alt.Scale(domain=['Stable Eigenspace'], range=['darkolivegreen'])).title(None)
)

eigendirection_arrows = alt.Chart(eigenspace).mark_point(shape="wedge",filled=True, fillOpacity=1,size=100,clip=True).encode(
        x='x2:Q',
        # x = alt.datum(0.874411),
        y='y2:Q',
        # y = alt.datum(-0.015899),
        angle=alt.Angle('angle_wedge').scale(domain=[0, 360]),
        # angle=alt.AngleValue(),
        color=alt.Color('label:N', scale=alt.Scale(domain=['Stable Eigenspace'], range=['darkolivegreen']), legend=None).title(None)
)

step_off_loc = alt.Chart(x_step_off).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Step-Off Point'], range=['red'])).title(None)
)

# scale = alt.Scale(domain=['L1','Moon'], range=['darkblue','gray'])
moon_loc = alt.Chart(moon).mark_point(filled=True,size=50,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Moon'], range=['gray'])).title(None)
)

L1_loc = alt.Chart(L1).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L1'], range=['darkblue'])).title(None)
)

earth_loc = alt.Chart(earth).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Earth'], range=['darkblue'])).title(None)
)

# Build plot
x_min = 0.78
x_max = 0.86
y_lim = (x_max-x_min)/2
new_chart = alt.Chart(stable_prop).mark_line(clip=True,strokeWidth=2).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='darkmulti')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Stable Negative Half-Manifold",f"of the orbit {xi_symbol}=0.01, {eta_symbol}=0 ({ps}, Lillian Shido)"]
)

new_chart_layer = alt.layer(orbit_plot, earth_loc, L1_loc, moon_loc, new_chart, eigendirections, fixed_point_loc, step_off_loc).resolve_scale(color='independent')
new_chart_layer.save(f'stable_negative_halfmanifold_backward_{ps}.png', ppi=200)
pdb.set_trace()