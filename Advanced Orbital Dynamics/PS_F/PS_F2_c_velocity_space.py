ps = "F2 part c"
# Author: Lillian Shido
# Date: 11/23/2025

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
# For plotting purposes:
full_period_prop = solve_ivp(spatial_ode, [0, period], IC, args=(mu,), rtol=1e-12,atol=1e-14)
orbit = pd.DataFrame({
    'name':'orbit',
    't':full_period_prop.t,
    'vx':full_period_prop.y[3],
    'vy':full_period_prop.y[4]
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

# # Propagate by t1
# t1 = 0.25*period
# print("\n".join([" ".join(f"{item:15.5f}" for item in row) for row in np.asarray(stm_t1)]))

# # Propagate a full period from t1 to get the monodromy @ t1
# t1_x = t1_prop.y[0,-1] 
# t1_y = t1_prop.y[1,-1] 
# t1_z = t1_prop.y[2,-1] 
# t1_vx = t1_prop.y[3,-1] 
# t1_vy = t1_prop.y[4,-1] 
# t1_vz = t1_prop.y[5,-1]
# t1 = pd.DataFrame({'name':["t1"],'x':[t1_x],'y':[t1_y]})

# IC_t1 = [
#     t1_x, t1_y, t1_z, t1_vx, t1_vy, t1_vz,
#     1,0,0,0,0,0, # Identity matrix for phi ICs
#     0,1,0,0,0,0,
#     0,0,1,0,0,0,
#     0,0,0,1,0,0,
#     0,0,0,0,1,0,
#     0,0,0,0,0,1
# ]
# t1_full_prop = solve_ivp(spatial_ode, [0, period], IC_t1, args=(mu,), rtol=1e-12,atol=1e-14)
# monodromy_t1 = t1_full_prop.y[6:42,-1].reshape(6,6)
# print("Monodromy @ t1")
# print("\n".join([" ".join(f"{item:15.5f}" for item in row) for row in np.asarray(monodromy_t1)]))

# monodromy_t1_calc = stm_t1 @ monodromy @ np.linalg.inv(stm_t1)
# print("Monodromy @ t1 from calcs")
# print("\n".join([" ".join(f"{item:15.5f}" for item in row) for row in np.asarray(monodromy_t1_calc)]))


fixed_points = dict()
for enum, f in enumerate(range(0, 100, 5)):
    fixed_points[enum] = dict()
    fixed_point_prop = solve_ivp(spatial_ode, [0, f/100*period], IC, args=(mu,), rtol=1e-12,atol=1e-14)
    stm_fixed_point = fixed_point_prop.y[6:42,-1].reshape(6,6)
    fixed_points[enum]['fp_x'] = fixed_point_prop.y[0,-1]
    fixed_points[enum]['fp_y'] = fixed_point_prop.y[1,-1]
    fixed_points[enum]['fp_vx'] = fixed_point_prop.y[3,-1]
    fixed_points[enum]['fp_vy'] = fixed_point_prop.y[4,-1]
    fixed_point_eigenvectors_list = []
    # Calculate the eigenvectors at the fixed point
    for i in range(0,6):
        fixed_point_eigenvectors_list.append((stm_fixed_point @ eigenvectors[:,i])/np.linalg.norm(stm_fixed_point @ eigenvectors[:,i]))
    fixed_point_eigenvectors = np.array(fixed_point_eigenvectors_list)
    fixed_points[enum]['eigenvectors'] = fixed_point_eigenvectors

# # Print the eigs
# for i in range(0,6):
#     print(f"t1 eigenvectors_{i+1}:")
#     print("\\begin{bmatrix}")
#     try:
#         for x in range(0,6): print(f"{t1_eigenvectors[i].flatten()[x]:15.5f}\\\\")
#     except:
#         pdb.set_trace()
#     print("\\end{bmatrix}")

# # Check their eigenvalues Av=λv
# t1_eigenvalues_list = []
# for i in range(0,6):
#     Av = monodromy_t1 @ t1_eigenvectors[i]
#     t1_eigenvalues_list.append(np.divide(Av,t1_eigenvectors[i]))
# t1_eigenvalues = np.array(t1_eigenvalues_list)
# # pdb.set_trace()
# # for i in range(0,6):
# #     print(f"Check eigenvalue match between t0 and t1")
# #     if np.isclose(eigenvalues[i],t1_eigenvalues[i],atol=1e-8).all():
# #         print("PASS!")
# #     else:
# #         print("FAIL.")
    
eigenspace = pd.DataFrame({})
velocity_eigenspace = pd.DataFrame({})
manifold = pd.DataFrame({})
for num, point in fixed_points.items():
    for i in range(0,6):
        if i==0 or i==1:
        # if i==0:
            if i==0:
                name="Stable Eigendirection"
            elif i==1:
                name="Unstable Eigendirection"
            scale = 0.05 # Extend eigenvector line out
            vscale = 0.1
            x_eig_max = point['fp_x'] + scale*point['eigenvectors'][i][0,0].real
            x_eig_min = point['fp_x'] + scale*-point['eigenvectors'][i][0,0].real
            y_eig_max = point['fp_y'] + scale*point['eigenvectors'][i][1,0].real
            y_eig_min = point['fp_y'] + scale*-point['eigenvectors'][i][1,0].real
            vx_eig_max = point['fp_vx'] + vscale*point['eigenvectors'][i][3,0].real
            vx_eig_min = point['fp_vx'] + vscale*-point['eigenvectors'][i][3,0].real
            vy_eig_max = point['fp_vy'] + vscale*point['eigenvectors'][i][4,0].real
            vy_eig_min = point['fp_vy'] + vscale*-point['eigenvectors'][i][4,0].real
            angle_max = np.rad2deg(atan2((y_eig_max - point['fp_y']),(x_eig_max - point['fp_x'])))
            angle_min = np.rad2deg(atan2((y_eig_min - point['fp_y']),(x_eig_min - point['fp_x'])))
            eigenspace_pos_data = pd.DataFrame({
                'name': f'{num}: {name} (+)',
                'label':f'{name}',
                'x':[point['fp_x']],
                'x2':[x_eig_max],
                'y':[point['fp_y']],
                'y2':[y_eig_max],
                'angle':[angle_max],
                'angle_wedge': [90-angle_max]
            })
            eigenspace = pd.concat([eigenspace, eigenspace_pos_data], ignore_index=True)
            eigenspace_neg_data = pd.DataFrame({
                'name': f'{num}: {name} (-)',
                'label': f'{name}',
                'x':[point['fp_x']],
                'x2':[x_eig_min],
                'y':[point['fp_y']],
                'y2':[y_eig_min],
                'angle':[angle_min],
                'angle_wedge': [90-angle_min]
            })
            eigenspace = pd.concat([eigenspace, eigenspace_neg_data], ignore_index=True)
            angle_max_v = np.rad2deg(atan2((vy_eig_max - point['fp_vy']),(vx_eig_max - point['fp_vx'])))
            angle_min_v = np.rad2deg(atan2((vy_eig_min - point['fp_vy']),(vx_eig_min - point['fp_vx'])))
            velocity_eigenspace_pos_data = pd.DataFrame({
                'name': f'{num}: {name} (+)',
                'label': f'{name}',
                'vx':[point['fp_vx']],
                'vx2':[vx_eig_max],
                'vy':[point['fp_vy']],
                'vy2':[vy_eig_max],
                'angle':[angle_max_v],
                'angle_wedge': [90-angle_max_v]
            })
            velocity_eigenspace = pd.concat([velocity_eigenspace, velocity_eigenspace_pos_data], ignore_index=True)
            velocity_eigenspace_neg_data = pd.DataFrame({
                'name': f'{num}: {name} (-)',
                'label': f'{name}',
                'vx':[point['fp_vx']],
                'vx2':[vx_eig_min],
                'vy':[point['fp_vy']],
                'vy2':[vy_eig_min],
                'angle':[angle_min_v],
                'angle_wedge': [90-angle_min_v]
            })
            velocity_eigenspace = pd.concat([velocity_eigenspace, velocity_eigenspace_neg_data], ignore_index=True)

# Build plot
x_min = -0.1
x_max = 0.1
y_lim = (x_max-x_min)/2
base = alt.Chart(orbit).mark_line(clip=True,strokeWidth=2).encode(
    x=alt.X('vx:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='vx [non-dim]')),
    # x=alt.X('vx:Q', axis=alt.Axis(title='vx [non-dim]')),
    y=alt.Y('vy:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='vy [non-dim]')),
    # y=alt.Y('vy:Q', axis=alt.Axis(title='vy [non-dim]')),
    color=alt.Color('name:N').title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Stable and Unstable Eigenvectors on the Velocity Space",f"of the orbit {xi_symbol}=0.01, {eta_symbol}=0 ({ps}, Lillian Shido)"]
)

velocity_eigendirections = alt.Chart(velocity_eigenspace).mark_rule(strokeWidth=0.5,clip=True).encode(
    x='vx:Q',
    y='vy:Q',
    x2='vx2:Q',
    y2='vy2:Q',
    color=alt.Color('label:N', scale=alt.Scale(domain=['Stable Eigendirection', 'Unstable Eigendirection'], range=['darkolivegreen', 'darkorange'])).title("Velocity Space Projections")
)

velocity_eigendirection_arrows = alt.Chart(velocity_eigenspace).mark_point(shape="wedge",filled=True, fillOpacity=1,size=100).encode(
        x='vx2:Q',
        # x = alt.datum(0.874411),
        y='vy2:Q',
        # y = alt.datum(-0.015899),
        angle=alt.Angle('angle_wedge').scale(domain=[0, 360]),
        # angle=alt.AngleValue(),
        color=alt.Color('label:N', scale=alt.Scale(domain=['Stable Eigendirection', 'Unstable Eigendirection'], range=['darkolivegreen', 'darkorange']), legend=None).title(None)
)

final_vel = alt.layer(base, velocity_eigendirections, velocity_eigendirection_arrows).resolve_scale(color='independent')

final_vel.save(f'eigenspaces_fixed_point_velocity_{ps}.png', ppi=200)

pdb.set_trace()