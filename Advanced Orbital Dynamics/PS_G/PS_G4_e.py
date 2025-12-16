ps = "G4 part e"
# Author: Lillian Shido
# Date: 12/15/2025

import pdb
import numpy as np
import pandas as pd
from math import atan2, sin, isclose, degrees
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import solve_ivp
from great_tables import GT, md, system_fonts,style, loc
import altair as alt
from copy import deepcopy
import plotly.express as px
import plotly.graph_objects as go

from symbols import xi_symbol, eta_symbol
from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, find_halfperiod_fixed_jacobi, spatial_ode, calc_spatial_Jacobi, calc_velocity_from_Jacobi, calc_L2, calc_initial_velocities, find_halfperiod, find_spatial_halfperiod
mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)
x_L2, y_L2 = calc_L2(mu, a_Moon)
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[y_L1],'z':[0]})
L2 = pd.DataFrame({'name':["L2"],'x':[x_L2],'y':[y_L2],'z':[0]})
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0],'z':[0]})
earth = pd.DataFrame({'name':["Earth"],'x':[x_Earth],'y':[0],'z':[0]})
mu = 0.01215058162343
desired_JC = 3.15000

# Calc L1 Lyapunov ICs
xi_L1 = -0.015
starting_x_L1 = x_L1 + xi_L1
xi_dot_L1, eta_dot_L1 = calc_initial_velocities(xi_L1, 0, x_L1, y_L1, mu)
ydot_guess_L1 = eta_dot_L1
IC_L1_Lyapunov = [
    starting_x_L1, 0, 0,
    0, ydot_guess_L1, 0
]
# Calc L2 Lyapunov ICs
xi_L2 = 0.01
starting_x_L2 = x_L2 + xi_L2
xi_dot_L2, eta_dot_L2 = calc_initial_velocities(xi_L2, 0, x_L2, y_L2, mu)
ydot_guess_L2 = eta_dot_L2
IC_L2_Lyapunov = [
    starting_x_L2, 0, 0,
    0, ydot_guess_L2, 0
]
IC_list = {'L1_Lyapunov':IC_L1_Lyapunov, 'L2_Lyapunov':IC_L2_Lyapunov}

# Initialize list to store orbit ICs
df_orbits = pd.DataFrame()
# Configuration
first_delta_x_value = 0.001

# Use the arrival states for xi (x0, y0, vx0, vy0) as the starting x
counter = 0
n_orbits = 200
orbit_desired_jc = pd.DataFrame()
for lyapunov, values in IC_list.items():
    orbit_x = values[0]
    orbit_ydot = values[5]
    for orbit in range(n_orbits):
        print(f"starting x0:{orbit_x}, starting ydot: {orbit_ydot}")
        iterations, tf, arrival_states, converged_initial_states = find_halfperiod(orbit_x, orbit_ydot, mu, tolerance=1e-12)
        jacobi = calc_spatial_Jacobi(mu, converged_initial_states[0], converged_initial_states[1], 0,converged_initial_states[2], converged_initial_states[3],0)
        delta_jc = desired_JC - jacobi
        if lyapunov == 'L1_Lyapunov':
            x_L = x_L1
            first_delta_x = -first_delta_x_value
        if lyapunov == 'L2_Lyapunov':
            x_L = x_L2
            first_delta_x = -first_delta_x_value
        if abs(delta_jc) > 1e-6:
            if orbit <= 1: # For the first two calculations, naively guess ydot is the same
                orbit_x = converged_initial_states[0] + first_delta_x # Step the x by delta_x
                orbit_ydot = converged_initial_states[3]
                print(f"found x0:{converged_initial_states[0]}, found ydot: {converged_initial_states[3]}")
                # Calc Jacobi
                # Compile data
                orbit_IC_data = pd.DataFrame({
                    "orbit":[orbit],
                    "iterations":[iterations],
                    "tf":[tf],
                    "xi":[converged_initial_states[0]-x_L],
                    "x":[converged_initial_states[0]],
                    "y":[converged_initial_states[1]],
                    "vx":[converged_initial_states[2]],
                    "vy":[converged_initial_states[3]],
                    "jacobi":[jacobi]
                })
                df_orbits = pd.concat([df_orbits, orbit_IC_data], ignore_index=True)
            else:
                # Calc the slope from the last pair of x0 and vy0
                x_jc_slope = (df_orbits.iloc[-1]['x']-df_orbits.iloc[-2]['x'])/(df_orbits.iloc[-1]['jacobi']-df_orbits.iloc[-2]['jacobi'])
                vy_x_slope = (df_orbits.iloc[-1]['vy']-df_orbits.iloc[-2]['vy'])/(df_orbits.iloc[-1]['x']-df_orbits.iloc[-2]['x'])
                delta_x = delta_jc*x_jc_slope
                orbit_x = converged_initial_states[0] + delta_x/5 # Step the x by delta_x
                orbit_ydot = converged_initial_states[3] + delta_x/5*vy_x_slope
                
                # Compile data
                orbit_IC_data = pd.DataFrame({
                    "orbit":[orbit],
                    "iterations":[iterations],
                    "tf":[tf],
                    "xi":[converged_initial_states[0]-x_L],
                    "x":[converged_initial_states[0]],
                    "y":[converged_initial_states[1]],
                    "vx":[converged_initial_states[2]],
                    "vy":[converged_initial_states[3]],
                    "jacobi":[jacobi]
                })
                df_orbits = pd.concat([df_orbits, orbit_IC_data], ignore_index=True)
        else:
            orbit_desired_jc_data = pd.DataFrame({
                "type":[lyapunov],
                "orbit":[orbit],
                "iterations":[iterations],
                "tf":[tf],
                "xi":[converged_initial_states[0]-x_L],
                "x":[converged_initial_states[0]],
                "y":[converged_initial_states[1]],
                "vx":[converged_initial_states[2]],
                "vy":[converged_initial_states[3]],
                "jacobi":[jacobi]
            })
            orbit_desired_jc = pd.concat([orbit_desired_jc, orbit_desired_jc_data], ignore_index=True)
            print("Found orbit with desired jacobi!")
            break

propagated_orbits = pd.DataFrame()
df_step_off = pd.DataFrame()
for row in orbit_desired_jc.iterrows():
    counter, tf, arrival_states, converged_IC = find_spatial_halfperiod(row[1]['x'],0,row[1]['vy'], mu)
    period = tf*2
    type_label = row[1]['type']
    df_IC = pd.DataFrame({
        "type": type_label,
        "period":[2*tf],
        "x":[converged_IC[0]],
        "y":[converged_IC[1]],
        "z":[converged_IC[2]],
        "vx":[converged_IC[3]],
        "vy":[converged_IC[4]],
        "vz":[converged_IC[5]],
    })

    # For plotting purposes:
    converged = [
        converged_IC[0], converged_IC[1], converged_IC[2], converged_IC[3], converged_IC[4], converged_IC[5],
        1,0,0,0,0,0, # Identity matrix for phi ICs
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,1
    ]
    
    halo_1_prop = solve_ivp(spatial_ode, [0, period], converged, args=(mu,), rtol=1e-14,atol=1e-16)
    orbit_1 = pd.DataFrame({
        'name': type_label,
        't':halo_1_prop.t,
        'x':halo_1_prop.y[0],
        'y':halo_1_prop.y[1],
        'z':halo_1_prop.y[2]
    })
    propagated_orbits = pd.concat([propagated_orbits, orbit_1], ignore_index=True)

heteroclinic_IC = [
    0.98785, -0.07835, 0,
    0.33414, -0.01246, 0
]
df_heteroclinic = pd.DataFrame()
# For plotting purposes:
converged = [
    heteroclinic_IC[0], heteroclinic_IC[1], heteroclinic_IC[2], heteroclinic_IC[3], heteroclinic_IC[4], heteroclinic_IC[5],
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]
L1_period = 1.89653*1.5
L2_period = 2.25323
L1_prop = solve_ivp(spatial_ode, [0, -L1_period], converged, args=(mu,), rtol=1e-14,atol=1e-16)
L1_heteroclinic = pd.DataFrame({
    'name': 'heteroclinic path',
    't':L1_prop.t,
    'x':L1_prop.y[0],
    'y':L1_prop.y[1],
    'z':L1_prop.y[2]
})
df_heteroclinic = pd.concat([df_heteroclinic, L1_heteroclinic], ignore_index=True)
L2_prop = solve_ivp(spatial_ode, [0, L1_period], converged, args=(mu,), rtol=1e-14,atol=1e-16)
L2_heteroclinic = pd.DataFrame({
    'name': 'heteroclinic path',
    't':L2_prop.t,
    'x':L2_prop.y[0],
    'y':L2_prop.y[1],
    'z':L2_prop.y[2]
})
df_heteroclinic = pd.concat([df_heteroclinic, L2_heteroclinic], ignore_index=True)

x_min = 0.8
x_max = 1.20
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2

moon_x_y = alt.Chart(moon).mark_point(filled=True,size=50,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Moon'], range=['gray'])).title(None)
)

L1_x_y = alt.Chart(L1).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L1'], range=['darkgreen'])).title(None)
)

L2_x_y = alt.Chart(L2).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L2'], range=['orangered'])).title(None)
)

earth_x_y = alt.Chart(earth).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Earth'], range=['darkblue'])).title(None)
)

orbit_x_y = alt.Chart(propagated_orbits).mark_line(strokeWidth=1,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N').title(None),
    order='t'
)

# 3-D chart
fig = px.line_3d(df_heteroclinic, x="x", y='y', z='z', color='name',
                 labels={
                     "x": "x [non-dim]",
                     "y": "y [non-dim]",
                     "z": "z [non-dim]"
                 })
fig.add_trace(go.Scatter3d(x=L1['x'],y=L1['y'],z=L1['z'],name='L1'))
fig.add_trace(go.Scatter3d(x=L2['x'],y=L2['y'],z=L2['z'],name='L2'))
fig.add_trace(go.Scatter3d(x=earth['x'],y=earth['y'],z=earth['z'],name='earth'))
fig.add_trace(go.Scatter3d(x=propagated_orbits['x'],y=propagated_orbits['y'],z=propagated_orbits['z'],mode='lines'))
fig.add_trace(go.Scatter3d(x=moon['x'],y=moon['y'],z=moon['z'],name='moon'))
fig.update_layout(
    title=dict(text=f"Stable and Unstable Manifolds in L1 and L2 Lyapunovs ({ps}, Lillian Shido)", font=dict(size=30), automargin=True, yref='paper'),
    legend=dict(title=dict(text=None)),
    scene = dict(
        xaxis=dict(range=[x_min, x_max]),
        yaxis=dict(range=[-y_lim,y_lim]),
        zaxis=dict(range=[-z_lim,z_lim]),
        aspectratio = dict(x=1,y=1,z=1),
        aspectmode='manual'  
    ),
)
# fig.show()

# Plot the heteroclinic transfer

x_min = 0.8
x_max = 1.20
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2

heteroclinic_transfer_chart = alt.Chart(df_heteroclinic).mark_line(strokeWidth=1,clip=True).encode(
# x_y_chart = alt.Chart(propagated_orbits).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('y:Q', axis=alt.Axis(title='y [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim, y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('vy:Q', axis=alt.Axis(title='vy [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(domain=['heteroclinic path'], range=['magenta'])).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Heteroclinic Path between L1 and L2",f"({ps}, Lillian Shido)"]
)

heteroclinic_transfer_layer = alt.layer(orbit_x_y, L1_x_y, L2_x_y, moon_x_y, heteroclinic_transfer_chart).resolve_scale(color='independent')
heteroclinic_transfer_layer.save(f'{ps}_heteroclinic.png', ppi=200)

pdb.set_trace()