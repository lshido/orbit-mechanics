ps = "G4 part a"
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
from methods import system_properties, calc_L1, find_halfperiod_fixed_jacobi, spatial_ode, calc_spatial_Jacobi, calc_velocity_from_Jacobi, calc_L2, calc_initial_velocities, find_halfperiod
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
# tf_1 = 2.77648121127569
orbit = pd.DataFrame()
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
                "period":[2*tf*t_char/3600/24],
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

lyapunov_specific_jc_table = (
    GT(orbit_desired_jc)
    .tab_header(
        title=md(f"Initial Conditions of L1 and L2 Lyapunov Orbits with JC=3.15000 ({ps}, Lillian Shido)")
    )
    .cols_label(
        period="{{Period}}<br>[days]",
        xi="{{:xi:}}<br>[nd]",
        x="{{x}}<br>[nd]",
        y="{{y}}<br>[nd]",
        vx="{{vx}}<br>[nd]",
        vy="{{vy}}<br>[nd]"
    )
    .fmt_number(
        columns=["period", "xi","x","y","vx","vy",'jacobi'],
        n_sigfig=5
    )
    .fmt_number(
        columns=['jacobi'],
        n_sigfig=6
    )
    .cols_align(
        align="center"
    )
    .cols_hide(
        columns=['orbit','iterations','tf']
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
lyapunov_specific_jc_table.show()

orbit = pd.DataFrame()
for row in orbit_desired_jc.iterrows():
    name = row[1]['type']
    period = 2*row[1]['tf']
    converged = [
        row[1]['x'], row[1]['y'], 0, row[1]['vx'], row[1]['vy'], 0,
        1,0,0,0,0,0, # Identity matrix for phi ICs
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,1
    ]
    
    prop = solve_ivp(spatial_ode, [0, period], converged, args=(mu,), rtol=1e-14,atol=1e-16)
    orbit_1 = pd.DataFrame({
        'name': name,
        't':prop.t,
        'x':prop.y[0],
        'y':prop.y[1],
        'z':prop.y[2]
    })
    orbit = pd.concat([orbit, orbit_1], ignore_index=True)

moon_x_y = alt.Chart(moon).mark_point(filled=True,size=50,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Moon'], range=['gray'])).title(None)
)
moon_y_z = alt.Chart(moon).mark_point(filled=True,size=50,clip=True).encode(
    x='y:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Moon'], range=['gray'])).title(None)
)
moon_x_z = alt.Chart(moon).mark_point(filled=True,size=50,clip=True).encode(
    x='x:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Moon'], range=['gray'])).title(None)
)

L1_x_y = alt.Chart(L1).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L1'], range=['darkgreen'])).title(None)
)
L1_y_z = alt.Chart(L1).mark_point(filled=True,size=30,clip=True).encode(
    x='y:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L1'], range=['darkgreen'])).title(None)
)
L1_x_z = alt.Chart(L1).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L1'], range=['darkgreen'])).title(None)
)


L2_x_y = alt.Chart(L2).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L2'], range=['orangered'])).title(None)
)
L2_y_z = alt.Chart(L2).mark_point(filled=True,size=30,clip=True).encode(
    x='y:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L2'], range=['orangered'])).title(None)
)
L2_x_z = alt.Chart(L2).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L2'], range=['orangered'])).title(None)
)

x_min = 0.8
x_max = 1.20
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2

# x-y plane chart
x_y_chart = alt.Chart(orbit).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='darkmulti')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["L1 and L2 Lyapunovs with JC=3.15000",f"({ps}, Lillian Shido)"]
)
x_y_chart_layer = alt.layer(x_y_chart, moon_x_y, L1_x_y,L2_x_y).resolve_scale(color='independent')
x_y_chart_layer.save(f'Lyapunovs {ps}.png', ppi=200)