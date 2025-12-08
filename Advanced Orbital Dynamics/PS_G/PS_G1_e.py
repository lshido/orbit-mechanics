ps = "G1 part e"
# Author: Lillian Shido
# Date: 12/8/2025

import pdb
import numpy as np
import pandas as pd
from math import atan2, sin, isclose, degrees
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import solve_ivp
from great_tables import GT, md, system_fonts
import altair as alt
from copy import deepcopy
import plotly.express as px

from symbols import xi_symbol, eta_symbol
from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, find_spatial_halfperiod_fixed_z, calc_spatial_Jacobi, spatial_ode, calc_poincare_exponents, calc_spatial_monodromy_half

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[y_L1]})
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0]})
earth = pd.DataFrame({'name':["Earth"],'x':[x_Earth],'y':[0]})
mu = 0.01215058162343
tf_1 = 2.77648121127569
tf_2 = 2.75330620148158
# Spatial ICs, but for z=zdot=0
# IC = [
#     1,2,3, 0, 0.19369724986446, 0,
#     1,0,0,0,0,0, # Identity matrix for phi ICs
#     0,1,0,0,0,0,
#     0,0,1,0,0,0,
#     0,0,0,1,0,0,
#     0,0,0,0,1,0,
#     0,0,0,0,0,1
# ]
# halo_prop_half = solve_ivp(spatial_ode, [0, tf_1/2], IC, args=(mu,), method='DOP853', rtol=1e-15,atol=1e-16)


IC_1 = [
    0.82575887090385, 0, 0.08, 0, 0.19369724986446, 0,
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]

IC_1_southern = [
    0.82575887090385, 0, -0.08, 0, 0.19369724986446, 0,
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]

IC_2 = [
    0.82356490862838, 0, 0.04, 0, 0.14924319723734, 0,
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]

IC_2_southern = [
    0.82356490862838, 0, -0.04, 0, 0.14924319723734, 0,
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]

df_orbits = pd.DataFrame()
my_ICs = {'halo_1_northern': IC_1, 'halo_1_southern': IC_1_southern, 'halo_2_northern': IC_2, 'halo_2_southern': IC_2_southern}
for key, value in my_ICs.items():
    orbit_x = value[0]
    orbit_z = value[2]
    orbit_ydot = value[4]
    iterations, tf, arrival_states, converged_initial_states = find_spatial_halfperiod_fixed_z(orbit_x, orbit_z, orbit_ydot, mu, tolerance=1e-12)
    jacobi = calc_spatial_Jacobi(mu, converged_initial_states[0], converged_initial_states[1], converged_initial_states[2], converged_initial_states[3], converged_initial_states[4], converged_initial_states[5])
    orbit_IC_data = pd.DataFrame({
        "label":[key],
        "name": [f'{key} (z={converged_initial_states[2]:.2f})'],
        "iterations":[iterations],
        "tf":[tf],
        "xi":[converged_initial_states[0]-x_L1],
        "x":[converged_initial_states[0]],
        "y":[converged_initial_states[1]],
        "z":[converged_initial_states[2]],
        "vx":[converged_initial_states[3]],
        "vy":[converged_initial_states[4]],
        "vz":[converged_initial_states[5]],
        "jacobi":[jacobi]
    })
    df_orbits = pd.concat([df_orbits, orbit_IC_data], ignore_index=True)
    
halo_table = (
    GT(df_orbits)
    .tab_header(
        title=md(f"Initial States and Properties of the Northern and Southern Halo Orbits <br>({ps}, Lillian Shido)")
    )
    .cols_label(
        tf="{{Half-Period}}<br>[non-dim]",
        x="{{x}}<br>[non-dim]",
        y="{{y}}<br>[non-dim]",
        z="{{z}}<br>[non-dim]",
        vx="{{v_x}}<br>[non-dim]",
        vy="{{v_y}}<br>[non-dim]",
        vz="{{v_z}}<br>[non-dim]",
        jacobi="JC"
    )
    .fmt_number(
        columns=["x","y","z","vx","vy","vz","tf","jacobi"],
        n_sigfig=6
    )
    .cols_align(
        align="center"
    )
    .cols_hide(
        columns=['xi','orbit','iterations','name']
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
halo_table.show()

orbit = pd.DataFrame()
for enum, halo in enumerate(df_orbits.iterrows()):
    x = halo[1]['x']
    IC = [
        halo[1]['x'], halo[1]['y'], halo[1]['z'], halo[1]['vx'], halo[1]['vy'], halo[1]['vz'],
        1,0,0,0,0,0, # Identity matrix for phi ICs
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,1
    ]
    halo_prop = solve_ivp(spatial_ode, [0, 2*halo[1]['tf']], IC, args=(mu,), method='DOP853', rtol=1e-14,atol=1e-16)
    halo_data = pd.DataFrame({
        'name':halo[1]['name'],
        't':halo_prop.t,
        'x':halo_prop.y[0],
        'y':halo_prop.y[1],
        'z':halo_prop.y[2]
    })
    try:
        orbit = pd.concat([orbit, halo_data], ignore_index=True)
    except:
        pdb.set_trace()

x_min = 0.76
x_max = 0.96
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2
title = "Northern and Southern Halo Orbits "
orbits_chart_xy = alt.Chart(orbit).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=[title + f"x-y plane ({ps}, Lillian Shido)"]
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
    color=alt.Color('name:N', scale=alt.Scale(domain=['L1'], range=['darkgreen'])).title(None)
)

earth_loc = alt.Chart(earth).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Earth'], range=['darkblue'])).title(None)
)

orbits_chart_xy_layer = alt.layer(orbits_chart_xy, L1_loc).resolve_scale(color='independent')
orbits_chart_xy_layer.save(f'halo_orbits_x-y_{ps}.png', ppi=200)

orbits_chart_xz = alt.Chart(orbit).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('z:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='z [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=[title + f"x-z plane ({ps}, Lillian Shido)"]
)

orbits_chart_xz_layer = alt.layer(orbits_chart_xz, L1_loc).resolve_scale(color='independent')
orbits_chart_xz_layer.save(f'halo_orbits_x-z_{ps}.png', ppi=200)

orbits_chart_yz = alt.Chart(orbit).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('y:Q', scale=alt.Scale(domain=[-z_lim,z_lim]), axis=alt.Axis(title='z [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('z:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=[title + f"z-y plane ({ps}, Lillian Shido)"]
)

orbits_chart_yz_layer = alt.layer(orbits_chart_yz, L1_loc).resolve_scale(color='independent')
orbits_chart_yz_layer.save(f'halo_orbits_y-z_{ps}.png', ppi=200)

fig = px.line_3d(orbit, x="x", y='y', z='z', color='name',
                 labels={
                     "x": "x [non-dim]",
                     "y": "y [non-dim]",
                     "z": "z [non-dim]"
                 })
fig.update_layout(
    title=dict(text=title + f"({ps}, Lillian Shido)", font=dict(size=30), automargin=True, yref='paper'),
    legend=dict(title=dict(text=None)),
    scene = dict(
        xaxis=dict(range=[x_min, x_max]),
        yaxis=dict(range=[-y_lim,y_lim]),
        zaxis=dict(range=[-z_lim,z_lim]),
        aspectratio = dict(x=1,y=1,z=1),
        aspectmode='manual'  
    ),
)
fig.show()



pdb.set_trace()