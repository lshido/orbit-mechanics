ps = "G1 part a"
# Author: Lillian Shido
# Date: 12/2/2025

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
from methods import system_properties, calc_L1, calc_initial_velocities, find_halfperiod, calc_Jacobi, spatial_ode, calc_poincare_exponents, calc_spatial_monodromy_half

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[y_L1]})
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0]})
earth = pd.DataFrame({'name':["Earth"],'x':[x_Earth],'y':[0]})

tf_1 = 2.77648121127569
tf_2 = 2.75330620148158
# Spatial ICs, but for z=zdot=0
IC_1 = [
    0.82575887090385, 0, 0.0800, 0, 0.19369724986446, 0,
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]

IC_2 = [
    0.82356490862838, 0, 0.0400, 0, 0.14924319723734, 0,
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]

# For plotting purposes:
halo_1_prop = solve_ivp(spatial_ode, [0, tf_1], IC_1, args=(mu,), rtol=1e-12,atol=1e-14)
halo_2_prop = solve_ivp(spatial_ode, [0, tf_2], IC_2, args=(mu,), rtol=1e-12,atol=1e-14)
orbit_1 = pd.DataFrame({
    'name':'halo_1',
    't':halo_1_prop.t,
    'x':halo_1_prop.y[0],
    'y':halo_1_prop.y[1],
    'z':halo_1_prop.y[2]
})
orbit_2 = pd.DataFrame({
    'name':'halo_2',
    't':halo_2_prop.t,
    'x':halo_2_prop.y[0],
    'y':halo_2_prop.y[1],
    'z':halo_2_prop.y[2]
})
orbit = pd.concat([orbit_1, orbit_2], ignore_index=True)

# Check that the orbits are periodic
# 1. Check det(monodromy)=1
monodromy_1 = halo_1_prop.y[6:42,-1].reshape(6,6)
monodromy_2 = halo_2_prop.y[6:42,-1].reshape(6,6)
error_monodromy_1 = abs(np.linalg.det(monodromy_1) - 1)
error_monodromy_2 = abs(np.linalg.det(monodromy_2) - 1)
df_monodromy_error = pd.DataFrame({
    "error_monodromy_1":[error_monodromy_1],
    "error_monodromy_2":[error_monodromy_2],
})
monodromy_error = (
    GT(df_monodromy_error)
    .tab_header(
        title=md(f"Monodromy Matrix Error of Halos<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        error_monodromy_1="{{Error Halo #1}}",
        error_monodromy_2="{{Error Halo #2}}",
    )
    .fmt_scientific(
        columns=["error_monodromy_1","error_monodromy_2"],
        n_sigfig=5
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
monodromy_error.show()

# 2. Check that there is a pair of trivial eigvals
eigvals_1, eigvecs_1 = np.linalg.eig(monodromy_1)
eigvals_2, eigvecs_2 = np.linalg.eig(monodromy_2)
df_eigenvalues = pd.DataFrame({
   "halo":['halo_1', 'halo_2'],
   "lambda_1":[eigvals_1[0],eigvals_2[0]], 
   "lambda_2":[eigvals_1[1],eigvals_2[1]], 
   "lambda_3":[eigvals_1[2],eigvals_2[2]], 
   "lambda_4":[eigvals_1[3],eigvals_2[3]], 
   "lambda_5":[eigvals_1[4],eigvals_2[4]], 
   "lambda_6":[eigvals_1[5],eigvals_2[5]], 
})
eigenvalue_table = (
    GT(df_eigenvalues)
    .tab_header(
        title=md(f"Eigenvalues of Halo Orbits<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        lambda_1="{{:lambda:_1}}",
        lambda_2="{{:lambda:_2}}",
        lambda_3="{{:lambda:_3}}",
        lambda_4="{{:lambda:_4}}",
        lambda_5="{{:lambda:_5}}",
        lambda_6="{{:lambda:_6}}",
    )
    .cols_align(
        align="center"
    )
    .fmt_number(
        columns=["lambda_1","lambda_2","lambda_3","lambda_4","lambda_5","lambda_6"],
        n_sigfig=6
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
eigenvalue_table.show()

# 3. Check that states return to the same points
error_halo_1_states = abs(halo_1_prop.y[0:6,-1] - halo_1_prop.y[0:6,0])
error_halo_2_states = abs(halo_2_prop.y[0:6,-1] - halo_2_prop.y[0:6,0])
df_error_states = pd.DataFrame({
    "halo":['halo_1','halo_2'],
    "error_x":[error_halo_1_states[0],error_halo_2_states[0]],
    "error_y":[error_halo_1_states[1],error_halo_2_states[1]],
    "error_z":[error_halo_1_states[2],error_halo_2_states[2]],
    "error_vx":[error_halo_1_states[3],error_halo_2_states[3]],
    "error_vy":[error_halo_1_states[4],error_halo_2_states[4]],
    "error_vz":[error_halo_1_states[5],error_halo_2_states[5]]
})
error_states_table = (
    GT(df_error_states)
    .tab_header(
        title=md(f"Error between Start and End States of Halos<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        error_x="{{Error x}}",
        error_y="{{Error y}}",
        error_z="{{Error z}}",
        error_vx="{{Error vx}}",
        error_vy="{{Error vy}}",
        error_vz="{{Error vz}}",
    )
    .cols_align(
        align="center"
    )
    .fmt_scientific(
        columns=["error_x","error_y","error_z","error_vx","error_vy","error_vz"],
        n_sigfig=5
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
error_states_table.show()

x_min = 0.75
x_max = 0.95
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2

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
    title=[f"Halo Orbits ({ps}, Lillian Shido)"]
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
    title=[f"Halo Orbits ({ps}, Lillian Shido)"]
)

orbits_chart_xz_layer = alt.layer(orbits_chart_xz, L1_loc).resolve_scale(color='independent')
orbits_chart_xz_layer.save(f'halo_orbits_x-z_{ps}.png', ppi=200)

fig = px.line_3d(orbit, x="x", y='y', z='z', color='name',
                 labels={
                     "x": "x [non-dim]",
                     "y": "y [non-dim]",
                     "z": "z [non-dim]"
                 })
fig.update_layout(
    title=dict(text=f"Halo Orbits #1 and #2 ({ps}, Lillian Shido)", font=dict(size=30), automargin=True, yref='paper'),
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