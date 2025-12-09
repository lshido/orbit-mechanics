ps = "G2 part a"
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
import plotly.graph_objects as go

from symbols import xi_symbol, eta_symbol
from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, find_spatial_halfperiod, find_spatial_halfperiod_fixed_z, calc_spatial_Jacobi, spatial_ode, calc_spatial_monodromy_half

# Properties of the system
mu_ignore, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
mu = 0.01215058162343
x_L1, y_L1 = calc_L1(mu, a_Moon)
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0]})
earth = pd.DataFrame({'name':["Earth"],'x':[x_Earth],'y':[0]})
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[y_L1],'z':[0]})
# tf_1 = 2.77648121127569
# tf_2 = 2.75330620148158

IC_1 = [
    0.82575887090385, 0, 0.08, 0, 0.19369724986446, 0,
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

counter_1, tf_1, arrival_states_1, converged_IC_1 = find_spatial_halfperiod(IC_1[0], IC_1[2], IC_1[4], mu)
counter_2, tf_2, arrival_states_2, converged_IC_2 = find_spatial_halfperiod(IC_2[0], IC_2[2], IC_2[4], mu)

vy_x_slope = (converged_IC_1[4]-converged_IC_2[4])/(converged_IC_1[0]-converged_IC_2[0])
z_x_slope = (converged_IC_1[2]-converged_IC_2[2])/(converged_IC_1[0]-converged_IC_2[0])
x_z_slope = (converged_IC_1[0]-converged_IC_2[0])/(converged_IC_1[2]-converged_IC_2[2])
vy_z_slope = (converged_IC_1[4]-converged_IC_2[4])/(converged_IC_1[2]-converged_IC_2[2])

# Configuratiby delta_z = -0.0001 # step in x
delta_z = -0.002 # z in x

# Start with halo_2 since it's the inside one 
orbit_x = converged_IC_1[0]
orbit_z = converged_IC_1[2]
orbit_ydot = converged_IC_1[4]
df_orbits = pd.DataFrame()
n_orbits=80
x_switch = 0.82345
for orbit in range(n_orbits):
    print(f"starting x0:{orbit_x}, starting z0: {orbit_z}, starting ydot: {orbit_ydot}")
    if orbit_z == 0:
        iterations, tf, arrival_states, converged_initial_states = find_spatial_halfperiod(orbit_x, orbit_z, orbit_ydot, mu, tolerance=1e-12)
    else:    
        iterations, tf, arrival_states, converged_initial_states = find_spatial_halfperiod_fixed_z(orbit_x, orbit_z, orbit_ydot, mu, tolerance=1e-12)
        print(f"found x0:{converged_initial_states[0]}, found z0: {converged_initial_states[2]}, found ydot: {converged_initial_states[4]}")
        if 1e-3 < orbit_z < 0.002:
            delta_z = -0.001
        if orbit_z < 1e-3:
            delta_z = -orbit_z
    orbit_x = converged_initial_states[0] + delta_z*x_z_slope # Step the x by delta_z
    orbit_z = converged_initial_states[2] + delta_z # Step the z by delta_z
    orbit_ydot = converged_initial_states[4] + delta_z*vy_z_slope
    # Calc Jacobi
    jacobi = calc_spatial_Jacobi(mu, converged_initial_states[0], converged_initial_states[1], converged_initial_states[2], converged_initial_states[3], converged_initial_states[4], converged_initial_states[5])
    # Compile data
    orbit_IC_data = pd.DataFrame({
        "name": [f'z={converged_initial_states[2]:.5f}'],
        "orbit":[orbit],
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
    if converged_initial_states[2] <= 0:
        break

halo_intersect = (
    GT(df_orbits[df_orbits['z']==0])
    .tab_header(
        title=md(f"Initial States and Properties of Halo that Intersects x-y Plane<br>({ps}, Lillian Shido)")
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
# halo_intersect.show()

orbit = pd.DataFrame()
halo = df_orbits.iloc[-1]
IC = [
    halo['x'], halo['y'], halo['z'], halo['vx'], halo['vy'], halo['vz'],
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]
try:
    halo_prop = solve_ivp(spatial_ode, [0, 2*halo['tf']], IC, args=(mu,), method='DOP853', rtol=1e-14,atol=1e-16)
except:
    pdb.set_trace()
monodromy = halo_prop.y[6:42,-1].reshape(6,6)
halo_data = pd.DataFrame({
    'name':halo['z'],
    't':halo_prop.t,
    'x':halo_prop.y[0],
    'y':halo_prop.y[1],
    'z':halo_prop.y[2]
})
try:
    orbit = pd.concat([orbit, halo_data], ignore_index=True)
except:
    pdb.set_trace()

# Check that the orbits are periodic
# 1. Check det(monodromy)=1
error_monodromy_1 = abs(np.linalg.det(monodromy) - 1)
df_monodromy_error = pd.DataFrame({
    "error_monodromy_1":[error_monodromy_1],
})
monodromy_error = (
    GT(df_monodromy_error)
    .tab_header(
        title=md(f"Monodromy Matrix Error of Halo Orbit at z=0<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        error_monodromy_1="{{Error}}",
    )
    .fmt_scientific(
        columns=["error_monodromy_1"],
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
eigvals_1, eigvecs_1 = np.linalg.eig(monodromy)
df_eigenvalues = pd.DataFrame({
   "halo":['halo_1'],
   "lambda_1":[f"{eigvals_1[0]:.7g}"], 
   "lambda_2":[f"{eigvals_1[1]:.7g}"], 
   "lambda_3":[f"{eigvals_1[2]:.7g}"], 
   "lambda_4":[f"{eigvals_1[3]:.7g}"], 
   "lambda_5":[f"{eigvals_1[4]:.7g}"], 
   "lambda_6":[f"{eigvals_1[5]:.7g}"], 
})
eigenvalue_table = (
    GT(df_eigenvalues)
    .tab_header(
        title=md(f"Eigenvalues of Halo Orbit at z=0<br>({ps}, Lillian Shido)")
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
    .cols_hide(
        columns=['halo']
    )
    # .fmt_number(
    #     columns=["lambda_1","lambda_2","lambda_3","lambda_4","lambda_5","lambda_6"],
    #     n_sigfig=6
    # )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
eigenvalue_table.show()

df_eigenvalues_abs = pd.DataFrame({
    "eig_1_abs":[abs(eigvals_1[0])],
    "eig_2_abs":[abs(eigvals_1[1])],
    "eig_3_abs":[abs(eigvals_1[2])],
    "eig_4_abs":[abs(eigvals_1[3])],
    "eig_5_abs":[abs(eigvals_1[4])],
    "eig_6_abs":[abs(eigvals_1[5])]
})
eig_abs_table = (
    GT(df_eigenvalues_abs)
    .tab_header(
        title=md(f"Absolute Values of the Eigenvalues of Halo Orbit at z=0<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        eig_1_abs="{{| :lambda:_1 |}}",
        eig_2_abs="{{| :lambda:_2 |}}",
        eig_3_abs="{{| :lambda:_3 |}}",
        eig_4_abs="{{| :lambda:_4 |}}",
        eig_5_abs="{{| :lambda:_5 |}}",
        eig_6_abs="{{| :lambda:_6 |}}"
    )
    .fmt_number(
        columns=["eig_1_abs","eig_2_abs","eig_3_abs","eig_4_abs","eig_5_abs","eig_6_abs"],
        n_sigfig=6
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
eig_abs_table.show()

# # 3. Check that states return to the same points
# error_halo_1_states = abs(halo_1_prop.y[0:6,-1] - halo_1_prop.y[0:6,0])
# error_halo_2_states = abs(halo_2_prop.y[0:6,-1] - halo_2_prop.y[0:6,0])
# df_error_states = pd.DataFrame({
#     "halo":['halo_1','halo_2'],
#     "error_x":[error_halo_1_states[0],error_halo_2_states[0]],
#     "error_y":[error_halo_1_states[1],error_halo_2_states[1]],
#     "error_z":[error_halo_1_states[2],error_halo_2_states[2]],
#     "error_vx":[error_halo_1_states[3],error_halo_2_states[3]],
#     "error_vy":[error_halo_1_states[4],error_halo_2_states[4]],
#     "error_vz":[error_halo_1_states[5],error_halo_2_states[5]]
# })
# error_states_table = (
#     GT(df_error_states)
#     .tab_header(
#         title=md(f"Error between Start and End States of Halos<br>({ps}, Lillian Shido)")
#     )
#     .cols_label(
#         error_x="{{Error x}}",
#         error_y="{{Error y}}",
#         error_z="{{Error z}}",
#         error_vx="{{Error vx}}",
#         error_vy="{{Error vy}}",
#         error_vz="{{Error vz}}",
#     )
#     .cols_align(
#         align="center"
#     )
#     .fmt_scientific(
#         columns=["error_x","error_y","error_z","error_vx","error_vy","error_vz"],
#         n_sigfig=5
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# # error_states_table.show()

x_min = 0.76
x_max = 0.96
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2

orbits_chart_xy = alt.Chart(orbit).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:Q', scale=alt.Scale(scheme='rainbow')).title("z"),
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
# orbits_chart_xy_layer.save(f'halo_orbits_x-y_{ps}.png', ppi=200)

orbits_chart_xz = alt.Chart(orbit).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('z:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='z [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:Q', scale=alt.Scale(scheme='rainbow')).title('z'),
    order='t'
).properties(
    width=400,
    height=400,
    title=[f"Halo Orbits x-z plane ({ps}, Lillian Shido)"]
)

orbits_chart_xz_layer = alt.layer(orbits_chart_xz, L1_loc).resolve_scale(color='independent')
# orbits_chart_xz_layer.save(f'halo_orbits_x-z_{ps}.png', ppi=200)

orbits_chart_yz = alt.Chart(orbit).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('z:Q', scale=alt.Scale(domain=[-z_lim,z_lim]), axis=alt.Axis(title='z [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:Q', scale=alt.Scale(scheme='rainbow')).title('z'),
    order='t'
).properties(
    width=400,
    height=400,
    title=[f"Halo Orbits z-y plane ({ps}, Lillian Shido)"]
)

orbits_chart_yz_layer = alt.layer(orbits_chart_yz, L1_loc).resolve_scale(color='independent')
# orbits_chart_yz_layer.save(f'halo_orbits_y-z_{ps}.png', ppi=200)

fig = px.line_3d(orbit, x="x", y='y', z='z', color='name',
                 labels={
                     "x": "x [non-dim]",
                     "y": "y [non-dim]",
                     "z": "z [non-dim]"
                 })
fig.add_trace(go.Scatter3d(x=L1['x'],y=L1['y'],z=L1['z'],name='L1'))
fig.update_layout(
    title=dict(text=f"3 Halo Orbits near L1 ({ps}, Lillian Shido)", font=dict(size=30), automargin=True, yref='paper'),
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

pdb.set_trace()