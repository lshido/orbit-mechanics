ps = "G2 part d"
# Author: Lillian Shido
# Date: 12/12/2025

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
from methods import system_properties, calc_L2, find_spatial_halfperiod, find_spatial_halfperiod_fixed_z, calc_spatial_Jacobi, spatial_ode, calc_spatial_monodromy_half, calc_stability_recip

# Properties of the system
mu_ignore, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
mu = 0.01215058162343
x_L2, y_L2 = calc_L2(mu, a_Moon)
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0]})
earth = pd.DataFrame({'name':["Earth"],'x':[x_Earth],'y':[0]})
L2 = pd.DataFrame({'name':["L2"],'x':[x_L2],'y':[y_L2],'z':[0]})
# tf_1 = 2.77648121127569
# tf_2 = 2.75330620148158

counter_1, tf_1, arrival_states_1, converged_IC_1 = find_spatial_halfperiod(1.181682, 0, -0.161623, mu)
counter_2, tf_2, arrival_states_2, converged_IC_2 = find_spatial_halfperiod(1.219682, 0, -0.42654, mu)
counter_3, tf_2, arrival_states_3, converged_IC_3 = find_spatial_halfperiod(1.275682, 0, -0.539256, mu)


# Configuratiby delta_z = -0.0001 # step in x
first_delta_z = 0.001
delta_z = 0.001 # z in x

# Start with halo_2 since it's the inside one 
orbit_x = converged_IC_1[0]
orbit_z = converged_IC_1[2]
orbit_ydot = converged_IC_1[4]
df_orbits = pd.DataFrame()
n_orbits=100
for orbit in range(n_orbits):
    print(f"starting x0:{orbit_x}, starting z0: {orbit_z}, starting ydot: {orbit_ydot}")
    if orbit_z == 0:
        iterations, tf, arrival_states, converged_initial_states = find_spatial_halfperiod(orbit_x, orbit_z, orbit_ydot, mu, tolerance=1e-12)
    else:    
        iterations, tf, arrival_states, converged_initial_states = find_spatial_halfperiod_fixed_z(orbit_x, orbit_z, orbit_ydot, mu, tolerance=1e-12)
        print(f"found x0:{converged_initial_states[0]}, found z0: {converged_initial_states[2]}, found ydot: {converged_initial_states[4]}")
        if 1e-3 < orbit_z < 0.002:
            delta_z = 0.001
        # if orbit_z < 1e-3:
        #     delta_z = -orbit_z
    if orbit <= 1: # For the first two calculations, naively guess x and ydot is the same
        orbit_x = converged_initial_states[0]
        orbit_ydot = converged_initial_states[4]
    else:
        # Calc the slope from the last pair of x0 and vy0
        vy_z_slope = (df_orbits.iloc[-1]['vy']-df_orbits.iloc[-2]['vy'])/(df_orbits.iloc[-1]['z']-df_orbits.iloc[-2]['z'])
        x_z_slope = (df_orbits.iloc[-1]['x']-df_orbits.iloc[-2]['x'])/(df_orbits.iloc[-1]['z']-df_orbits.iloc[-2]['z'])
        orbit_x = converged_initial_states[0] + delta_z*x_z_slope # Step the x by delta_z
        orbit_ydot = converged_initial_states[4] + delta_z*vy_z_slope
    # Calc Jacobi
    orbit_z = converged_initial_states[2] + delta_z # Step the z by delta_z
    jacobi = calc_spatial_Jacobi(mu, converged_initial_states[0], converged_initial_states[1], converged_initial_states[2], converged_initial_states[3], converged_initial_states[4], converged_initial_states[5])
    # Compile data
    orbit_IC_data = pd.DataFrame({
        "name": [f'z={converged_initial_states[2]:.5f}'],
        "orbit":[orbit],
        "iterations":[iterations],
        "tf":[tf],
        "xi":[converged_initial_states[0]-x_L2],
        "x":[converged_initial_states[0]],
        "y":[converged_initial_states[1]],
        "z":[converged_initial_states[2]],
        "vx":[converged_initial_states[3]],
        "vy":[converged_initial_states[4]],
        "vz":[converged_initial_states[5]],
        "jacobi":[jacobi]
    })
    df_orbits = pd.concat([df_orbits, orbit_IC_data], ignore_index=True)
    if converged_initial_states[2] < 0:
        break

df_eigenvalues_string = pd.DataFrame()
df_eigenvalues_abs = pd.DataFrame()
df_stability = pd.DataFrame()
orbit = pd.DataFrame()
for enum, halo in enumerate(df_orbits.iterrows()):
    tf = halo[1]['tf']
    x = halo[1]['x']
    z = halo[1]['z']
    jacobi = halo[1]['jacobi']
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
        'name':halo[1]['z'],
        't':halo_prop.t,
        'x':halo_prop.y[0],
        'y':halo_prop.y[1],
        'z':halo_prop.y[2]
    })
    try:
        orbit = pd.concat([orbit, halo_data], ignore_index=True)
    except:
        pdb.set_trace()
    monodromy = halo_prop.y[6:42,-1].reshape(6,6)
    eigenvalues = np.linalg.eigvals(monodromy)
    print(f"Orbit {enum} monodromy det: {np.linalg.det(monodromy)}")
    eigenvalues_string_data = pd.DataFrame({
        "orbit":[enum],
        "x":[x],
        "z":[z],
        "period":[2*tf],
        "jacobi":[jacobi],
        "eig_1_string":[f"{eigenvalues[0]:.5g}"],
        "eig_2_string":[f"{eigenvalues[1]:.5g}"],
        "eig_3_string":[f"{eigenvalues[2]:.5g}"],
        "eig_4_string":[f"{eigenvalues[3]:.5g}"],
        "eig_5_string":[f"{eigenvalues[4]:.5g}"],
        "eig_6_string":[f"{eigenvalues[5]:.5g}"]
    })
    df_eigenvalues_string = pd.concat([df_eigenvalues_string, eigenvalues_string_data], ignore_index=True)
    eigenvalues_abs_data = pd.DataFrame({
        "orbit":[enum],
        "x":[x],
        "z":[z],
        "period":[2*tf],
        "jacobi":[jacobi],
        "eig_1_abs":[abs(eigenvalues[0])],
        "eig_2_abs":[abs(eigenvalues[1])],
        "eig_3_abs":[abs(eigenvalues[2])],
        "eig_4_abs":[abs(eigenvalues[3])],
        "eig_5_abs":[abs(eigenvalues[4])],
        "eig_6_abs":[abs(eigenvalues[5])]
    })
    df_eigenvalues_abs = pd.concat([df_eigenvalues_abs, eigenvalues_abs_data], ignore_index=True)
    for n in range(0,6):
        if n==0:
            label = 'in-plane'
        elif n==5:
            label = 'out-of-plane'
        else:
            continue
        stability_data = pd.DataFrame({
            "orbit":[enum],
            "z":[z],
            "period":[2*tf],
            "jacobi":[jacobi],
            "label":[f'{label}'],
            f"stability":[calc_stability_recip(eigenvalues[n]).real]
        })
        df_stability = pd.concat([df_stability, stability_data], ignore_index=True)

eigenvalue_table = (
    GT(df_eigenvalues_string)
    .tab_header(
        title=md(f"Eigenvalues of Halo Orbits<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        x='x',
        z="z",
        period="Period",
        jacobi="Jacobi",
        eig_1_string="{{:lambda:_1}}",
        eig_2_string="{{:lambda:_2}}",
        eig_3_string="{{:lambda:_3}}",
        eig_4_string="{{:lambda:_4}}",
        eig_5_string="{{:lambda:_5}}",
        eig_6_string="{{:lambda:_6}}",
    )
    .cols_align(
        align="center"
    )
    .fmt_number(
        columns=["jacobi", "period"],
        n_sigfig=6
    )
    .fmt_number(
        columns=['z'],
        decimals=3
    )
    .fmt_number(
        columns=['x'],
        decimals=6
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
eigenvalue_table.show()

eigenvalue_abs_table = (
    GT(df_eigenvalues_abs)
    .tab_header(
        title=md(f"Normalized Eigenvalues of Halo Orbits<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        x='x',
        z="z",
        period="Period",
        jacobi="Jacobi",
        eig_1_abs="{{| :lambda:_1 |}}",
        eig_2_abs="{{| :lambda:_2 |}}",
        eig_3_abs="{{| :lambda:_3 |}}",
        eig_4_abs="{{| :lambda:_4 |}}",
        eig_5_abs="{{| :lambda:_5 |}}",
        eig_6_abs="{{| :lambda:_6 |}}",
    )
    .cols_align(
        align="center"
    )
    .fmt_number(
        columns=["jacobi", "period", "eig_1_abs","eig_2_abs","eig_3_abs","eig_4_abs","eig_5_abs","eig_6_abs"],
        n_sigfig=6
    )
    .fmt_number(
        columns=['z'],
        decimals=3
    )
    .fmt_number(
        columns=['x'],
        decimals=6
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
eigenvalue_abs_table.show()


x_min = 1.00
x_max = 1.25
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2

orbits_chart_xy = alt.Chart(orbit[~np.isclose(orbit['name'],0)]).mark_line(clip=True,strokeWidth=1).encode(
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

# scale = alt.Scale(domain=['L2','Moon'], range=['darkblue','gray'])
moon_loc = alt.Chart(moon).mark_point(filled=True,size=50,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Moon'], range=['gray'])).title(None)
)

L2_loc = alt.Chart(L2).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L2'], range=['darkgreen'])).title(None)
)

earth_loc = alt.Chart(earth).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Earth'], range=['darkblue'])).title(None)
)

orbits_chart_xy_layer = alt.layer(orbits_chart_xy, L2_loc).resolve_scale(color='independent')
orbits_chart_xy_layer.save(f'halo_orbits_x-y_{ps}.png', ppi=200)

orbits_chart_xz = alt.Chart(orbit[~np.isclose(orbit['name'],0)]).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('z:Q', scale=alt.Scale(domain=[-z_lim,z_lim]), axis=alt.Axis(title='z [non-dim]')),
    # y=alt.Y('z:Q', axis=alt.Axis(title='z [non-dim]')),
    color=alt.Color('name:Q', scale=alt.Scale(scheme='rainbow')).title('z'),
    order='t'
).properties(
    width=400,
    height=400,
    title=[f"Halo Orbits x-z plane ({ps}, Lillian Shido)"]
)

orbits_chart_xz_layer = alt.layer(orbits_chart_xz, L2_loc).resolve_scale(color='independent')
orbits_chart_xz_layer.save(f'halo_orbits_x-z_{ps}.png', ppi=200)

orbits_chart_yz = alt.Chart(orbit[~np.isclose(orbit['name'],0)]).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # x=alt.X('y:Q', axis=alt.Axis(title='y [non-dim]')),
    y=alt.Y('z:Q', scale=alt.Scale(domain=[-z_lim,z_lim]), axis=alt.Axis(title='z [non-dim]')),
    # y=alt.Y('z:Q', axis=alt.Axis(title='z [non-dim]')),
    color=alt.Color('name:Q', scale=alt.Scale(scheme='rainbow')).title('z'),
    order='t'
).properties(
    width=400,
    height=400,
    title=[f"Halo Orbits y-z plane ({ps}, Lillian Shido)"]
)

# orbits_chart_yz_layer = alt.layer(orbits_chart_yz, L2_loc).resolve_scale(color='independent')
orbits_chart_yz.save(f'halo_orbits_y-z_{ps}.png', ppi=200)

fig = px.line_3d(orbit[~np.isclose(orbit['name'],0)], x="x", y='y', z='z', color='name',
                 labels={
                     "x": "x [non-dim]",
                     "y": "y [non-dim]",
                     "z": "z [non-dim]"
                 })
fig.add_trace(go.Scatter3d(x=L2['x'],y=L2['y'],z=L2['z'],name='L2'))
fig.update_layout(
    title=dict(text=f"Halo Orbits near L2 ({ps}, Lillian Shido)", font=dict(size=30), automargin=True, yref='paper'),
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