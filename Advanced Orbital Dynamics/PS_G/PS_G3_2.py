ps = "G3 part 2"
# Author: Lillian Shido
# Date: 12/14/2025

import pdb
import numpy as np
import pandas as pd
from math import atan2, sin, isclose, degrees
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import solve_ivp
from great_tables import GT, md, html, style, loc, system_fonts
import altair as alt
from copy import deepcopy
import plotly.express as px
import plotly.graph_objects as go

from symbols import xi_symbol, eta_symbol
from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, find_spatial_halfperiod, calc_Jacobi, spatial_ode, calc_initial_velocities

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[y_L1],'z':[0]})
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0],'z':[0]})
earth = pd.DataFrame({'name':["Earth"],'x':[x_Earth],'y':[0],'z':[0]})
mu = 0.01215058162343
# tf = 2.77648121127569
orbit = pd.DataFrame()
# Choose a halo orbit to work with
IC_halo_1 = [
    0.82575887090385, 0, 0.08, 
    0, 0.19369724986446, 0
]
# Choose a Lyapunov to work with
# Calc starting guess values
xi = 0.034
eta = 0
xi_dot_0, eta_dot_0 = calc_initial_velocities(xi, eta, x_L1, y_L1, mu)
starting_x = x_L1 + xi
starting_y = y_L1 + eta
starting_xdot = xi_dot_0
ydot_guess = eta_dot_0
IC_lyapunov = [
    starting_x, starting_y, 0,
    starting_xdot, ydot_guess, 0
]
IC_list = {'halo':IC_halo_1, 'lyapunov':IC_lyapunov}
df_step_off = pd.DataFrame()
for IC_key, IC_values in IC_list.items():
    counter, tf, arrival_states, converged_IC = find_spatial_halfperiod(IC_values[0], IC_values[2], IC_values[4], mu)
    period = tf*2
    df_IC = pd.DataFrame({
        "type":IC_key,
        "period":[2*tf],
        "x":[converged_IC[0]],
        "y":[converged_IC[1]],
        "z":[converged_IC[2]],
        "vx":[converged_IC[3]],
        "vy":[converged_IC[4]],
        "vz":[converged_IC[5]],
    })

    # For plotting purposes:
    converged_1 = [
        converged_IC[0], converged_IC[1], converged_IC[2], converged_IC[3], converged_IC[4], converged_IC[5],
        1,0,0,0,0,0, # Identity matrix for phi ICs
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,1
    ]
    
    halo_1_prop = solve_ivp(spatial_ode, [0, period], converged_1, args=(mu,), rtol=1e-14,atol=1e-16)
    orbit_1 = pd.DataFrame({
        'name': IC_key,
        't':halo_1_prop.t,
        'x':halo_1_prop.y[0],
        'y':halo_1_prop.y[1],
        'z':halo_1_prop.y[2]
    })
    orbit = pd.concat([orbit, orbit_1], ignore_index=True)
    monodromy = halo_1_prop.y[6:42,-1].reshape(6,6)
    eigenvalues, eigenvectors = np.linalg.eig(monodromy)

    # Print the eigs
    for i in range(0,6):
        print(f"lambda_{i+1}: {eigenvalues[i]:.5f}")
        print(f"eigenvectors_{i+1}:")
        print("\\begin{bmatrix}")
        try:
            for x in range(0,6): print(f"{eigenvectors[:,i].flatten()[x]:15.5f}\\\\")
        except:
            pdb.set_trace()
        print("\\end{bmatrix}")

    # Fixed point for stepping-off
    x_fixed = converged_IC[0]
    y_fixed = converged_IC[1]
    z_fixed = converged_IC[2]
    vx_fixed = converged_IC[3]
    vy_fixed = converged_IC[4]
    vz_fixed = converged_IC[5]
    fixed_point = pd.DataFrame({
        'label': [f"Fixed Point @ {eta_symbol}=0"],'name': ['Fixed Point'],'x':[x_fixed],'y':[y_fixed]})

    # Calc the step-off point
    d_val = eigenvalues[1].real #[nondim]
    d_dim = d_val*l_char #[km]
    print(f'step off: {d_dim}[km]')
    print(f'step off: {d_val}[non-dim]')
    for stability_label in ['stable','unstable']:
        if stability_label == 'unstable':
            s = 0
        if stability_label == 'stable':
            s = 1
        nu_W = eigenvectors[:,s]/np.linalg.norm(eigenvectors[:,s][0:3])
        print(nu_W)
        scaled_nu_W = d_val*nu_W
        for manifold_label in ['positive','negative']:
            if manifold_label == 'positive':
                x0 = x_fixed + scaled_nu_W[0].real
                y0 = y_fixed + scaled_nu_W[1].real
                z0 = z_fixed + scaled_nu_W[2].real
                vx0 = vx_fixed + scaled_nu_W[3].real
                vy0 = vy_fixed + scaled_nu_W[4].real
                vz0 = vz_fixed + scaled_nu_W[5].real
            if manifold_label == 'negative':
                x0 = x_fixed - scaled_nu_W[0].real
                y0 = y_fixed - scaled_nu_W[1].real
                z0 = z_fixed - scaled_nu_W[2].real
                vx0 = vx_fixed - scaled_nu_W[3].real
                vy0 = vy_fixed - scaled_nu_W[4].real
                vz0 = vz_fixed - scaled_nu_W[5].real
            step_off_data = pd.DataFrame({
                'type': [IC_key],
                'name':[f'{IC_key}: {stability_label} - {manifold_label}'],
                'stability':stability_label,
                'manifold':manifold_label,
                'x':[x0],
                'y':[y0],
                'z':[z0],
                'vx':[vx0],
                'vy':[vy0],
                'vz':[vz0]
            })
            df_step_off = pd.concat([df_step_off, step_off_data], ignore_index=True)

# Propagate the negative and positive manifolds for stable and unstable 
def crossxEventMoon(t, sv, mu):
    return sv[1] if sv[0] > x_Moon else np.nan
crossxEventMoon.terminal = True

def crossxEventEarth(t, sv, mu):
    return sv[1] if sv[0] < x_Earth else np.nan
crossxEventEarth.terminal = True
events_list = [crossxEventEarth, crossxEventMoon]
df_prop = pd.DataFrame()
prop_length = 2*period
for row in df_step_off.iterrows():
    IC = [
        row[1]['x'], row[1]['y'], row[1]['z'], row[1]['vx'], row[1]['vy'], row[1]['vz'],
        1,0,0,0,0,0, # Identity matrix for phi ICs
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,1
    ]
    if row[1]['stability']=='unstable':
        propagation_time = prop_length
    if row[1]['stability']=='stable':
        propagation_time = -prop_length
    prop = solve_ivp(spatial_ode, [0, propagation_time], IC, events=events_list, args=(mu,), rtol=1e-12,atol=1e-14)
    prop_data = pd.DataFrame({
        'type': row[1]['type'],
        'name': row[1]['name'],
        'stability':row[1]['stability'],
        'manifold':row[1]['manifold'],
        't':prop.t,
        'x':prop.y[0],
        'y':prop.y[1],
        'z':prop.y[2]
    })
    df_prop = pd.concat([df_prop, prop_data], ignore_index=True)

# Calculate the distance between the trajectory and earth or moon
def calc_distance(x,y,z,specified_point):
    point_x = specified_point.iloc[0]['x']
    point_y = specified_point.iloc[0]['y']
    point_z = specified_point.iloc[0]['z']
    distance = np.linalg.norm([x-point_x, y-point_y, z-point_z])
    return distance 

df_prop['moon_distance_nd'] = df_prop.apply(lambda x: calc_distance(x['x'],x['y'],x['z'],moon), axis=1)
df_prop['moon_distance_km'] = df_prop.apply(lambda x: calc_distance(x['x'],x['y'],x['z'],moon)*l_char, axis=1)
df_prop['earth_distance_nd'] = df_prop.apply(lambda x: calc_distance(x['x'],x['y'],x['z'],earth), axis=1)
df_prop['earth_distance_km'] = df_prop.apply(lambda x: calc_distance(x['x'],x['y'],x['z'],earth)*l_char, axis=1)
moon_closest = df_prop[df_prop['moon_distance_nd']==df_prop['moon_distance_nd'].min()]
moon_closest.insert(loc=0,column='closest',value='moon')
earth_closest = df_prop[df_prop['earth_distance_nd']==df_prop['earth_distance_nd'].min()]
earth_closest.insert(loc=0,column='closest',value='earth')
df_closest = pd.concat([earth_closest, moon_closest], ignore_index=True)
df_closest['time_to_approach'] = df_closest['t'].apply(lambda x: abs(x*t_char/3600/24))
closest_table = (
    GT(df_closest)
    .tab_header(
        title=md(f"Closest Approach to Earth and Moon ({ps}, Lillian Shido)")
    )
    .cols_label(
        moon_distance_nd="{{Moon Distance}}<br>[nd]",
        moon_distance_km="{{Moon Distance}}<br>[km]",
        earth_distance_nd="{{Earth Distance}}<br>[nd]",
        earth_distance_km="{{Earth Distance}}<br>[km]",
        time_to_approach="{{Time to Approach}}<br>[days]",
    )
    .fmt_number(
        columns=["earth_distance_km","moon_distance_km","time_to_approach"],
        decimals=3
    )
    .tab_style(
        style=style.fill(color="yellow"),
        locations=loc.body(columns=["moon_distance_km"], rows=[1])
    )
    .tab_style(
        style=style.fill(color="yellow"),
        locations=loc.body(columns=["earth_distance_km"], rows=[0])
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
    .cols_hide(columns=['x','y','z','t','moon_distance_nd','earth_distance_nd'])
)
# closest_table.show()

x_min = -0.8
x_max = 1.2
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2

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

earth_x_y = alt.Chart(earth).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Earth'], range=['darkblue'])).title(None)
)
earth_y_z = alt.Chart(earth).mark_point(filled=True,size=30,clip=True).encode(
    x='y:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Earth'], range=['darkblue'])).title(None)
)
earth_x_z = alt.Chart(earth).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Earth'], range=['darkblue'])).title(None)
)

orbit_x_y = alt.Chart(orbit).mark_line(strokeWidth=2,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['halo','lyapunov'], range=['red','purple'])).title(None),
    order='t'
)
orbit_y_z = alt.Chart(orbit).mark_line(strokeWidth=2,clip=True).encode(
    x='y:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['halo','lyapunov'], range=['red','purple'])).title(None),
    order='t'
)
orbit_x_z = alt.Chart(orbit).mark_line(strokeWidth=2,clip=True).encode(
    x='x:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['halo','lyapunov'], range=['red','purple'])).title(None),
    order='t'
)

# x-y plane chart
x_y_chart = alt.Chart(df_prop).mark_line(clip=True,strokeWidth=2).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    strokeDash=alt.condition(
        alt.datum.type == 'lyapunov',
        alt.value([5, 5]),  # dashed line: 5 pixels  dash + 5 pixels space
        alt.value([0])),
    color=alt.Color('name:N', scale=alt.Scale(scheme='darkmulti')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Stable and Unstable Manifold Trajectories for Halo and Lyapunovs (x-y plane)",f"({ps}, Lillian Shido)"]
)
x_y_chart_layer = alt.layer(x_y_chart, orbit_x_y, moon_x_y, earth_x_y, L1_x_y).resolve_scale(color='independent')
x_y_chart_layer.save(f'stable and unstable manifolds_x_y{ps}.png', ppi=200)

# y-z plane chart
y_z_chart = alt.Chart(df_prop).mark_line(clip=True,strokeWidth=2).encode(
    x=alt.X('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # x=alt.X('y:Q', axis=alt.Axis(title='y [non-dim]')),
    y=alt.Y('z:Q', scale=alt.Scale(domain=[-z_lim,z_lim]), axis=alt.Axis(title='z [non-dim]')),
    # y=alt.Y('z:Q', axis=alt.Axis(title='z [non-dim]')),
    strokeDash=alt.condition(
        alt.datum.type == 'lyapunov',
        alt.value([5, 5]),  # dashed line: 5 pixels  dash + 5 pixels space
        alt.value([0])),
    color=alt.Color('name:N', scale=alt.Scale(scheme='darkmulti')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Stable and Unstable Manifold Trajectories for Halo and Lyapunovs (y-z plane)",f"({ps}, Lillian Shido)"]
)
y_z_chart_layer = alt.layer(y_z_chart, orbit_y_z, moon_y_z, earth_y_z, L1_y_z).resolve_scale(color='independent')
y_z_chart_layer.save(f'stable and unstable manifolds_y_z{ps}.png', ppi=200)

# x-z plane chart
x_z_chart = alt.Chart(df_prop).mark_line(clip=True,strokeWidth=2).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('z:Q', scale=alt.Scale(domain=[-z_lim,z_lim]), axis=alt.Axis(title='z [non-dim]')),
    # y=alt.Y('z:Q', axis=alt.Axis(title='z [non-dim]')),
    strokeDash=alt.condition(
        alt.datum.type == 'lyapunov',
        alt.value([5, 5]),  # dashed line: 5 pixels  dash + 5 pixels space
        alt.value([0])),
    color=alt.Color('name:N', scale=alt.Scale(scheme='darkmulti')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Stable and Unstable Manifold Trajectories for Halo and Lyapunovs (x-z plane)",f"({ps}, Lillian Shido)"]
)
x_z_chart_layer = alt.layer(x_z_chart, orbit_x_z, moon_x_z, earth_x_z, L1_x_z).resolve_scale(color='independent')
x_z_chart_layer.save(f'stable and unstable manifolds_x_z{ps}.png', ppi=200)

# 3-D chart
fig = px.line_3d(df_prop, x="x", y='y', z='z', color='name',
                 labels={
                     "x": "x [non-dim]",
                     "y": "y [non-dim]",
                     "z": "z [non-dim]"
                 })
fig.add_trace(go.Scatter3d(x=L1['x'],y=L1['y'],z=L1['z'],name='L1'))
fig.add_trace(go.Scatter3d(x=earth['x'],y=earth['y'],z=earth['z'],name='earth'))
fig.add_trace(go.Scatter3d(x=moon['x'],y=moon['y'],z=moon['z'],name='moon'))
fig.add_trace(go.Scatter3d(x=orbit[orbit['name']=='halo']['x'],y=orbit[orbit['name']=='halo']['y'],z=orbit[orbit['name']=='halo']['z'],mode='lines',name='halo'))
fig.add_trace(go.Scatter3d(x=orbit[orbit['name']=='lyapunov']['x'],y=orbit[orbit['name']=='lyapunov']['y'],z=orbit[orbit['name']=='lyapunov']['z'],mode='lines',name='lyapunov'))
fig.update_layout(
    title=dict(text=f"Stable and Unstable Manifolds in L1 Halo ({ps}, Lillian Shido)", font=dict(size=30), automargin=True, yref='paper'),
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