ps = "G3 part 4 closest"
# Author: Lillian Shido
# Date: 12/14/2025

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
from methods import system_properties, calc_L1, find_spatial_halfperiod, spatial_ode, calc_spatial_Jacobi

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[y_L1],'z':[0]})
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0],'z':[0]})
earth = pd.DataFrame({'name':["Earth"],'x':[x_Earth],'y':[0],'z':[0]})
mu = 0.01215058162343
# tf_1 = 2.77648121127569
orbit = pd.DataFrame()
# Choose a halo orbit to work with
IC_1 = [
    0.82575887090385, 0, 0.08, 0, 0.19369724986446, 0,
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]
counter_1, tf_1, arrival_states_1, converged_IC_1 = find_spatial_halfperiod(IC_1[0], IC_1[2], IC_1[4], mu)
period = tf_1*2
df_IC = pd.DataFrame({
    "period":[2*tf_1],
    "x":[converged_IC_1[0]],
    "y":[converged_IC_1[1]],
    "z":[converged_IC_1[2]],
    "vx":[converged_IC_1[3]],
    "vy":[converged_IC_1[4]],
    "vz":[converged_IC_1[5]],
})
ic_table = (
    GT(df_IC)
    .tab_header(
        title=md(f"Initial Conditions of Selected L1 Halo Orbit ({ps}, Lillian Shido)")
    )
    .cols_label(
        period="Period",
    )
    .fmt_number(
        columns=["period","x","y","z","vx","vy","vz"],
        n_sigfig=6
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
    .cols_move_to_start(columns=["period"])
)
# ic_table.show()

# For plotting purposes:
converged_1 = [
    converged_IC_1[0], converged_IC_1[1], converged_IC_1[2], converged_IC_1[3], converged_IC_1[4], converged_IC_1[5],
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]

halo_1_prop = solve_ivp(spatial_ode, [0, period], converged_1, args=(mu,), rtol=1e-14,atol=1e-16)
orbit_1 = pd.DataFrame({
    'name':'halo_1',
    't':halo_1_prop.t,
    'x':halo_1_prop.y[0],
    'y':halo_1_prop.y[1],
    'z':halo_1_prop.y[2]
})
orbit = pd.concat([orbit, orbit_1], ignore_index=True)
monodromy = halo_1_prop.y[6:42,-1].reshape(6,6)
eigenvalues, eigenvectors = np.linalg.eig(monodromy)

# Identify 20 fixed points around the orbit
fixed_points = dict()
# Generating each fixed point with reduced error!
t0 = 0
monodromy_error = pd.DataFrame()
last_monodromy = deepcopy(monodromy)
last_IC = deepcopy(converged_1)
# Build an array of eigenvectors that is easier to use
my_eigenvectors_list = []
for i in range(0,6):
    my_eigenvectors_list.append(eigenvectors[:,i])
my_eigenvectors = np.array(my_eigenvectors_list)
last_eigenvectors = deepcopy(my_eigenvectors)
for enum, f in enumerate(range(0, 100, 5)):
    fixed_points[enum] = dict()
    fixed_point_prop = solve_ivp(spatial_ode, [t0, f/100*period], last_IC, args=(mu,), rtol=1e-12,atol=1e-14)
    # Set the new t0 
    t0 = fixed_point_prop.t[-1]
    stm_fixed_point = fixed_point_prop.y[6:42,-1].reshape(6,6)
    fixed_points[enum]['fp_x'] = fixed_point_prop.y[0,-1]
    fixed_points[enum]['fp_y'] = fixed_point_prop.y[1,-1]
    fixed_points[enum]['fp_z'] = fixed_point_prop.y[2,-1]
    fixed_points[enum]['fp_vx'] = fixed_point_prop.y[3,-1]
    fixed_points[enum]['fp_vy'] = fixed_point_prop.y[4,-1]
    fixed_points[enum]['fp_vz'] = fixed_point_prop.y[5,-1]
    # Update the ICs for the next round
    last_IC = [
        fixed_point_prop.y[0,-1], fixed_point_prop.y[1,-1], fixed_point_prop.y[2,-1], fixed_point_prop.y[3,-1], fixed_point_prop.y[4,-1], fixed_point_prop.y[5,-1],
        1,0,0,0,0,0, # Identity matrix for phi ICs
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,1
    ]
    fixed_point_eigenvectors_list = []
    fixed_point_eigenvalues_list = []
    # Get the monodromy matrix at fixed point
    fixed_point_monodromy = stm_fixed_point @ last_monodromy @ np.linalg.inv(stm_fixed_point)
    eigvals, eigvecs = np.linalg.eig(fixed_point_monodromy)
    fixed_points[enum]['eigvals'] = eigvals
    error = abs(np.linalg.det(fixed_point_monodromy)-1)
    print(f"Fixed Point: {enum}, Monodromy error: {error}")
    monodromy_error_data = pd.DataFrame({
        "Fixed Point": [enum+1],
        "Error": [error]
    })
    monodromy_error = pd.concat([monodromy_error, monodromy_error_data], ignore_index=True)
    # Calculate the eigenvectors at the fixed point
    for i in range(0,6):
        fixed_point_eigenvectors = (stm_fixed_point @ last_eigenvectors[i])/np.linalg.norm(stm_fixed_point @ last_eigenvectors[i])
        fixed_point_eigenvectors_list.append(fixed_point_eigenvectors)
        # Use the monodromy matrix at fixed point to calc eigenvalues for each eigenvector
        fixed_point_eigenvalues = np.divide(fixed_point_monodromy @ fixed_point_eigenvectors, fixed_point_eigenvectors)
        fixed_point_eigenvalues_list.append(fixed_point_eigenvalues)
    last_monodromy = fixed_point_monodromy
    fixed_point_eigenvalues = np.array(fixed_point_eigenvalues_list)
    fixed_point_eigenvectors = np.array(fixed_point_eigenvectors_list)
    last_eigenvectors = fixed_point_eigenvectors
    fixed_points[enum]['eigenvectors'] = fixed_point_eigenvectors
    fixed_points[enum]['eigenvalues'] = fixed_point_eigenvalues

# Calc the step-off point
df_step_off = pd.DataFrame()
d_val = eigenvalues[1].real #[nondim]
d_dim = d_val*l_char #[km]
print(f'step off: {d_dim}[km]')
print(f'step off: {d_val}[non-dim]')

for num, point in fixed_points.items():
    # Fixed point for stepping-off
    x_fixed = point['fp_x']
    y_fixed = point['fp_y']
    z_fixed = point['fp_z']
    vx_fixed = point['fp_vx']
    vy_fixed = point['fp_vy']
    vz_fixed = point['fp_vz']
    for stability_label in ['stable','unstable']:
        if stability_label == 'unstable':
            s = 0
        if stability_label == 'stable':
            s = 1
        nu_W = point['eigenvectors'][s]/np.linalg.norm(point['eigenvectors'][s][0:3])
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
            jacobi = calc_spatial_Jacobi(mu, x0, y0, z0, vx0, vy0, vz0)
            step_off_data = pd.DataFrame({
                'orbit': [num],
                'name':[f'{num}: {stability_label} - {manifold_label}'],
                'stability':stability_label,
                'manifold':manifold_label,
                'x':[x0],
                'y':[y0],
                'z':[z0],
                'vx':[vx0],
                'vy':[vy0],
                'vz':[vz0],
                'jacobi':[jacobi]
            })
            df_step_off = pd.concat([df_step_off, step_off_data], ignore_index=True)
step_off_table = (
    GT(df_step_off)
    .tab_header(
        title=md(f"Initial Conditions of Step-Off Points to Propagate Manifold Trajectories of Selected L1 Halo ({ps}, Lillian Shido)")
    )
    .fmt_number(
        columns=["x","y","z","vx","vy","vz",'jacobi'],
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
step_off_table.show()


# Propagate the negative and positive manifolds for stable and unstable 
def crossxEventMoon(t, sv, mu):
    return sv[1] if sv[0] > x_Moon else np.nan
crossxEventMoon.terminal = True

def crossxEventEarth(t, sv, mu):
    return sv[1] if sv[0] < x_Earth else np.nan
crossxEventEarth.terminal = True
events_list = [crossxEventEarth, crossxEventMoon]
df_prop = pd.DataFrame()
prop_length = 5*period
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
closest_table.show()

x_min = 0.75
x_max = 0.95
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

earth_closest_fp = pd.DataFrame({
    'name':["Delivers Closest Earth Pass"],
    'x':[fixed_points[18]['fp_x']],
    'y':[fixed_points[18]['fp_y']],
    'z':[fixed_points[18]['fp_z']]
})
earth_closest_x_y = alt.Chart(earth_closest_fp).mark_point(shape='cross',filled=True,size=100,clip=True, angle=45).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Delivers Closest Earth Pass'], range=['magenta'])).title(None)
)
earth_closest_y_z = alt.Chart(earth_closest_fp).mark_point(shape='cross',filled=True,size=100,clip=True, angle=45).encode(
    x='y:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Delivers Closest Earth Pass'], range=['magenta'])).title(None)
)
earth_closest_x_z = alt.Chart(earth_closest_fp).mark_point(shape='cross',filled=True,size=100,clip=True, angle=45).encode(
    x='x:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Delivers Closest Earth Pass'], range=['magenta'])).title(None)
)

moon_closest_fp = pd.DataFrame({
    'name':["Delivers Closest Moon Pass"],
    'x':[fixed_points[8]['fp_x']],
    'y':[fixed_points[8]['fp_y']],
    'z':[fixed_points[8]['fp_z']]
})
moon_closest_x_y = alt.Chart(moon_closest_fp).mark_point(shape='cross',filled=True,size=100,clip=True, angle=45).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Delivers Closest Moon Pass'], range=['red'])).title(None)
)
moon_closest_y_z = alt.Chart(moon_closest_fp).mark_point(shape='cross',filled=True,size=100,clip=True, angle=45).encode(
    x='y:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Delivers Closest Moon Pass'], range=['red'])).title(None)
)
moon_closest_x_z = alt.Chart(moon_closest_fp).mark_point(shape='cross',filled=True,size=100,clip=True, angle=45).encode(
    x='x:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Delivers Closest Moon Pass'], range=['red'])).title(None)
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

orbit_x_y = alt.Chart(orbit).mark_line(strokeWidth=1,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['halo_1'], range=['black'])).title(None),
    order='t'
)
orbit_y_z = alt.Chart(orbit).mark_line(strokeWidth=1,clip=True).encode(
    x='y:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['halo_1'], range=['black'])).title(None),
    order='t'
)
orbit_x_z = alt.Chart(orbit).mark_line(strokeWidth=1,clip=True).encode(
    x='x:Q',
    y='z:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['halo_1'], range=['black'])).title(None),
    order='t'
)

# x-y plane chart
x_y_chart = alt.Chart(df_prop).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    strokeDash=alt.condition(
        alt.datum.stability == 'unstable',
        alt.value([4, 4]),  # dashed line: 5 pixels  dash + 5 pixels space
        alt.value([0])),
    color=alt.Color('name:N', scale=alt.Scale(scheme='darkmulti')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Stable and Unstable Manifold Trajectories for Halo Orbit (x-y plane)",f"({ps}, Lillian Shido)"]
)
x_y_chart_layer = alt.layer(x_y_chart, orbit_x_y, moon_x_y, earth_x_y, L1_x_y, earth_closest_x_y, moon_closest_x_y).resolve_scale(color='independent')
x_y_chart_layer.save(f'stable and unstable manifolds_x_y{ps}.png', ppi=200)

# y-z plane chart
y_z_chart = alt.Chart(df_prop).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # x=alt.X('y:Q', axis=alt.Axis(title='y [non-dim]')),
    y=alt.Y('z:Q', scale=alt.Scale(domain=[-z_lim,z_lim]), axis=alt.Axis(title='z [non-dim]')),
    # y=alt.Y('z:Q', axis=alt.Axis(title='z [non-dim]')),
    strokeDash=alt.condition(
        alt.datum.stability == 'unstable',
        alt.value([4, 4]),  # dashed line: 5 pixels  dash + 5 pixels space
        alt.value([0])),
    color=alt.Color('name:N', scale=alt.Scale(scheme='darkmulti')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Stable and Unstable Manifold Trajectories for Halo Orbit (y-z plane)",f"({ps}, Lillian Shido)"]
)
y_z_chart_layer = alt.layer(y_z_chart, orbit_y_z, moon_y_z, earth_y_z, L1_y_z, earth_closest_y_z, moon_closest_y_z).resolve_scale(color='independent')
y_z_chart_layer.save(f'stable and unstable manifolds_y_z{ps}.png', ppi=200)

# x-z plane chart
x_z_chart = alt.Chart(df_prop).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('z:Q', scale=alt.Scale(domain=[-z_lim,z_lim]), axis=alt.Axis(title='z [non-dim]')),
    # y=alt.Y('z:Q', axis=alt.Axis(title='z [non-dim]')),
    strokeDash=alt.condition(
        alt.datum.stability == 'unstable',
        alt.value([4, 4]),  # dashed line: 5 pixels  dash + 5 pixels space
        alt.value([0])),
    color=alt.Color('name:N', scale=alt.Scale(scheme='darkmulti')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Stable and Unstable Manifold Trajectories for Halo Orbit (x-z plane)",f"({ps}, Lillian Shido)"]
)
x_z_chart_layer = alt.layer(x_z_chart, orbit_x_z, moon_x_z, earth_x_z, L1_x_z, earth_closest_x_z, moon_closest_x_z).resolve_scale(color='independent')
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
fig.add_trace(go.Scatter3d(x=orbit[orbit['name']=='halo_1']['x'],y=orbit[orbit['name']=='halo_1']['y'],z=orbit[orbit['name']=='halo_1']['z'],mode='lines',name='halo_1'))
fig.add_trace(go.Scatter3d(x=moon['x'],y=moon['y'],z=moon['z'],name='moon'))
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
closest_table.show()

pdb.set_trace()