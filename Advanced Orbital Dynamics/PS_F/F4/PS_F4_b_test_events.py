ps = "F4 part b"
# Author: Lillian Shido
# Date: 11/25/2025

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
xi = 0.036
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

# -0.174029
# iterations, tf, arrival_states, converged_initial_states = find_halfperiod(starting_x, ydot_guess, mu, tolerance=1e-12)
iterations, tf, arrival_states, converged_initial_states = find_halfperiod(starting_x, -0.248204, mu, tolerance=1e-12)
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
    'name':'Periodic Orbit',
    't':full_period_prop.t,
    'x':full_period_prop.y[0],
    'y':full_period_prop.y[1]
})

x_min = 0.7
x_max = 0.9
y_lim = (x_max-x_min)/2

orbits_chart = alt.Chart(orbit).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None).legend(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Periodic Orbit",f"for the orbit {xi_symbol}={ps}, {eta_symbol}=0 ({ps}, Lillian Shido)"]
)
orbits_chart.save(f'orbit_{ps}.png', ppi=200)


# Get eigs and check they match
eigenvalues, eigenvectors = np.linalg.eig(monodromy)
for i in range(0,6):
    print(f"Check eigenset {i+1}")
    if np.isclose(monodromy*eigenvectors[:,i],eigenvalues[i]*eigenvectors[:,i],atol=1e-8).all():
        print("PASS!")
    else:
        print("FAIL.")
# Classify Eigenvalues and put them into a dataframe
eigenvalue_df = pd.DataFrame({})
eigenvalue_string_df = pd.DataFrame({})
for enum, eigenvalue in enumerate(eigenvalues):
    if not isclose(eigenvalue.imag, 0):
        data = eigenvalue
        stability = 'Center'
    elif isclose(eigenvalue.real, 1, rel_tol=1e-6):
        data = eigenvalue.real
        stability = 'Center'
    elif eigenvalue > 1:
        data = eigenvalue.real
        stability = 'Unstable'
    elif eigenvalue < 1:
        data = eigenvalue.real
        stability = 'Stable'
    eigval_data = pd.DataFrame({f"lambda_{enum+1}":[data, stability]})
    eigval_string_data = pd.DataFrame({f"lambda_{enum+1}":[f'{data:.5g}', stability]})
    eigenvalue_df = pd.concat([eigenvalue_df, eigval_data], axis=1)
    eigenvalue_string_df = pd.concat([eigenvalue_string_df, eigval_string_data], axis=1)

eigenvalue_table = (
    GT(eigenvalue_string_df)
    .tab_header(
        title=md(f"Eigenvalues of Periodic Orbit at xi={xi}, eta=0<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        lambda_1="{{:lambda:_1}}<br>[non-dim]",
        lambda_2="{{:lambda:_2}}<br>[non-dim]",
        lambda_3="{{:lambda:_3}}<br>[non-dim]",
        lambda_4="{{:lambda:_4}}<br>[non-dim]",
        lambda_5="{{:lambda:_5}}<br>[non-dim]",
        lambda_6="{{:lambda:_6}}<br>[non-dim]",
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
eigenvalue_table.show()

# Identify 20 fixed points around the orbit
fixed_points = dict()
# Generating each fixed point with reduced error!
t0 = 0
monodromy_error = pd.DataFrame()
last_monodromy = deepcopy(monodromy)
last_IC = deepcopy(IC)
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
    fixed_points[enum]['fp_vx'] = fixed_point_prop.y[3,-1]
    fixed_points[enum]['fp_vy'] = fixed_point_prop.y[4,-1]
    # Update the ICs for the next round
    last_IC = [
        fixed_point_prop.y[0,-1], fixed_point_prop.y[1,-1], 0, fixed_point_prop.y[3,-1], fixed_point_prop.y[4,-1], 0,
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

# Choose the stable eigenvalue for the step-off distance
d_val = eigenvalues[1].real #[nondim]
d_dim = d_val*l_char #[km]
print(f'step off: {d_dim}[km]')
print(f'step off: {d_val}[non-dim]')

# Generate the ICs for stable/unstable positive/negative half-manifold for each fixed point 
for num, point in fixed_points.items():
    x_fixed = point['fp_x']
    y_fixed = point['fp_y']
    vx_fixed = point['fp_vx']
    vy_fixed = point['fp_vy']
    unstable_eigenvector = point['eigenvectors'][0]
    stable_eigenvector = point['eigenvectors'][1]
    # Normalize against position components
    nu_W_unstable = unstable_eigenvector/np.linalg.norm(unstable_eigenvector[0:3])
    nu_W_stable = stable_eigenvector/np.linalg.norm(stable_eigenvector[0:3])
    # Scale by the step-off distance
    scaled_nu_W_unstable = d_val*nu_W_unstable
    scaled_nu_W_stable = d_val*nu_W_stable
    # Stable positive half-manifold (towards the Moon?)
    x0_stable_pos = x_fixed + scaled_nu_W_stable[0].real[0]
    y0_stable_pos = y_fixed + scaled_nu_W_stable[1].real[0]
    vx0_stable_pos = vx_fixed + scaled_nu_W_stable[3].real[0]
    vy0_stable_pos = vy_fixed + scaled_nu_W_stable[4].real[0]
    fixed_points[num]['IC_stable_pos'] = [x0_stable_pos, y0_stable_pos, 0, vx0_stable_pos, vy0_stable_pos, 0]
    # Stable negative half-manifold (towards the Earth?)
    x0_stable_neg = x_fixed - scaled_nu_W_stable[0].real[0]
    y0_stable_neg = y_fixed - scaled_nu_W_stable[1].real[0]
    vx0_stable_neg = vx_fixed - scaled_nu_W_stable[3].real[0]
    vy0_stable_neg = vy_fixed - scaled_nu_W_stable[4].real[0]
    fixed_points[num]['IC_stable_neg'] = [x0_stable_neg, y0_stable_neg, 0, vx0_stable_neg, vy0_stable_neg, 0]
    # Unstable positive half-manifold (towards the Moon?)
    x0_unstable_pos = x_fixed + scaled_nu_W_unstable[0].real[0]
    y0_unstable_pos = y_fixed + scaled_nu_W_unstable[1].real[0]
    vx0_unstable_pos = vx_fixed + scaled_nu_W_unstable[3].real[0]
    vy0_unstable_pos = vy_fixed + scaled_nu_W_unstable[4].real[0]
    fixed_points[num]['IC_unstable_pos'] = [x0_unstable_pos, y0_unstable_pos, 0, vx0_unstable_pos, vy0_unstable_pos, 0]
    # Unstable negative half-manifold (towards the Earth?)
    x0_unstable_neg = x_fixed - scaled_nu_W_unstable[0].real[0]
    y0_unstable_neg = y_fixed - scaled_nu_W_unstable[1].real[0]
    vx0_unstable_neg = vx_fixed - scaled_nu_W_unstable[3].real[0]
    vy0_unstable_neg = vy_fixed - scaled_nu_W_unstable[4].real[0]
    fixed_points[num]['IC_unstable_neg'] = [x0_unstable_neg, y0_unstable_neg, 0, vx0_unstable_neg, vy0_unstable_neg, 0]

# Now construct the half-manifolds
# Need to propagate for each point and put the t,x,y data in a dataframe
# Label each propagation (stable/unstable, neg/pos) so it plots right (we can get rid of the legend later)
# Separate the dataframes so we can plot and layer separately.
def crossxEventMoon(t, sv, mu):
    return sv[1] if sv[0] > x_Moon else np.nan
crossxEventMoon.terminal = True

def crossxEventEarth(t, sv, mu):
    return sv[1] if sv[0] < x_Earth else np.nan
crossxEventEarth.terminal = True

def crossxEvent(t, sv, mu):
    return sv[1]
crossxEvent.terminal = 3

phi_IC = [
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]

events_list = [crossxEventEarth, crossxEventMoon]
tscale=2
negative_half_manifolds = pd.DataFrame()
positive_half_manifolds = pd.DataFrame()
for num, point in fixed_points.items():
    for IC in ['IC_stable_pos', 'IC_stable_neg', 'IC_unstable_pos', 'IC_unstable_neg']:
    # for IC in ['IC_stable_pos', 'IC_stable_neg']:
        # Adjust forward or backward time based on stable(backward)/unstable(forward)
        if '_stable_' in IC:
            tf = -tscale*period
        elif '_unstable_' in IC:
            tf = tscale*period

        complete_IC = point[IC] + phi_IC
        prop = solve_ivp(spatial_ode, [0, tf], complete_IC, events=events_list, args=(mu,), rtol=1e-12,atol=1e-14)
        # Construct the manifold trajectories that flow toward the Earth (neg half-manifolds, I think)
        if 'neg' in IC:
            negative_half_manifolds_data = pd.DataFrame({
                'name': f'{num}: {IC}',
                'num': num,
                'halfmanifold': IC,
                't':prop.t,
                'x':prop.y[0],
                'y':prop.y[1],
                'vx':prop.y[3],
                'vy':prop.y[4],
            })
            negative_half_manifolds = pd.concat([negative_half_manifolds, negative_half_manifolds_data], ignore_index=True)
        # Construct the manifold trajectories that flow toward the Moon (pos half-manifolds, I think)
        elif 'pos' in IC:
            positive_half_manifolds_data = pd.DataFrame({
                'name': f'{num}: {IC}',
                'num': num,
                'halfmanifold': IC,
                't':prop.t,
                'x':prop.y[0],
                'y':prop.y[1],
                'vx':prop.y[3],
                'vy':prop.y[4]
            })
            positive_half_manifolds = pd.concat([positive_half_manifolds, positive_half_manifolds_data], ignore_index=True)


# Build Plots
orbit_plot = alt.Chart(orbit).mark_line(clip=True,strokeWidth=2).encode(
    x=alt.X('x:Q'),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q'),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(domain=['Periodic Orbit'], range=['red'])).title(None),
    order='t'
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

# Build plot

x_min = -0.7
x_max = 1.2
# x_min = 0.8
# x_max = 0.88
y_lim = (x_max-x_min)/2

pos_manifolds_chart = alt.Chart(positive_half_manifolds).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None).legend(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Stable and Unstable Manifold Trajectories",f"for the orbit {xi_symbol}={xi}, {eta_symbol}=0 ({ps}, Lillian Shido)"]
)

neg_manifolds_chart = alt.Chart(negative_half_manifolds).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q'),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q'),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None).legend(None),
    order='t'
)
all_manifolds_chart_layer = alt.layer(pos_manifolds_chart, neg_manifolds_chart, L1_loc, earth_loc, moon_loc, orbit_plot).resolve_scale(color='independent')
all_manifolds_chart_layer.save(f'test_events_all_manifolds_{ps}.png', ppi=200)

# pos_manifolds_chart_layer = alt.layer(pos_manifolds_chart, L1_loc, earth_loc, moon_loc, orbit_plot, ).resolve_scale(color='independent')
# pos_manifolds_chart_layer.save(f'pos_manifolds_{ps}.png', ppi=200)



# neg_manifolds_chart = alt.Chart(negative_half_manifolds).mark_line(clip=True,strokeWidth=1).encode(
#     x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
#     # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
#     y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
#     # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
#     color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None).legend(None),
#     order='t'
# ).properties(
#     width=400,
#     height=400,
#     title=["Negative Manifold Trajectories",f"for the orbit {xi_symbol}={xi}, {eta_symbol}=0 ({ps}, Lillian Shido)"]
# )
# neg_manifolds_chart_layer = alt.layer(neg_manifolds_chart, L1_loc, earth_loc, moon_loc, orbit_plot, ).resolve_scale(color='independent')
# neg_manifolds_chart_layer.save(f'neg_manifolds_{ps}.png', ppi=200)

# all_manifolds_chart = alt.Chart(massive_df).mark_line(clip=True,strokeWidth=1).encode(
#     x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
#     # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
#     y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
#     # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
#     color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None).legend(None),
#     order='t'
# ).properties(
#     width=400,
#     height=400,
#     title=["All Manifold Trajectories",f"for the orbit {xi_symbol}=0.01, {eta_symbol}=0 ({ps}, Lillian Shido)"]
# )
# all_manifolds_chart_layer = alt.layer(all_manifolds_chart, L1_loc, earth_loc, moon_loc, orbit_plot, moon_closest_fp_loc, earth_closest_fp_loc).resolve_scale(color='independent')
# all_manifolds_chart_layer.save(f'all_manifolds_{ps}.png', ppi=200)

# unstable_manifolds_chart = alt.Chart(massive_df[massive_df['halfmanifold'].str.contains('_unstable_')]).mark_line(clip=True,strokeWidth=1).encode(
#     x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
#     # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
#     y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
#     # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
#     color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None),
#     order='t'
# ).properties(
#     width=400,
#     height=400,
#     title=["Unstable Manifold Trajectories",f"for the orbit {xi_symbol}=0.01, {eta_symbol}=0 ({ps}, Lillian Shido)"]
# )
# unstable_manifolds_chart_layer = alt.layer(unstable_manifolds_chart, L1_loc, earth_loc, moon_loc, orbit_plot).resolve_scale(color='independent')
# unstable_manifolds_chart_layer.save(f'unstable_manifolds_{ps}.png', ppi=200)

pdb.set_trace()