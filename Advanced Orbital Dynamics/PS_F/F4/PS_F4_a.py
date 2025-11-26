ps = "F4 part a"
# Author: Lillian Shido
# Date: 11/24/2025

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
    'name':'Periodic Orbit',
    't':full_period_prop.t,
    'x':full_period_prop.y[0],
    'y':full_period_prop.y[1]
})

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

# eigenvalue_table = (
#     GT(eigenvalue_string_df)
#     .tab_header(
#         title=md(f"Eigenvalues of Periodic Orbit at xi=0.01, eta=0<br>({ps}, Lillian Shido)")
#     )
#     .cols_label(
#         lambda_1="{{:lambda:_1}}<br>[non-dim]",
#         lambda_2="{{:lambda:_2}}<br>[non-dim]",
#         lambda_3="{{:lambda:_3}}<br>[non-dim]",
#         lambda_4="{{:lambda:_4}}<br>[non-dim]",
#         lambda_5="{{:lambda:_5}}<br>[non-dim]",
#         lambda_6="{{:lambda:_6}}<br>[non-dim]",
#     )
#     .cols_align(
#         align="center"
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# eigenvalue_table.show()

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

# monodromy_error_table = (
#     GT(monodromy_error)
#     .tab_header(
#         title=md(f"Reduced Error of Monodromy Matrix<br>via optimization of transformations<br>(F2 part c, Lillian Shido)")
#     )
#     .fmt_scientific(
#         columns = ['Error'],
#         n_sigfig = 5
#     )
#     .cols_align(
#         align="center"
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# monodromy_error_table.show()

# For each fixed point, write down which eigenvector is the stable one and unstable one
# verify_stability = pd.DataFrame()
# for num, point in fixed_points.items():
#     lambda_1 = point['eigvals'][0].real
#     lambda_2 = point['eigvals'][1].real
#     if lambda_1 > 1 and lambda_2 < 1:
#         stability_1 = 'Unstable'
#         stability_2 = 'Stable'
#     else:
#         stability_1 = 'NA'
#         stability_2 = 'NA'
#     verify_stability_data = pd.DataFrame({
#         'Fixed Point': [num+1],
#         'lambda_1': [lambda_1],
#         'lambda_2': [lambda_2],
#         'stability_1': [stability_1],
#         'stability_2': [stability_2]
#     })
#     verify_stability = pd.concat([verify_stability, verify_stability_data], ignore_index=True)

# verify_stability_table = (
#     GT(verify_stability)
#     .tab_header(
#         title=md(f"Stability Verification of Eigenvalues 1 and 2<br>{ps}, Lillian Shido)")
#     )
#     .cols_label(
#         lambda_1 = "{{:lambda:_1}}",
#         lambda_2 = "{{:lambda:_2}}",
#         stability_1 = "{{:lambda:_1}} Stability",
#         stability_2 = "{{:lambda:_2}} Stability"
#     )
#     .fmt_number(
#         columns = ['lambda_1', 'lambda_2'],
#         n_sigfig = 5
#     )
#     .cols_align(
#         align="center"
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# verify_stability_table.show()

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
    if sv[0] > x_Moon:
        return sv[1]
    else:
        return 1.0 # Return value of y
crossxEventMoon.terminal = True

def crossxEventEarth(t, sv, mu):
    if sv[0] < x_Earth:
        crossxEventEarth.terminal = True
    else:
        crossxEventEarth.terminal = False
    return sv[1]

def crossxEvent(t, sv, mu):
    return sv[1]

phi_IC = [
    1,0,0,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1
]
tscale=1.4
df_xcross = pd.DataFrame()
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

        # if num == 10:
        #     crossxEvent.terminal = 1
        # if 11 <= num <= 19:
        #     crossxEvent.terminal = 3
        # else:
        #     crossxEvent.terminal = 3

        complete_IC = point[IC] + phi_IC
        prop = solve_ivp(spatial_ode, [0, tf], complete_IC, events=crossxEvent, args=(mu,), rtol=1e-12,atol=1e-14)
        # Construct the manifold trajectories that flow toward the Earth (neg half-manifolds, I think)
        if 'neg' in IC:
            # Build x-axis cross df
            for n, row in enumerate(prop.y_events[0]):
                if row[0] < x_Earth:
                    df_xcross_data = pd.DataFrame({
                        'name': [f'{num}: {IC}'],
                        'num': [num],
                        'type': [IC],
                        't': [prop.t_events[0][n]],
                        'x': [row[0]],
                        'y': [row[1]]
                    })
            df_xcross = pd.concat([df_xcross, df_xcross_data], ignore_index=True)
            negative_half_manifolds_data = pd.DataFrame({
                'name': f'{num}: {IC}',
                'num': num,
                'type': IC,
                't':prop.t,
                'x':prop.y[0],
                'y':prop.y[1],
                'vx':prop.y[3],
                'vy':prop.y[4],
            })
            negative_half_manifolds = pd.concat([negative_half_manifolds, negative_half_manifolds_data], ignore_index=True)
        # Construct the manifold trajectories that flow toward the Moon (pos half-manifolds, I think)
        elif 'pos' in IC:
            # Build x-axis cross df
            for n, row in enumerate(prop.y_events[0]):
                if row[0] > x_Moon:
                    df_xcross_data = pd.DataFrame({
                        'name': [f'{num}: {IC}'],
                        'num': [num],
                        'type': [IC],
                        't': [prop.t_events[0][n]],
                        'x': [row[0]],
                        'y': [row[1]]
                    })
            df_xcross = pd.concat([df_xcross, df_xcross_data], ignore_index=True)
            positive_half_manifolds_data = pd.DataFrame({
                'name': f'{num}: {IC}',
                'num': num,
                'type': IC,
                't':prop.t,
                'x':prop.y[0],
                'y':prop.y[1],
                'vx':prop.y[3],
                'vy':prop.y[4]
            })
            positive_half_manifolds = pd.concat([positive_half_manifolds, positive_half_manifolds_data], ignore_index=True)

# Propagate to plot
# For every n, row in df_xcross
# propagate
# Need ICs (from fixed_points), need tf (from df_xcross)
# Put into the same dataframe
massive_df = pd.DataFrame()
for row in df_xcross.iterrows():
    num = row[1]['num']
    type = row[1]['type']
    tf = row[1]['t']
    point_IC = fixed_points[num][type]
    complete_IC = point_IC + phi_IC
    xaxis_prop = solve_ivp(spatial_ode, [0, tf], complete_IC, args=(mu,), rtol=1e-12,atol=1e-14)
    # Want this data
    massive_df_data = pd.DataFrame({
        'num': num,
        'type': type,
        't':xaxis_prop.t,
        'x':xaxis_prop.y[0],
        'y':xaxis_prop.y[1],
        'vx':xaxis_prop.y[3],
        'vy':xaxis_prop.y[4]
    })
    massive_df = pd.concat([massive_df, massive_df_data], ignore_index=True)

# Calculate minimum and maximum closest pass to Moon
moon_farthest = df_xcross.iloc[df_xcross[df_xcross['name'].str.contains('pos')]['x'].idxmax()]
moon_closest = df_xcross.iloc[df_xcross[df_xcross['name'].str.contains('pos')]['x'].idxmin()]
time_moon_closest = abs(moon_closest['t'])*t_char/3600/24
# For table
moon_pass = pd.DataFrame({
    'fixed_point_closest_pass': [moon_closest['name']],
    'closest_point_from_moon': [(moon_closest['x'] - x_Moon)*l_char], # km
    'farthest_point_from_moon': [(moon_farthest['x'] - x_Moon)*l_char], # km
    'interval_to_closest_pass': [time_moon_closest] # days  
})
# For plot
moon_closest_fp_x = fixed_points[moon_closest['num']][moon_closest['name'].split()[1]][0]
moon_closest_fp_y = fixed_points[moon_closest['num']][moon_closest['name'].split()[1]][1]
moon_closest_fp = pd.DataFrame({
    'name':["Delivers Closest Moon Pass"],
    'fixed point': [moon_closest['num']],
    'subspace': [moon_closest['name'].split('_')[1]],
    'half-manifold': [moon_closest['name'].split('_')[2]],
    'x':[moon_closest_fp_x],
    'y':[moon_closest_fp_y],
})

moon_pass_table = (
    GT(moon_pass)
    .tab_header(
        title=md(f"Manifold Trajectories Moon Pass Overview<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        closest_point_from_moon = '{{Closest Point from Moon}}<br>[km]',
        farthest_point_from_moon = '{{Farthest Point from Moon}}<br>[km]',
        interval_to_closest_pass = '{{Interval to Closest Pass}}<br>[days]',
        fixed_point_closest_pass = '{{Fixed Point of Closest Pass}}'
    )
    .fmt_scientific(
        columns = ['closest_point_from_moon','farthest_point_from_moon'],
        n_sigfig = 5
    )
    .fmt_number(
        columns = ['interval_to_closest_pass'],
        n_sigfig = 5
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
moon_pass_table.show()

moon_fp_table = (
    GT(moon_closest_fp)
    .tab_header(
        title=md(f"Fixed Point that Delivers Closest Moon Pass<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        x = '{{x}}<br>[non-dim]',
        y = '{{y}}<br>[non-dim]',
    )
    .fmt_number(
        columns = ['x','y'],
        n_sigfig = 6
    )
    .cols_align(
        align="center"
    )
    .cols_hide(
        columns = ['name']
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
moon_fp_table.show()

# Calculate minimum and maximum closest pass to Earth
earth_closest = df_xcross.iloc[df_xcross[df_xcross['name'].str.contains('neg')]['x'].idxmax()]
earth_farthest = df_xcross.iloc[df_xcross[df_xcross['name'].str.contains('neg')]['x'].idxmin()]
time_earth_closest = abs(earth_closest['t'])*t_char/3600/24
# For table
earth_pass = pd.DataFrame({
    'fixed_point_closest_pass': [earth_closest['name']],
    'closest_point_from_earth': [(x_Earth - earth_closest['x'])*l_char], # km
    'farthest_point_from_earth': [(x_Earth - earth_farthest['x'])*l_char], # km
    'interval_to_closest_pass': [time_earth_closest] # days  
})

# Calculate minimum and maximum closest pass to Earth
earth_closest_fp_x = fixed_points[earth_closest['num']][earth_closest['name'].split()[1]][0]
earth_closest_fp_y = fixed_points[earth_closest['num']][earth_closest['name'].split()[1]][1]
earth_closest_fp = pd.DataFrame({
    'name':["Delivers Closest Earth Pass"],
    'fixed point': [earth_closest['num']],
    'subspace': [earth_closest['name'].split('_')[1]],
    'half-manifold': [earth_closest['name'].split('_')[2]],
    'x':[earth_closest_fp_x],
    'y':[earth_closest_fp_y],
})

earth_pass_table = (
    GT(earth_pass)
    .tab_header(
        title=md(f"Manifold Trajectories Earth Pass Overview<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        closest_point_from_earth = '{{Closest Point from Earth}}<br>[km]',
        farthest_point_from_earth = '{{Farthest Point from Earth}}<br>[km]',
        interval_to_closest_pass = '{{Interval to Closest Pass}}<br>[days]',
        fixed_point_closest_pass = '{{Fixed Point of Closest Pass}}'
    )
    .fmt_scientific(
        columns = ['closest_point_from_earth','farthest_point_from_earth'],
        n_sigfig = 5
    )
    .fmt_number(
        columns = ['interval_to_closest_pass'],
        n_sigfig = 5
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
earth_pass_table.show()

earth_fp_table = (
    GT(earth_closest_fp)
    .tab_header(
        title=md(f"Fixed Point that Delivers Closest Earth Pass<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        x = '{{x}}<br>[non-dim]',
        y = '{{y}}<br>[non-dim]',
    )
    .fmt_number(
        columns = ['x','y'],
        n_sigfig = 6
    )
    .cols_align(
        align="center"
    )
    .cols_hide(
        columns = ['name']
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
earth_fp_table.show()

pdb.set_trace()

eigenspace = pd.DataFrame({})
velocity_eigenspace = pd.DataFrame({})
manifold = pd.DataFrame({})
for num, point in fixed_points.items():
    for i in range(0,6):
        if i==0 or i==1:
        # if i==0:
            if i==0:
                name="Unstable Eigendirection"
            elif i==1:
                name="Stable Eigendirection"
            scale = 0.1 # Extend eigenvector line out
            vscale = 0.01
            x_eig_max = point['fp_x'] + scale*point['eigenvectors'][i][0,0].real
            x_eig_min = point['fp_x'] + scale*-point['eigenvectors'][i][0,0].real
            y_eig_max = point['fp_y'] + scale*point['eigenvectors'][i][1,0].real
            y_eig_min = point['fp_y'] + scale*-point['eigenvectors'][i][1,0].real
            vx_eig_max = point['fp_x'] + vscale*point['eigenvectors'][i][3,0].real
            vx_eig_min = point['fp_x'] + vscale*-point['eigenvectors'][i][3,0].real
            vy_eig_max = point['fp_y'] + vscale*point['eigenvectors'][i][4,0].real
            vy_eig_min = point['fp_y'] + vscale*-point['eigenvectors'][i][4,0].real
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
            velocity_eigenspace_pos_data = pd.DataFrame({
                'name': f'{num}: {name} (+)',
                'vx':[point['fp_x']],
                'vx2':[vx_eig_max],
                'vy':[point['fp_y']],
                'vy2':[vy_eig_max],
                'angle':[np.rad2deg(atan2((vy_eig_max - point['fp_y']),(vx_eig_max - point['fp_x'])))]
            })
            velocity_eigenspace = pd.concat([velocity_eigenspace, velocity_eigenspace_pos_data], ignore_index=True)
            velocity_eigenspace_neg_data = pd.DataFrame({
                'name': f'{num}: {name} (-)',
                'vx':[point['fp_x']],
                'vx2':[vx_eig_min],
                'vy':[point['fp_y']],
                'vy2':[vy_eig_min],
                'angle':[np.rad2deg(atan2((vy_eig_min - point['fp_y']),(vx_eig_min - point['fp_x'])))]
            })
            velocity_eigenspace = pd.concat([velocity_eigenspace, velocity_eigenspace_neg_data], ignore_index=True)

# Build Plots
orbit_plot = alt.Chart(orbit).mark_line(clip=True,strokeWidth=2).encode(
    x=alt.X('x:Q'),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q'),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(domain=['Periodic Orbit'], range=['red'])).title(None),
    order='t'
)

# fixed_point_loc = alt.Chart(fixed_point).mark_point(filled=True,size=30,clip=True).encode(
#     x='x:Q',
#     y='y:Q',
#     color=alt.Color('name:N').title(None)
# )

# eigendirections = alt.Chart(eigenspace).mark_rule(strokeWidth=0.5,clip=True).encode(
#     x='x:Q',
#     y='y:Q',
#     x2='x2:Q',
#     y2='y2:Q',
#     color=alt.Color('label:N', scale=alt.Scale(domain=['Unstable Eigenspace'], range=['darkolivegreen'])).title(None)
# )

# eigendirection_arrows = alt.Chart(eigenspace).mark_point(shape="wedge",filled=True, fillOpacity=1,size=100,clip=True).encode(
#         x='x2:Q',
#         # x = alt.datum(0.874411),
#         y='y2:Q',
#         # y = alt.datum(-0.015899),
#         angle=alt.Angle('angle_wedge').scale(domain=[0, 360]),
#         # angle=alt.AngleValue(),
#         color=alt.Color('label:N', scale=alt.Scale(domain=['Unstable Eigenspace'], range=['darkolivegreen']), legend=None).title(None)
# )

# step_off_loc = alt.Chart(x_step_off).mark_point(filled=True,size=30,clip=True).encode(
#     x='x:Q',
#     y='y:Q',
#     color=alt.Color('name:N', scale=alt.Scale(domain=['Step-Off Point'], range=['red'])).title(None)
# )

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
# zoom in
# x_min = 0.800
# x_max = 0.870
# Zoom out
x_min = 0.82
x_max = 1.03
y_lim = (x_max-x_min)/2
stable_positive_manifold_chart = alt.Chart(positive_half_manifolds[positive_half_manifolds['name'].str.contains('_stable_')]).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='viridis')).legend(None).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Manifold trajectories that flow toward the moon",f"for the orbit {xi_symbol}=0.01, {eta_symbol}=0 ({ps}, Lillian Shido)"]
)

positive_manifold_chart_layer = alt.layer(stable_positive_manifold_chart, L1_loc, moon_loc, orbit_plot).resolve_scale(color='independent')
# positive_manifold_chart_layer.save(f'flow_toward_moon_{ps}.png', ppi=200)

x_min = -0.5
x_max = 1.2
# x_min = 0.8
# x_max = 1.05
y_lim = (x_max-x_min)/2
stable_negative_manifold_chart = alt.Chart(negative_half_manifolds[negative_half_manifolds['name'].str.contains('_stable_')]).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Manifold trajectories that flow toward the Earth",f"for the orbit {xi_symbol}=0.01, {eta_symbol}=0 ({ps}, Lillian Shido)"]
)

negative_manifold_chart = alt.Chart(negative_half_manifolds).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None).legend(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Manifold trajectories that flow toward the Earth",f"for the orbit {xi_symbol}=0.01, {eta_symbol}=0 ({ps}, Lillian Shido)"]
)

negative_manifold_chart_layer = alt.layer(negative_manifold_chart, L1_loc, earth_loc, orbit_plot).resolve_scale(color='independent')
# negative_manifold_chart_layer.save(f'flow_toward_earth_{ps}.png', ppi=200)

stable_negative_manifold_chart = alt.Chart(negative_half_manifolds[negative_half_manifolds['name'].str.contains("_stable_")]).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Stable Manifold Trajectories",f"for the orbit {xi_symbol}=0.01, {eta_symbol}=0 ({ps}, Lillian Shido)"]
)

stable_positive_manifold_chart = alt.Chart(positive_half_manifolds[positive_half_manifolds['name'].str.contains("_stable_")]).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q'),
    y=alt.Y('y:Q'),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None).legend(None),
    order='t'
)

stable_manifold_chart_layer = alt.layer(stable_negative_manifold_chart, stable_positive_manifold_chart, L1_loc, earth_loc, moon_loc, orbit_plot).resolve_scale(color='independent')
stable_manifold_chart_layer.save(f'stable_flow_{ps}.png', ppi=200)


unstable_negative_manifold_chart = alt.Chart(negative_half_manifolds[negative_half_manifolds['name'].str.contains("_unstable_")]).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Unstable Manifold Trajectories",f"for the orbit {xi_symbol}=0.01, {eta_symbol}=0 ({ps}, Lillian Shido)"]
)

unstable_positive_manifold_chart = alt.Chart(positive_half_manifolds[positive_half_manifolds['name'].str.contains("_unstable_")]).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q'),
    y=alt.Y('y:Q'),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None).legend(None),
    order='t'
)

unstable_manifold_chart_layer = alt.layer(unstable_negative_manifold_chart, unstable_positive_manifold_chart, L1_loc, earth_loc, moon_loc, orbit_plot).resolve_scale(color='independent')
unstable_manifold_chart_layer.save(f'unstable_flow_{ps}.png', ppi=200)

pdb.set_trace()