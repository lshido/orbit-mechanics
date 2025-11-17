ps = "F1 part d"
# Author: Lillian Shido
# Date: 11/16/2025

import pdb
import numpy as np
from math import pi, sqrt

import pandas as pd
from scipy.integrate import solve_ivp
import altair as alt
from great_tables import GT, md, system_fonts

from symbols import tau_symbol, pi_symbol
from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, build_A_matrix_collinear, spatial_ode, calc_spatial_Jacobi, calc_ZVC_Jacobi,calc_velocity_from_Jacobi, calc_spatial_Jacobi_array

import warnings
warnings.filterwarnings("ignore",category=PendingDeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)

A = build_A_matrix_collinear(mu, x_L1, y_L1, 0)
eigenvalues, eigenvectors = np.linalg.eig(A)

# Generate eigenspaces
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[y_L1]})
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0]})
eigenspace = pd.DataFrame({})
manifold = pd.DataFrame({})
for i in range(0,6):
    if i==0 or i==1:
        print(f"vector {i+1} i-comp: {eigenvectors[:,i][0].real[0,0]}")
        print(f"vector {i+1} y-comp: {eigenvectors[:,i][1].real[0,0]}")
        if i==0:
            name="Stable Eigenspace"
            label_manifold="Stable Half-Manifold"
        elif i==1:
            name="Unstable Eigenspace"
            label_manifold="Unstable Half-Manifold"
        scale = 1 # Extend eigenvector line out
        x_eig_max = x_L1 + scale*eigenvectors[:,i][0].real[0,0]
        x_eig_min = x_L1 + scale*-eigenvectors[:,i][0].real[0,0]
        y_eig_max = y_L1 + scale*eigenvectors[:,i][1].real[0,0]
        y_eig_min = y_L1 + scale*-eigenvectors[:,i][1].real[0,0]
        eigenspace_data = pd.DataFrame({
            'name': name,
            'x':[x_eig_min,x_L1,x_eig_max],
            'y':[y_eig_min,y_L1,y_eig_max]
        })
        eigenspace = pd.concat([eigenspace, eigenspace_data], ignore_index=True)

#  Generate the initial velocities
jc = calc_ZVC_Jacobi(mu, x_L1, y_L1)
IC_list = [(0.8370,0.0001),(0.8371,0.0002),(0.8373,0.0004)]
initial_velocities = pd.DataFrame()
for n,(x0,y0) in enumerate(IC_list):
    x0_km = (x0-x_L1)*l_char
    y0_km = y0*l_char
    vx0 = calc_velocity_from_Jacobi(jc, mu, x0, y0, 0)
    initial_velocities_data = pd.DataFrame({
        "Case":[f"Case {n+1}"],
        "x0_km":[x0_km],
        "y0_km":[y0_km],
        "x0":[x0],
        "y0":[y0],
        "vx0":[vx0],
        "vy0":[0],
        "jacobi":[calc_spatial_Jacobi(mu,x0,y0,0,vx0,0,0)]
    })
    initial_velocities = pd.concat([initial_velocities, initial_velocities_data], ignore_index=True)

# Configure the positive and negative time for propagation of the nonlinear equations
tf = 0.5
pos_tspan = [0,tf*pi]
neg_tspan = [0,-tf*pi]
# Specify time constant
tau = 1/eigenvalues[1].real

nonlinear_result = pd.DataFrame()
tau_result = pd.DataFrame()
for case, row in enumerate(initial_velocities.iterrows()):
        x0 = row[1]['x0']
        y0 = row[1]['y0']
        vx0 = row[1]['vx0']
        IC = [
            x0,y0,0,vx0,0,0,# IC states
            1,0,0,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,0,0,
            0,0,1,0,0,0,
            0,0,0,1,0,0,
            0,0,0,0,1,0,
            0,0,0,0,0,1
        ]
        pos_t_prop = solve_ivp(spatial_ode, pos_tspan, IC, args=(mu,), rtol=1e-12,atol=1e-14)
        nonlinear_result_data = pd.DataFrame({
            'name': f"Case {case+1}",
            't':pos_t_prop.t,
            'x':pos_t_prop.y[0],
            'y':pos_t_prop.y[1],
            'vx':pos_t_prop.y[3],
            'vy':pos_t_prop.y[4]
        })
        nonlinear_result = pd.concat([nonlinear_result, nonlinear_result_data], ignore_index=True)
        neg_t_prop = solve_ivp(spatial_ode, neg_tspan, IC, args=(mu,), rtol=1e-12,atol=1e-14)
        nonlinear_result_data = pd.DataFrame({
            'name': f"Case {case+1}",
            't':neg_t_prop.t,
            'x':neg_t_prop.y[0],
            'y':neg_t_prop.y[1],
            'vx':neg_t_prop.y[3],
            'vy':neg_t_prop.y[4]
        })
        nonlinear_result = pd.concat([nonlinear_result, nonlinear_result_data], ignore_index=True)
        pos_tau_prop = solve_ivp(spatial_ode, [0,tau], IC, args=(mu,), rtol=1e-12,atol=1e-14)
        tau_result_data = pd.DataFrame({
            'name': f"Case {case+1}",
            't':pos_tau_prop.t,
            'x':pos_tau_prop.y[0],
            'y':pos_tau_prop.y[1]
        })
        tau_result = pd.concat([tau_result, tau_result_data], ignore_index=True)

# Calc JC for each step
nonlinear_result['jacobi'] = nonlinear_result.apply(lambda row: calc_spatial_Jacobi(mu, row['x'],row['y'],0,row['vx'],row['vy'],0),axis=1)
nonlinear_result['error_jc'] = nonlinear_result.apply(lambda row: abs(row['jacobi']-jc),axis=1)

# Get 0, min, max for tau results
max_tau = tau_result.loc[tau_result['t']==tau_result['t'].max()]
min_tau = tau_result.loc[tau_result['t']==tau_result['t'].min()]
zero_tau = tau_result.loc[tau_result['t']==0].drop_duplicates()

#  Calc ZVCs
df_zvc = pd.DataFrame()
tolerance = 1e-12
max_iterations = 100
C = calc_ZVC_Jacobi(mu, x_L1, y_L1)
# Give an initial guess for y that is around the Earth
for x in np.arange(0.8365,0.8378,1e-5): # Find the curve for x
    for y in np.arange(0,5e-4,1e-5): # These are my guesses
        counter = 0
        while True:
            counter = counter + 1
            d = sqrt((x+mu)**2 + y**2)
            r = sqrt((x-1+mu)**2 + y**2)
            f = x**2 + y**2 + (2*(1-mu)/d) + (2*mu/r) - C
            f_prime_y = 2*y*( 1 - (1-mu)/d**3 - mu/r**3)
            delta = f/f_prime_y
            if abs(f) > tolerance:
                if counter > max_iterations: # check for max iterations
                    new_C = calc_ZVC_Jacobi(mu, x, y)
                    C_error = abs(new_C - C)
                    if C_error < 1e-12:
                        zvc_data = pd.DataFrame({
                            "name":'ZVC',
                            "x":[x],
                            "y":[y],
                            "new_C":[new_C],
                            "error":[error_C],
                        })
                        df_zvc = pd.concat([df_zvc,zvc_data],ignore_index=True)
                    else:
                        break
                elif d == 0 or r == 0: # check for dividing by 0
                    break
                elif abs(f_prime_y) <= tolerance: # check for zero derivative
                    break
                else:
                    y = y - delta
                    continue
            else:
                new_C = calc_ZVC_Jacobi(mu, x, y)
                error_C = abs(new_C - C)
                zvc_data = pd.DataFrame({
                    "name":'ZVC',
                    "x":[x],
                    "y":[y],
                    "new_C":[new_C],
                    "error":[error_C],
                })
                df_zvc = pd.concat([df_zvc,zvc_data],ignore_index=True)
                break
for y in np.arange(1e-6,1.2e-3,1e-5): # Find the curve for y
    for x in np.arange(0.83,0.838,1e-5): # These are my guesses for x
        counter = 0
        while True:
            counter = counter + 1
            d = sqrt((x+mu)**2 + y**2)
            r = sqrt((x-1+mu)**2 + y**2)
            f = x**2 + y**2 + (2*(1-mu)/d) + (2*mu/r) - C
            f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d**3 - 2*mu*(x-1+mu)/r**3
            delta = f/f_prime_x
            if abs(f) > tolerance:
                if counter > max_iterations: # check for max iterations
                    new_C = calc_ZVC_Jacobi(mu, x, y)
                    C_error = abs(new_C - C)
                    if C_error < 1e-12:
                        zvc_data = pd.DataFrame({
                            "name":'ZVC',
                            "x":[x],
                            "y":[y],
                            "new_C":[new_C],
                            "error":[error_C],
                        })
                        df_zvc = pd.concat([df_zvc,zvc_data],ignore_index=True)
                    else:
                        break
                elif d == 0 or r == 0: #check for dividing by 0
                    break
                elif abs(f_prime_x) <= tolerance: #check for zero derivative
                    break
                else:
                    x = x - delta
                    continue
            else:
                new_C = calc_ZVC_Jacobi(mu, x, y)
                error_C = abs(new_C - C)
                zvc_data = pd.DataFrame({
                    "name":'ZVC',
                    "x":[x],
                    "y":[y],
                    "new_C":[new_C],
                    "error":[error_C],
                })
                df_zvc = pd.concat([df_zvc,zvc_data],ignore_index=True)
                break

neg_y_zvc_data = pd.DataFrame({
    'name':'ZVC',
    'x':df_zvc['x'],
    'y':-df_zvc['y'],
    'new_C':df_zvc['new_C'],
    'error':df_zvc['error']
})
df_zvc = pd.concat([df_zvc,neg_y_zvc_data],ignore_index=True)

# Build plot
x_min = 0.8365
x_max = 0.8385
y_lim = (x_max-x_min)/2
base = alt.Chart(nonlinear_result).mark_line(strokeWidth=1,clip=True).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='plasma')).title(f"t = {tf}{pi_symbol}"),
    order='t',
).properties(
    width=400,
    height=400,
    title=f"Nonlinear Results near L1 ({ps}, Lillian Shido)"
)

zvc = alt.Chart(df_zvc).mark_point(size=5,filled=True,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['ZVC'], range=['forestgreen'])).title(None)
)

eigspace = alt.Chart(eigenspace).mark_line(strokeWidth=1.5,strokeDash=(4,4),clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(scheme='greenblue')).title(None)
)

time_constant = alt.Chart(tau_result).mark_line(strokeWidth=1.5,clip=True).encode(
    x=alt.X('x:Q'),
    y=alt.Y('y:Q'),
    color=alt.Color('name:N', scale=alt.Scale(scheme='darkmulti')).title(f"t = {tau_symbol}"),
    order='t',
)

t_tau = alt.Chart(max_tau).mark_point(clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(scheme='darkmulti')).title(f"Location @ t = {tau_symbol}")
)

t_zero = alt.Chart(zero_tau).mark_point(shape='diamond', clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(scheme='darkmulti')).title(f"Location @ t = 0")
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
    color=alt.Color('name:N', scale=alt.Scale(domain=['L1'], range=['darkblue'])).title(None)
)

final = alt.layer(base, time_constant, t_tau, t_zero, eigspace, L1_loc, zvc).resolve_scale(color='independent')
final.save(f'nonlinear_L1 {ps}.png', ppi=200)

error = alt.Chart(nonlinear_result).mark_point(size=10,clip=True,filled=True).encode(
     x=alt.X('t:Q', axis=alt.Axis(title='Time [non-dim]')),
     y=alt.Y('error_jc:Q', axis=alt.Axis(format='.1e',title='JC Error')),
     color=alt.Color('name:N').title(None)
).properties(
    width=400,
    height=400,
    title=f"Jacobi Constant Error ({ps}, Lillian Shido)"
)
error.save(f'error_jc {ps}.png',ppi=200)

# initial_velocities_table = (
#     GT(initial_velocities)
#     .tab_header(
#         title=md(f"Initial Conditions for each pair of (x0,y0)<br>({ps}, Lillian Shido)")
#     )
#     .cols_label(
#         x0_km="{{x_0}} [km]",
#         y0_km="{{y_0}} [km]",
#         x0="x [non-dim]",
#         y0="y [non-dim]",
#         vx0="{{v_x0}} [non-dim]",
#         vy0="{{v_y0}} [non-dim]",
#         jacobi="JC"
#     )
#     .opt_horizontal_padding(scale=2)
#     .fmt_number(
#         columns=["x0_km", "y0_km", "vx0","jacobi"],
#         n_sigfig=6
#     )
#     .fmt_number(
#         columns=["x0"],
#         n_sigfig=4
#     )
#     .fmt_number(
#         columns=["y0"],
#         n_sigfig=1
#     )
#     .cols_align(align="center")
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
#     .opt_table_outline()
# )
# initial_velocities_table.show()