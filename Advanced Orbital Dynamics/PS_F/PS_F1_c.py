ps = "F1 part c"
# Author: Lillian Shido
# Date: 11/15/2025

import pdb
import numpy as np
from math import pi, sqrt

import pandas as pd
from scipy.integrate import solve_ivp
import altair as alt
from great_tables import GT, md, system_fonts

from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, build_A_matrix_collinear, spatial_ode, calc_spatial_Jacobi, calc_ZVC_Jacobi,calc_velocity_from_Jacobi

import warnings
warnings.filterwarnings("ignore",category=PendingDeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)

# Configure the positive and negative time 
tf = 0.5*pi
pos_tspan = [0,tf]
neg_tspan = [0,-tf]

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
tolerance = 1e-12
max_iterations = 100
jc = calc_ZVC_Jacobi(mu, x_L1, y_L1)
IC_list = [30,100,250]
initial_velocities = pd.DataFrame()
for x0_km in IC_list:
    x0 = x0_km/l_char + x_L1
    vy0 = calc_velocity_from_Jacobi(jc, mu, x0, 0, 0)
    initial_velocities_data = pd.DataFrame({
        "x0_km":[x0_km],
        "x0":[x0],
        "vy0":[vy0]
    })
    initial_velocities = pd.concat([initial_velocities, initial_velocities_data], ignore_index=True)

nonlinear_result = pd.DataFrame()
for row in initial_velocities.iterrows():
        x0_km = row[1]['x0_km']
        x0 = row[1]['x0']
        vy0 = row[1]['vy0']
        IC = [
            x0,0,0,0,vy0,0,# IC states
            1,0,0,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,0,0,
            0,0,1,0,0,0,
            0,0,0,1,0,0,
            0,0,0,0,1,0,
            0,0,0,0,0,1
        ]
        pos_t_prop = solve_ivp(spatial_ode, pos_tspan, IC, args=(mu,), rtol=1e-12,atol=1e-14)
        nonlinear_result_data = pd.DataFrame({
            'name': f"x0 = {x0_km} km",
            't':pos_t_prop.t,
            'x':pos_t_prop.y[0],
            'y':pos_t_prop.y[1]
        })
        nonlinear_result = pd.concat([nonlinear_result, nonlinear_result_data], ignore_index=True)
        # neg_t_prop = solve_ivp(spatial_ode, neg_tspan, IC, args=(mu,), rtol=1e-12,atol=1e-14)
        # nonlinear_result_data = pd.DataFrame({
        #     'name': f"x0 = {x0_km} km",
        #     't':neg_t_prop.t,
        #     'x':neg_t_prop.y[0],
        #     'y':neg_t_prop.y[1]
        # })
        # nonlinear_result = pd.concat([nonlinear_result, nonlinear_result_data], ignore_index=True)

# Build plot
x_min = 0.834
x_max = 0.840
y_lim = (x_max-x_min)/2
base = alt.Chart(nonlinear_result).mark_line(strokeWidth=1,clip=True).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='plasma')).title(None),
    order='t',
).properties(
    width=400,
    height=400,
    title=f"Nonlinear Results near L1 ({ps}, Lillian Shido)"
)

eigspace = alt.Chart(eigenspace).mark_line(strokeWidth=1.5,strokeDash=(4,4),clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(scheme='greenblue')).title(None)
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

final = alt.layer(base, eigspace, L1_loc).resolve_scale(color='independent')
final.save(f'nonlinear_L1 {ps}.png', ppi=200)

# initial_velocities_table = (
#     GT(initial_velocities)
#     .tab_header(
#         title=md(f"Initial Conditions for each value of x0<br>({ps}, Lillian Shido)")
#     )
#     .cols_label(
#         x0_km="{{x_0}} [km]",
#         x0="x [non-dim]",
#         vy0="{{v_y0}} [non-dim]"
#     )
#     .opt_horizontal_padding(scale=3)
#     .fmt_number(
#         columns=["x0", "vy0"],
#         n_sigfig=6
#     )
#     .cols_align(align="center")
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
#     .opt_table_outline()
# )
# initial_velocities_table.show()