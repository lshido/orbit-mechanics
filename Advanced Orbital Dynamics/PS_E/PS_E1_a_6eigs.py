ps = "E1 part a"
# Problem D4 Part b iii
# Author: Lillian Shido
# Date: 10/26/2025

import pdb
import copy
import numpy as np
import pandas as pd
from math import pi, sqrt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
from scipy import stats
from great_tables import GT, md, html, style, loc, system_fonts
from pypalettes import load_cmap

from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, calc_initial_velocities, find_halfperiod, planar_ode, calc_poincare_exponents

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)

# Initial Conditions
xi = 0.01
eta = 0

# Distance from L1
xi_from_L1_dim = xi*l_char
eta_from_L1_dim = eta*l_char
df_distance_from_L1 = pd.DataFrame({    
    "xi_from_L1":[xi],
    "eta_from_L1":[eta],
    "xi_from_L1_dim":[xi_from_L1_dim],
    "eta_from_L1_dim":[eta_from_L1_dim]
})

# Calc starting guess values
xi_dot_0, eta_dot_0 = calc_initial_velocities(xi, eta, x_L1, y_L1, mu)
starting_x = x_L1 + xi
starting_y = y_L1 + eta
starting_xdot = xi_dot_0
ydot_guess = eta_dot_0


iterations, tf, arrival_states, converged_initial_states = find_halfperiod(starting_x, ydot_guess, mu, tolerance=1e-12)

IC = [
    converged_initial_states[0], converged_initial_states[1], converged_initial_states[2], converged_initial_states[3],
    1,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
]
# Monodromy matrix from full period
full_period_prop = solve_ivp(planar_ode, [0, 2*tf], IC, args=(mu,), rtol=1e-12,atol=1e-14)
monodromy_full = full_period_prop.y[4:20,-1].reshape(4,4)

# Monodromy matrix from half period
half_period_prop = solve_ivp(planar_ode, [0, tf], IC, args=(mu,), rtol=1e-12,atol=1e-14)
stm_half = half_period_prop.y[4:20,-1].reshape(4,4)
G = np.array([
    [1,0,0,0],
    [0,-1,0,0],
    [0,0,-1,0],
    [0,0,0,1]
    ])
I = np.identity(2)
omega = np.array([
    [0,1],
    [-1,0]
    ])
term1 = np.bmat([
    [np.zeros((2,2)) , -I],
    [I , -2*omega] 
    ])
term2 = np.transpose(stm_half)
term3 = np.bmat([
    [-2*omega , I],
    [-I , np.zeros((2,2))]
])
monodromy_half = G*term1*term2*term3*G*stm_half

df_det_error = pd.DataFrame({
    'full_det':[np.linalg.det(monodromy_full)],
    'full_error':[abs(np.linalg.det(monodromy_full)-1)],
    'half_det':[np.linalg.det(monodromy_half)],
    'half_error':[abs(np.linalg.det(monodromy_half)-1)]
})

eigenvalues = np.linalg.eigvals(monodromy_half)
df_eigenvalues = pd.DataFrame({
    "eig_1":[eigenvalues[0]],
    "eig_2":[eigenvalues[1]],
    "eig_3":[eigenvalues[2]],
    "eig_4":[eigenvalues[3]],
})
pdb.set_trace()

df_poincare_exponents = pd.DataFrame({
    "pe_1":[calc_poincare_exponents(eigenvalues[0],2*tf)],
    "pe_2":[calc_poincare_exponents(eigenvalues[1],2*tf)],
    "pe_3":[calc_poincare_exponents(eigenvalues[2],2*tf)],
    "pe_4":[calc_poincare_exponents(eigenvalues[3],2*tf)],
})

pdb.set_trace()

# # Configure det error table
# columns = [
# "eig_1", "eig_2", "eig_3", "eig_4"]
pe_table = (
    GT(df_poincare_exponents)
    .tab_header(
        title=md(f"Poincaré Exponents<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        pe_1="{{:omega:_1}}",
        pe_2="{{:omega:_2}}",
        pe_3="{{:omega:_3}}",
        pe_4="{{:omega:_4}}"
    )
    .fmt_number(
        columns=["pe_1", "pe_2", "pe_3", "pe_4"],
        decimals=5
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
pe_table.show()


# # # Configure det error table
# # columns = [
# # "eig_1", "eig_2", "eig_3", "eig_4"]
# eig_table = (
#     GT(df_eigenvalues)
#     .tab_header(
#         title=md(f"Eigenvalues of the Monodromy Matrix<br>({ps}, Lillian Shido)")
#     )
#     .cols_label(
#         eig_1="{{:lambda:_1}}",
#         eig_2="{{:lambda:_2}}",
#         eig_3="{{:lambda:_3}}",
#         eig_4="{{:lambda:_4}}"
#     )
#     .fmt_scientific(
#         columns=["eig_1", "eig_2", "eig_3", "eig_4"],
#         decimals=9
#     )
#     .cols_align(
#         align="center"
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# eig_table.show()


# # Configure tables
# # Configure det error table
# # columns = ['full_period','full_error','half_period','half_error']
# det_error_table = (
#     GT(df_det_error)
#     .tab_header(
#         title=md(f"Accuracy of the Monodromy Matrix<br>({ps}, Lillian Shido)")
#     )
#     .tab_spanner(
#         label='Full Period',
#         columns=['full_det','full_error']
#     )
#     .tab_spanner(
#         label='Half Period',
#         columns=['half_det','half_error']
#     )
#     .cols_label(
#         full_det='Determinant',
#         full_error='Error',
#         half_det='Determinant',
#         half_error='Error'
#     )
#     .fmt_scientific(
#         columns=["full_error","half_error"],
#         decimals=5
#     )
#     .cols_align(
#         align="center"
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# det_error_table.show()
