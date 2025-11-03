ps = "E2 part a"
# Author: Lillian Shido
# Date: 11/2/2025

import pdb
import numpy as np
import pandas as pd
from math import pi
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import solve_ivp
from great_tables import GT, md, system_fonts

from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, calc_initial_velocities, find_halfperiod, calc_Jacobi, spatial_ode, calc_poincare_exponents, calc_spatial_monodromy_half

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)

# Initial Conditions
xi = 0.025
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

iterations, tf, arrival_states, converged_initial_states = find_halfperiod(starting_x, -0.174029, mu, tolerance=1e-12)

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
half_period_prop = solve_ivp(spatial_ode, [0, tf], IC, args=(mu,), rtol=1e-12,atol=1e-14)
stm_half = half_period_prop.y[6:42,-1].reshape(6,6)
monodromy = calc_spatial_monodromy_half(stm_half)
df_monodromy = pd.DataFrame({
    'Determinant':[np.linalg.det(monodromy)],
    'Accuracy':[abs(np.linalg.det(monodromy))-1]
})

# Plot the full orbit
full_period_prop = solve_ivp(spatial_ode, [0, 2*tf], IC, args=(mu,), rtol=1e-12,atol=1e-14)
period = full_period_prop.t[-1]

df_period = pd.DataFrame({
    'Period_nondim':[period],
    'Period_dim': [period*t_char/3600/24],
    'primary_period':[2*pi*t_char/3600/24]
})

df_state_error = pd.DataFrame({
    'error_x_nondim':[full_period_prop.y[0,0]-full_period_prop.y[0,-1]],
    'error_y_nondim':[full_period_prop.y[1,0]-full_period_prop.y[1,-1]],
    'error_vx_nondim':[full_period_prop.y[3,0]-full_period_prop.y[3,-1]],
    'error_vy_nondim':[full_period_prop.y[4,0]-full_period_prop.y[4,-1]],
    'error_x_dim':[(full_period_prop.y[0,0]-full_period_prop.y[0,-1])*l_char],
    'error_y_dim':[(full_period_prop.y[1,0]-full_period_prop.y[1,-1])*l_char],
    'error_vx_dim':[(full_period_prop.y[3,0]-full_period_prop.y[3,-1])*l_char],
    'error_vy_dim':[(full_period_prop.y[4,0]-full_period_prop.y[4,-1])*l_char]
})

# Initialize plot
fig1, ax1 = plt.subplots()
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
# ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
# ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
# ax1.scatter(x_Earth, 0, label='Earth', s=20, color='blue')
ax1.scatter(x_Moon, 0, label='Moon', s=20, color='gray')
ax1.scatter(x_L1, y_L1, label='L1', s=10, color='purple')
ax1.plot(full_period_prop.y[0], full_period_prop.y[1], linewidth=1, label="Periodic Orbit")
ax1.axhline(y=0, color='r', linestyle='--', linewidth=0.5)
ax1.set_aspect('equal', 'box')
ax1.autoscale()
ax1.legend(fontsize=6)
plt.grid()
plt.xlabel("x [non-dim]")
plt.ylabel("y [non-dim]")
ax1.tick_params(axis='both', which='major', labelsize=6)
ax1.set_title(f'Periodic Orbit near L1 with $\\xi$={xi:.2f}, $\\eta$={eta}\n({ps}, Lillian Shido)')
plt.savefig(f'periodic_orbit_{ps}.png', dpi=300, bbox_inches='tight')

monodromy_error_table = (
    GT(df_monodromy)
    .tab_header(
        title=md(f"Accuracy of Monodromy Matrix<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        Determinant="Determinant",
        Accuracy="Accuracy"
    )
    .fmt_scientific(
        columns=["Accuracy"],
        n_sigfig=6
    )
    .fmt_number(
        columns=["Determinant"],
        n_sigfig=8
    )
    .cols_align(align="center")
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
monodromy_error_table.show()

# state_error_table = (
#     GT(df_state_error)
#     .tab_header(
#         title=md(f"State Error after one period<br>({ps}, Lillian Shido)")
#     )
#     .tab_spanner(
#         label="Non-Dimensional",
#         columns=["error_x_nondim","error_y_nondim","error_vx_nondim","error_vy_nondim"]
#     )
#     .tab_spanner(
#         label="Dimensional",
#         columns=["error_x_dim","error_y_dim","error_vx_dim","error_vy_dim"]
#     )
#     .cols_label(
#         error_x_nondim="{{x}}<br>[non-dim]",
#         error_y_nondim="{{y}}<br>[non-dim]",
#         error_vx_nondim="{{v_x}}<br>[non-dim]",
#         error_vy_nondim="{{v_y}}<br>[non-dim]",
#         error_x_dim="{{x}}<br>[km]",
#         error_y_dim="{{y}}<br>[km]",
#         error_vx_dim="{{v_x}}<br>[km/s]",
#         error_vy_dim="{{v_y}}<br>[km/s]"
#     )
#     .fmt_scientific(
#         columns=["error_x_nondim","error_y_nondim","error_vx_nondim","error_vy_nondim","error_x_dim","error_y_dim","error_vx_dim","error_vy_dim"],
#         n_sigfig=6
#     )
#     .cols_align(align="center")
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# state_error_table.show()

# period_table = (
#     GT(df_period)
#     .tab_header(
#         title=md(f"Period of Orbit<br>({ps}, Lillian Shido)")
#     )
#     .cols_label(
#         Period_nondim="{{Orbit Period}}<br>[non-dim]",
#         Period_dim="{{Orbit Period}}<br>[days]",
#         primary_period="{{Period of}}<br>Primaries [days]"
#     )
#     .fmt_number(
#         columns=["Period_nondim", "Period_dim", "primary_period"],
#         decimals=4
#     )
#     .cols_align(align="center")
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# period_table.show()

pdb.set_trace()