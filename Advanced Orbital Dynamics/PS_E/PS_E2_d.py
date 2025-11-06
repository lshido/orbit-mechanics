ps = "E2 part d"
# Continuation Algorithm: Natural Parameter Process with Dynamic Step Sizes
# Problem D4 Part d
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
import warnings

from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L2, calc_initial_velocities, find_halfperiod, calc_Jacobi, spatial_ode, calc_spatial_monodromy_half, calc_spatial_Jacobi, calc_stability_index

warnings.filterwarnings("ignore",category=PendingDeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L2, y_L2 = calc_L2(mu, a_Moon)

# Initial Conditions
xi = 0.01
eta = 0

# Distance from L2
xi_from_L2_dim = xi*l_char
eta_from_L2_dim = eta*l_char
df_distance_from_L2 = pd.DataFrame({    
    "xi_from_L2":[xi],
    "eta_from_L2":[eta],
    "xi_from_L2_dim":[xi_from_L2_dim],
    "eta_from_L2_dim":[eta_from_L2_dim]
})

# Calc starting guess values
xi_dot_0, eta_dot_0 = calc_initial_velocities(xi, eta, x_L2, y_L2, mu)
starting_x = x_L2 + xi
starting_y = y_L2 + eta
starting_xdot = xi_dot_0
ydot_guess = eta_dot_0

# Initialize list to store orbit ICs
df_orbits = pd.DataFrame(
    columns = [
        "orbit",
        "iterations",
        "tf",
        "xi",
        "x",
        "y",
        "vx",
        "vy",
        "jacobi"
    ]
)

# Configuration
first_delta_x = 0.001
delta_x = 0.007 # step in x

orbit_x = starting_x
orbit_ydot = ydot_guess
# Use the arrival states for xi (x0, y0, vx0, vy0) as the starting x
for orbit in range(21):
    print(f"starting x0:{orbit_x}, starting ydot: {orbit_ydot}")
    iterations, tf, arrival_states, converged_initial_states = find_halfperiod(orbit_x, orbit_ydot, mu, tolerance=1e-12)
    if orbit <= 1: # For the first two calculations, naively guess ydot is the same
        orbit_x = converged_initial_states[0] + first_delta_x # Step the x by delta_x
        orbit_ydot = converged_initial_states[3]
        print(f"found x0:{converged_initial_states[0]}, found ydot: {converged_initial_states[3]}")
        # Calc Jacobi
        jacobi = calc_spatial_Jacobi(mu, converged_initial_states[0], converged_initial_states[1], 0, converged_initial_states[2], converged_initial_states[3],0)
        # Compile data
        orbit_IC_data = pd.DataFrame({
            "orbit":[orbit],
            "iterations":[iterations],
            "tf":[tf],
            "xi":[converged_initial_states[0]-x_L2],
            "x":[converged_initial_states[0]],
            "y":[converged_initial_states[1]],
            "vx":[converged_initial_states[2]],
            "vy":[converged_initial_states[3]],
            "jacobi":[jacobi]
        })
        df_orbits = pd.concat([df_orbits, orbit_IC_data], ignore_index=True)
    else:
        # Calc the slope from the last pair of x0 and vy0
        my_slope = (df_orbits.iloc[-1]['vy']-df_orbits.iloc[-2]['vy'])/(df_orbits.iloc[-1]['x']-df_orbits.iloc[-2]['x'])
        if orbit_x > 1.24:
            delta_x = 0.008
        if orbit_x > 1.25:
            delta_x = 0.009
        orbit_x = converged_initial_states[0] + delta_x # Step the x by delta_x
        orbit_ydot = converged_initial_states[3] + delta_x*my_slope
        print(f"found x0:{converged_initial_states[0]}, found ydot: {converged_initial_states[3]}")
        # Calc Jacobi
        jacobi = calc_spatial_Jacobi(mu, converged_initial_states[0], converged_initial_states[1], 0, converged_initial_states[2], converged_initial_states[3], 0) 
        # Compile data
        orbit_IC_data = pd.DataFrame({
            "orbit":[orbit],
            "iterations":[iterations],
            "tf":[tf],
            "xi":[converged_initial_states[0]-x_L2],
            "x":[converged_initial_states[0]],
            "y":[converged_initial_states[1]],
            "vx":[converged_initial_states[2]],
            "vy":[converged_initial_states[3]],
            "jacobi":[jacobi]
        })
        df_orbits = pd.concat([df_orbits, orbit_IC_data], ignore_index=True)

# Table iterations each orbit takes to converge
df_orbit_iterations = df_orbits[['orbit','xi','iterations']]

# Set up colors for plot and their labels
cmap = load_cmap(name="Classic_Cyclic")
colors = cmap(np.linspace(0,1,21))
# colors = load_cmap(name="Classic_Cyclic", cmap_type='discrete').colors
labels = [fr'$\xi$={i}' if i!=0 else f'Baseline' for i, c in enumerate(colors)]

df_eigenvalues = pd.DataFrame()

# Initialize family plot
fig1, ax1 = plt.subplots()
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax1.scatter(x_Earth, 0, label='Earth', s=20, color='blue')
ax1.scatter(x_Moon, 0, label='Moon', s=20, color='gray')
ax1.scatter(x_L2, y_L2, label='L2', s=10, color='purple')
ax1.axhline(y=0, color='r', linestyle='--', linewidth=0.5)
for enum, orbit in enumerate(df_orbits.iterrows()):
    tf = orbit[1]["tf"]
    x0 = orbit[1]["x"]
    vy0 = orbit[1]["vy"]
    xi0 = orbit[1]["xi"]
    IC = [
        x0,0,0,0,vy0,0,# IC states
        1,0,0,0,0,0, # Identity matrix for phi ICs
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,1
    ]
    full_period_prop = solve_ivp(spatial_ode, [0, 2*tf], IC, args=(mu,), rtol=1e-13,atol=1e-14)
    half_period_prop = solve_ivp(spatial_ode, [0, tf], IC, args=(mu,), rtol=1e-13,atol=1e-14)
    stm_half = half_period_prop.y[6:42,-1].reshape(6,6)
    monodromy_half = calc_spatial_monodromy_half(stm_half)
    eigenvalues = np.linalg.eigvals(monodromy_half)
    # Collect eigenvalues
    eigenvalues_data = pd.DataFrame({
        "orbit":[enum],
        "xi":[xi0],
        "period":[2*tf],
        "jacobi":[jacobi],
        "eig_1":[eigenvalues[0]],
        "eig_1_abs":[abs(eigenvalues[0])],
        "eig_1_string":[f"{eigenvalues[0]:.5g}"],
        "eig_2":[eigenvalues[1]],
        "eig_2_abs":[abs(eigenvalues[1])],
        "eig_2_string":[f"{eigenvalues[1]:.5g}"],
        "eig_3":[eigenvalues[2]],
        "eig_3_abs":[abs(eigenvalues[2])],
        "eig_3_string":[f"{eigenvalues[2]:.5g}"],
        "eig_4":[eigenvalues[3]],
        "eig_4_abs":[abs(eigenvalues[3])],
        "eig_4_string":[f"{eigenvalues[3]:.5g}"],
        "eig_5":[eigenvalues[4]],
        "eig_5_abs":[abs(eigenvalues[4])],
        "eig_5_string":[f"{eigenvalues[4]:.5g}"],
        "eig_6":[eigenvalues[5]],
        "eig_6_abs":[abs(eigenvalues[5])],
        "eig_6_string":[f"{eigenvalues[5]:.5g}"]
    })
    df_eigenvalues = pd.concat([df_eigenvalues, eigenvalues_data], ignore_index=True)
    
    
    arc = [np.column_stack([full_period_prop.y[0], full_period_prop.y[1]])]
    try:
        if enum==4:
            line_collection = LineCollection(arc, colors="red", linewidth=1.5, label="Stable => Unstable")
        elif enum==10:
            line_collection = LineCollection(arc, colors="green", linewidth=1.5, label="Unstable => Stable")
        elif enum==17:
            line_collection = LineCollection(arc, colors="blue", linewidth=1.5, label="Stable => Unstable")
        else:
            line_collection = LineCollection(arc, colors="orange", linewidth=0.5)
    except:
        pdb.set_trace()
    ax1.add_collection(line_collection)

ax1.set_aspect('equal', 'box')
ax1.autoscale()
ax1.legend(loc="upper center", fontsize=6)
plt.grid()
plt.xlabel("x [non-dim]")
plt.ylabel("y [non-dim]")
ax1.tick_params(axis='both', which='major', labelsize=6)
ax1.set_title(f'Lyapunovs near L2 starting from $\\xi$={xi:.2f}, $\\eta$={eta}\n({ps}, Lillian Shido)')
plt.savefig(f'Lyapunov_family_{ps}.png', dpi=300, bbox_inches='tight')

df_stability_index = pd.DataFrame()
for enum, row in enumerate(df_eigenvalues.iterrows()):
    eig_1 = row[1]['eig_1']
    eig_2 = row[1]['eig_2']
    eig_3 = row[1]['eig_3']
    eig_4 = row[1]['eig_4']
    eig_5 = row[1]['eig_5']
    eig_6 = row[1]['eig_6']
    if enum <= 9:
        stability_data = pd.DataFrame({
            "stability_1": [calc_stability_index(eig_1, eig_2)],
            "stability_2": [calc_stability_index(eig_3, eig_4)],
            "stability_3": [calc_stability_index(eig_5, eig_6)]
        })
        df_stability_index = pd.concat([df_stability_index, stability_data], ignore_index=True)
    else:
        stability_data = pd.DataFrame({
            "stability_1": [calc_stability_index(eig_1, eig_4)],
            "stability_2": [calc_stability_index(eig_2, eig_3)],
            "stability_3": [calc_stability_index(eig_5, eig_6)]
        })
        df_stability_index = pd.concat([df_stability_index, stability_data], ignore_index=True)

df_family_stability = df_eigenvalues[['orbit','xi','period']].join(df_stability_index)

eig_table = (
    GT(df_eigenvalues)
    .tab_header(
        title=md(f"All 6 Eigenvalues in L1 Lyapunov Family<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        orbit="Orbit",
        xi="{{:xi:}}",
        period="{{Period}}<br>[non-dim]",
        eig_1_string="{{:lambda:_1}}",
        eig_2_string="{{:lambda:_2}}",
        eig_3_string="{{:lambda:_3}}",
        eig_4_string="{{:lambda:_4}}",
        eig_5_string="{{:lambda:_5}}",
        eig_6_string="{{:lambda:_6}}"
    )
    .fmt_number(
        columns=["xi","period"],
        decimals=3
    )
    # .fmt_number(
    #     columns=["eig_1","eig_2","eig_3","eig_4","eig_5","eig_6"],
    #     n_sigfig=6
    # )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
    .cols_hide(columns=["jacobi","eig_1_abs","eig_2_abs","eig_3_abs","eig_4_abs","eig_5_abs","eig_6_abs","eig_1","eig_2","eig_3","eig_4","eig_5","eig_6"])
)
eig_table.show()

eig_abs_table = (
    GT(df_eigenvalues)
    .tab_header(
        title=md(f"All 6 Eigenvalues in L1 Lyapunov Family<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        orbit="Orbit",
        xi="{{:xi:}}",
        period="{{Period}}<br>[non-dim]",
        eig_1_abs="{{| :lambda:_1 |}}",
        eig_2_abs="{{| :lambda:_2 |}}",
        eig_3_abs="{{| :lambda:_3 |}}",
        eig_4_abs="{{| :lambda:_4 |}}",
        eig_5_abs="{{| :lambda:_5 |}}",
        eig_6_abs="{{| :lambda:_6 |}}"
    )
    .fmt_number(
        columns=["xi","period"],
        decimals=3
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
    .cols_hide(columns=["jacobi","eig_1_string","eig_2_string","eig_3_string","eig_4_string","eig_5_string","eig_6_string","eig_1","eig_2","eig_3","eig_4","eig_5","eig_6"])
)
eig_abs_table.show()


stability_table = (
    GT(df_family_stability)
    .tab_header(
        title=md(f"Stability indices in L1 Lyapunov Family<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        orbit="Orbit",
        xi="{{:xi:}}",
        period="{{Period}}<br>[non-dim]",
        stability_1="Stability Index 1<br>{{:nu:_1}}",
        stability_2="Stability Index 2<br>{{:nu:_2}}",
        stability_3="Stability Index 3<br>{{:nu:_3}}",
    )
    .fmt_number(
        columns=["xi","period"],
        decimals=3
    )
    .fmt_number(
        columns=["stability_1","stability_2","stability_3",],
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
stability_table.show()
