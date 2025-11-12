ps = "E1 part b"
# Continuation Algorithm: Natural Parameter Process with Dynamic Step Sizes
# Problem E1 part b
# Author: Lillian Shido
# Date: 10/26/2025

import pdb
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from great_tables import GT, md, system_fonts
from pypalettes import load_cmap
from math import isclose

from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, calc_initial_velocities, find_halfperiod, calc_Jacobi, planar_ode, calc_poincare_exponents, calc_planar_monodromy_half

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
first_delta_x = 0.002
delta_x = 0.010 # step in x

orbit_x = starting_x
orbit_ydot = ydot_guess
# Use the arrival states for xi (x0, y0, vx0, vy0) as the starting x
for orbit in range(11):
    print(f"starting x0:{orbit_x}, starting ydot: {orbit_ydot}")
    iterations, tf, arrival_states, converged_initial_states = find_halfperiod(orbit_x, orbit_ydot, mu, tolerance=1e-12)
    if orbit <= 1: # For the first two calculations, naively guess ydot is the same
        orbit_x = converged_initial_states[0] + first_delta_x # Step the x by delta_x
        orbit_ydot = converged_initial_states[3]
        print(f"found x0:{converged_initial_states[0]}, found ydot: {converged_initial_states[3]}")
        # Calc Jacobi
        jacobi = calc_Jacobi(mu, converged_initial_states[0], converged_initial_states[1], converged_initial_states[2], converged_initial_states[3])
        # Compile data
        orbit_IC_data = pd.DataFrame({
            "orbit":[orbit],
            "iterations":[iterations],
            "tf":[tf],
            "xi":[converged_initial_states[0]-x_L1],
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
        if orbit_x > 0.92:
            delta_x = 0.005
        if orbit_x > 0.95:
            delta_x = 0.001
        orbit_x = converged_initial_states[0] + delta_x # Step the x by delta_x
        orbit_ydot = converged_initial_states[3] + delta_x*my_slope
        print(f"found x0:{converged_initial_states[0]}, found ydot: {converged_initial_states[3]}")
        # Calc Jacobi
        jacobi = calc_Jacobi(mu, converged_initial_states[0], converged_initial_states[1], converged_initial_states[2], converged_initial_states[3])
        # Compile data
        orbit_IC_data = pd.DataFrame({
            "orbit":[orbit],
            "iterations":[iterations],
            "tf":[tf],
            "xi":[converged_initial_states[0]-x_L1],
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
colors = cmap(np.linspace(0,1,20))
# colors = load_cmap(name="Classic_Cyclic", cmap_type='discrete').colors
labels = [fr'$\xi$={i}' if i!=0 else f'Baseline' for i, c in enumerate(colors)]

# Initialize eigenvalue dataframe
df_eigenvalues_comp = pd.DataFrame(
    columns=[
        "orbit",
        "xi",
        "period",
        "eig_1_real","eig_1_imag","eig_1_plot",
        "eig_2_real","eig_2_imag","eig_2_plot",
        "eig_3_real","eig_3_imag","eig_3_plot",
        "eig_4_real","eig_4_imag","eig_4_plot"
    ]
)
df_eigenvalues = pd.DataFrame(
    columns=[
        "orbit",
        "xi",
        "period",
        "eig_1",
        "eig_1_abs",
        "eig_2",
        "eig_2_abs",
        "eig_3",
        "eig_3_abs",
        "eig_4",
        "eig_4_abs"
    ]
)
df_eigenvalues_float = pd.DataFrame(
    columns=[
        "orbit",
        "xi",
        "period",
        "eig_1",
        "eig_1_abs",
        "eig_2",
        "eig_2_abs",
        "eig_3",
        "eig_3_abs",
        "eig_4",
        "eig_4_abs"
    ]
)
df_eigenvalues_abs = pd.DataFrame(
    columns=[
        "orbit",
        "xi",
        "eig_1_abs",
        "eig_2_abs",
        "eig_3_abs",
        "eig_4_abs"
    ]
)

# Initialize Poincaré exponents
df_poincare_exponents = pd.DataFrame(
    columns=[
        "orbit",
        "xi",
        "omega_1_real","omega_1_imag",
        "omega_2_real","omega_2_imag",
        "omega_3_real","omega_3_imag",
        "omega_4_real","omega_4_imag"
    ]
)
# Initialize family plot
fig1, ax1 = plt.subplots()
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax1.scatter(x_Earth, 0, label='Earth', s=20, color='blue')
ax1.scatter(x_Moon, 0, label='Moon', s=20, color='gray')
ax1.scatter(x_L1, y_L1, label='L1', s=10, color='purple')
ax1.axhline(y=0, color='r', linestyle='--', linewidth=0.5)
for enum, orbit in enumerate(df_orbits.iterrows()):
    tf = orbit[1]["tf"]
    x0 = orbit[1]["x"]
    vy0 = orbit[1]["vy"]
    xi0 = orbit[1]["xi"]
    IC = [
        x0, 0, 0, vy0, # IC states
        1,0,0,0, # Identity matrix for phi ICs
        0,1,0,0,
        0,0,1,0,
        0,0,0,1
    ]
    full_period_prop = solve_ivp(planar_ode, [0, 2*tf], IC, args=(mu,), rtol=1e-13,atol=1e-14)
    half_period_prop = solve_ivp(planar_ode, [0, tf], IC, args=(mu,), rtol=1e-13,atol=1e-14)
    stm_half = half_period_prop.y[4:20,-1].reshape(4,4)
    monodromy_half = calc_planar_monodromy_half(stm_half)
    # monodromy_full = full_period_prop.y[4:20,-1].reshape(4,4)
    eigenvalues = np.linalg.eigvals(monodromy_half)
    # eigenvalues = np.linalg.eigvals(monodromy_full)
    # Collect eigenvalues
    eigenvalues_comp_data = pd.DataFrame({
        "orbit":[enum],
        "xi":[xi0],
        "period":[2*tf],
        "eig_1_real":[eigenvalues[0].real],
        "eig_1_imag":[eigenvalues[0].imag],
        "eig_1_plot":[{"x":[eigenvalues[0].real],"y":[eigenvalues[0].imag]}],
        "eig_2_real":[eigenvalues[1].real],
        "eig_2_imag":[eigenvalues[1].imag],
        "eig_2_plot":[{"x":[eigenvalues[1].real],"y":[eigenvalues[1].imag]}],
        "eig_3_real":[eigenvalues[2].real],
        "eig_3_imag":[eigenvalues[2].imag],
        "eig_3_plot":[{"x":[eigenvalues[2].real],"y":[eigenvalues[2].imag]}],
        "eig_4_real":[eigenvalues[3].real],
        "eig_4_imag":[eigenvalues[3].imag],
        "eig_4_plot":[{"x":[eigenvalues[3].real],"y":[eigenvalues[3].imag]}]
    })
    df_eigenvalues_comp = pd.concat([df_eigenvalues_comp, eigenvalues_comp_data], ignore_index=True)
    eigenvalues_string_data = pd.DataFrame({
        "orbit":[enum],
        "xi":[xi0],
        "period":[2*tf],
        "eig_1":[f"{eigenvalues[0]:.6g}"],
        "eig_1_abs":[abs(eigenvalues[0])],
        "eig_2":[f"{eigenvalues[1]:.6g}"],
        "eig_2_abs":[abs(eigenvalues[1])],
        "eig_3":[f"{eigenvalues[2]:.6g}"],
        "eig_3_abs":[abs(eigenvalues[2])],
        "eig_4":[f"{eigenvalues[3]:.6g}"],
        "eig_4_abs":[abs(eigenvalues[3])]
    })
    df_eigenvalues = pd.concat([df_eigenvalues, eigenvalues_string_data], ignore_index=True)

    eigenvalues_float_data = pd.DataFrame({
        "orbit":[enum],
        "xi":[xi0],
        "period":[2*tf],
        "eig_1":[eigenvalues[0]],
        "eig_1_abs":[abs(eigenvalues[0])],
        "eig_2":[eigenvalues[1]],
        "eig_2_abs":[abs(eigenvalues[1])],
        "eig_3":[eigenvalues[2]],
        "eig_3_abs":[abs(eigenvalues[2])],
        "eig_4":[eigenvalues[3]],
        "eig_4_abs":[abs(eigenvalues[3])]
    })
    df_eigenvalues_float = pd.concat([df_eigenvalues_float, eigenvalues_float_data], ignore_index=True)


    eigenvalues_abs_data = pd.DataFrame({
        "orbit":[enum],
        "xi":[xi0],
        "eig_1_abs":[abs(eigenvalues[0])],
        "eig_2_abs":[abs(eigenvalues[1])],
        "eig_3_abs":[abs(eigenvalues[2])],
        "eig_4_abs":[abs(eigenvalues[3])]
    })
    df_eigenvalues_abs = pd.concat([df_eigenvalues_abs, eigenvalues_abs_data], ignore_index=True)


    # Calc Poincaré exponents
    poincare_exponents_data = pd.DataFrame({
        "orbit":[enum],
        "xi":[xi0],
        "omega_1_real":[calc_poincare_exponents(eigenvalues[0],2*tf).real],       "omega_1_imag":[calc_poincare_exponents(eigenvalues[0],2*tf).imag],
        "omega_2_real":[calc_poincare_exponents(eigenvalues[1],2*tf).real],       "omega_2_imag":[calc_poincare_exponents(eigenvalues[1],2*tf).imag],
        "omega_3_real":[calc_poincare_exponents(eigenvalues[2],2*tf).real],       "omega_3_imag":[calc_poincare_exponents(eigenvalues[2],2*tf).imag],
        "omega_4_real":[calc_poincare_exponents(eigenvalues[3],2*tf).real],        "omega_4_imag":[calc_poincare_exponents(eigenvalues[3],2*tf).imag]
    })
    df_poincare_exponents = pd.concat([df_poincare_exponents, poincare_exponents_data], ignore_index=True)
    # Plot arcs
    arc = [np.column_stack([full_period_prop.y[0], full_period_prop.y[1]])]
    try:
        line_collection = LineCollection(arc, colors=colors[enum], linewidth=0.5, label=fr"$\xi$={xi0:.3f}")
    except:
        pdb.set_trace()
    ax1.add_collection(line_collection)

ax1.set_aspect('equal', 'box')
ax1.autoscale()
ax1.legend(fontsize=6)
plt.grid()
plt.xlabel("x [non-dim]")
plt.ylabel("y [non-dim]")
ax1.tick_params(axis='both', which='major', labelsize=6)
ax1.set_title(f'Lyapunovs near L1 starting from $\\xi$={xi:.2f}, $\\eta$={eta}\nwith predicted $\dot{{y_0}}$, Dynamic $\\Delta{{x}}$ ({ps}, Lillian Shido)')
plt.savefig(f'Lyapunov_family_{ps}.png', dpi=300, bbox_inches='tight')


# Make a table with just the non-unity pairs
df_stability_index = pd.DataFrame(columns=[
    "lambda_i",
    "reciprocal_lambda_i",
    "lambda_j",
    "stability"
]
)
for enum, row in enumerate(df_eigenvalues_float.iterrows()):
    nonunity = []
    if not isclose(abs(row[1]['eig_1']), 1, rel_tol=1e-5): nonunity.append(row[1]['eig_1'])
    if not isclose(abs(row[1]['eig_2']), 1, rel_tol=1e-5): nonunity.append(row[1]['eig_2'])
    if not isclose(abs(row[1]['eig_3']), 1, rel_tol=1e-5): nonunity.append(row[1]['eig_3'])
    if not isclose(abs(row[1]['eig_4']), 1, rel_tol=1e-5): nonunity.append(row[1]['eig_4'])
    if len(nonunity) != 2:
        print("Not a pair!")
        pdb.set_trace()
    stability = 1/2*(nonunity[0] + nonunity[1])
    stability_data = pd.DataFrame({
        "orbit": [row[1]['orbit']],
        "xi": [row[1]['xi']],
        "lambda_i": [f"{nonunity[0]:.6g}"],
        "reciprocal_lambda_i":[f"{1/nonunity[0]:.6g}"],
        "lambda_j": [f"{nonunity[1]:.6g}"],
        "stability": [f"{stability:.6g}"]
    })
    df_stability_index = pd.concat([df_stability_index, stability_data], ignore_index=True)

# Plot eigs
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,2,1)
ax2 = fig1.add_subplot(1,2,2)
ax1.xaxis.set_major_locator(ticker.MultipleLocator(500))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(500))
ax2.yaxis.set_major_locator(ticker.MultipleLocator(1e1))
circle_outline = plt.Circle((0, 0), 1, fill=False, edgecolor='black', linewidth=1, linestyle='dashed',zorder=1.5)
for row in df_eigenvalues_float.iterrows():
    for i in range(1,5):
        if np.isclose(np.linalg.norm(row[1][f'eig_{i}']),1,atol=1e-6):
            continue
        else:
            # if np.linalg.norm(row[1][f'eig_{i}']) < 1.01 and np.linalg.norm(row[1][f'eig_{i}']) > 1/1.01:
            # if not np.isreal(np.linalg.norm(row[1][f'eig_{i}'])):
            # if abs(row[1][f'eig_{i}'].imag)<1e-8:
            ax1.scatter(row[1][f'eig_{i}'].real, row[1][f'eig_{i}'].imag, color=colors[row[1]['orbit']])
            ax1.annotate(row[1]['orbit'], [row[1][f'eig_{i}'].real, row[1][f'eig_{i}'].imag])
            # else:
            ax2.scatter(row[1]['xi'], row[1][f'eig_{i}'].real, color=colors[row[1]['orbit']])
            ax2.annotate(row[1]['orbit'],[row[1]['xi'], row[1][f'eig_{i}'].real])
plt.legend(ncol=2,framealpha=1)
lim=4e-4
ax1.axis('equal')
legend_elements = [Line2D([0],[0], color=colors[row[1]['orbit']], label=rf"$\xi$={row[1]['xi']:.3f}") for row in df_eigenvalues.iterrows()]
fig1.legend(handles=legend_elements, loc='outside right lower')
# ax2.legend(handles=legend_elements, loc='outside left upper')
# ax1.set_aspect(aspect=1, adjustable="box")
# ax1.set(xlim=(x_L1-lim, x_L1+lim),ylim=(-lim,lim))
ax1.set_title('Eigenvalues on the Complex Plane')
ax1.set_xlabel("Real")
ax1.set_ylabel("Imaginary")
ax1.add_artist(circle_outline)
ax2.set_title('Eigenvalues on a log scale')
ax2.set_yscale('log')
ax2.set_xlabel(r"$\xi$")
ax2.set_ylabel(r'$\lambda$')
# ax2.set(yticks=[2e2,4e2,6e2,8e2,1e3,2e3], yticklabels=[r"$2x10^2$",r"$4x10^2$",r"$6x10^2$",r"$8x10^2$",r"$1x10^3$",r"$2x10^3$"]) 
# plt.colorbar(sc)
fig1.suptitle(rf"Eigenvalues for each $\xi$ ({ps}, Lillian Shido)")
ax1.grid()
ax2.grid()
plt.savefig(f'eigs_complex_{ps}.png', dpi=300, bbox_inches='tight')
plt.show()

# eig_abs_table = (
#     GT(df_eigenvalues_abs)
#     .tab_header(
#         title=md(f"Absolute Values of the Eigenvalues of Orbits in L1 Lyapunov Family<br>1/2 period-derived monodromy matrix ({ps}, Lillian Shido)")
#     )
#     .cols_label(
#         orbit="Orbit",
#         xi="{{:xi:}}",
#         eig_1_abs="{{| :lambda:_1 |}}",
#         eig_2_abs="{{| :lambda:_2 |}}",
#         eig_3_abs="{{| :lambda:_3 |}}",
#         eig_4_abs="{{| :lambda:_4 |}}"
#     )
#     .fmt_number(
#         columns=["xi"],
#         decimals=3
#     )
#     .fmt_number(
#         columns=["eig_1_abs","eig_2_abs","eig_3_abs","eig_4_abs",],
#         n_sigfig=6
#     )
#     .cols_align(
#         align="center"
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# eig_abs_table.show()

stability_table = (
    GT(df_stability_index)
    .tab_header(
        title=md(f"Stability Indices of Non-Unity Eigenvalues of Orbits in L1 Lyapunov Family<br>1/2 period-derived monodromy matrix ({ps}, Lillian Shido)")
    )
    .cols_label(
        orbit="Orbit",
        xi="{{:xi:}}",
        lambda_i="{{:lambda:_i}}",
        reciprocal_lambda_i="{{Reciprocal of :lambda:_i}}",
        lambda_j="{{:lambda:_j}}",
        stability="{{:nu:}}"
    )
    .fmt_number(
        columns=["orbit"],
        decimals=0
    )
    .fmt_number(
        columns=["xi","period"],
        decimals=3
    )
    .fmt_number(
        columns=["nu"],
        n_sigfig=6
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
    .cols_move_to_start(columns=["orbit","xi"])
)
stability_table.show()


# Configure tables
# eig_table = (
#     GT(df_eigenvalues)
#     .tab_header(
#         title=md(f"Eigenvalues of Orbits in L1 Lyapunov Family<br>1/2 period-derived monodromy matrix ({ps}, Lillian Shido)")
#     )
#     .cols_label(
#         orbit="Orbit",
#         xi="{{:xi:}}",
#         period="Period [non-dim]",
#         eig_1="{{:lambda:_1}}",
#         eig_1_abs="{{| :lambda:_1 |}}",
#         eig_2="{{:lambda:_2}}",
#         eig_2_abs="{{| :lambda:_2 |}}",
#         eig_3="{{:lambda:_3}}",
#         eig_3_abs="{{| :lambda:_3 |}}",
#         eig_4="{{:lambda:_4}}",
#         eig_4_abs="{{| :lambda:_4 |}}"
#     )
#     .fmt_number(
#         columns=["xi","period"],
#         decimals=3
#     )
#     .fmt_number(
#         columns=["eig_1_abs","eig_2_abs","eig_3_abs","eig_4_abs"],
#         n_sigfig=6
#     )
#     .tab_style(
#         style=style.fill(color="yellow"),
#         locations=loc.body(columns=["eig_1","eig_2"], rows=[0,1,2,3,4])
#     )
#     .tab_style(
#         style=style.fill(color="yellow"),
#         locations=loc.body(columns=["eig_1","eig_4"], rows=[5,6,7,8,9,10])
#     )
#     .cols_align(
#         align="center"
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
#     .cols_hide(columns=["eig_1_abs","eig_2_abs","eig_3_abs","eig_4_abs"])
# )
# eig_table.show()


# Configure det error table
# columns = [
# # "orbit","xi","eig_1", "eig_2", "eig_3", "eig_4"]
# eig_comp_table = (
#     GT(df_eigenvalues_comp)
#     .tab_header(
#         title=md(f"Eigenvalues of Orbits in L1 Lyapunov Family<br>1/2 period-derived monodromy matrix ({ps}, Lillian Shido)")
#     )
#     .tab_spanner(
#         label="{{:lambda:_1}}",
#         columns=["eig_1_real","eig_1_imag","eig_1_plot"]
#     )
#     .tab_spanner(
#         label="{{:lambda:_2}}",
#         columns=["eig_2_real","eig_2_imag","eig_2_plot"]
#     )
#     .tab_spanner(
#         label="{{:lambda:_3}}",
#         columns=["eig_3_real","eig_3_imag","eig_3_plot"]
#     )
#     .tab_spanner(
#         label="{{:lambda:_4}}",
#         columns=["eig_4_real","eig_4_imag","eig_4_plot"]
#     )
#     .cols_label(
#         orbit="Orbit",
#         xi="{{:xi:}}",
#         eig_1_real="Real",
#         eig_1_imag="Imag",
#         eig_2_real="Real",
#         eig_2_imag="Imag",
#         eig_3_real="Real",
#         eig_3_imag="Imag",
#         eig_4_real="Real",
#         eig_4_imag="Imag"
#     )
#     .fmt_number(
#         columns=["xi"],
#         decimals=3
#     )
#     .fmt_number(
#         columns=[
#             "eig_1_real",
#             "eig_2_real",
#             "eig_3_real",
#             "eig_4_real",
#         ],
#         n_sigfig=7
#     )
#     .fmt_number(
#         columns=[
#             "eig_1_imag",
#             "eig_2_imag",
#             "eig_3_imag",
#             "eig_4_imag"
#         ],
#         n_sigfig=7,
#         pattern="{x} i"
#     )
#     .tab_style(
#         style=style.borders(
#             sides="right",
#             color="lightgray",
#             style="solid",
#             weight="1px"
#         ),
#         locations=loc.body(columns=[1, 3, 5, 7])
#     )
#     .cols_align(
#         align="center"
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
#     .cols_hide(
#         columns=["eig_1_plot","eig_2_plot","eig_3_plot","eig_4_plot"]
#     )
#     # .fmt_nanoplot(columns = "eig_1_plot", reference_area = [-1,1], reference_line = 0)
#     # .fmt_nanoplot(columns = "eig_2_plot", reference_area = [-1,1], reference_line = 0)
#     # .fmt_nanoplot(columns = "eig_3_plot", reference_area = [-1,1], reference_line = 0)
#     # .fmt_nanoplot(columns = "eig_4_plot", reference_area = [-1,1], reference_line = 0)
# )
# eig_comp_table.show()

# Configure tables
# Configure det error table
# columns = [
# "orbit","xi","eig_1", "eig_2", "eig_3", "eig_4"]
# pe_table = (
#     GT(df_poincare_exponents)
#     .tab_header(
#         title=md(f"Poincaré Exponents of Orbits in L1 Lyapunov Family<br>1/2 period-derived monodromy matrix ({ps}, Lillian Shido)")
#     )
#     .tab_spanner(
#         label="{{:omega:_1}}",
#         columns=["omega_1_real","omega_1_imag"]
#     )
#     .tab_spanner(
#         label="{{:omega:_2}}",
#         columns=["omega_2_real","omega_2_imag"]
#     )
#     .tab_spanner(
#         label="{{:omega:_3}}",
#         columns=["omega_3_real","omega_3_imag"]
#     )
#     .tab_spanner(
#         label="{{:omega:_4}}",
#         columns=["omega_4_real","omega_4_imag"]
#     )
#     .cols_label(
#         orbit="Orbit",
#         xi="{{:xi:}}",
#         omega_1_real="Real",
#         omega_1_imag="Imag",
#         omega_2_real="Real",
#         omega_2_imag="Imag",
#         omega_3_real="Real",
#         omega_3_imag="Imag",
#         omega_4_real="Real",
#         omega_4_imag="Imag"
#     )
#     .fmt_number(
#         columns=["xi"],
#         decimals=3
#     )
#     .fmt_number(
#         columns=[
#             "omega_1_real",
#             "omega_2_real",
#             "omega_3_real",
#             "omega_4_real",
#         ],
#         n_sigfig=4
#     )
#     .fmt_number(
#         columns=[
#             "omega_1_imag",
#             "omega_2_imag",
#             "omega_3_imag",
#             "omega_4_imag"
#         ],
#         n_sigfig=4,
#         pattern="{x} i"
#     )
#     .tab_style(
#         style=style.borders(
#             sides="right",
#             color="lightgray",
#             style="solid",
#             weight="1px"
#         ),
#         locations=loc.body(columns=[1, 3, 5, 7])
#     )
#     .cols_align(
#         align="center"
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# pe_table.show()

# # Configure table for question 7
# #     columns = ['Orbit','xi','iterations']
# iterations_table = (
#     GT(df_orbit_iterations)
#     .tab_header(
#         title=md(f"Iterations to Produce Periodic Orbit<br>({ps}, Lillian Shido)")
#     )
#     .cols_label(
#         orbit="Orbit",
#         xi="{{:xi:}}",
#         iterations="Iterations",
#     )
#     .fmt_number(
#         columns=["xi"],
#         decimals=3
#     )
#     .cols_align(
#         align="center"
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# iterations_table.show()
