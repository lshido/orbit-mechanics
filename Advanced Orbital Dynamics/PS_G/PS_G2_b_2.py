ps = "G2 part b_2"
# Author: Lillian Shido
# Date: 12/8/2025

import pdb
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from great_tables import GT, md, html, style, loc, system_fonts
from pypalettes import load_cmap
from math import isclose
import warnings
import altair as alt

from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, calc_initial_velocities, find_halfperiod, calc_spatial_Jacobi, spatial_ode, calc_poincare_exponents, calc_spatial_monodromy_half, calc_stability_index

warnings.filterwarnings("ignore",category=PendingDeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0]})
earth = pd.DataFrame({'name':["Earth"],'x':[x_Earth],'y':[0]})
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[y_L1],'z':[0]})
# Initial Conditions
xi = -0.01
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
df_orbits = pd.DataFrame()

# Configuration
first_delta_x = -0.0005
delta_x = -0.001 # step in x

starting_x_1 = 0.822915
ydot_guess_1 =  0.131296

starting_x_2 = 0.780915
ydot_guess_2 =  0.445412

starting_x_3 = 0.712415
ydot_guess_3 =  0.611142

orbit_x = starting_x_1
orbit_ydot = ydot_guess_1
# Use the arrival states for xi (x0, y0, vx0, vy0) as the starting x
n_orbits = 2
for orbit in range(n_orbits):
    print(f"starting x0:{orbit_x}, starting ydot: {orbit_ydot}")
    iterations, tf, arrival_states, converged_initial_states = find_halfperiod(orbit_x, orbit_ydot, mu, tolerance=1e-12)
    if orbit <= 1: # For the first two calculations, naively guess ydot is the same
        orbit_x = converged_initial_states[0] + first_delta_x # Step the x by delta_x
        orbit_ydot = converged_initial_states[3]
        print(f"found x0:{converged_initial_states[0]}, found ydot: {converged_initial_states[3]}")
        # Calc Jacobi
        jacobi = calc_spatial_Jacobi(mu, converged_initial_states[0], converged_initial_states[1], 0,converged_initial_states[2], converged_initial_states[3],0)
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
        # if 0.75 < orbit_x < 0.81:
        #     delta_x = first_delta_x*8
        # if 0.713 < orbit_x < 0.77:
        #     delta_x = first_delta_x*3
        # if orbit_x < 0.713:
        #     delta_x = first_delta_x/2 
        orbit_x = converged_initial_states[0] + delta_x # Step the x by delta_x
        orbit_ydot = converged_initial_states[3] + delta_x*my_slope
        print(f"found x0:{converged_initial_states[0]}, found ydot: {converged_initial_states[3]}")
        # Calc Jacobi
        jacobi = calc_spatial_Jacobi(mu, converged_initial_states[0], converged_initial_states[1], 0,converged_initial_states[2], converged_initial_states[3],0)
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
colors = cmap(np.linspace(0,1,n_orbits))
# colors = load_cmap(name="Classic_Cyclic", cmap_type='discrete').colors
labels = [fr'$\xi$={i}' if i!=0 else f'Baseline' for i, c in enumerate(colors)]

# Initialize eigenvalue dataframe
df_eigenvalues = pd.DataFrame()
df_abs_eigenvalues = pd.DataFrame()
df_eigenvalues_float = pd.DataFrame()
df_jc_period = pd.DataFrame()
df_oop_stability = pd.DataFrame()
df_oop_eigs = pd.DataFrame()

# Initialize Poincaré exponents
df_poincare_exponents = pd.DataFrame()
# Initialize family plot
plot_orbits = pd.DataFrame()
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
    jacobi = calc_spatial_Jacobi(mu, x0, 0, 0, 0, vy0, 0)
    full_period_prop = solve_ivp(spatial_ode, [0, 2*tf], IC, args=(mu,), rtol=1e-13,atol=1e-14)
    half_period_prop = solve_ivp(spatial_ode, [0, tf], IC, args=(mu,), rtol=1e-13,atol=1e-14)
    stm_half = half_period_prop.y[6:42,-1].reshape(6,6)
    monodromy_half = calc_spatial_monodromy_half(stm_half)
    eigenvalues = np.linalg.eigvals(monodromy_half)
    # Collect eigenvalues
    eigenvalues_data = pd.DataFrame({
        "orbit":[enum],
        "xi":[xi0],
        "x0":[x0],
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
    
    oop_stability_data = pd.DataFrame({
        "orbit":[enum],
        "xi":[xi0],
        "x0":[x0],
        "period":[2*tf],
        "jacobi":[jacobi],
        "label":['nu_3'],
        f"stability":[calc_stability_index(eigenvalues[4],eigenvalues[5]).real]
    })
    df_oop_stability = pd.concat([df_oop_stability, oop_stability_data], ignore_index=True)

    oop_eigs_data = pd.DataFrame({
        "orbit":[enum],
        "xi":[xi0],
        "x0":[x0],
        "period":[2*tf],
        "jacobi":[jacobi],
        "eig_5_string":[f"{eigenvalues[4]:.5g}"],
        "eig_6_string":[f"{eigenvalues[5]:.5g}"],
        "eig_5_abs":[abs(eigenvalues[4])],
        "eig_6_abs":[abs(eigenvalues[5])]
    })
    df_oop_eigs = pd.concat([df_oop_eigs, oop_eigs_data], ignore_index=True)

    for e in range(0,6):
        abs_eigenvalues_data = pd.DataFrame({
            "orbit":[enum],
            "xi":[xi0],
            "x0":[x0],
            "period":[2*tf],
            "jacobi":[jacobi],
            "label":[f'lambda_{e+1}'],
            f"lambda":[abs(eigenvalues[e])]
        })
        df_abs_eigenvalues = pd.concat([df_abs_eigenvalues, abs_eigenvalues_data], ignore_index=True)
    eigenvalues_float_data = pd.DataFrame({
        "orbit":[enum],
        "xi":[xi0],
        "x0":[x0],
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

    jc_period_data = pd.DataFrame({
        "orbit":[enum],
        "xi":[xi0],
        "x0":[x0],
        "period":[2*tf],
        "jacobi":[jacobi]
    })
    df_jc_period = pd.concat([df_jc_period, jc_period_data], ignore_index=True)

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
    lyapunov_orbit_data = pd.DataFrame({
        "name":f'x={x0:.4f}',
        "orbit": enum,
        "t":full_period_prop.t,
        "x":full_period_prop.y[0],
        "y":full_period_prop.y[1]
    })
    plot_orbits = pd.concat([plot_orbits, lyapunov_orbit_data], ignore_index=True)


# Plot eigs
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,2,1)
ax2 = fig1.add_subplot(1,2,2)
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax2.yaxis.set_major_locator(ticker.MultipleLocator(1e1))
circle_outline = plt.Circle((0, 0), 1, fill=False, edgecolor='black', linewidth=1, linestyle='dashed',zorder=1.5)
for row in df_eigenvalues.iterrows():
    for i in range(1,7):
        if np.linalg.norm(row[1][f'eig_{i}']) < 1.01 and np.linalg.norm(row[1][f'eig_{i}']) > 1/1.01:
        # if not np.isreal(np.linalg.norm(row[1][f'eig_{i}'])):
        # if abs(row[1][f'eig_{i}'].imag)<1e-8:
            sc = ax1.scatter(row[1][f'eig_{i}'].real, row[1][f'eig_{i}'].imag, color=colors[row[1]['orbit']])
            ax1.annotate(row[1]['orbit'], [row[1][f'eig_{i}'].real, row[1][f'eig_{i}'].imag])
        else:
            ax2.scatter(row[1]['xi'], row[1][f'eig_{i}'].real, color=colors[row[1]['orbit']])
            ax2.annotate(row[1]['orbit'],[row[1]['xi'], row[1][f'eig_{i}'].real])
plt.legend(ncol=2,framealpha=1)
lim=4e-4
ax1.axis('equal')
legend_elements = [Line2D([0],[0], color=colors[row[1]['orbit']], label=rf"$\xi$={row[1]['xi']:.3f}") for row in df_eigenvalues.iterrows()]
fig1.legend(handles=legend_elements, loc='outside right lower')

ax1.set_title('Complex Eigenvalues')
ax1.set_xlabel("Real")
ax1.set_ylabel("Imaginary")
ax1.add_artist(circle_outline)
ax2.set_title('Real Eigenvalues')
ax2.set_yscale('log')
ax2.set_xlabel(r"$\xi$")
ax2.set_ylabel(r'$\lambda$')

fig1.suptitle(rf"Eigenvalues for each $\xi$ ({ps}, Lillian Shido)")
ax1.grid()
ax2.grid()
plt.savefig(f'eigs_complex_{ps}.png', dpi=300, bbox_inches='tight')
# plt.show()

# Make a table with just the non-unity pairs
df_stability_index = pd.DataFrame()
for enum, row in enumerate(df_eigenvalues.iterrows()):
    eig_1 = row[1]['eig_1']
    eig_2 = row[1]['eig_2']
    eig_3 = row[1]['eig_3']
    eig_4 = row[1]['eig_4']
    eig_5 = row[1]['eig_5']
    eig_6 = row[1]['eig_6']
    if enum <= 10:
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

df_family_stability = df_eigenvalues[['x0','orbit','xi','period']].join(df_stability_index)

eig_table = (
    GT(df_eigenvalues)
    .tab_header(
        title=md(f"All 6 Eigenvalues in L1 Lyapunov Family<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        orbit="Orbit",
        xi="{{:xi:}}",
        x0="{{x_0}}<br>[non-dim]",
        period="{{Period}}<br>[non-dim]",
        eig_1_string="{{:lambda:_1}}",
        eig_2_string="{{:lambda:_2}}",
        eig_3_string="{{:lambda:_3}}",
        eig_4_string="{{:lambda:_4}}",
        eig_5_string="{{:lambda:_5}}",
        eig_6_string="{{:lambda:_6}}"
    )
    .fmt_number(
        columns=["xi",'x0', "period"],
        decimals=4
    )
    .cols_align(
        align="center"
    )
    # .tab_style(
    #     style=style.fill(color="yellow"),
    #     locations=loc.body(columns=["stability_3", "x0"], rows=[3, 4, 20, 21, 41, 42])
    # )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
    .cols_hide(columns=["orbit", "jacobi","eig_1_abs","eig_2_abs","eig_3_abs","eig_4_abs","eig_5_abs","eig_6_abs","eig_1","eig_2","eig_3","eig_4","eig_5","eig_6"])
)
eig_table.show()


oop_eigs_table = (
    GT(df_oop_eigs)
    .tab_header(
        title=md(f"Out-Of-Plane Eigenvalues in L1 Lyapunov Family<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        orbit="Orbit",
        xi="{{:xi:}}",
        x0="{{x_0}}<br>[non-dim]",
        period="{{Period}}<br>[non-dim]",
        jacobi="JC",
        eig_5_string="{{:lambda:_5}}",
        eig_6_string="{{:lambda:_6}}",
        eig_5_abs="{{| :lambda:_5 |}}",
        eig_6_abs="{{| :lambda:_6 |}}"
    )
    .fmt_number(
        columns=["xi",'x0', "period", "jacobi"],
        decimals=4
    )
    .fmt_number(
        columns=["eig_5_abs","eig_6_abs"],
        n_sigfig=5
    )
    .cols_hide(
        columns=['orbit']
    )
    .cols_align(
        align="center"
    )
    # .tab_style(
    #     style=style.fill(color="yellow"),
    #     locations=loc.body(columns=["eig_5_string", "eig_6_string", "eig_5_abs", "eig_6_abs", "x0"], rows=[3, 4, 20, 21, 41, 42])
    # )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
oop_eigs_table.show()

eig_abs_table = (
    GT(df_eigenvalues)
    .tab_header(
        title=md(f"All 6 Eigenvalues in L1 Lyapunov Family<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        orbit="Orbit",
        xi="{{:xi:}}",
        x0="{{x_0}}<br>[non-dim]",
        period="{{Period}}<br>[non-dim]",
        eig_1_abs="{{| :lambda:_1 |}}",
        eig_2_abs="{{| :lambda:_2 |}}",
        eig_3_abs="{{| :lambda:_3 |}}",
        eig_4_abs="{{| :lambda:_4 |}}",
        eig_5_abs="{{| :lambda:_5 |}}",
        eig_6_abs="{{| :lambda:_6 |}}"
    )
    .fmt_number(
        columns=["x0","xi","period"],
        decimals=4
    )
    .fmt_number(
        columns=["eig_1_abs","eig_2_abs","eig_3_abs","eig_4_abs","eig_5_abs","eig_6_abs"],
        n_sigfig=6
    )
    .cols_align(
        align="center"
    )
    # .tab_style(
    #     style=style.fill(color="yellow"),
    #     locations=loc.body(columns=["eig_5_abs", "eig_6_abs", "x0"], rows=[3, 4, 20, 21, 41, 42])
    # )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
    .cols_hide(columns=["jacobi","eig_1_string","eig_2_string","eig_3_string","eig_4_string","eig_5_string","eig_6_string","eig_1","eig_2","eig_3","eig_4","eig_5","eig_6"])
)
# eig_abs_table.show()

# stability_table = (
#     GT(df_family_stability)
#     .tab_header(
#         title=md(f"Stability indices in L1 Lyapunov Family<br>({ps}, Lillian Shido)")
#     )
#     .cols_label(
#         orbit="Orbit",
#         xi="{{:xi:}}",
#         x0="{{x_0}}<br>[non-dim]",
#         period="{{Period}}<br>[non-dim]",
#         stability_1="Stability Index 1<br>{{:nu:_1}}",
#         stability_2="Stability Index 2<br>{{:nu:_2}}",
#         stability_3="Stability Index 3<br>{{:nu:_3}}",
#     )
#     .fmt_number(
#         columns=["x0","xi","period"],
#         decimals=4
#     )
#     .fmt_number(
#         columns=["stability_1","stability_2","stability_3",],
#         n_sigfig=6
#     )
#     .cols_align(
#         align="center"
#     )
#     .tab_style(
#         style=style.fill(color="yellow"),
#         locations=loc.body(columns=["stability_3", "x0"], rows=[3, 4, 20, 21, 41, 42])
#     )
#     .opt_table_outline()
#     .opt_stylize()
#     .opt_table_font(font=system_fonts(name="industrial"))
#     .opt_horizontal_padding(scale=2)
# )
# stability_table.show()

x_min = 0.35
x_max = 1.35
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2

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

orbit_3_df = plot_orbits[plot_orbits['orbit']==3]
orbit_3_df['label']='Bifurcating Orbit'

orbit_3_chart = alt.Chart(orbit_3_df).mark_line(clip=True, strokeWidth=1.5).encode(
    x='x:Q',
    y='y:Q',
    order='t',
    color=alt.Color('label:N', scale=alt.Scale(domain=['Bifurcating Orbit'], range=['darkblue'])).title(None)
)

lyapunovs_chart = alt.Chart(plot_orbits).mark_line(clip=True,strokeWidth=1).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='rainbow')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=[f"L1 Lyapunovs ({ps}, Lillian Shido)"]
)

lyapunovs_chart_layer = alt.layer(lyapunovs_chart, L1_loc, orbit_3_chart).resolve_scale(color='independent')
# lyapunovs_chart_layer.save(f'L1_lyapunovs_{ps}.png', ppi=200)

stability_chart = alt.Chart(df_oop_stability).mark_point(size=20,clip=True,filled=True).encode(
    x=alt.X('x0:Q', scale=alt.Scale(domain=[0.7,0.84]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x0:Q', axis=alt.Axis(title='x [non-dim]')),
    # y=alt.Y('stability:Q', scale=alt.Scale(domain=[0.5,1.8]), axis=alt.Axis(title='normalized lambda')),
    y=alt.Y('stability:Q', axis=alt.Axis(title='stability')),
    color=alt.Color('label:N').title(None),
).properties(
    width=400,
    # height=100,
    title=[f"Out-of-Plane Stability ({ps}, Lillian Shido)"]
)

stability_chart.save(f'out_of_plane_stability_{ps}.png', ppi=200)

pdb.set_trace()
