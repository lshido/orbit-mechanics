ps = "E3"
# Author: Lillian Shido
# Date: 11/6/2025

import sys
import pdb
import numpy as np
import pandas as pd
from math import pi, sqrt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from great_tables import GT, md, html, style, loc, system_fonts

from constants import r_Earth, mu_Earth
from methods import spatial_2bp_ode, calc_spatial_monodromy_half, calc_stability_index

# Circular orbit about Earth
e = 0
mu_2BP = mu_Earth

df_eigenvalues = pd.DataFrame()
df_stability = pd.DataFrame()
df_initial_conditions = pd.DataFrame()
for altitude in [100, 500, 2000, 10000]:
    a = r_Earth + altitude
    r = a
    v0 = (mu_2BP/r)**(1/2) 
    period = 2*pi/(mu_2BP/a**3)**(1/2) # sec
    pdb.set_trace()
    IC = [
        r,0,0,0,v0,0,# IC states
        1,0,0,0,0,0, # Identity matrix for phi ICs
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,1
    ]
    initial_conditions = pd.DataFrame({
        "x":[r],
        "y":[0],
        "z":[0],
        "vx":[0],
        "vy":[v0],
        "vz":[0]
    })
    df_initial_conditions = pd.concat([df_initial_conditions, initial_conditions], ignore_index=True)
    full_period_prop = solve_ivp(spatial_2bp_ode, [0, period], IC, args=(mu_2BP,r), rtol=1e-13,atol=1e-14)
    half_period_prop = solve_ivp(spatial_2bp_ode, [0, period/2], IC, args=(mu_2BP,r), rtol=1e-13,atol=1e-14)
    stm_half = half_period_prop.y[6:42,-1].reshape(6,6)
    monodromy_full = full_period_prop.y[6:42,-1].reshape(6,6)
    # monodromy_half = calc_spatial_monodromy_half(stm_half)

    np.set_printoptions(threshold=sys.maxsize)
    # eigenvalues = np.linalg.eigvals(monodromy_full)
    eigenvalues = np.linalg.eigvals(monodromy_full)
    print("\n".join([" ".join(f"{item:15.5f}" for item in row) for row in monodromy_full]))
    eigenvalues_data = pd.DataFrame({
        "altitude":[altitude],
        "r":[r],
        "ydot":[v0],
        "period":[period/3600],
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

    stability_data = pd.DataFrame({
        "altitude":[altitude],
        "r":[r],
        "ydot":[v0],
        "period":[period/3600],
        "stability_1":[calc_stability_index(eigenvalues[0],eigenvalues[1])],
        "stability_2":[calc_stability_index(eigenvalues[2],eigenvalues[3])],
        "stability_3":[calc_stability_index(eigenvalues[4],eigenvalues[5])]
    })
    df_stability = pd.concat([df_stability, stability_data], ignore_index=True)
# monodromy = calc_spatial_monodromy_half(stm_half)

fig, ax = plt.subplots()
# ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
# ax.yaxis.set_major_locator(ticker.MultipleLocator(50))
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
# ax.yaxis.set_minor_locator(ticker.MultipleLocator(25))
ax.plot(full_period_prop.y[0], full_period_prop.y[1], label='orbit')
ax.axis('equal')
ax.set(xlim=(-10000, 10000), ylim=(-10000,10000))
ax.legend(loc="upper center", fontsize=6)
plt.grid()
plt.xlabel("x [km]")
plt.ylabel("y [km]")
ax.tick_params(axis='both', which='major', labelsize=6)
ax.set_title(f'Circular orbit\n({ps}, Lillian Shido)')
plt.savefig(f'circular_orbit_{ps}.png', dpi=300, bbox_inches='tight')

eig_table = (
    GT(df_eigenvalues)
    .tab_header(
        title=md(f"All 6 Eigenvalues, Earth-Orbit 2BP<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        r="{{radius}}<br>[km]",
        altitude="{{altitude}}<br>[km]",
        ydot="{{velocity}}<br>[km/s]",
        period="{{Period}}<br>[hours]",
        eig_1="{{:lambda:_1}}",
        eig_2="{{:lambda:_2}}",
        eig_3="{{:lambda:_3}}",
        eig_4="{{:lambda:_4}}",
        eig_5="{{:lambda:_5}}",
        eig_6="{{:lambda:_6}}"
    )
    .fmt_number(
        columns=["r","ydot", "eig_1","eig_2","eig_3","eig_4","eig_5","eig_6"],
        n_sigfig=6
    )
    .fmt_number(
        columns=["period"],
        n_sigfig=3
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
    .cols_hide(columns=["jacobi","eig_1_abs","eig_2_abs","eig_3_abs","eig_4_abs","eig_5_abs","eig_6_abs","eig_1_string","eig_2_string","eig_3_string","eig_4_string","eig_5_string","eig_6_string"])
)
eig_table.show()

stability_table = (
    GT(df_stability)
    .tab_header(
        title=md(f"Stability indices, Earth-Orbit 2BP<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        r="{{radius}}<br>[km]",
        altitude="{{altitude}}<br>[km]",
        ydot="{{velocity}}<br>[km/s]",
        period="{{Period}}<br>[hours]",
        stability_1="Stability Index 1<br>{{:nu:_1}}",
        stability_2="Stability Index 2<br>{{:nu:_2}}",
        stability_3="Stability Index 3<br>{{:nu:_3}}",
    )
    .fmt_number(
        columns=["r","ydot","stability_1","stability_2","stability_3",],
        n_sigfig=6
    )
    .fmt_number(
        columns=["period"],
        n_sigfig=3
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

IC_table = (
    GT(df_initial_conditions.iloc[0].to_frame().transpose())
    .tab_header(
        title=md(f"Initial Conditions, Earth-Orbit 2BP<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        x="{{x}}<br>[km]",
        y="{{y}}<br>[km]",
        z="{{z}}<br>[km]",
        vx="{{v_x}}<br>[km/s]",
        vy="{{v_y}}<br>[km/s]",
        vz="{{v_z}}<br>[km/s]"
    )
    .fmt_number(
        columns=["x","y","z","vx","vy","vz"],
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
IC_table.show()