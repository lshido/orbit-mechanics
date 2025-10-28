# Problem D2 part(d)
ps = "D2 part d"
# Varying-Time Position Targeter with max allowable delta_t
# Author: Lillian Shido
# Date: 10/21/2025

import pdb
import copy
from constants import mu_Earth, mu_Moon, a_Moon
import numpy as np
import pandas as pd
from math import pi, sqrt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
from great_tables import GT, md, html, style, loc, system_fonts
from pypalettes import load_cmap

mu = mu_Moon/(mu_Earth + mu_Moon)

def system_properties(mu_major,mu_minor, a):
    mu = mu_minor/(mu_major+mu_minor)
    x_major = -mu
    x_minor = 1-mu
    l_char = a
    t_char = sqrt(l_char**3/(mu_major+mu_minor))
    return mu,l_char, t_char, x_major, x_minor

def calc_error(actual, ideal):
    return abs(actual-ideal)/abs(ideal)

def ode(t,sv,mu):
    # Set up the EOM ODEs
    eoms = [
        sv[2],
        sv[3],
        2*sv[3] + sv[0] - (1 - mu) * (sv[0] + mu) / ((sv[0] + mu)**2 + sv[1]**2)**(3/2)-\
        mu * (sv[0] - 1 + mu) / ((sv[0] - 1 + mu)**2 + sv[1]**2)**(3/2),
        -2*sv[2] + sv[1] - (1 - mu) * sv[1] / ((sv[0] + mu)**2 + sv[1]**2)**(3/2) -\
        mu * sv[1]/((sv[0] - 1 + mu)**2 + sv[1]**2)**(3/2)
    ]
    # Calc the partials using the current x and y values
    d = sqrt((sv[0]+mu)**2 + sv[1]**2);
    r = sqrt((sv[0]-1+mu)**2 + sv[1]**2);
    U_xx = 1 - (1-mu)/d**3 - mu/r**3 + 3*(1-mu)*(sv[0]+mu)**2/d**5 + 3*mu*(sv[0]-1+mu)**2/r**5;
    U_yy = 1 - (1-mu)/d**3 - mu/r**3 + 3*(1-mu)*sv[1]**2/d**5 + 3*mu*sv[1]**2/r**5;
    U_xy = 3*(1-mu)*(sv[0]+mu)*sv[1]/d**5 + 3*mu*(sv[0]-1+mu)*sv[1]/r**5;
    # Set up the STM ODEs
    stm = [
    sv[12],
    sv[13],
    sv[14],
    sv[15],
    sv[16],
    sv[17],
    sv[18],
    sv[19],
    U_xx*sv[4] + U_xy*sv[8] + 2*sv[16],
    U_xx*sv[5] + U_xy*sv[9] + 2*sv[17],
    U_xx*sv[6] + U_xy*sv[10] + 2*sv[18],
    U_xx*sv[7] + U_xy*sv[11] + 2*sv[19],
    U_xy*sv[4] + U_yy*sv[8] - 2*sv[12],
    U_xy*sv[5] + U_yy*sv[9] - 2*sv[13],
    U_xy*sv[6] + U_yy*sv[10] - 2*sv[14],
    U_xy*sv[7] + U_yy*sv[11] - 2*sv[15]
    ]
    # Combine them into one big matrix
    combined = eoms + stm
    return combined

def crossxEvent(t, sv, mu):
    return sv[1] # Return value of y
crossxEvent.terminal = True

mu,l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)

# Trjaectory Initial Conditions
original_IC = [
    0.488, 0.200, # x and y
    -0.880, 0.200 # x_dot and y_dot
    ]

# ODE Initial Conditions
sv0_IC = [
    original_IC[0], # x
    original_IC[1], # y
    original_IC[2], # x_dot
    original_IC[3], # y_dot
    1,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
    ]


# Set the span of the integrator
t_final = 1.5*pi;
tspan = [0, t_final];

sol = solve_ivp(ode, tspan, sv0_IC, events=crossxEvent, args=(mu,), rtol=1e-12,atol=1e-14)
tf_IC = sol.t_events[0][0]

# Set initial conditions for targeter reference arc
r0_IC = np.array([
    [original_IC[0]],
    [original_IC[1]]
])
v0_IC = np.array([
    [original_IC[2]],
    [original_IC[3]]
])

# Set up dataframes for report
# Question 4
df_delta_v0 = pd.DataFrame(
    columns = ["Case", "delta_v0_mag", "delta_v0_x", "delta_v0_y", "delta_v0_mag_dim", "delta_v0_x_dim", "delta_v0_y_dim"]
)
# Question 3: For each iteration in each case, list the time interval_dim, delta_x_dim, delta_y_dim, delta_xdot_dim, delta_ydot_dim, delta_v_mag_dim
df_iterations_per_case = pd.DataFrame(
    columns = ["Case", "Iteration", "tf_dim", "delta_tf_dim", "delta_x_dim", "delta_y_dim", "delta_v0_x_dim", "delta_v0_y_dim", "delta_v0_mag_dim"]
)
# Set up colors for plot and their labels
cmap = load_cmap(name="Classic_Cyclic")
colors = cmap(np.linspace(0,1,50))
labels = [f'Iteration {i}' if i!=1 else f'Reference' for i, c in enumerate(colors, start=1) ]
# Run the targeter for a set of targets. Try scaling.
scale = 15
rf_target_list = [
    np.array([[-0.3*scale], [0.05*scale]])
    # np.array([[-0.1*scale],[0*scale]])
    # np.array([[-0.35*scale],[-.1*scale]])
]
tolerance = 1e-12 # Set tolerance
max_error = 0.1
for case, rf_target in enumerate(rf_target_list):
    r0 = copy.deepcopy(r0_IC)
    v0 = copy.deepcopy(v0_IC)
    tf = copy.deepcopy(tf_IC)
    sv0 = copy.deepcopy(sv0_IC)
    arc = []
    # Set up plots for report
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    lim = (np.linalg.norm(rf_target-r0))
    plt.xlim(-lim,lim)
    plt.ylim(-lim,lim)
    plt.axis('square')
    ax.scatter(x_Earth,0, label='Earth')
    ax.scatter(x_Moon,0, label='Moon')
    ax.scatter(original_IC[0], original_IC[1], label='Start', s=20, marker="<")
    ax.scatter(rf_target[0,0], rf_target[1,0], label='Target', s=20, marker="x")
    counter = 0
    while True:
        counter = counter + 1
        # Step 1: Using ICs, propagate EOM+STM until tf, get STM from tf and rf
        prop = solve_ivp(ode, [0, tf], sv0, args=(mu,), rtol=1e-13,atol=1e-14)
        stm = prop.y[4:20,-1].reshape(4,4) # turn into 4x4 phi matrix
        # if n==1:
        #     pdb.set_trace()
        # Pull out phi_14 and phi_24 into a vector
        phi_14_24 = np.array([
            [stm[0,3]],
            [stm[1,3]]
        ])
        rf = prop.y[0:2,-1].reshape(2,1) # position at tf, turn into 2x1 vector
        vf = prop.y[2:4,-1].reshape(2,1) # velocity at tf, turn into 2x1 vector
        # Step 2: Compare rf with rf_target
        error = rf_target - rf
        # Add the arc for this iteration
        arc = [np.column_stack([prop.y[0], prop.y[1]])]
        try:
            line_collection = LineCollection(arc, colors=colors[counter-1], label=labels[counter-1], linewidth=0.8)
        except:
            pdb.set_trace()
        ax.add_collection(line_collection)
        # Check if the error is within acceptable margins
        if (abs(error) > tolerance).any() and counter < 50: # If not, recalculate the delta_v0 and try again
            print(f"Iteration: {counter}, tf: {tf}, error: {error}")
            # Step 3: Calc new delta_vy0 and delta_t
            phi_vf = np.c_[phi_14_24, vf] # form phi and vf matrix
            deltas_vy0_tf = np.linalg.inv(phi_vf) @ error # multiply the matrices together (dot product)
            delta_vy0 = deltas_vy0_tf[0,0] # Pull out delta_vy0
            delta_tf = deltas_vy0_tf[1,0] # Pull out delta_tf
            # Check error. If error is too large, limit delta t to a portion of existing delta_t
            if (abs(error) > max_error).any():
                delta_tf = delta_tf/3
            # print(f"delta_v0: {delta_v0[0,0]}, {delta_v0[1,0]}, Iteration: {counter}")
            # Step 4: Use new deltas_ydot_t to calc new ICs
            v0 = np.array([
                [v0[0,0]], # vx0 remains unchanged
                [v0[1,0] + delta_vy0] # new vy0
            ])
            r0 = r0 # Initial position is constant; we only vary velocity in this problem
            tf = tf + delta_tf # new tf
            # Rebuild the ICs for combined EOM+STM ODEs (sv0)
            sv0 = [
                r0[0,0], # x
                r0[1,0], # y
                v0[0,0], # x_dot
                v0[1,0], # y_dot
                1,0,0,0, # Identity matrix for phi ICs
                0,1,0,0,
                0,0,1,0,
                0,0,0,1
            ]
            # Build data for report Question 3
            iteration_data = pd.DataFrame({
                "Case":[f"Case {case+1} | Non-dimensional Targets x: {rf_target[0,0]}, y: {rf_target[1,0]}"],
                "Iteration":[counter],
                "tf_dim":[tf*t_char/3600/24],
                "delta_tf_dim":[delta_tf*t_char/3600/24],
                "delta_x_dim":[error[0,0]*l_char],
                "delta_y_dim":[error[1,0]*l_char],
                "delta_v0_x_dim":[0],
                "delta_v0_y_dim":[delta_vy0*l_char/t_char],
                "delta_v0_mag_dim":[np.linalg.norm(delta_vy0)*l_char/t_char], # delta_vx0=0
            })
            df_iterations_per_case = pd.concat([df_iterations_per_case, iteration_data], ignore_index=True)
            continue
        else: # If error is within acceptable margins, break out of iterative loop
            delta_v0 = v0 - v0_IC
            delta_v0_x = delta_v0[0,0]
            delta_v0_y = delta_v0[1,0]
            delta_v0_mag = np.linalg.norm(delta_v0)

            # Build data for report Question 4
            delta_v0_data = pd.DataFrame({
                "Case":[f"{case+1}"],
                'delta_v0_mag':[delta_v0_mag],
                'delta_v0_x':[delta_v0_x],
                'delta_v0_y':[delta_v0_y],
                'delta_v0_mag_dim':[delta_v0_mag*l_char/t_char],
                'delta_v0_x_dim':[delta_v0_x*l_char/t_char],
                'delta_v0_y_dim':[delta_v0_y*l_char/t_char]
            })
            df_delta_v0 = pd.concat([df_delta_v0, delta_v0_data], ignore_index=True)
            
            # Build data for report Question 3
            iteration_data = pd.DataFrame({
                "Case":[f"Case {case+1} | Non-dimensional Targets x: {rf_target[0,0]}, y: {rf_target[1,0]}"],
                "Iteration":[counter],
                "tf_dim":[tf*t_char/3600/24],
                "delta_tf_dim":[delta_tf*t_char/3600/24],
                "delta_x_dim":[error[0,0]*l_char],
                "delta_y_dim":[error[1,0]*l_char],
                "delta_v0_x_dim":[delta_v0[0,0]*l_char/t_char],
                "delta_v0_y_dim":[delta_v0[1,0]*l_char/t_char],
                "delta_v0_mag_dim":[np.linalg.norm(delta_v0)*l_char/t_char],
            })
            df_iterations_per_case = pd.concat([df_iterations_per_case, iteration_data], ignore_index=True)
            
            # Finish building the plot and save
            ax.set_title(f'Case {case+1}: target_x={rf_target[0,0]}, target_y={rf_target[1,0]}\nIterations: {counter} ({ps}, Lillian Shido)')
            ax.legend(fontsize=8)
            plt.savefig(f'{ps}_Case_{case+1}.png', dpi=300, bbox_inches='tight')
            break
    continue

# Configure the table for Question #5:
    # columns = ["Case", "delta_v0_mag", "delta_v0_x", "delta_v0_y", "delta_v0_mag_dim", "delta_v0_x_dim", "delta_v0_y_dim"]
delta_v0_table = (
    GT(df_delta_v0)
    .tab_header(
        title=md(f"Final Total Change in Initial Velocity -<br>Time-Varying Targeter ({ps}, Lillian Shido)")
    )
    .tab_stub(rowname_col="Case")
    .tab_stubhead(label="Case")
    .tab_spanner(
        label="Non-Dimensional",
        columns=["delta_v0_mag","delta_v0_x","delta_v0_y"]
    )
    .tab_spanner(
        label="Dimensional",
        columns=["delta_v0_mag_dim","delta_v0_x_dim","delta_v0_y_dim"]
    )
    .cols_label(
        delta_v0_mag="{{:Delta:v_mag_f}}",
        delta_v0_x="{{:Delta:v_x_f}}",
        delta_v0_y="{{:Delta:v_y_f}}",
        delta_v0_mag_dim="{{:Delta:v_mag_f}}<br>[km/s]",
        delta_v0_x_dim="{{:Delta:v_x_f}}<br>[km/s]",
        delta_v0_y_dim="{{:Delta:v_y_f}}<br>[km/s]"
    )
    # .tab_style(
    #     style=style.borders(
    #         sides="right",
    #         color="lightgray",
    #         style="solid",
    #         weight="1px"
    #     ),
    #     locations=loc.body(columns=[2, 4, 6, 8])
    # )
    .fmt_number(
        columns=["delta_v0_mag", "delta_v0_x", "delta_v0_y", "delta_v0_mag_dim", "delta_v0_x_dim", "delta_v0_y_dim"],
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
delta_v0_table.show()

# Configure table for question 7
    # columns = ["Case", "Iteration", "tf_dim", "delta_x_dim", "delta_y_dim", "delta_v0_x_dim", "delta_v0_y_dim", "delta_v0_mag_dim"]
iteration_per_case_table = (
    GT(df_iterations_per_case)
    .tab_header(
        title=md(f"States at each Iteration - Time-Varying Targeter<br>({ps}, Lillian Shido)")
    )
    .tab_stub(rowname_col="Iteration", groupname_col="Case")
    .tab_stubhead(label="Iteration") 
    .tab_spanner(
        label="Change in Position", 
        columns=["delta_x_dim","delta_y_dim"]
    )
    .tab_spanner(
        label="Change in Velocity",
        columns=["delta_v0_x_dim", "delta_v0_y_dim", "delta_v0_mag_dim"]
    )
    .cols_label(
        tf_dim=html("Time Interval<br>[days]"),
        delta_tf_dim="{{:delta:t_f}}<br>[days]",
        delta_x_dim="{{:Delta:x}}<br>[km]",
        delta_y_dim="{{:Delta:y}}<br>[km]",
        delta_v0_mag_dim="{{:Delta:v_mag}}<br>[km/s]",
        delta_v0_x_dim="{{:Delta:v_x}}<br>[km/s]",
        delta_v0_y_dim="{{:Delta:v_y}}<br>[km/s]"
    )
    # .tab_style(
    #     style=style.borders(
    #         sides="right",
    #         color="lightgray",
    #         style="solid",
    #         weight="1px"
    #     ),
    #     locations=loc.body(columns=[2, 4, 6, 8])
    # )
    .fmt_number(
        columns=["delta_v0_mag_dim", "delta_v0_x_dim", "delta_v0_y_dim"],
        decimals=5
    )
    .fmt_number(
        columns=["tf_dim", "delta_tf_dim", "delta_x_dim", "delta_y_dim"],
        decimals=3
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
iteration_per_case_table.show()


# pdb.set_trace()