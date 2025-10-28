ps = "D4 part a"
# Problem D4 Part a
# Author: Lillian Shido
# Date: 10/26/2025

import pdb
import copy
from constants import mu_Earth, mu_Moon, a_Moon
import numpy as np
import pandas as pd
from math import pi, sqrt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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

def calc_L1(mu, a):
    gamma = 0.2
    tolerance = 1e-12
    counter = 0
    while True:
        counter = counter + 1
        f = 1-mu-gamma-((1-mu)/(1-gamma)**2)+(mu/gamma**2)
        f_prime = -1-(2*(1-mu)/(1-gamma)**3)-(2*mu/gamma**3)
        if abs(f) > tolerance:
            gamma = gamma - f/f_prime;
            continue
        else:
            x = 1 - mu - gamma
            # Add check for accleration
            d_x = x*a+mu
            r_x = x*a-1+mu
            d = ((x*a+mu)**2)**(1/2)
            r = ((x*a-1+mu)**2)**(1/2)
            accel = -(1-mu)/d**3*d_x - mu/r**3*r_x
            # Add check for partial wrt gamma
            partial = (1-mu)/(1-gamma)**2 - mu/(gamma)**2 - (1-mu-gamma)
            x_L1 = x
            y_L1 = 0
            counter = 0
            break
    return x_L1, y_L1

def calc_initial_velocities(xi_0, eta_0, x_L, y_L, mu):
    """
    Calculates the initial velocities around the collinear libration points.
    
    Args:
        xi_0: x-direction distance offset from the libration point
        eta_0: y-direction distance offset from the libration point
        x_L: The x value of the libration point
        y_L: The y value of the libration point
        mu: The mu value for the system

    Returns:
        xi_dot: The initial velocity in the xi direction
        eta_dot: The initial velocity in the eta direction
    """
    d = sqrt((x_L+mu)**2 + y_L**2)
    r = sqrt((x_L-1+mu)**2 + y_L**2)
    U_xx = 1 - (1-mu)/d**3 - mu/r**3 + 3*(1-mu)*(x_L+mu)**2/d**5 + 3*mu*(x_L-1+mu)**2/r**5
    U_yy = 1 - (1-mu)/d**3 - mu/r**3
    B_1 = 2 - (U_xx + U_yy)/2
    B_2_squared = -U_xx*U_yy
    s = (B_1 + (B_1**2 + B_2_squared)**(1/2))**(1/2)
    B_3 = (s**2 + U_xx)/(2*s)
    xi_dot_0 = eta_0*s/B_3
    eta_dot_0 = -B_3*xi_0*s
    return xi_dot_0, eta_dot_0


def eoms(t,sv,mu):
    
    # Set up the EOM ODEs
    eoms = [
        sv[2],
        sv[3],
        2*sv[3] + sv[0] - (1 - mu) * (sv[0] + mu) / ((sv[0] + mu)**2 + sv[1]**2)**(3/2)-\
        mu * (sv[0] - 1 + mu) / ((sv[0] - 1 + mu)**2 + sv[1]**2)**(3/2),
        -2*sv[2] + sv[1] - (1 - mu) * sv[1] / ((sv[0] + mu)**2 + sv[1]**2)**(3/2) -\
        mu * sv[1]/((sv[0] - 1 + mu)**2 + sv[1]**2)**(3/2)
    ]
    return eoms

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
crossxEvent.terminal = 2

def eval_acceleration(x,y,xdot,ydot,mu):
    # Simply evaluate the EOMs with the position and velocity
    a_x = 2*ydot + x - (1 - mu) * (x + mu) / ((x + mu)**2 + y**2)**(3/2)-\
        mu * (x - 1 + mu) / ((x - 1 + mu)**2 + y**2)**(3/2)
    a_y = -2*xdot + y - (1 - mu) * y / ((x + mu)**2 + y**2)**(3/2) -\
        mu * y/((x - 1 + mu)**2 + y**2)**(3/2)
    return a_x, a_y

def eval_kinematic(x0,y0,xdot_0,ydot_0,a_x, a_y,t_span):
    # Evaluate the kinematic equations with the time
    xdot = xdot_0 + a_x*t_span
    ydot = ydot_0 + a_y*t_span
    x = x0 + xdot_0*t_span
    y = y0 + ydot_0*t_span
    return x,y,xdot,ydot

def calc_Jacobi(mu, x, y, vx, vy):
    d = sqrt((x+mu)**2 + y**2)
    r = sqrt((x-1+mu)**2 + y**2)
    x_y_sq = (x**2+y**2)/2
    term_1 = (1-mu)/d
    term_2 = mu/r
    pseudo_U = term_1 + term_2 + x_y_sq
    v_squared = vx**2 + vy**2
    C = 2*pseudo_U - v_squared
    return C

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

xi_dot_0, eta_dot_0 = calc_initial_velocities(xi, eta, x_L1, y_L1, mu)

x0_IC = x_L1 + xi
y0_IC = y_L1 + eta
vx0_IC = xi_dot_0
vy0_IC = eta_dot_0

# Trjaectory Initial Conditions
original_IC = [
    x0_IC, y0_IC, # x and y
    vx0_IC, vy0_IC # x_dot and y_dot
    ]
df_ICs = pd.DataFrame({
    'x0_IC':[x0_IC],
    'y0_IC':[y0_IC],
    'vx0_IC':[vx0_IC],
    'vy0_IC':[vy0_IC],
    'x0_IC_dim':[x0_IC*l_char],
    'y0_IC_dim':[y0_IC*l_char],
    'vx0_IC_dim':[vx0_IC*l_char/t_char],
    'vy0_IC_dim':[vy0_IC*l_char/t_char]
})

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

# Set initial conditions for targeter reference arc
r0_IC = np.array([
    [original_IC[0]],
    [original_IC[1]]
])
v0_IC = np.array([
    [original_IC[2]],
    [original_IC[3]]
])
print(f"Initial Conditions: x: {r0_IC[0,0]}, y:{r0_IC[1,0]}, vx:{v0_IC[0,0]}, vy:{v0_IC[1,0]}")
# Set the span of the integrator
t_final = 15*pi;
tspan = [0, t_final];

sol = solve_ivp(ode, tspan, sv0_IC, events=crossxEvent, args=(mu,), rtol=1e-12,atol=1e-14)
try:
    tf_IC = sol.t_events[0][1]
except:
    pdb.set_trace()
initial_x_cross_vx = sol.y[2,-1] # xdot column
print(f"xdot at first crossing: {initial_x_cross_vx}")

# Set up colors for plot and their labels
# colors = load_cmap(name="Classic_Cyclic", cmap_type='discrete').colors
cmap = load_cmap(name="Classic_Cyclic")
colors = cmap(np.linspace(0,1,50))
labels = [f'Iteration {i}' if i!=1 else f'Reference' for i, c in enumerate(colors, start=1) ]

tolerance = 1e-12 # Set tolerance

r0 = copy.deepcopy(r0_IC)
v0 = copy.deepcopy(v0_IC)
sv0 = copy.deepcopy(sv0_IC)
tf = copy.deepcopy(tf_IC)
vy0 = v0[1,0]
arc = []

# Question 3: For each iteration in each case, list the time interval_dim, delta_x_dim, delta_y_dim, delta_xdot_dim, delta_ydot_dim, delta_v_mag_dim
df_iterations_per_case = pd.DataFrame(
    columns = ["Iteration", "x_dim", "y_dim", "v_x_dim", "v_y_dim"]
)

# Set up plots for report
fig1, ax1 = plt.subplots(figsize=(6.5, 6.5))
ax1.set_aspect('equal', 'box')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
# ax1.scatter(x_Earth, 0, label='Earth', s=20, color='blue')
ax1.scatter(x_Moon, 0, label='Moon', s=20, color='gray')
plt.grid()
ax1.scatter(original_IC[0], original_IC[1], label='Start', s=20, marker="<", color='green')
ax1.scatter(x_L1, y_L1, label='L1', s=10, color='purple')
ax1.axhline(y=0, color='r', linestyle='--', linewidth=0.5)
counter = 0
while True:
    counter = counter + 1
    # Step 1: Using ICs, propagate EOM+STM until second x-axis cross, get STM, get xdot, ydot
    prop = solve_ivp(ode, tspan, sv0, events=crossxEvent, args=(mu,), rtol=1e-12,atol=1e-14)
    stm = prop.y[4:20,-1].reshape(4,4) # turn into 4x4 phi matrix
    tf = prop.t_events[0][1]
    phi_34 = stm[2,3]
    phi_24 = stm[1,3]
    rf = prop.y[0:2,-1].reshape(2,1) # position at tf, turn into 2x1 vector
    vf = prop.y[2:4,-1].reshape(2,1) # velocity at tf, turn into 2x1 vector
    vxf = prop.y[2,-1] # get xdot to compare against tolerance
    vyf = prop.y[3,-1] # get ydot to calculate the acceleration
    # Add the arc for this iteration
    arc = [np.column_stack([prop.y[0], prop.y[1]])]
    try:
        line_collection = LineCollection(arc, colors=colors[counter-1], label=labels[counter-1], linewidth=0.8)
    except:
        pdb.set_trace()
    ax1.add_collection(line_collection)
    # Check if the vxf is close to 0 within acceptable margins
    if abs(vxf) > tolerance and counter < 20: # If not, recalculate the delta_v0 and try again
        # Calc the acceleration in x
        a_x, a_y = eval_acceleration(rf[0,0], rf[1,0], vxf, vyf, mu)
        # Step 3: Calc new delta_vy0
        delta_vy0 = -vxf / (phi_34 - phi_24*(a_x/vyf))
        # Step 4: Use new deltas_ydot_t to calc new ICs
        vy0 = vy0 + delta_vy0
        r0 = r0 # Initial position is constant; we only vary velocity in this problem
        # Rebuild the ICs for combined EOM+STM ODEs (sv0)
        sv0 = [
            r0[0,0], # x
            r0[1,0], # y
            v0[0,0], # x_dot
            vy0, # y_dot
            1,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,
            0,0,1,0,
            0,0,0,1
        ]
        # Build data for report Question 3
        iteration_data = pd.DataFrame({
            "Iteration":[counter],
            "x_dim":[rf[0,0]*l_char],
            "y_dim":[rf[1,0]*l_char],
            "v_x_dim":[vxf*l_char/t_char],
            "v_y_dim":[vyf*l_char/t_char]
        })
        df_iterations_per_case = pd.concat([df_iterations_per_case, iteration_data], ignore_index=True)
        continue
    else: # If error is within acceptable margins, break out of iterative loop
        # Build the dataframe for the states
        # columns = "arrival_time", "x", "y", "xdot", "ydot", "arrival_time_dim", "x_dim", "y_dim", "xdot_dim", "ydot_dim"
        df_arrival_states = pd.DataFrame({
            "arrival_time":[tf],
            "x":[rf[0,0]],
            "y":[rf[1,0]],
            "xdot":[vf[0,0]],
            "ydot":[vf[1,0]],
            "arrival_time_dim":[tf*t_char/3600/24],
            "x_dim":[rf[0,0]*l_char],
            "y_dim":[rf[1,0]*l_char],
            "xdot_dim":[vf[0,0]*l_char/t_char],
            "ydot_dim":[vf[1,0]*l_char/t_char]
        })
        print(f"final xdot:{vf[0,0]}")

        # Build the dataframe for the final initial states
        df_final_initial_states = pd.DataFrame({
            "half_period":[tf],
            "initial_x":[sv0[0]],
            "initial_y":[sv0[1]],
            "initial_xdot":[sv0[2]],
            "initial_ydot":[sv0[3]],
            "half_period_dim":[tf*t_char/3600/24],
            "initial_x_dim":[sv0[0]*l_char],
            "initial_y_dim":[sv0[1]*l_char],
            "initial_xdot_dim":[sv0[2]*l_char/t_char],
            "initial_ydot_dim":[sv0[3]*l_char/t_char]
        })

        # Build data for report Question 3
        iteration_data = pd.DataFrame({
            "Iteration":[counter],
            "x_dim":[rf[0,0]*l_char],
            "y_dim":[rf[1,0]*l_char],
            "v_x_dim":[vxf*l_char/t_char],
            "v_y_dim":[vyf*l_char/t_char]
        })
        df_iterations_per_case = pd.concat([df_iterations_per_case, iteration_data], ignore_index=True)
        # Finish building the plot and save
        # ax1.scatter(rf[0,0], 0, marker='x', label='Perpendicular Cross', color="green")
        ax1.set_title(f'Perpendicular Crossing near L1 with $\\xi$={xi}, $\\eta$={eta}\nIterations: {counter} ({ps}, Lillian Shido)')
        # ax1.legend(loc='lower left', fontsize=6)
        plt.savefig(f'Iterations_{ps}.png', dpi=300, bbox_inches='tight')
        break
    

full_period_prop = solve_ivp(ode, [0, 2*tf], sv0, args=(mu,), rtol=1e-12,atol=1e-14)
monodromy = full_period_prop.y[4:20,-1].reshape(4,4)
fig2, ax2 = plt.subplots(figsize=(6.5, 6.5))
ax2.set_aspect('equal', 'box')
ax2.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
# ax2.scatter(x_Earth, 0, label='Earth', s=20, color='blue')
ax2.scatter(x_Moon, 0, label='Moon', s=20, color='gray')
plt.grid()
ax2.scatter(original_IC[0], original_IC[1], label='Start', s=20, marker="v", color='green')
ax2.scatter(x_L1, y_L1, label='L1', s=10, color='purple')
ax2.plot(full_period_prop.y[0], full_period_prop.y[1], label='Full Period Arc')
plt.quiver(full_period_prop.y[0,100], full_period_prop.y[1,100],
    full_period_prop.y[2,100], full_period_prop.y[3,100],
    color='red', label='Direction', zorder=2.5, width=0.005
)
ax2.axhline(y=0, color='r', linestyle='--', linewidth=0.5)
ax2.set_title(f'Full Period Propagation near L1 with $\\xi$={xi}, $\\eta$={eta}\n$x_0$={sv0[0]:.3f},$y_0$={sv0[1]:.3f},$v_{{x0}}$={sv0[2]:.3f},$v_{{y0}}$={sv0[3]:.3f} ({ps}, Lillian Shido)')
ax2.legend(loc='upper right',fontsize=8)
plt.savefig(f'Full_period_{ps}.png', dpi=300, bbox_inches='tight')

# pdb.set_trace()

df_check_return = pd.DataFrame({
    "location":['Start','End'],
    "time":[0,2*tf],
    "jacobi":[calc_Jacobi(mu, full_period_prop.y[0,0], full_period_prop.y[1,0], full_period_prop.y[2,0],full_period_prop.y[3,0]), calc_Jacobi(mu, full_period_prop.y[0,-1],full_period_prop.y[1,-1],full_period_prop.y[2,-1],full_period_prop.y[3,-1])],
    "x":[full_period_prop.y[0,0], full_period_prop.y[0,-1]],
    "y":[full_period_prop.y[1,0], full_period_prop.y[1,-1]],
    "xdot":[full_period_prop.y[2,0], full_period_prop.y[2,-1]],
    "ydot":[full_period_prop.y[3,0], full_period_prop.y[3,-1]],
})

df_state_error = pd.DataFrame({
    "jacobi":[calc_Jacobi(mu, full_period_prop.y[0,-1], full_period_prop.y[1,-1], full_period_prop.y[2,-1], full_period_prop.y[3,-1]) - calc_Jacobi(mu, full_period_prop.y[0,0],full_period_prop.y[1,0], full_period_prop.y[2,0], full_period_prop.y[3,0])],
    "x_error":[full_period_prop.y[0,-1]-full_period_prop.y[0,0]],
    "y_error":[full_period_prop.y[1,-1]- full_period_prop.y[1,0]],
    "xdot_error":[full_period_prop.y[2,-1]- full_period_prop.y[2,0]],
    "ydot_error":[full_period_prop.y[3,-1]- full_period_prop.y[3,0]],
    "x_error_dim":[(full_period_prop.y[0,-1]-full_period_prop.y[0,0])*l_char],
    "y_error_dim":[(full_period_prop.y[1,-1]- full_period_prop.y[1,0])*l_char],
    "xdot_error_dim":[(full_period_prop.y[2,-1]- full_period_prop.y[2,0])*l_char/t_char],
    "ydot_error_dim":[(full_period_prop.y[3,-1]- full_period_prop.y[3,0])*l_char/t_char]
})

# Configure table for question 7
    # columns = ["Case", "Iteration", "tf_dim", "delta_x_dim", "delta_y_dim", "delta_v0_x_dim", "delta_v0_y_dim", "delta_v0_mag_dim"]
iteration_per_case_table = (
    GT(df_iterations_per_case)
    .tab_header(
        title=md(f"States at each Iteration - 1D Periodic Targeter<br>({ps}, Lillian Shido)")
    )
    .tab_stub(rowname_col="Iteration")
    .tab_stubhead(label="Iteration") 
    .tab_spanner(
        label="Position", 
        columns=["x_dim","y_dim"]
    )
    .tab_spanner(
        label="Velocity",
        columns=["v_x_dim", "v_y_dim"]
    )
    .cols_label(
        Iteration="Iteration",
        x_dim="{{x}}<br>[km]",
        y_dim="{{y}}<br>[km]",
        v_x_dim="{{v_x}}<br>[km/s]",
        v_y_dim="{{v_y}}<br>[km/s]"
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
        columns=["v_x_dim", "v_y_dim"],
        decimals=5
    )
    .fmt_number(
        columns=["x_dim","y_dim"],
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

# Configure table for question 1 
    # columns = "half_period", "initial_x", "initial_y", "initial_xdot", "initial_ydot",
    # "half_period_dim", "initial_x_dim", "initial_y_dim", "initial_xdot_dim", "initial_ydot_dim"
state_error_table = (
    GT(df_state_error)
    .tab_header(
        title=md(f"State Error<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        jacobi="Jacobi Constant",
        x_error="{{x}}<br>[non-dim]",
        y_error="{{y}}<br>[non-dim]",
        xdot_error="{{v_x}}<br>[non-dim]",
        ydot_error="{{v_y}}<br>[non-dim]",
        x_error_dim="{{x}}<br>[km]",
        y_error_dim="{{y}}<br>[km]",
        xdot_error_dim="{{v_x}}<br>[km/s]",
        ydot_error_dim="{{v_y}}<br>[km/s]"
    )
    .fmt_scientific(
        columns=["jacobi", "x_error","y_error","xdot_error","ydot_error","x_error_dim","y_error_dim","xdot_error_dim","ydot_error_dim"],
        decimals=4
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
state_error_table.show()

check_return_table = (
    GT(df_check_return)
    .tab_header(
        title=md(f"Start and End States for Full Period Propagation<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        location="Location",
        jacobi="Jacobi Constant",
        time="{{Time}}<br>[non-dim]",
        x="{{x}}<br>[non-dim]",
        y="{{y}}<br>[non-dim]",
        xdot="{{v_x}}<br>[non-dim]",
        ydot="{{v_y}}<br>[non-dim]"
    )
    .fmt_number(
        columns=["time","jacobi", "x","y","xdot","ydot"],
        decimals=5
    )
    .fmt_number(
        columns=[],
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
check_return_table.show()

# # Configure table for question 1 
#     # columns = "half_period", "initial_x", "initial_y", "initial_xdot", "initial_ydot",
#     # "half_period_dim", "initial_x_dim", "initial_y_dim", "initial_xdot_dim", "initial_ydot_dim"
# final_initial_states_table = (
#     GT(df_final_initial_states)
#     .tab_header(
#         title=md(f"Initial State Vector<br>({ps}, Lillian Shido)")
#     )
#     .tab_spanner(
#         label="Non-Dimensional", 
#         columns=["half_period", "initial_x", "initial_y", "initial_xdot", "initial_ydot"]
#     )
#     .tab_spanner(
#         label="Dimensional",
#         columns=["half_period_dim", "initial_x_dim", "initial_y_dim", "initial_xdot_dim", "initial_ydot_dim"]
#     )
#     .cols_label(
#         half_period="{{Half-Period}}<br>[non-dim]",
#         initial_x="{{x_0}}<br>[non-dim]",
#         initial_y="{{y_0}}<br>[non-dim]",
#         initial_xdot="{{v_x0}}<br>[non-dim]",
#         initial_ydot="{{v_y0}}<br>[non-dim]",
#         half_period_dim="{{Half-Period}}<br>[days]",
#         initial_x_dim="{{x_0}}<br>[km]",
#         initial_y_dim="{{y_0}}<br>[km]",
#         initial_xdot_dim="{{v_x0}}<br>[km/s]",
#         initial_ydot_dim="{{v_y0}}<br>[km/s]"
#     )
#     .fmt_number(
#         columns=["half_period", "initial_x", "initial_y", "initial_xdot", "initial_ydot", "half_period_dim", "initial_x_dim", "initial_y_dim", "initial_xdot_dim", "initial_ydot_dim"],
#         decimals=4
#     )
#     .fmt_number(
#         columns=[],
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
# final_initial_states_table.show()

# Configure Tables
# Configure table for question 1 
# columns = "xi_from_L1", "eta_from_L1", "xi_from_L1_dim", "eta_from_L1_dim"
distance_from_L1_table = (
    GT(df_distance_from_L1)
    .tab_header(
        title=md(f"Distance from L1 Libration Point<br>({ps}, Lillian Shido)")
    )
    .tab_spanner(
        label="Non-Dimensional", 
        columns=["xi_from_L1", "eta_from_L1"]
    )
    .tab_spanner(
        label="Dimensional",
        columns=["xi_from_L1_dim", "eta_from_L1_dim"]
    )
    .cols_label(
        xi_from_L1="{{:xi:}}<br>[non-dim]",
        eta_from_L1="{{:eta:}}<br>[non-dim]",
        xi_from_L1_dim="{{:xi:}}<br>[km]",
        eta_from_L1_dim="{{:eta:}}<br>[km]"
    )
    .fmt_number(
        columns=["xi_from_L1", "eta_from_L1"],
        decimals=2
    )
    .fmt_number(
        columns=["xi_from_L1_dim", "eta_from_L1_dim"],
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
distance_from_L1_table.show()

# pdb.set_trace()