ps = "D4 part b iii"
# Problem D4 Part b iii
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
from scipy import stats
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

def find_halfperiod(starting_x, ydot_guess, mu, tolerance=1e-12, max_iterations=50):
    """
    Calculates the resulting half period states for a perpindicular x-axis crossing
    
    args:
        starting_x: starting x location that you want to propagate half period for 
        ydot_guess: the guessed ydot to start the propagation
        mu: the mu value for the system
        tolerance: convergence tolerance for the targeter, default: 1e-12
        max_iterations: maximum number of iterations before terminating, default=50

    returns:
        iterations: number of iterations targeter took to converge. 50 is the maximum
        tf: half period
        xf: x at arrival
        yf: y at arrival
        vxf: xdot at arrival
        vyf: ydot at arrival
    """
    x0 = starting_x # Initial x guess 
    vy0 = ydot_guess
    
    counter = 0
    while True:
        IC_guess = [
            x0, 0, 0, vy0, # Initial states
            1,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,
            0,0,1,0,
            0,0,0,1
        ]
        counter = counter + 1
        # Step 1: Using ICs, propagate EOM+STM until second x-axis cross, get STM, get xdot, ydot
        prop = solve_ivp(ode, [0, 50*pi], IC_guess, events=crossxEvent, args=(mu,), rtol=1e-12,atol=1e-14)
        stm = prop.y[4:20,-1].reshape(4,4) # turn into 4x4 phi matrix
        tf = prop.t_events[0][1]
        phi_34 = stm[2,3]
        phi_24 = stm[1,3]
        rf = prop.y[0:2,-1].reshape(2,1) # position at tf, turn into 2x1 vector
        vf = prop.y[2:4,-1].reshape(2,1) # velocity at tf, turn into 2x1 vector
        vxf = prop.y[2,-1] # get xdot to compare against tolerance
        vyf = prop.y[3,-1] # get ydot to calculate the acceleration
        # print for info
        # print(f"{counter}, tf: {tf}, x:{IC_guess[0]}, y:{IC_guess[1]}, vx:{IC_guess[2]}, vy0:{IC_guess[3]}")
        # Check if the vxf is close to 0 within acceptable margins
        if abs(vxf) > tolerance and counter <= max_iterations: # If not, recalculate the delta_v0 and try again
            # Calc the acceleration in x
            a_x, a_y = eval_acceleration(rf[0,0], rf[1,0], vxf, vyf, mu)
            # Step 3: Calc new delta_vy0
            delta_vy0 = -vxf / (phi_34 - phi_24*(a_x/vyf))
            # Step 4: Use new deltas_ydot_t to calc new ICs
            if delta_vy0 > 0.1:
                delta_vy0 = delta_vy0/3
            vy0 = vy0 + delta_vy0
            continue
        else: # If error is within acceptable margins, break out of iterative loop
            arrival_states = [rf[0,0], rf[1,0], vf[0,0], vf[1,0]]
            converged_initial_states = [IC_guess[0], IC_guess[1], IC_guess[2], IC_guess[3]]            
            break
    return counter, tf, arrival_states, converged_initial_states

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

# Find arrival states for starting state using xi as starting x
iterations, tf, arrival_states, _ = find_halfperiod(starting_x, ydot_guess, mu, tolerance=1e-12, max_iterations=50)
print(f"First guess:\nIterations:{iterations},half-period:{tf}, x:{arrival_states[0]}, y:{arrival_states[1]}, vx:{arrival_states[2]}, vy:{arrival_states[3]}")
# pdb.set_trace()
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
first_delta_x = -0.001
delta_x = -0.008 # step in x

orbit_x = arrival_states[0]
orbit_ydot = arrival_states[3]
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
            "xi":[arrival_states[0]-x_L1],
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
            "xi":[arrival_states[0]-x_L1],
            "x":[converged_initial_states[0]],
            "y":[converged_initial_states[1]],
            "vx":[converged_initial_states[2]],
            "vy":[converged_initial_states[3]],
            "jacobi":[jacobi]
        })
        df_orbits = pd.concat([df_orbits, orbit_IC_data], ignore_index=True)

# Table iterations each orbit takes to converge
df_orbit_iterations = df_orbits[['orbit','xi','iterations']]
# Calculate the ydot vs x slope of the previous set
# result = stats.linregress(df_orbits['x'].astype(float), df_orbits['vy'].astype(float))
# print(f"linregress slope: {result.slope}")
# # Check slope calc
# my_slope = (df_orbits.iloc[-1]['vy']-df_orbits.iloc[-2]['vy'])/(df_orbits.iloc[-1]['x']-df_orbits.iloc[-2]['x'])
# # print(f"my_slope: {my_slope}")

# # Plot vy0 as a function of x0
# fig2, ax2 = plt.subplots()
# ax2.plot(df_orbits['x'],df_orbits['vy'])
# plt.grid()
# plt.xlabel(r"${{x_0}}$""\n[non-dim]")
# plt.ylabel(r"$\dot{{y_0}}$""\n[non-dim]")
# ax2.set_title(rf'$\dot{{y_0}}$ as a function of ${{x_0}}$'f'\n({ps}, Lillian Shido)')
# plt.savefig(f'ydot0_wrt_x0_{ps}.png', dpi=300, bbox_inches='tight')

# Set up colors for plot and their labels
colors = load_cmap(name="Classic_Cyclic", cmap_type='discrete').colors
labels = [fr'$\xi$={i}' if i!=0 else f'Baseline' for i, c in enumerate(colors)]

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
    full_period_prop = solve_ivp(ode, [0, 2*tf], IC, args=(mu,), rtol=1e-12,atol=1e-14)
    arc = [np.column_stack([full_period_prop.y[0], full_period_prop.y[1]])]
    try:
        line_collection = LineCollection(arc, colors=colors[enum], linewidth=0.5, label=fr"$\xi$={xi0:.3f}")
    except:
        pdb.set_trace()
    ax1.add_collection(line_collection)
    # plt.quiver(full_period_prop.y[0,100], full_period_prop.y[1,100],
    # full_period_prop.y[2,100], full_period_prop.y[3,100],
    # color='red', zorder=2.5, width=0.005
# )
ax1.set_aspect('equal', 'box')
ax1.autoscale()
ax1.legend(fontsize=6)
plt.grid()
plt.xlabel("x [non-dim]")
plt.ylabel("y [non-dim]")
ax1.tick_params(axis='both', which='major', labelsize=6)
ax1.set_title(f'Lyapunovs near L1 starting from $\\xi$={xi:.2f}, $\\eta$={eta}\nwith predicted $\dot{{y_0}}$, $\\Delta{{x}}$={delta_x:.3f} ({ps}, Lillian Shido)')
plt.savefig(f'Lyapunov_family_{ps}.png', dpi=300, bbox_inches='tight')

# Configure tables
# Configure table for question 7
#     columns = ['Orbit','xi','iterations']
iterations_table = (
    GT(df_orbit_iterations)
    .tab_header(
        title=md(f"Iterations to Produce Periodic Orbit<br>({ps}, Lillian Shido)")
    )
    .cols_label(
        orbit="Orbit",
        xi="{{:xi:}}",
        iterations="Iterations",
    )
    .fmt_number(
        columns=["xi"],
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
iterations_table.show()
