ps = "D4 part b"
# Problem D4 Part b
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

def calc_Jacobi(mu, x, y):
    d = sqrt((x+mu)**2 + y**2)
    r = sqrt((x-1+mu)**2 + y**2)
    C = x**2 + y**2 + (2*(1-mu)/d) + (2*mu/r)
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

# Initialize the dataframe to hold the set of initial conditions for each orbit we find
df_initial_conditions = pd.DataFrame(
    columns=['Orbit','Iterations','x0','y0','vx0','vy0']
)

r0 = copy.deepcopy(r0_IC)
v0 = copy.deepcopy(v0_IC)
IC_guess = copy.deepcopy(sv0_IC)
tf = copy.deepcopy(tf_IC)
x0 = r0[0,0]
vy0 = v0[1,0]
arc = []

# Configuration
delta_x = 0.001 # Step size for x
tolerance = 1e-12 # Set convergence tolerance
step_size = 0.01

# Find 10 orbits
for orbit in range(11):
    counter = 0
    while True:
        counter = counter + 1
        # Step 1: Using ICs, propagate EOM+STM until second x-axis cross, get STM, get xdot, ydot
        prop = solve_ivp(ode, tspan, IC_guess, events=crossxEvent, args=(mu,), rtol=1e-12,atol=1e-14)
        stm = prop.y[4:20,-1].reshape(4,4) # turn into 4x4 phi matrix
        tf = prop.t_events[0][1]
        phi_34 = stm[2,3]
        phi_24 = stm[1,3]
        rf = prop.y[0:2,-1].reshape(2,1) # position at tf, turn into 2x1 vector
        vf = prop.y[2:4,-1].reshape(2,1) # velocity at tf, turn into 2x1 vector
        vxf = prop.y[2,-1] # get xdot to compare against tolerance
        vyf = prop.y[3,-1] # get ydot to calculate the acceleration
        # Check if the vxf is close to 0 within acceptable margins
        if abs(vxf) > tolerance: # If not, recalculate the delta_v0 and try again
            # Calc the acceleration in x
            a_x, a_y = eval_acceleration(rf[0,0], rf[1,0], vxf, vyf, mu)
            # Step 3: Calc new delta_vy0
            delta_vy0 = -vxf / (phi_34 - phi_24*(a_x/vyf))
            # Step 4: Use new deltas_ydot_t to calc new ICs
            if delta_vy0 > step_size:
                delta_vy0 = delta_vy0/3
            vy0 = vy0 + delta_vy0
            # Rebuild the ICs for combined EOM+STM ODEs (IC_guess)
            IC_guess = [
                x0, # x
                r0[1,0], # y
                v0[0,0], # x_dot
                vy0, # y_dot
                1,0,0,0, # Identity matrix for phi ICs
                0,1,0,0,
                0,0,1,0,
                0,0,0,1
            ]
            print(f"{counter}, x:{x0}, y:{r0[1,0]}, vx:{v0[0,0]}, vy0:{vy0}, dvy0:{delta_vy0}")
            continue
        else: # If error is within acceptable margins, break out of iterative loop
            # Store these ICs so we can propagate them in bulk later
            orbit_initial_conditions = pd.DataFrame({
            'Orbit':[orbit],
            'Iterations':[counter],
            'half-period':[tf],
            'x0':[x0],
            'y0':[0],
            'vx0':[0],
            'vy0':[vy0]
            })
            df_initial_conditions = pd.concat([df_initial_conditions, orbit_initial_conditions], ignore_index=True)
            # Set the new x0 and vy0 values
            x0 = x0 + delta_x
            # Build new IC guess for next orbit
            IC_guess = [
                x0, # Stepped x
                0,  # y=0, collinear point
                0,  # vx=0, perpindicular cross
                vy0, # Use the same vy0 as before
                1,0,0,0, # Identity matrix for phi ICs
                0,1,0,0,
                0,0,1,0,
                0,0,0,1
            ]
            break
pdb.set_trace()

# Set up colors for plot and their labels
colors = load_cmap(name="Classic_Cyclic", cmap_type='discrete').colors
# Initialize family plot
fig1, ax1 = plt.subplots(figsize=(6.5, 6.5))
ax1.set_aspect('equal', 'box')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax1.scatter(x_Earth, 0, label='Earth', s=20, color='blue')
ax1.scatter(x_Moon, 0, label='Moon', s=20, color='gray')
plt.grid()
ax1.scatter(original_IC[0], original_IC[1], label='Start', s=20, marker="<", color='green')
ax1.scatter(x_L1, y_L1, label='L1', s=10, color='purple')
ax1.axhline(y=0, color='r', linestyle='--', linewidth=0.5)
for row in df_initial_conditions.iterrows():
    tf = row[1]['half-period']
    x0 = row[1]['x0']
    y0 = row[1]['y0']
    vx0 = row[1]['vx0']
    vy0 = row[1]['vy0']
    IC = [
        x0, y0, vx0, vy0, # IC states
        1,0,0,0, # Identity matrix for phi ICs
        0,1,0,0,
        0,0,1,0,
        0,0,0,1
    ]
    full_period_prop = solve_ivp(ode, [0, 2*tf], IC, args=(mu,), rtol=1e-12,atol=1e-14)
    arc = [np.column_stack([full_period_prop.y[0], full_period_prop.y[1]])]
    try:
        line_collection = LineCollection(arc, colors=colors[counter-1], linewidth=0.8)
    except:
        pdb.set_trace()
    ax1.add_collection(line_collection)
    # plt.quiver(full_period_prop.y[0,100], full_period_prop.y[1,100],
    # full_period_prop.y[2,100], full_period_prop.y[3,100],
    # color='red', zorder=2.5, width=0.005
)
ax1.set_title(f'Perpendicular Crossing near L1 with $\\xi$={xi}, $\\eta$={eta}\nIterations: {counter} ({ps}, Lillian Shido)')
ax1.legend(loc='lower left', fontsize=6)
plt.savefig(f'Lyapunov_family_{ps}.png', dpi=300, bbox_inches='tight')