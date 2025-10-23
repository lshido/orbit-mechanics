# Problem D2 part(a)
# Code for Numerical Integrator
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
sv0 = [
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

sol = solve_ivp(ode, tspan, sv0, events=crossxEvent, args=(mu,), rtol=1e-12,atol=1e-14)
tf = sol.t_events[0][0]
stm = sol.y[4:20,-1].reshape(4,4) # turn into 4x4 phi matrix
phi_rv = np.array([
    [stm[0,2], stm[0,3]],
    [stm[1,2], stm[1,3]]
])

# Set initial conditions for targeter reference arc
r0_IC = np.array([
    [original_IC[0]],
    [original_IC[1]]
])
v0_IC = np.array([
    [original_IC[2]],
    [original_IC[3]]
])
# Run the targeter for a set of targets
rf_target_list = [
    np.array([[-0.3], [0.05]]),
    np.array([[-0.1],[0]]),
    np.array([[-0.35],[-.1]])
]
tolerance = 1e-8 # Set tolerance
for rf_target in rf_target_list:
    r0 = copy.deepcopy(r0_IC)
    v0 = copy.deepcopy(v0_IC)
    counter = 0
    while True:
        counter = counter + 1
        # Step 1: Using ICs, propagate EOM+STM until tf, get STM from tf and rf
        prop = solve_ivp(ode, [0, tf], sv0, args=(mu,), rtol=1e-13,atol=1e-14)
        rf = prop.y[0:2,-1].reshape(2,1) # position at tf, turn into 2x1 vector
        vf = prop.y[2:4,-1].reshape(2,1) # velocity at tf, turn into 2x1 vector
        # Step 2: Compare rf with rf_target
        error = rf_target - rf
        # Check if the error is within acceptable margins
        if (abs(error) > tolerance).all(): # If not, recalculate the delta_v0 and try again
            # Step 3: Calc new delta_v0
            delta_v0 = np.linalg.inv(phi_rv) @ error # multiply the matrices together (dot product)
            print(f"delta_v0: {delta_v0[0,0]}, {delta_v0[1,0]}, Iteration: {counter}")
            # Step 4: Use new delta_v0 to calc new ICs
            v0 = v0_IC + delta_v0
            r0 = r0 # Initial position is constant; we only vary velocity in this problem
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
            continue
        else: # If error is within acceptable margins, break out of iterative loop
            delta_vf = vf - v0_IC
            print(f"Position: {rf}\nVelocity: {vf}\n Iterations: {counter}")
            break
    continue

pdb.set_trace()