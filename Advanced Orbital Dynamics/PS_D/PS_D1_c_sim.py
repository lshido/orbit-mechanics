
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

stm_tf = sol.y[4:20,-1].reshape(4,4)

mu,l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)

# def predict_w_STM(adjustment):
state = ['x','y','x_dot','y_dot']
df = pd.DataFrame(
    columns=['% Change','Changed State','Changed IC','Predicted State','delta[non-dim]','delta[dim]','final[non-dim]','final[dim]']
    )
for change in [0.01, 0.1]:
    # Change [x, vx, vy]
    for n in [0, 2, 3]:
        delta_IC = sv0[n]*change
        changed_IC = sv0[n] + delta_IC
        for m in [0, 1, 2, 3]:
            delta_f = stm_tf[m,n]*delta_IC 
            # Calc the final values in x,y,vx,vy
            final = sv0[m] + delta_f
            if m<2:
                dim_delta_f = delta_f*l_char
                dim_final = final*l_char
            else:
                dim_delta_f = delta_f*l_char/t_char
                dim_final = final*l_char/t_char;
            new_data = pd.DataFrame({
                "% Change":[change*100],
                "Changed State":[state[n]],
                "Changed IC":[changed_IC],
                "Predicted State":[state[m]],
                "delta[non-dim]":[delta_f],
                "delta[dim]":[dim_delta_f],
                "final[non-dim]":[final],
                "final[dim]":[dim_final]
            })
            df = pd.concat([df, new_data], ignore_index=True)

case_1_data = df.loc[df['% Change']==1.0].loc[df['Changed State']=='x']
case_2a_data = df.loc[df['% Change']==1.0].loc[df['Changed State']=='y_dot']
case_2b_data = df.loc[df['% Change']==10.0].loc[df['Changed State']=='y_dot']
case_3a_data = df.loc[df['% Change']==1.0].loc[df['Changed State']=='x_dot']
case_3b_data = df.loc[df['% Change']==10.0].loc[df['Changed State']=='x_dot']

def build_new_ICs(case_data, original_IC, states):
    IC = copy.deepcopy(original_IC)
    state_to_change = states.index(case_data['Changed State'].to_list()[0])
    update_value = case_data['Changed IC'].to_list()[0]
    IC[state_to_change] = update_value
    return IC

case_1_IC = build_new_ICs(case_1_data, original_IC, state)
case_2a_IC = build_new_ICs(case_2a_data, original_IC, state)
case_2b_IC = build_new_ICs(case_2b_data, original_IC, state)
case_3a_IC = build_new_ICs(case_3a_data, original_IC, state)
case_3b_IC = build_new_ICs(case_3b_data, original_IC, state)

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

fig, ax = plt.subplots()
ax.plot(sol.y[0],sol.y[1], label='nominal')
df_actual = pd.DataFrame(
    columns=['IC','delta_x','delta_y','delta_xdot','delta_ydot']
    )
nominal_final = sol.y[0:4,-1]
for enum, IC in enumerate([case_1_IC, case_2a_IC, case_2b_IC, case_3a_IC, case_3b_IC]):
    actual = solve_ivp(eoms, [0, 0.47753], IC, args=(mu,), rtol=1e-12,atol=1e-14)
    ax.plot(actual.y[0],actual.y[1], label=f'actual_[{enum}]')
    final = actual.y[0:4,-1]
    delta_x = final[0] - nominal_final[0]
    delta_y = final[1] - nominal_final[1]
    delta_xdot = final[2] - nominal_final[2]
    delta_ydot = final[3] - nominal_final[3]
    delta_x_dim = delta_x*l_char
    delta_y_dim = delta_y*l_char
    delta_xdot_dim = delta_xdot*l_char/t_char
    delta_ydot_dim = delta_ydot*l_char/t_char
    new_actual_data = pd.DataFrame({
        'IC':[IC],
        'delta_x':[delta_x_dim],
        'delta_y':[delta_y_dim],
        'delta_xdot':[delta_xdot_dim],
        'delta_ydot':[delta_ydot_dim]
    })
    df_actual = pd.concat([df_actual,new_actual_data],ignore_index=True)

ax.scatter(x_Earth,0)
ax.scatter(x_Moon,0)
ax.legend()
plt.axis('square')
plt.xlim(-0.3,0.5)
plt.ylim(-0.3,0.5)
plt.show()


pdb.set_trace()
