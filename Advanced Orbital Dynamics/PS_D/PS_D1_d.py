
import pdb
import copy
from constants import mu_Earth, mu_Moon, a_Moon
import numpy as np
import pandas as pd
from math import pi, sqrt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from great_tables import GT, md, html, style, loc


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

def combined_ode(t,sv,mu):
    
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

def crossxEvent(t, sv, mu):
    return sv[1] # Return value of y
crossxEvent.terminal = True

def eval_acceleration(x,y,xdot,ydot):
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

# Trjaectory Initial Conditions
original_IC = [
    0.488, 0.200, # x and y
    -0.880, 0.200 # xdot and ydot
    ]

# ODE Initial Conditions
sv0 = [
    original_IC[0], # x
    original_IC[1], # y
    original_IC[2], # xdot
    original_IC[3], # ydot
    1,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
    ]

# Set the span of the integrator
t_final = 1.5*pi;
tspan = [0, t_final];

sol = solve_ivp(combined_ode, tspan, sv0, events=crossxEvent, args=(mu,), rtol=1e-12,atol=1e-14)

stm_tf = sol.y[4:20,-1].reshape(4,4)
sv_tf = sol.y[0:4,-1]
mu,l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)

original_time = sol.t_events[0][0]
# "Predict" the states given the change in tf
states_crossX = sol.y[0:4,-1] # The states at tf in the original problem
df_predicted_states = pd.DataFrame(
    columns = ['% change','delta_t','t_final','x_predict','y_predict','xdot_predict','ydot_predict',\
    'x_predict_dim','y_predict_dim','xdot_predict_dim','ydot_predict_dim']
) 
for change in [0.01, 0.1]:
    t_span = original_time*change
    x_f,y_f,xdot_f,ydot_f = states_crossX[0], states_crossX[1], states_crossX[2], states_crossX[3]
    a_x, a_y = eval_acceleration(x_f,y_f,xdot_f,ydot_f)
    x_new, y_new, xdot_new, ydot_new = eval_kinematic(x_f,y_f,xdot_f,ydot_f,a_x,a_y,t_span)
    predicted_data = pd.DataFrame({
        '% change':[change*100],
        'delta_t':[t_span],
        't_final':[original_time + t_span],
        'x_predict':[x_new],
        'y_predict':[y_new],
        'xdot_predict':[xdot_new],
        'ydot_predict':[ydot_new],
        'x_predict_dim':[x_new*l_char],
        'y_predict_dim':[y_new*l_char],
        'xdot_predict_dim':[xdot_new*l_char/t_char],
        'ydot_predict_dim':[ydot_new*l_char/t_char]
    })
    df_predicted_states = pd.concat([df_predicted_states,predicted_data],ignore_index=True)

df_changed_time = pd.DataFrame(
    columns=['% change','delta_t','t_final','x','y','xdot','ydot','x_dim','y_dim','xdot_dim','ydot_dim']
)
for change in [0.01, 0.1]:
    new_time = original_time + original_time*change
    actual = solve_ivp(eoms, [0, new_time], original_IC, args=(mu,), rtol=1e-12,atol=1e-14)
    x = actual.y[0,-1]
    y = actual.y[1,-1]
    xdot = actual.y[2,-1]
    ydot = actual.y[3,-1]
    new_data = pd.DataFrame({
        '% change':[change*100],
        'delta_t':[original_time*change],
        't_final':[new_time],
        'x':[x],
        'y':[y],
        'xdot':[xdot],
        'ydot':[ydot],
        'x_dim':[x*l_char],
        'y_dim':[y*l_char],
        'xdot_dim':[xdot*l_char/t_char],
        'ydot_dim':[ydot*l_char/t_char]
    })
    df_changed_time = pd.concat([df_changed_time, new_data], ignore_index=True)

df_time_error = pd.DataFrame(
    columns=['% change','delta_t','t_final','error_x','error_y','error_xdot','error_ydot']
)
for row in [0,1]:
    error_x = calc_error(df_predicted_states['x_predict'][row],df_changed_time['x'][row])*100
    error_y = calc_error(df_predicted_states['y_predict'][row],df_changed_time['y'][row])*100
    error_xdot = calc_error(df_predicted_states['xdot_predict'][row],df_changed_time['xdot'][row])*100
    error_ydot = calc_error(df_predicted_states['ydot_predict'][row],df_changed_time['ydot'][row])*100
    error_data = pd.DataFrame({
        '% change':[df_predicted_states['% change'][row]],
        'delta_t':[df_predicted_states['delta_t'][row]],
        't_final':[df_predicted_states['t_final'][row]],
        'error_x':[error_x],
        'error_y':[error_y],
        'error_xdot':[error_xdot],
        'error_ydot':[error_ydot]
    })
    df_time_error = pd.concat([df_time_error,error_data],ignore_index=True)

# Build dim comparison table
df_compare_dim = pd.DataFrame(
    columns=['% change', 'delta_t', 't_final', 'x_predict_dim', 'y_predict_dim', 'xdot_predict_dim', 'ydot_predict_dim',\
    'x_dim','y_dim','xdot_dim','ydot_dim']
)
df_compare_dim['% change'] = df_predicted_states['% change']
df_compare_dim['delta_t'] = df_predicted_states['delta_t']
df_compare_dim['t_final'] = df_predicted_states['t_final']
df_compare_dim['x_predict_dim'] = df_predicted_states['x_predict_dim']
df_compare_dim['y_predict_dim'] = df_predicted_states['y_predict_dim']
df_compare_dim['xdot_predict_dim'] = df_predicted_states['xdot_predict_dim']
df_compare_dim['ydot_predict_dim'] = df_predicted_states['ydot_predict_dim']
df_compare_dim['x_dim'] = df_changed_time['x_dim']
df_compare_dim['y_dim'] = df_changed_time['y_dim']
df_compare_dim['xdot_dim'] = df_changed_time['xdot_dim']
df_compare_dim['ydot_dim'] = df_changed_time['ydot_dim']

# Build final state table
state_table

# Configure the table showing the predicted values
compare_table = (
    GT(df_compare_dim)
    .tab_header(
        title=md("Predicted vs Numerical Results (Dimensional)")
    )
    .tab_stub(rowname_col="% change")
    .tab_stubhead(label="% change")
    .tab_spanner(
        label="Predicted",
        columns=["x_predict_dim", "y_predict_dim", "xdot_predict_dim", "ydot_predict_dim"]
    )
    .tab_spanner(
        label="Numerical",
        columns=["x_dim","y_dim","xdot_dim","ydot_dim"]
    )
    .tab_spanner(
        label="Time",
        columns=["% change", "delta_t", "t_final"]
    )
    .cols_label(
        x_predict_dim=html("x<br>[km]"),
        y_predict_dim=html("y<br>[km]"),
        xdot_predict_dim=html("x_dot<br>[km/s]"),
        ydot_predict_dim=html("y_dot<br>[km/s]"),
        x_dim=html("x<br>[km]"),
        y_dim=html("y<br>[km]"),
        xdot_dim=html("x_dot<br>[km/s]"),
        ydot_dim=html("y_dot<br>[km/s]"),
    )
    .tab_style(
        style=style.borders(
            sides="right",
            color="lightgray",
            style="solid",
            weight="1px"
        ),
        locations=loc.body(columns=[2, 4, 6, 8])
    )
    .fmt_number(
        columns=["delta_t", "t_final", "xdot_dim", "ydot_dim", \
        "xdot_predict_dim", "ydot_predict_dim"],
        decimals=5
    )
    .fmt_number(
        columns=["x_dim", "y_dim","x_predict_dim", "y_predict_dim"],
        decimals=3
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
)
compare_table.show()

# Configure Error table
error_table = (
    GT(df_time_error)
    .tab_header(
        title=md("Percentage Error: Predicted vs Numerical Results")
    )
    .tab_stub(rowname_col="% change")
    .tab_stubhead(label="% change")
    .tab_spanner(
        label="Position",
        columns=["error_x", "error_y"]
    )
    .tab_spanner(
        label="Velocity",
        columns=["error_xdot", "error_ydot"]
    )
    .tab_spanner(
        label="Time",
        columns=["% change", "delta_t", "t_final"]
    )
    .cols_label(
        error_x=html("x<br>[%]"),
        error_y=html("y<br>[%]"),
        error_xdot=html("x_dot<br>[%]"),
        error_ydot=html("y_dot<br>[%]"),
    )
    .opt_horizontal_padding(scale=3)
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
        columns=["delta_t", "t_final"],
        decimals=5
    )
    .fmt_number(
        columns=["error_x", "error_y", "error_xdot", "error_ydot"],
        decimals=4
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
)
error_table.show()

pdb.set_trace()
