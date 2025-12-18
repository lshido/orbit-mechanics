ps = "H2"
# Author: Lillian Shido
# Date: 12/16/2025

import time
import pdb
import numpy as np
import pandas as pd
from math import pi
from scipy.integrate import solve_ivp
import altair as alt
import plotly.express as px
import plotly.graph_objects as go

from methods import calc_velocity_from_Jacobi, planar_ode, system_properties, calc_L1, calc_L2
from constants import mu_Earth, mu_Moon, a_Moon


mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu)
x_L2, y_L2 = calc_L2(mu)
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[0],'z':[0],'vx':[0]})
L2 = pd.DataFrame({'name':["L2"],'x':[x_L2],'y':[0],'z':[0],'vx':[0]})
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0],'z':[0],'vx':[0]})
earth = pd.DataFrame({'name':["Earth"],'x':[x_Earth],'y':[0],'z':[0],'vx':[0]})

# x (will be chosen)
y = 0 # Fixed
vx = 0 # Perpindicular crossing
# vy (will be calculated from JC)
period = 50*2*pi

def crossxEvent(t, sv, mu):
    return sv[1] # Return value of y
crossxEvent.terminal = 300
crossxEvent.direction = 1
events_list = [crossxEvent]
df_map = pd.DataFrame()

JC_list = [
    3.5, 3.2, 3.1883, 3.18, 3.1722, 3.0121, 2.988, 2.8, 2.5
]
total_start_time = time.perf_counter()

for JC in JC_list:
    print(f"Building map for JC={JC}")
    if JC == 3.1883:
        label = f"JC={JC}\n(L1 Gateway Opens)"
    if JC == 3.1722:
        label = f"JC={JC}\n(L2 Gateway Opens)"
    if JC == 3.0121:
        label = f"JC={JC}\n(L3 Gateway Opens)"
    if JC == 2.988:
        label = f"JC={JC}\n(L4 and L5 Gateways Open)"
    elif JC not in [3.1883, 3.1722, 3.0121, 2.988]:
        label = f"JC={JC}"
    for x in [-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.1,1.20]:
        start_time = time.perf_counter()
        print(f"\tCalculating crossings for x0={x}")
        try:
            vy = calc_velocity_from_Jacobi(JC, mu, x, y, 0)
        except ZeroDivisionError:
            print(f"Zero Division Error! Skipping x={x}") 
            continue
        except ValueError:
            print(f"Value Error! Skipping x={x}") 
            continue
        IC = [
            x, y, vx, vy,
            1,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,
            0,0,1,0,
            0,0,0,1
        ]
        prop = solve_ivp(planar_ode, [0, period], IC, events=events_list, args=(mu,), rtol=1e-12,atol=1e-14)
        map_data = pd.DataFrame({
            "label": label,
            "JC":JC,
            "x0": x,
            "number": np.arange(len(prop.t_events[0])),
            "time": prop.t_events[0],
            "x": prop.y_events[0][:,0],
            "y": prop.y_events[0][:,1],
            "vx": prop.y_events[0][:,2],
            "vy": prop.y_events[0][:,3],
        })
        df_map = pd.concat([df_map, map_data], ignore_index=True)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print(f"\tComplete! Elapsed time: {elapsed_time:.4f} seconds")

earth_x_y = alt.Chart(earth).mark_point(shape="diamond", strokeWidth=0.5, filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='vx:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Earth'], range=['darkblue'])).title(None)
)

L1_x_y = alt.Chart(L1).mark_point(shape="diamond", strokeWidth=0.5, filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='vx:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L1'], range=['darkgreen'])).title(None)
)

moon_x_y = alt.Chart(moon).mark_point(shape="diamond", strokeWidth=0.5, filled=True,size=50,clip=True).encode(
    x='x:Q',
    y='vx:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Moon'], range=['gray'])).title(None)
)

L2_x_y = alt.Chart(L2).mark_point(shape="diamond", strokeWidth=0.5, filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='vx:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L2'], range=['darkgreen'])).title(None)
)

# vy-y map plot
map_chart = alt.Chart(df_map).mark_point(filled=True,size=5,clip=True,color='darkblue').encode(
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    x=alt.X('x:Q', scale=alt.Scale(domain=[-0.8,1.2]), axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('vx:Q', scale=alt.Scale(domain=[-4,4]), axis=alt.Axis(title='vx [non-dim]')),
).properties(
    width=200,
    height=200,
    title=[f"Poincaré Map for Earth-Moon for various JC",f"({ps}, Lillian Shido)"]
).facet(
    facet=alt.Facet('label:N', title='Poincaré Maps in Earth-Moon System for varying Jacobi Constant', sort='descending'),
    columns=3
)

# map_chart_layer = alt.layer(map_chart, earth_x_y, L1_x_y, moon_x_y, L2_x_y)
map_chart.save(f'{ps}_map_JC.png', ppi=200)

total_end_time = time.perf_counter()
total_elapsed_time = total_end_time - total_start_time
print(f"Total Elapsed time: {total_elapsed_time/60:.4f} minutes")


pdb.set_trace()