ps = "H1 part a"
# Author: Lillian Shido
# Date: 12/16/2025

import pdb
import numpy as np
import pandas as pd
from math import pi
from scipy.integrate import solve_ivp
import altair as alt
import plotly.express as px
import plotly.graph_objects as go

from methods import calc_velocity_from_Jacobi, planar_ode

mu = 0.5
JC = 4.5
# x (will be chosen)
y = 0 # Fixed
vx = 0 # Perpindicular crossing
# vy (will be calculated from JC)
period = 50*2*pi

def crossxEvent(t, sv, mu):
    return sv[1] # Return value of y
crossxEvent.terminal = 100
crossxEvent.direction = 1
events_list = [crossxEvent]
df_map = pd.DataFrame()

range1 = np.linspace(-0.2, -0.4, 3)
range2 = np.linspace(-0.41, -0.5, 5)
range3 = np.linspace(-0.51, -0.7, 3)
for x in [-0.2,-0.25,-0.3,-0.4,-0.45,-0.49,-0.51,-0.55,-0.6,-0.63,-0.67,-0.7,-0.75,-0.8]:
    try:
        vy = calc_velocity_from_Jacobi(JC, mu, x, y, 0)
    except ZeroDivisionError:
        print(f"Zero Division Error! Skipping x={x}") 
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
        "x0": x,
        "number": np.arange(len(prop.t_events[0])),
        "time": prop.t_events[0],
        "x": prop.y_events[0][:,0],
        "y": prop.y_events[0][:,1],
        "vx": prop.y_events[0][:,2],
        "vy": prop.y_events[0][:,3],
    })
    df_map = pd.concat([df_map, map_data], ignore_index=True)

qpo = df_map[df_map['x0']==-0.55].iloc[0]
df_qpo = pd.DataFrame()
IC_po = [
    qpo['x'], qpo['y'], qpo['vx'], qpo['vy'],
    1,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
]
period_qpo = 2*2*pi
qpo_prop = solve_ivp(planar_ode, [0, period_qpo], IC, args=(mu,), rtol=1e-12,atol=1e-14)
df_periodic = pd.DataFrame({
    'name': 'Periodic Orbit @ x=-0.3',
    't':qpo_prop.t,
    'x':qpo_prop.y[0],
    'y':qpo_prop.y[1],
    'z':qpo_prop.y[2]
})

# vy-y map plot
map_chart = alt.Chart(df_map).mark_point(filled=True,size=10,clip=True).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[-0.9,0]), axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('vx:Q', scale=alt.Scale(domain=[-5,5]), axis=alt.Axis(title='vx [non-dim]')),
    color=alt.Color('x0:N', scale=alt.Scale(scheme='sinebow')).title(None)
).properties(
    width=200,
    height=400,
    title=["Hénon's Map",f"({ps}, Lillian Shido)"]
)

map_line = alt.Chart(df_map).mark_line(clip=True,strokeWidth=0.5,interpolate="cardinal").encode(
    x='x:Q',
    y='vx:Q',
    color=alt.Color('x0:N').title(None)
)

annotation = alt.Chart(df_map).mark_text(
    align='right',
    baseline='middle',
    fontSize = 8,
    dx = -2
).encode(
    x='x:Q',
    y='vx:Q',
    text='number:Q'
).transform_filter(
    (alt.datum.x0==-0.55) & (alt.datum.number < 16)
)

map_chart_layer = alt.layer(map_chart, annotation)
map_chart_layer.save(f'{ps}_map.png', ppi=200)

pdb.set_trace()