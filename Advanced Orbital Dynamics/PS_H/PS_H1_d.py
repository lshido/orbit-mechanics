ps = "H1 part d"
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

from methods import calc_velocity_from_Jacobi, planar_ode, system_properties
from constants import mu_Earth, mu_Moon, a_Moon


mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
JC = 3.18
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

for x in [-0.2,-0.25,-0.3,-0.4,-0.45,-0.49,-0.51,-0.55,-0.6,-0.63,-0.67,-0.7,-0.75,-0.8]:
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
        "x0": x,
        "number": np.arange(len(prop.t_events[0])),
        "time": prop.t_events[0],
        "x": prop.y_events[0][:,0],
        "y": prop.y_events[0][:,1],
        "vx": prop.y_events[0][:,2],
        "vy": prop.y_events[0][:,3],
    })
    df_map = pd.concat([df_map, map_data], ignore_index=True)

# vy-y map plot
map_chart = alt.Chart(df_map).mark_point(filled=True,size=10,clip=True).encode(
    x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', scale=alt.Scale(domain=[-0.9,0]), axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('vx:Q', scale=alt.Scale(domain=[-5,5]), axis=alt.Axis(title='vx [non-dim]')),
    color=alt.Color('x0:N', scale=alt.Scale(scheme='sinebow')).title(None)
).properties(
    width=200,
    height=400,
    title=[f"Poincaré Map for JC={JC}",f"({ps}, Lillian Shido)"]
)

map_line = alt.Chart(df_map).mark_line(clip=True,strokeWidth=0.5,interpolate="cardinal").encode(
    x='x:Q',
    y='vx:Q',
    color=alt.Color('x0:N').title(None)
)

map_chart_layer = alt.layer(map_chart)
map_chart_layer.save(f'{ps}_map.png', ppi=200)

# 3-D chart
fig = px.scatter(df_map, x="x", y='vx', color='x0',
                 labels={
                     "x": "x [non-dim]",
                     "y": "vx [non-dim]"
                 })
fig.update_layout(
    title=dict(text=f"Poincaré Map JC={JC} ({ps}, Lillian Shido)", font=dict(size=30), automargin=True, yref='paper'),
    legend=dict(title=dict(text=None)),
    scene = dict(
        xaxis=dict(range=[-0.9, 0]),
        yaxis=dict(range=[-5,5]),
        aspectratio = dict(x=1,y=1,z=1),
        aspectmode='manual'  
    ),
)
fig.show()

pdb.set_trace()