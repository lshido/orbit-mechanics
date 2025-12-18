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

qpo = pd.DataFrame({
    'x0':[-0.550000],
    'x':[-0.550000],
    'y':[0.000000],
    'vx':[0.000000],
    'vy':[4.093273],
})
IC_qpo = [
    qpo['x'][0], qpo['y'][0], qpo['vx'][0], qpo['vy'][0],
    1,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
]
period_qpo = 3*2*pi
qpo_prop = solve_ivp(planar_ode, [0, period_qpo], IC_qpo, args=(mu,), rtol=1e-12,atol=1e-14)
df_qpo = pd.DataFrame({
    'name': 'QPO @ x=-0.55',
    't':qpo_prop.t,
    'x':qpo_prop.y[0],
    'y':qpo_prop.y[1],
})

# qpo plot
x_min = 0.8
x_max = 1.20
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2
qpo_chart = alt.Chart(df_qpo).mark_line(strokeWidth=1,clip=True).encode(
    # x=alt.X('x:Q', scale=alt.Scale(domain=[-0.9,0]), axis=alt.Axis(title='x [non-dim]')),
    x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    # y=alt.Y('vx:Q', scale=alt.Scale(domain=[-5,5]), axis=alt.Axis(title='vx [non-dim]')),
    y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='sinebow')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=["Quasi-Periodic Orbit at x=-0.55 [nd]",f"({ps}, Lillian Shido)"]
)

qpo_chart_layer = alt.layer(qpo_chart)
qpo_chart_layer.save(f'{ps}_qpo.png', ppi=200)

pdb.set_trace()