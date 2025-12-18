ps = "H1 part c"
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
from great_tables import GT, md, system_fonts,style, loc


from methods import planar_ode, calc_Jacobi

mu = 0.5
JC = 3.5

qpo = pd.DataFrame({
    'x0':[-0.750000],
    'x':[-0.750000],
    'y':[0.000000],
    'vx':[0.000000],
    'vy':[1.364734],
})
IC_qpo = [
    qpo['x'][0], qpo['y'][0], qpo['vx'][0], qpo['vy'][0],
    1,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
]
jacobi = calc_Jacobi(mu, qpo['x'][0], qpo['y'][0], qpo['vx'][0], qpo['vy'][0])
x0 = qpo['x'][0]
period_qpo = 3*2*pi
qpo_prop = solve_ivp(planar_ode, [0, period_qpo], IC_qpo, args=(mu,), rtol=1e-12,atol=1e-14)
df_qpo = pd.DataFrame({
    'name': f'qpo @ x={x0:0.3f}',
    't':qpo_prop.t,
    'x':qpo_prop.y[0],
    'y':qpo_prop.y[1],
})

df_dpo_table = pd.DataFrame({
    'x':[qpo['x'][0]],
    'y':[qpo['y'][0]],
    'vx':[qpo['vx'][0]],
    'vy':[qpo['vy'][0]],
    'jacobi':[jacobi]
})

ic_table = (
    GT(df_dpo_table)
    .tab_header(
        title=md(f"Initial Conditions of the Quasi-Perodic Orbit Found on Poincaré Map ({ps}, Lillian Shido)")
    )
    .cols_label(
        x="{{x}}<br>[nd]",
        y="{{y}}<br>[nd]",
        vx="{{vx}}<br>[nd]",
        vy="{{vy}}<br>[nd]"
    )
    .fmt_number(
        columns=["x","y","vx","vy"],
        decimals=4
    )
    .fmt_number(
        columns=['jacobi'],
        decimals=5
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
ic_table.show()

# qpo plot
x_min = -0.8
x_max = -0.2
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2
qpo_chart = alt.Chart(df_qpo).mark_line(strokeWidth=1,clip=True).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='sinebow')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=[f"Quasi-Periodic Orbit at x={x0:0.3f}[nd] for JC={JC}",f"({ps}, Lillian Shido)"]
)

qpo_chart_layer = alt.layer(qpo_chart)
qpo_chart_layer.save(f'{ps}_qpo.png', ppi=200)

pdb.set_trace()