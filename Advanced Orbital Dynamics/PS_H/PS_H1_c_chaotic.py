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


from methods import calc_Jacobi, planar_ode

mu = 0.5
JC = 3.5

chaotic = pd.DataFrame({
    'x0':[-0.400000],
    'x':[-0.400000],
    'y':[0.000000],
    'vx':[0.000000],
    'vy':[2.787671],
})
IC_chaotic = [
    chaotic['x'][0], chaotic['y'][0], chaotic['vx'][0], chaotic['vy'][0],
    1,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
]
jacobi = calc_Jacobi(mu, chaotic['x'][0], chaotic['y'][0], chaotic['vx'][0], chaotic['vy'][0])
x0 = chaotic['x'][0]
period_chaotic = 3*2*pi
chaotic_prop = solve_ivp(planar_ode, [0, period_chaotic], IC_chaotic, args=(mu,), rtol=1e-12,atol=1e-14)
df_chaotic = pd.DataFrame({
    'name': f'chaotic @ x={x0:0.3f}',
    't':chaotic_prop.t,
    'x':chaotic_prop.y[0],
    'y':chaotic_prop.y[1],
})

df_chaotic_table = pd.DataFrame({
    'x':[chaotic['x'][0]],
    'y':[chaotic['y'][0]],
    'vx':[chaotic['vx'][0]],
    'vy':[chaotic['vy'][0]],
    'jacobi':[jacobi]
})

ic_table = (
    GT(df_chaotic_table)
    .tab_header(
        title=md(f"Initial Conditions of the Chaotic Orbit Found on Poincaré Map ({ps}, Lillian Shido)")
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

# chaotic plot
x_min = -1.2
x_max = 1.2
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2
chaotic_chart = alt.Chart(df_chaotic).mark_line(strokeWidth=1,clip=True).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    # x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    # y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N', scale=alt.Scale(scheme='sinebow')).title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=[f"Chaotic Orbit at x={x0:0.3f}[nd] for JC={JC}",f"({ps}, Lillian Shido)"]
)

chaotic_chart_layer = alt.layer(chaotic_chart)
chaotic_chart_layer.save(f'{ps}_chaotic.png', ppi=200)

pdb.set_trace()