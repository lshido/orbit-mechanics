ps = "H1 part b"
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

from methods import planar_ode, find_halfperiod_with_jacobi, calc_Jacobi

mu = 0.5
JC = 4.5

periodic_orbit = pd.DataFrame({
    'x0':[-0.300000],
    'x':[-0.300000],
    'y':[0.000000],
    'vx':[0.000000],
    'vy':[1.356466]
})
# Find the exact periodic orbit
tf_po, converged_IC_po = find_halfperiod_with_jacobi(JC, periodic_orbit['x'][0], periodic_orbit['vy'][0], mu, -0.01)
period_po = 2*tf_po
x0 = converged_IC_po[0]
jacobi = calc_Jacobi(mu, converged_IC_po[0], converged_IC_po[1], converged_IC_po[2], converged_IC_po[3])
# Now propagate to plot the orbit
IC_po = [
    converged_IC_po[0], converged_IC_po[1], converged_IC_po[2], converged_IC_po[3],
    1,0,0,0, # Identity matrix for phi ICs
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
]
po_prop = solve_ivp(planar_ode, [0, period_po], IC_po, args=(mu,), rtol=1e-12,atol=1e-14)
df_po = pd.DataFrame({
    'name': f'Periodic Orbit @ x={x0:.3f}',
    't':po_prop.t,
    'x':po_prop.y[0],
    'y':po_prop.y[1],
})

df_error = pd.DataFrame({
    'x':[abs(po_prop.y[0,0]-po_prop.y[0,-1])],
    'y':[abs(po_prop.y[1,0]-po_prop.y[1,-1])],
    'vx':[abs(po_prop.y[2,0]-po_prop.y[2,-1])],
    'vy':[abs(po_prop.y[3,0]-po_prop.y[3,-1])]
})

error_table = (
    GT(df_error)
    .tab_header(
        title=md(f"Error between Initial and Final States ({ps}, Lillian Shido)")
    )
    .fmt_scientific(
        columns=["x","y","vx","vy"],
        decimals=4
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
error_table.show()

df_po_table = pd.DataFrame({
    'period':[period_po],
    'x':[converged_IC_po[0]],
    'y':[converged_IC_po[1]],
    'vx':[converged_IC_po[2]],
    'vy':[converged_IC_po[3]],
    'jacobi':[jacobi]
})

ic_table = (
    GT(df_po_table)
    .tab_header(
        title=md(f"Initial Conditions of the Periodic Orbit Found on Hénon's Map ({ps}, Lillian Shido)")
    )
    .cols_label(
        x="{{x}}<br>[nd]",
        y="{{y}}<br>[nd]",
        vx="{{vx}}<br>[nd]",
        vy="{{vy}}<br>[nd]"
    )
    .fmt_number(
        columns=["x","y","vx","vy",'period'],
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

# periodic orbit plot
x_min = 0.8
x_max = 1.20
y_lim = (x_max-x_min)/2
z_lim = (x_max-x_min)/2
po_chart = alt.Chart(df_po).mark_line(strokeWidth=1,clip=True).encode(
    # x=alt.X('x:Q', scale=alt.Scale(domain=[-0.9,0]), axis=alt.Axis(title='x [non-dim]')),
    x=alt.X('x:Q', axis=alt.Axis(title='x [non-dim]')),
    # y=alt.Y('vx:Q', scale=alt.Scale(domain=[-5,5]), axis=alt.Axis(title='vx [non-dim]')),
    y=alt.Y('y:Q', axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N').title(None),
    order='t'
).properties(
    width=400,
    height=400,
    title=[f"Periodic Orbit @ x={x0:.3f}[nd]",f"({ps}, Lillian Shido)"]
)

qpo_chart_layer = alt.layer(po_chart)
# qpo_chart_layer.save(f'{ps}_po.png', ppi=200)

pdb.set_trace()