ps = "H2"
# Author: Lillian Shido
# Date: 12/17/2025

import pandas as pd
from great_tables import GT, md, system_fonts

from methods import system_properties, calc_L1, calc_L2
from constants import mu_Earth, mu_Moon, a_Moon


mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu)
x_L2, y_L2 = calc_L2(mu)

df_object_locations = pd.DataFrame({
    'Object':['Earth', 'Moon', 'L1', 'L2', 'L3', 'L4', 'L5'],
    'x':[-0.01215, 0.98785, 0.83692, 1.1557, -1.0051, 0.48785, 0.48785],
    'y':[0,0,0,0,0,0.86603,-0.86603],
    'JC':['--','--',3.1883,3.1722,3.0121,2.988,2.988]
})

objects_table = (
    GT(df_object_locations)
    .tab_header(
        title=md(f"Locations and Jacobi Constants of Objects in Earth-Moon System<br>({ps}, Lillian Shido)")
    )
    .cols_align(
        align="center"
    )
    .opt_table_outline()
    .opt_stylize()
    .opt_table_font(font=system_fonts(name="industrial"))
    .opt_horizontal_padding(scale=2)
)
objects_table.show()