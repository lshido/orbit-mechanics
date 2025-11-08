ps = "E1 part a"
# Problem D4 Part b iii
# Author: Lillian Shido
# Date: 10/26/2025

import pdb
import copy
import numpy as np
import pandas as pd
from math import pi, sqrt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
from scipy import stats
from great_tables import GT, md, html, style, loc, system_fonts
from pypalettes import load_cmap

from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, build_A_matrix_collinear

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)

A = build_A_matrix_collinear(mu, x_L1, y_L1, 0)
eigenvalues, eigenvectors = np.linalg.eig(A)
try:
    print("\n".join([" ".join(f"{item:15.5f}" for item in row) for row in eigenvectors.toarray()]))
except:
    pdb.set_trace()
print("\n".join([" ".join(f"{item:15.5f}" for item in row) for row in eigenvalues]))