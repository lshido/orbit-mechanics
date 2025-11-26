ps = "F1 part b"
# Author: Lillian Shido
# Date: 11/14/2025

import pdb
import numpy as np
from math import pi, sqrt

import pandas as pd
from scipy.integrate import solve_ivp
import altair as alt


from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, build_A_matrix_collinear, spatial_ode, calc_spatial_Jacobi, calc_ZVC_Jacobi

import warnings
warnings.filterwarnings("ignore",category=PendingDeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)

# Determine the stable and unstable manifolds 
d_dim = 30 #[km]
d_val = d_dim/l_char #[non-dim]
print(d_val)
tf = 1*pi
pos_tspan = [0,tf]
neg_tspan = [0,-tf]

A = build_A_matrix_collinear(mu, x_L1, y_L1, 0)
eigenvalues, eigenvectors = np.linalg.eig(A)
# Plot eigenvectors
L1 = pd.DataFrame({'name':["L1"],'x':[x_L1],'y':[y_L1]})
moon = pd.DataFrame({'name':["Moon"],'x':[x_Moon],'y':[0]})
eigenspace = pd.DataFrame({})
manifold = pd.DataFrame({})
for i in range(0,6):
    if i==0 or i==1:
        print(f"vector {i+1} i-comp: {eigenvectors[:,i][0].real[0,0]}")
        print(f"vector {i+1} y-comp: {eigenvectors[:,i][1].real[0,0]}")
        if i==0:
            name="Stable Eigenspace"
            label=f"$E^S$"
            linestyle="solid"
            label_point="$x_S^+$"
            label_manifold="Stable Half-Manifold"
            color="blue"
        elif i==1:
            name="Unstable Eigenspace"
            label=f"$E^U$"
            linestyle="dashed"
            label_point="$x_U^+$"
            label_manifold="Unstable Half-Manifold"
            color="red"
        scale = 1 # Extend eigenvector line out
        x_eig_max = x_L1 + scale*eigenvectors[:,i][0].real[0,0]
        x_eig_min = x_L1 + scale*-eigenvectors[:,i][0].real[0,0]
        y_eig_max = y_L1 + scale*eigenvectors[:,i][1].real[0,0]
        y_eig_min = y_L1 + scale*-eigenvectors[:,i][1].real[0,0]
        nu_W = eigenvectors[:,i]/np.linalg.norm(eigenvectors[:,i][0:3])
        print(nu_W)
        scaled_nu_W = d_val*nu_W
        x0_pos = x_L1 + scaled_nu_W[0].real[0,0]
        y0_pos = y_L1 + scaled_nu_W[1].real[0,0]
        vx0_pos = scaled_nu_W[3].real[0,0]
        vy0_pos = scaled_nu_W[4].real[0,0]
        print(f"{x0_pos:.4f}")
        print(f"{y0_pos:.4f}")
        print(f"{vx0_pos:.4f}")
        print(f"{vy0_pos:.4f}")
        # check the distance from eq
        check_distance = np.linalg.norm(np.array([[x0_pos-x_L1],[y0_pos-y_L1]]))*l_char
        print(f"Check distance!: {check_distance}")
        eigenspace_data = pd.DataFrame({
            'name': name,
            'x':[x_eig_min,x_L1,x_eig_max],
            'y':[y_eig_min,y_L1,y_eig_max]
        })
        eigenspace = pd.concat([eigenspace, eigenspace_data], ignore_index=True)
        IC_pos = [
            x0_pos,y0_pos,0,vx0_pos,vy0_pos,0,# IC states
            1,0,0,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,0,0,
            0,0,1,0,0,0,
            0,0,0,1,0,0,
            0,0,0,0,1,0,
            0,0,0,0,0,1
        ]
        pos_prop = solve_ivp(spatial_ode, pos_tspan, IC_pos, args=(mu,), rtol=1e-12,atol=1e-14)
        manifold_data = pd.DataFrame({
            'name': f"{label_manifold} +",
            't':pos_prop.t,
            'x':pos_prop.y[0],
            'y':pos_prop.y[1]
        })
        manifold = pd.concat([manifold, manifold_data], ignore_index=True)
        neg_prop = solve_ivp(spatial_ode, neg_tspan, IC_pos, args=(mu,), rtol=1e-12,atol=1e-14)
        manifold_data = pd.DataFrame({
            'name': f"{label_manifold} +",
            't':neg_prop.t,
            'x':neg_prop.y[0],
            'y':neg_prop.y[1]
        })
        manifold = pd.concat([manifold, manifold_data], ignore_index=True)
        # Now do the other sides of the manifolds
        x0_neg = x_L1 - scaled_nu_W[0].real[0,0]
        y0_neg = y_L1 - scaled_nu_W[1].real[0,0]
        vx0_neg = -scaled_nu_W[3].real[0,0]
        vy0_neg = -scaled_nu_W[4].real[0,0]
        IC_neg = [
            x0_neg,y0_neg,0,vx0_neg,vy0_neg,0,# IC states
            1,0,0,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,0,0,
            0,0,1,0,0,0,
            0,0,0,1,0,0,
            0,0,0,0,1,0,
            0,0,0,0,0,1
        ]
        pos_prop = solve_ivp(spatial_ode, pos_tspan, IC_neg, args=(mu,), rtol=1e-12,atol=1e-14)
        manifold_data = pd.DataFrame({
            'name':f"{label_manifold} -",
            't':pos_prop.t,
            'x':pos_prop.y[0],
            'y':pos_prop.y[1]
        })
        manifold = pd.concat([manifold, manifold_data], ignore_index=True)
        neg_prop = solve_ivp(spatial_ode, neg_tspan, IC_neg, args=(mu,), rtol=1e-12,atol=1e-14)
        manifold_data = pd.DataFrame({
            'name': f"{label_manifold} -",
            't':neg_prop.t,
            'x':neg_prop.y[0],
            'y':neg_prop.y[1]
        })
        manifold = pd.concat([manifold, manifold_data], ignore_index=True)



#  Calc ZVCs
df_zvc = pd.DataFrame()
tolerance = 1e-12
max_iterations = 100
C = calc_ZVC_Jacobi(mu, x_L1, y_L1)
# Give an initial guess for y that is around the Earth
for x in np.arange(0,1.5,1e-3): # Find the curve for -1.5 < x < 1.5
    for y in np.arange(0,1.5,0.1): # These are my guesses
        counter = 0
        while True:
            counter = counter + 1
            d = sqrt((x+mu)**2 + y**2)
            r = sqrt((x-1+mu)**2 + y**2)
            f = x**2 + y**2 + (2*(1-mu)/d) + (2*mu/r) - C
            f_prime_y = 2*y*( 1 - (1-mu)/d**3 - mu/r**3)
            delta = f/f_prime_y
            if abs(f) > tolerance:
                if counter > max_iterations: # check for max iterations
                    new_C = calc_ZVC_Jacobi(mu, x, y)
                    C_error = abs(new_C - C)
                    if C_error < 1e-12:
                        zvc_data = pd.DataFrame({
                            "name":'ZVC',
                            "x":[x],
                            "y":[y],
                            "new_C":[new_C],
                            "error":[error_C],
                        })
                        df_zvc = pd.concat([df_zvc,zvc_data],ignore_index=True)
                    else:
                        break
                elif d == 0 or r == 0: # check for dividing by 0
                    break
                elif abs(f_prime_y) <= tolerance: # check for zero derivative
                    break
                else:
                    y = y - delta
                    continue
            else:
                new_C = calc_ZVC_Jacobi(mu, x, y)
                error_C = abs(new_C - C)
                zvc_data = pd.DataFrame({
                    "name":'ZVC',
                    "x":[x],
                    "y":[y],
                    "new_C":[new_C],
                    "error":[error_C],
                })
                df_zvc = pd.concat([df_zvc,zvc_data],ignore_index=True)
                break
for y in np.arange(0.5,1.5,1e-3): # Find the curve for 0 < y < 1.21
    for x in np.arange(-1.5,1.5,0.1): # These are my guesses for x
        counter = 0
        while True:
            counter = counter + 1
            d = sqrt((x+mu)**2 + y**2)
            r = sqrt((x-1+mu)**2 + y**2)
            f = x**2 + y**2 + (2*(1-mu)/d) + (2*mu/r) - C
            f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d**3 - 2*mu*(x-1+mu)/r**3
            delta = f/f_prime_x
            if abs(f) > tolerance:
                if counter > max_iterations: # check for max iterations
                    new_C = calc_ZVC_Jacobi(mu, x, y)
                    C_error = abs(new_C - C)
                    if C_error < 1e-12:
                        zvc_data = pd.DataFrame({
                            "name":'ZVC',
                            "x":[x],
                            "y":[y],
                            "new_C":[new_C],
                            "error":[error_C],
                        })
                        df_zvc = pd.concat([df_zvc,zvc_data],ignore_index=True)
                    else:
                        break
                elif d == 0 or r == 0: #check for dividing by 0
                    break
                elif abs(f_prime_x) <= tolerance: #check for zero derivative
                    break
                else:
                    x = x - delta
                    continue
            else:
                new_C = calc_ZVC_Jacobi(mu, x, y)
                error_C = abs(new_C - C)
                zvc_data = pd.DataFrame({
                    "name":'ZVC',
                    "x":[x],
                    "y":[y],
                    "new_C":[new_C],
                    "error":[error_C],
                })
                df_zvc = pd.concat([df_zvc,zvc_data],ignore_index=True)
                break

neg_y_zvc_data = pd.DataFrame({
    'name':'ZVC',
    'x':df_zvc['x'],
    'y':-df_zvc['y'],
    'new_C':df_zvc['new_C'],
    'error':df_zvc['error']
})
df_zvc = pd.concat([df_zvc,neg_y_zvc_data],ignore_index=True)

# Build plot
x_min = 0.65
x_max = 1.05
y_lim = (x_max-x_min)/2
base = alt.Chart(manifold).mark_line(clip=True).encode(
    x=alt.X('x:Q', scale=alt.Scale(domain=[x_min,x_max]), axis=alt.Axis(title='x [non-dim]')),
    y=alt.Y('y:Q', scale=alt.Scale(domain=[-y_lim,y_lim]), axis=alt.Axis(title='y [non-dim]')),
    color=alt.Color('name:N').title(None),
    order='t',
).properties(
    width=400,
    height=400,
    title=f"Stable and Unstable Manifolds with ZVC at L1 ({ps}, Lillian Shido)"
)

zvc = alt.Chart(df_zvc).mark_point(size=5,filled=True,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['ZVC'], range=['forestgreen'])).title(None)
)

eigspace = alt.Chart(eigenspace).mark_line(strokeDash=(4,4),clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N').title(None)
)

# scale = alt.Scale(domain=['L1','Moon'], range=['darkblue','gray'])
moon_loc = alt.Chart(moon).mark_point(filled=True,size=50,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['Moon'], range=['gray'])).title(None)
)

L1_loc = alt.Chart(L1).mark_point(filled=True,size=30,clip=True).encode(
    x='x:Q',
    y='y:Q',
    color=alt.Color('name:N', scale=alt.Scale(domain=['L1'], range=['darkblue'])).title(None)
)

final = alt.layer(base, zvc, eigspace, moon_loc, L1_loc).resolve_scale(color='independent')
final.save(f'manifolds_zvc {ps}.png', ppi=200)