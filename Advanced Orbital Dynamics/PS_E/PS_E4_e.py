ps = "E4 part e"
# Author: Lillian Shido
# Date: 11/9/2025

import pdb
import numpy as np
from math import pi
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


from constants import mu_Earth, mu_Moon, a_Moon
from methods import system_properties, calc_L1, build_A_matrix_collinear, spatial_ode, calc_spatial_Jacobi

import warnings
warnings.filterwarnings("ignore",category=PendingDeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)

mu = mu_Moon/(mu_Earth + mu_Moon)

# Properties of the system
mu, l_char, t_char, x_Earth, x_Moon = system_properties(mu_Earth, mu_Moon, a_Moon)
x_L1, y_L1 = calc_L1(mu, a_Moon)
jacobi_L1 = calc_spatial_Jacobi(mu, x_L1, y_L1, 0, 0, 0, 0)

# Determine the stable and unstable manifolds 
d_dim = 0.1 #[km]
d_val = d_dim/l_char #[non-dim]
print(d_val)
tf = 2.5
pos_tspan = [0,tf*pi]
neg_tspan = [0,-tf*pi]

A = build_A_matrix_collinear(mu, x_L1, y_L1, 0)
eigenvalues, eigenvectors = np.linalg.eig(A)

tau = 1/eigenvalues[1]
print(f"time constant [non-dim]:{tau}")
print(f"time constant [dim]:{tau*t_char/3600/24}")
tau_eval = np.arange(0,tf,tau.real)

# Plot eigenvectors
fig1 = plt.figure()
ax1 = fig1.add_subplot()
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax1.scatter(x_L1, y_L1, s=25, label="L1", color="green",zorder=3)
ax1.scatter(x_Moon, 0, s=30, label="Moon", color="gray")
ax1.scatter(x_Earth,0,s=30, label="Earth", color="blue")
for i in range(0,6):
    if i==1 or i==0:
        print(f"vector {i+1} i-comp: {eigenvectors[:,i][0].real[0,0]}")
        print(f"vector {i+1} y-comp: {eigenvectors[:,i][1].real[0,0]}")
        if i==0:
            label=f"$E^S$"
            linestyle="solid"
            label_point="$x_S^+$"
            label_manifold="$W_{{loc}}^S$"
            color="blue"
        elif i==1:
            label=f"$E^U$"
            linestyle="dashed"
            label_point="$x_U^+$"
            label_manifold="$W_{{loc}}^U$"
            color="red"
        x_eig = x_L1 + eigenvectors[:,i][0].real[0,0]
        y_eig = y_L1 + eigenvectors[:,i][1].real[0,0]
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
        ax1.axline((x_L1,y_L1), (x_eig,y_eig), label=label,linestyle=linestyle)
        # ax1.scatter(x0_pos, y0_pos, label="$x_U^+$", s=10,zorder=3)
        x0_neg = x_L1 - scaled_nu_W[0].real[0,0]
        y0_neg = y_L1 - scaled_nu_W[1].real[0,0]
        vx0_neg = -scaled_nu_W[3].real[0,0]
        vy0_neg = -scaled_nu_W[4].real[0,0]
        check_neg_distance = np.linalg.norm(np.array([[x0_neg-x_L1],[y0_neg-y_L1]]))*l_char
        print(f"Check neg distance!: {check_neg_distance}")
        # ax1.scatter(x0_neg, y0_neg, label="$x_U^-$", s=10,zorder=3)
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
        ax1.plot(pos_prop.y[0], pos_prop.y[1],label=label_manifold+" +", color=color, linestyle=linestyle, zorder=2.5)
        neg_prop = solve_ivp(spatial_ode, neg_tspan, IC_pos, args=(mu,), rtol=1e-12,atol=1e-14)
        ax1.plot(neg_prop.y[0], neg_prop.y[1], color=color, linestyle=linestyle, zorder=2.5)
        # Now do the other sides of the manifolds
        IC_neg = [
            x0_neg,y0_neg,0,vx0_neg,vy0_neg,0,# IC states
            1,0,0,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,0,0,
            0,0,1,0,0,0,
            0,0,0,1,0,0,
            0,0,0,0,1,0,
            0,0,0,0,0,1
        ]
        minus_prop = solve_ivp(spatial_ode, pos_tspan, IC_neg, args=(mu,), rtol=1e-12,atol=1e-14)
        ax1.plot(minus_prop.y[0], minus_prop.y[1],label=label_manifold+" -", color=color, linestyle=linestyle, zorder=2.5)
        neg_prop = solve_ivp(spatial_ode, neg_tspan, IC_neg, args=(mu,), rtol=1e-12,atol=1e-14)
        ax1.plot(neg_prop.y[0], neg_prop.y[1], color=color, linestyle=linestyle, zorder=2.5)
        tau_pos_prop = solve_ivp(spatial_ode, [0, tau.real], IC_pos, args=(mu,), rtol=1e-12,atol=1e-14)
        # ax1.scatter(tau_pos_prop.y[0], tau_pos_prop.y[1], marker="*", color="blue", zorder=3, label=r'$\tau$')
        tau_neg_prop = solve_ivp(spatial_ode, [0, tau.real], IC_neg, args=(mu,), rtol=1e-12,atol=1e-14)
        # ax1.scatter(tau_neg_prop.y[0], tau_neg_prop.y[1], marker="*", color="blue", zorder=3)
        
        plus_jacobi = calc_spatial_Jacobi(mu, pos_prop.y[0,-1], pos_prop.y[1,-1], pos_prop.y[2,-1], pos_prop.y[3,-1],pos_prop.y[4,-1],pos_prop.y[5,-1])
        minus_jacobi = calc_spatial_Jacobi(mu, minus_prop.y[0,-1], minus_prop.y[1,-1], minus_prop.y[2,-1], minus_prop.y[3,-1],minus_prop.y[4,-1],minus_prop.y[5,-1])

        distance_pos = 0
        distance_neg = 0
        for m in range(len(tau_pos_prop.y[0])):
            if m < len(tau_pos_prop.y[0])-1:
                try:
                    delta_x_pos = tau_pos_prop.y[0,m+1] - tau_pos_prop.y[0,m]
                except:
                    pdb.set_trace()
                delta_y_pos = tau_pos_prop.y[1,m+1] - tau_pos_prop.y[1,m]
                s_pos = np.linalg.norm([[delta_x_pos], [delta_y_pos]])
                distance_pos = distance_pos + s_pos
                delta_x_neg = tau_neg_prop.y[0,m+1] - tau_neg_prop.y[0,m]
                delta_y_neg = tau_neg_prop.y[1,m+1] - tau_neg_prop.y[1,m]
                s_neg = np.linalg.norm([[delta_x_neg], [delta_y_neg]])
                distance_neg = distance_neg + s_neg
            else:
                break
        print(f"Distance of positive trajectory: {distance_pos} [non-dim]")
        print(f"Distance of positive trajectory: {distance_pos*l_char} [km]")
        print(f"Distance of negative trajectory: {distance_neg} [non-dim]")
        print(f"Distance of negative trajectory: {distance_neg*l_char} [km]")

plt.legend(ncol=2,framealpha=1)
lim=4e-4
ax1.axis('equal')
# ax1.set_aspect(aspect=1, adjustable="box")
# ax1.set(xlim=(x_L1-lim, x_L1+lim),ylim=(-lim,lim))
ax1.set_xlabel("x [non-dim]")
ax1.set_ylabel("y [non-dim]")
plt.title("Stable and Unstable Half-Manifolds Propagated from L1"+"\n"+rf"t={tf}$\pi$ ({ps}, Lillian Shido)")
plt.grid()
plt.savefig(f'manifolds_{ps}.png', dpi=300, bbox_inches='tight')
plt.show()