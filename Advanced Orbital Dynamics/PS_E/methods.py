from math import pi, sqrt
import numpy as np
from scipy.integrate import solve_ivp
import pdb


def system_properties(mu_major,mu_minor, a):
    mu = mu_minor/(mu_major+mu_minor)
    x_major = -mu
    x_minor = 1-mu
    l_char = a
    t_char = sqrt(l_char**3/(mu_major+mu_minor))
    return mu,l_char, t_char, x_major, x_minor

def calc_error(actual, ideal):
    return abs(actual-ideal)/abs(ideal)

def calc_L1(mu, a):
    gamma = 0.2
    tolerance = 1e-12
    counter = 0
    while True:
        counter = counter + 1
        f = 1-mu-gamma-((1-mu)/(1-gamma)**2)+(mu/gamma**2)
        f_prime = -1-(2*(1-mu)/(1-gamma)**3)-(2*mu/gamma**3)
        if abs(f) > tolerance:
            gamma = gamma - f/f_prime;
            continue
        else:
            x = 1 - mu - gamma
            # Add check for accleration
            d_x = x*a+mu
            r_x = x*a-1+mu
            d = ((x*a+mu)**2)**(1/2)
            r = ((x*a-1+mu)**2)**(1/2)
            accel = -(1-mu)/d**3*d_x - mu/r**3*r_x
            # Add check for partial wrt gamma
            partial = (1-mu)/(1-gamma)**2 - mu/(gamma)**2 - (1-mu-gamma)
            x_L1 = x
            y_L1 = 0
            counter = 0
            break
    return x_L1, y_L1

def calc_initial_velocities(xi_0, eta_0, x_L, y_L, mu):
    """
    Calculates the initial velocities around the collinear libration points.
    
    Args:
        xi_0: x-direction distance offset from the libration point
        eta_0: y-direction distance offset from the libration point
        x_L: The x value of the libration point
        y_L: The y value of the libration point
        mu: The mu value for the system

    Returns:
        xi_dot: The initial velocity in the xi direction
        eta_dot: The initial velocity in the eta direction
    """
    d = sqrt((x_L+mu)**2 + y_L**2)
    r = sqrt((x_L-1+mu)**2 + y_L**2)
    U_xx = 1 - (1-mu)/d**3 - mu/r**3 + 3*(1-mu)*(x_L+mu)**2/d**5 + 3*mu*(x_L-1+mu)**2/r**5
    U_yy = 1 - (1-mu)/d**3 - mu/r**3
    B_1 = 2 - (U_xx + U_yy)/2
    B_2_squared = -U_xx*U_yy
    s = (B_1 + (B_1**2 + B_2_squared)**(1/2))**(1/2)
    B_3 = (s**2 + U_xx)/(2*s)
    xi_dot_0 = eta_0*s/B_3
    eta_dot_0 = -B_3*xi_0*s
    return xi_dot_0, eta_dot_0

def planar_eoms(t,sv,mu):
    
    # Set up the EOM ODEs
    eoms = [
        sv[2],
        sv[3],
        2*sv[3] + sv[0] - (1 - mu) * (sv[0] + mu) / ((sv[0] + mu)**2 + sv[1]**2)**(3/2)-\
        mu * (sv[0] - 1 + mu) / ((sv[0] - 1 + mu)**2 + sv[1]**2)**(3/2),
        -2*sv[2] + sv[1] - (1 - mu) * sv[1] / ((sv[0] + mu)**2 + sv[1]**2)**(3/2) -\
        mu * sv[1]/((sv[0] - 1 + mu)**2 + sv[1]**2)**(3/2)
    ]
    return eoms

def planar_ode(t,sv,mu):
    # Set up the EOM ODEs
    eoms = [
        sv[2],
        sv[3],
        2*sv[3] + sv[0] - (1 - mu) * (sv[0] + mu) / ((sv[0] + mu)**2 + sv[1]**2)**(3/2)-\
        mu * (sv[0] - 1 + mu) / ((sv[0] - 1 + mu)**2 + sv[1]**2)**(3/2),
        -2*sv[2] + sv[1] - (1 - mu) * sv[1] / ((sv[0] + mu)**2 + sv[1]**2)**(3/2) -\
        mu * sv[1]/((sv[0] - 1 + mu)**2 + sv[1]**2)**(3/2)
    ]
    # Calc the partials using the current x and y values
    d = sqrt((sv[0]+mu)**2 + sv[1]**2);
    r = sqrt((sv[0]-1+mu)**2 + sv[1]**2);
    U_xx = 1 - (1-mu)/d**3 - mu/r**3 + 3*(1-mu)*(sv[0]+mu)**2/d**5 + 3*mu*(sv[0]-1+mu)**2/r**5;
    U_yy = 1 - (1-mu)/d**3 - mu/r**3 + 3*(1-mu)*sv[1]**2/d**5 + 3*mu*sv[1]**2/r**5;
    U_xy = 3*(1-mu)*(sv[0]+mu)*sv[1]/d**5 + 3*mu*(sv[0]-1+mu)*sv[1]/r**5;

    # Set up the STM ODEs
    stm = [
    sv[12],
    sv[13],
    sv[14],
    sv[15],
    sv[16],
    sv[17],
    sv[18],
    sv[19],
    U_xx*sv[4] + U_xy*sv[8] + 2*sv[16],
    U_xx*sv[5] + U_xy*sv[9] + 2*sv[17],
    U_xx*sv[6] + U_xy*sv[10] + 2*sv[18],
    U_xx*sv[7] + U_xy*sv[11] + 2*sv[19],
    U_xy*sv[4] + U_yy*sv[8] - 2*sv[12],
    U_xy*sv[5] + U_yy*sv[9] - 2*sv[13],
    U_xy*sv[6] + U_yy*sv[10] - 2*sv[14],
    U_xy*sv[7] + U_yy*sv[11] - 2*sv[15]
    ]
    # Combine them into one big matrix
    combined = eoms + stm
    return combined

def spatial_eoms(t,sv, mu):
    # Set up the EOM ODEs
    d = ((sv[0] + mu)**2 + sv[1]**2 + sv[2]**2)**(1/2)
    r = ((sv[0] - 1 + mu)**2 + sv[1]**2 + sv[2]**2)**(1/2)
    eoms = [
        sv[3],
        sv[4],
        sv[5],
        2*sv[4] + sv[0] - (1 - mu) * (sv[0] + mu) / d**3 -\
        mu * (sv[0] - 1 + mu) / r**3,
        -2*sv[3] + sv[1] - (1 - mu) * sv[1] / d**3 -\
        mu * sv[1]/r**3,
        -(1-mu)*sv[2]/d**3 - mu*sv[2]/r**3
    ]
    return eoms

def spatial_ode(t,sv,mu):
    # Set up the EOM ODEs
    eoms = spatial_eoms(t,sv,mu)

    # Calc the partials using the current x and y values
    d = ((sv[0] + mu)**2 + sv[1]**2 + sv[2]**2)**(1/2)
    r = ((sv[0] - 1 + mu)**2 + sv[1]**2 + sv[2]**2)**(1/2)
    U_xx = 1 - (1-mu)/d**3 - mu/r**3 + 3*(1-mu)*(sv[0]+mu)**2/d**5 + 3*mu*(sv[0]-1+mu)**2/r**5
    U_yy = 1 - (1-mu)/d**3 - mu/r**3 + 3*(1-mu)*sv[1]**2/d**5 + 3*mu*sv[1]**2/r**5
    U_zz = -(1-mu)/d**3 - mu/r**3 + 3*(1-mu)*sv[2]**2/d**5 + 3*mu*sv[2]**2/r**5
    U_xy = 3*(1-mu)*(sv[0]+mu)*sv[1]/d**5 + 3*mu*(sv[0]-1+mu)*sv[1]/r**5
    U_xz = 3*(1-mu)*(sv[0]+mu)*sv[2]/d**5 + 3*mu*(sv[0]-1+mu)*sv[2]/r**5
    U_yz = 3*(1-mu)*sv[1]*sv[2]/d**5 + 3*mu*sv[1]*sv[2]/r**5
    
    # Build A matrix
    quad_1 = np.zeros((3,3))
    quad_2 = np.identity(3)
    quad_3 = np.array([
        [U_xx,U_xy,U_xz],
        [U_xy,U_yy,U_yz],
        [U_xz,U_yz,U_zz]
    ])
    quad_4 = np.array([
        [0,2,0],
        [-2,0,0],
        [0,0,0]
    ])
    A = np.bmat([
        [quad_1, quad_2],
        [quad_3, quad_4]
    ])
    # Set up the STM ODEs
    phi = sv[6:42].reshape(6,6)
    stm = A*phi
    # Combine them into one big matrix
    combined = eoms + np.squeeze(stm.reshape(-1)).tolist()[0]
    return combined

def crossxEvent(t, sv, mu):
    return sv[1] # Return value of y
crossxEvent.terminal = 2

def eval_acceleration(x,y,xdot,ydot,mu):
    # Simply evaluate the EOMs with the position and velocity
    a_x = 2*ydot + x - (1 - mu) * (x + mu) / ((x + mu)**2 + y**2)**(3/2)-\
        mu * (x - 1 + mu) / ((x - 1 + mu)**2 + y**2)**(3/2)
    a_y = -2*xdot + y - (1 - mu) * y / ((x + mu)**2 + y**2)**(3/2) -\
        mu * y/((x - 1 + mu)**2 + y**2)**(3/2)
    return a_x, a_y

def eval_kinematic(x0,y0,xdot_0,ydot_0,a_x, a_y,t_span):
    # Evaluate the kinematic equations with the time
    xdot = xdot_0 + a_x*t_span
    ydot = ydot_0 + a_y*t_span
    x = x0 + xdot_0*t_span
    y = y0 + ydot_0*t_span
    return x,y,xdot,ydot

def calc_Jacobi(mu, x, y, vx, vy):
    d = sqrt((x+mu)**2 + y**2)
    r = sqrt((x-1+mu)**2 + y**2)
    x_y_sq = (x**2+y**2)/2
    term_1 = (1-mu)/d
    term_2 = mu/r
    pseudo_U = term_1 + term_2 + x_y_sq
    v_squared = vx**2 + vy**2
    C = 2*pseudo_U - v_squared
    return C

def calc_spatial_Jacobi(mu, x, y, z, vx, vy, vz):
    d = sqrt((x+mu)**2 + y**2 + z**2)
    r = sqrt((x-1+mu)**2 + y**2 + z**2)
    x_y_sq = (x**2+y**2)/2
    term_1 = (1-mu)/d
    term_2 = mu/r
    pseudo_U = term_1 + term_2 + x_y_sq
    v_squared = vx**2 + vy**2 + vz**2
    C = 2*pseudo_U - v_squared
    return C

def find_halfperiod(starting_x, ydot_guess, mu, tolerance=1e-12, max_iterations=50):
    """
    Calculates the resulting half period states for a perpindicular x-axis crossing
    
    args:
        starting_x: starting x location that you want to propagate half period for 
        ydot_guess: the guessed ydot to start the propagation
        mu: the mu value for the system
        tolerance: convergence tolerance for the targeter, default: 1e-12
        max_iterations: maximum number of iterations before terminating, default=50

    returns:
        iterations: number of iterations targeter took to converge. 50 is the maximum
        tf: half period
        xf: x at arrival
        yf: y at arrival
        vxf: xdot at arrival
        vyf: ydot at arrival
    """
    x0 = starting_x # Initial x guess 
    vy0 = ydot_guess
    
    counter = 0
    while True:
        IC_guess = [
            x0, 0, 0, vy0, # Initial states
            1,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,
            0,0,1,0,
            0,0,0,1
        ]
        counter = counter + 1
        # Step 1: Using ICs, propagate EOM+STM until second x-axis cross, get STM, get xdot, ydot
        prop = solve_ivp(planar_ode, [0, 50*pi], IC_guess, events=crossxEvent, args=(mu,), rtol=1e-12,atol=1e-14)
        stm = prop.y[4:20,-1].reshape(4,4) # turn into 4x4 phi matrix
        tf = prop.t_events[0][1]
        phi_34 = stm[2,3]
        phi_24 = stm[1,3]
        rf = prop.y[0:2,-1].reshape(2,1) # position at tf, turn into 2x1 vector
        vf = prop.y[2:4,-1].reshape(2,1) # velocity at tf, turn into 2x1 vector
        vxf = prop.y[2,-1] # get xdot to compare against tolerance
        vyf = prop.y[3,-1] # get ydot to calculate the acceleration
        # print for info
        # print(f"{counter}, tf: {tf}, x:{IC_guess[0]}, y:{IC_guess[1]}, vx:{IC_guess[2]}, vy0:{IC_guess[3]}")
        # Check if the vxf is close to 0 within acceptable margins
        if abs(vxf) > tolerance and counter <= max_iterations: # If not, recalculate the delta_v0 and try again
            # Calc the acceleration in x
            a_x, a_y = eval_acceleration(rf[0,0], rf[1,0], vxf, vyf, mu)
            # Step 3: Calc new delta_vy0
            delta_vy0 = -vxf / (phi_34 - phi_24*(a_x/vyf))
            # Step 4: Use new deltas_ydot_t to calc new ICs
            if delta_vy0 > 0.1:
                delta_vy0 = delta_vy0/3
            vy0 = vy0 + delta_vy0
            continue
        else: # If error is within acceptable margins, break out of iterative loop
            arrival_states = [rf[0,0], rf[1,0], vf[0,0], vf[1,0]]
            converged_initial_states = [IC_guess[0], IC_guess[1], IC_guess[2], IC_guess[3]]            
            break
    return counter, tf, arrival_states, converged_initial_states

def calc_poincare_exponents(eig, period):
    omega = 1/period*np.log(eig)
    return omega

def calc_planar_monodromy_half(stm_half):
    G = np.array([
    [1,0,0,0],
    [0,-1,0,0],
    [0,0,-1,0],
    [0,0,0,1]
    ])
    I = np.identity(2)
    omega = np.array([
        [0,1],
        [-1,0]
        ])
    term1 = np.bmat([
        [np.zeros((2,2)) , -I],
        [I , -2*omega] 
        ])
    term2 = np.transpose(stm_half)
    term3 = np.bmat([
        [-2*omega , I],
        [-I , np.zeros((2,2))]
    ])
    monodromy_half = G*term1*term2*term3*G*stm_half
    return monodromy_half

def calc_spatial_monodromy_half(stm_half):
    G = np.array([
        [1,0,0,0,0,0],
        [0,-1,0,0,0,0],
        [0,0,1,0,0,0],
        [0,0,0,-1,0,0],
        [0,0,0,0,1,0],
        [0,0,0,0,0,-1]
        ])
    I = np.identity(3)
    omega = np.array([
        [0,1,0],
        [-1,0,0],
        [0,0,0]
        ])
    term1 = np.bmat([
        [np.zeros((3,3)) , -I],
        [I , -2*omega] 
        ])
    term2 = np.transpose(stm_half)
    term3 = np.bmat([
        [-2*omega , I],
        [-I , np.zeros((3,3))]
    ])
    monodromy_half = G*term1*term2*term3*G*stm_half
    return monodromy_half

def calc_stability_index(lambda_i, lambda_j):
    return 1/2*(lambda_i + lambda_j)