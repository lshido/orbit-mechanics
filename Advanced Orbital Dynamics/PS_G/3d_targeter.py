def find_spatial_halfperiod(starting_x, z_guess, ydot_guess, mu, tolerance=1e-12, max_iterations=50):
    """
    Calculates the resulting half period states for a 3-dimensional perpindicular x-axis crossing
    
    args:
        starting_x: starting x location that you want to propagate half period for 
        z_guess: the guessed z to start the propagation
        ydot_guess: the guessed ydot to start the propagation
        mu: the mu value for the system
        tolerance: convergence tolerance for the targeter, default: 1e-12
        max_iterations: maximum number of iterations before terminating, default=50

    returns:
        counter: number of iterations targeter took to converge. 50 is the maximum
        tf: half period
        arrival_states: list of values of the final states
        converged_initial_states: list of values of the converged initial states
    """
    print(f"Looking for the initial conditions for half-period with initial guesses starting x = {starting_x}, z= {z_guess}, and vy guess = {ydot_guess}")
    x0 = starting_x # Starting x
    z0 = z_guess # Initial z guess
    vy0 = ydot_guess # Initial ydot guess
    
    counter = 0
    while True:
        IC_guess = [
            x0, 0, z0, 0, vy0, 0,
            1,0,0,0,0,0, # Identity matrix for phi ICs
            0,1,0,0,0,0,
            0,0,1,0,0,0,
            0,0,0,1,0,0,
            0,0,0,0,1,0,
            0,0,0,0,0,1
        ]
        counter = counter + 1
        # Step 1: Using ICs, propagate EOM+STM until second x-axis cross, get STM, get xdot, ydot
        prop = solve_ivp(spatial_ode, [0, 50*pi], IC_guess, events=crossxEvent, args=(mu,), rtol=1e-12,atol=1e-14)
        stm = prop.y[6:42,-1].reshape(6,6) # turn into 6x6 phi matrix
        tf = prop.t_events[0][1]
        phi_43 = stm[3,2]
        phi_45 = stm[3,4]
        phi_63 = stm[5,2]
        phi_65 = stm[5,4]
        phi_23 = stm[1,2]
        phi_25 = stm[1,4]
        rf = prop.y[0:3,-1].reshape(3,1) # position at tf, turn into 3x1 vector
        vf = prop.y[3:6,-1].reshape(3,1) # velocity at tf, turn into 3x1 vector
        xf = prop.y[0,-1]
        yf = prop.y[1,-1]
        zf = prop.y[2,-1]
        vxf = prop.y[3,-1] # get xdot to compare against tolerance
        vyf = prop.y[4,-1] # get ydot to calculate the acceleration
        vzf = prop.y[5,-1] 
        # Check if the vxf is close to 0 within acceptable margins
        if abs(vxf) > tolerance and abs(vzf) > tolerance and counter <= max_iterations: # If not, recalculate the delta_v0 and try again
            # Calc the acceleration in x
            a_x, a_y, a_z = eval_spatial_acceleration(xf, yf, zf, vxf, vyf, mu)
            # Step 3: Calc new delta_vy0
            term_1 = phi_43 - a_x*phi_23/vyf
            term_2 = phi_45 - a_x*phi_25/vyf
            term_3 = phi_63 - a_z*phi_23/vyf
            term_4 = phi_65 - a_z*phi_25/vyf
            a = np.array([[term_1, term_2],[term_3, term_4]])
            b = np.array([-vxf,-vzf])
            sol = np.linalg.solve(a,b)
            delta_vy0 = sol[0]
            delta_z0 = sol[1]
            # pdb.set_trace()
            # Step 4: Use new deltas_ydot_t to calc new ICs
            if delta_vy0 > 0.1:
                delta_vy0 = delta_vy0/3
            if delta_z0 > 0.001:
                delta_z0 = delta_z0/10
            vy0 = vy0 + delta_vy0
            z0 = z0 + delta_z0
            continue
        else: # If error is within acceptable margins, break out of iterative loop
            arrival_states = [rf[0,0], rf[1,0], rf[2,0], vf[0,0], vf[1,0], vf[2,0]]
            converged_initial_states = [IC_guess[0], IC_guess[1], IC_guess[2], IC_guess[3], IC_guess[4], IC_guess[5]]            
            break
    return counter, tf, arrival_states, converged_initial_states

def eval_spatial_acceleration(x,y,z,xdot,ydot,mu):
    # Simply evaluate the EOMs with the position and velocity
    d = ((x + mu)**2 + y**2 + z**2)**(1/2)
    r = ((x - 1 + mu)**2 + y**2 + z**2)**(1/2)
    a_x = 2*ydot + x - (1 - mu) * (x + mu) / d**3 -\
        mu * (x - 1 + mu) / r**3
    a_y = -2*xdot + y - (1 - mu) * y / d**3 -\
        mu * y/r**3
    a_z = -(1-mu)*z/d**3 - mu*z/r**3
    return a_x, a_y, a_z

def crossxEvent(t, sv, mu):
    return sv[1] # Return value of y
crossxEvent.terminal = 2

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