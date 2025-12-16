# Fixed Jacobi
def find_halfperiod_fixed_jacobi(JC, x_guess, ydot_guess, mu, tolerance=1e-12, max_iterations=50):
    """
    Calculates the resulting half period states for a perpindicular x-axis crossing
    
    args:
        JC: desired JC
        x_guess: starting x location that you want to propagate half period for 
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
    print(f"Looking for the initial conditions for half-period with initial guesses starting x = {x_guess} and vy guess = {ydot_guess}")
    x0 = x_guess # Initial x guess 
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
        try:
            tf = prop.t_events[0][1]
        except:
            pdb.set_trace()
        phi_31 = stm[2,0]
        phi_34 = stm[2,3]
        phi_21 = stm[1,0]
        phi_24 = stm[1,3]
        phi_42 = stm[3,1]
        phi_44 = stm[3,3]
        rf = prop.y[0:2,-1].reshape(2,1) # position at tf, turn into 2x1 vector
        vf = prop.y[2:4,-1].reshape(2,1) # velocity at tf, turn into 2x1 vector
        xf = prop.y[0,-1] # get xdot to compare against tolerance
        yf = prop.y[1,-1] # get xdot to compare against tolerance
        vxf = prop.y[2,-1] # get xdot to compare against tolerance
        vyf = prop.y[3,-1] # get ydot to calculate the acceleration
        calced_vyf = calc_velocity_from_Jacobi(JC, mu, xf, yf, 0)
        delta_vyf = vyf - calced_vyf
        # print for info
        # Check if the vxf is close to 0 within acceptable margins
        if (abs(vxf) > tolerance or abs(calced_vyf) > tolerance)and counter <= max_iterations: # If not, recalculate the delta_v0 and try again
            # Calc the acceleration in x
            a_x, a_y = eval_acceleration(rf[0,0], rf[1,0], vxf, vyf, mu)
            # Step 3: Calc new delta_vy0
            term_1 = phi_31 - a_x*phi_21/vyf
            term_2 = phi_34 - a_x*phi_24/vyf
            term_3 = phi_42 - a_y*phi_21/vyf
            term_4 = phi_44 - a_y*phi_24/vyf
            a = np.array([[term_1, term_2],[term_3, term_4]])
            b = np.array([-vxf, delta_vyf])
            sol = np.linalg.solve(a,b) # [delta_x0;delta_vy0]
            delta_x0 = sol[0]
            delta_vy0 = sol[1]
            # Step 4: Use new deltas_ydot_t to calc new ICs
            if delta_vy0 > 0.001:
                delta_vy0 = delta_vy0/10
            if delta_x0 > 0.001:
                delta_x0 = delta_x0/10
            vy0 = vy0 + delta_vy0
            x0 = x0 + delta_x0
            continue
        else: # If error is within acceptable margins, break out of iterative loop
            arrival_states = [rf[0,0], rf[1,0], vf[0,0], vf[1,0]]
            converged_initial_states = [IC_guess[0], IC_guess[1], IC_guess[2], IC_guess[3]]            
            break
    return counter, tf, arrival_states, converged_initial_states