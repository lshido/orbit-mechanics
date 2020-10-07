R_earth = 6378.1363; % [km]

% Solve for a and e
r_p = 1.5*R_earth;
r_a = 200*R_earth;
fprintf("radius at periapsis: %.4e km\n", r_p);
fprintf("radius at apoapsis: %.4e km\n", r_a);

syms e a
equations = [a*(1-e) == r_p, a*(1+e) == r_a];
variables = [a e];
[solved_a, solved_e] = solve(equations, variables);
fprintf("semi-major axis: %.4e km\n", solved_a);
fprintf("eccentricity: %.4e\n", solved_e);