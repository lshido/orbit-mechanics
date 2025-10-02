clear all

% Calculate C for the given position and velocity.
% Gravitational Parameters [km^3/s^2]
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

mu = mu_Moon/(mu_Earth + mu_Moon);
fprintf("mu %d\n", mu)

% mu_Earth_Moon = 1.2151e-02;

% position and velocity in NON-DIMENSIONAL units
r_vector = [-0.270 -0.420];
v_vector = [0.300 -1.000];

x_0 = r_vector(1);
y_0 = r_vector(2);
d_0 = sqrt((x_0+mu)^2 + y_0^2);
r_0 = sqrt((x_0-1+mu)^2 + y_0^2);
pseudo_U = (1-mu)/d_0 + mu/r_0 + (x_0^2+y_0^2)/2;
v_squared = norm(v_vector)^2;
C = 2*pseudo_U - v_squared;
fprintf("Jacobi Constant C: %d\n", C)

% To plot the ZVC for z = 0:
% Step 1: For a given x
% Step 2: Guess a y
% Step 3: Calculate Jacobi constant result
% Step 4: If Jacobi constant result is less than tolerance, break
% Step 5: If not, update y with Newton Raphson
% Step 6: Repeat


% Identify where it hits the curve when y=0
x_result = zeros(0,4);
tolerance = 1e-12;
y_0 = 0;
x = 1;
counter = 0;
while 1
    counter = counter + 1;
    d = sqrt((x+mu)^2 + y_0^2);
    r = sqrt((x-1+mu)^2 + y_0^2);
    f = x^2 + y_0^2 + (2*(1-mu)/d) + (2*mu/r) - C;
    f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d^3 - 2*mu*(x-1+mu)/r^3;
    delta = f/f_prime_x;
    % Step 4:
    if abs(f) > tolerance
        x = x - delta;
        continue
    else
        fprintf("current x: %d, current y_0: %d\n", x, y_0)
        new_C = x^2 + y_0^2 + (2*(1-mu)/d) + (2*mu/r);
        error_C = abs(new_C - C)/C*100;
        x_result(end+1,:) = [x y_0 error_C counter];
        counter = 0;
        break
    end
end
x_table = array2table(x_result, 'VariableNames', {'x','y', '% Error of C','Iterations'});
disp(x_table)

% Identify where it hits the curve when x=0
y_result = zeros(0,4);
tolerance = 1e-12;
% Step 1
x_0 = 0;
y=1;
% Step 3
counter = 0;
while 1
    counter = counter + 1;
    d = sqrt((x_0+mu)^2 + y^2);
    r = sqrt((x_0-1+mu)^2 + y^2);
    f = x_0^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
    f_prime_y = 2*y*( 1 - (1-mu)/d^3 - mu/r^3);
    delta = f/f_prime_y;
    % Step 4:
    if abs(f) > tolerance
        y = y - delta;
        continue
    else
        fprintf("current x_0: %d, current y: %d\n", x_0, y)
        new_C = x_0^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
        error_C = abs(new_C - C)/C*100;
        y_result(end+1,:) = [x_0 y error_C counter];
        counter = 0;
        break
    end
end

y_table = array2table(y_result, 'VariableNames', {'x','y', '% Error of C','Iterations'});
disp(y_table)

