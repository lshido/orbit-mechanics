% Problem C1 part(b)
% Matlab Code: ZVC Plotting Algorithm
% Author: Lillian Shido
% Date: 10/6/2025

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

% Position of primary bodies
x_Earth = -mu;
x_Moon = 1-mu;

% Locations of L1 and L2
x_L1 = 8.3692e-01;
y_L1 = 0;

x_L2 = 1.1557e+00;
y_L2 = 0;

x_L3 = -1.0051;
y_L3 = 0;

x_L4 = cosd(60) - mu;
y_L4 = sind(60);

x_L5 = cosd(60) - mu;
y_L5 = -sind(60);

% Find Jacobi Constant C
function C = JC_check(mu, x, y)
    d = sqrt((x+mu)^2 + y^2);
    r = sqrt((x-1+mu)^2 + y^2);
    C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
end

% Calc error
function error = calc_error(actual, ideal)
    error =  abs(actual - ideal)/ideal;
end

% To plot the ZVC for z = 0:
% Step 1: For a given x
% Step 2: Guess a y
% Step 3: Calculate Jacobi constant result
% Step 4: If Jacobi constant result is less than tolerance, break
% Step 5: If not, update y with Newton Raphson
% Step 6: Repeat

% Calc ZVCs
zvc_result = zeros(0,6);
tolerance = 1e-12;
max_iterations = 100;
% Give an initial guess for y that is around the Earth
% for x = linspace(-1.5,1.5,1e3) %  Find the curve for -1.5 < x < 1.5
for x = linspace(-1.5,1.5,1e3) %  Find the curve for -1.5 < x < 1.5
    fprintf("x sweep: %f\n", x)
    for y = 0:0.1:1.5 % These are my guesses
        counter = 0;
        while 1
            counter = counter + 1;
            d = sqrt((x+mu)^2 + y^2);
            r = sqrt((x-1+mu)^2 + y^2);
            f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
            f_prime_y = 2*y*( 1 - (1-mu)/d^3 - mu/r^3);
            delta = f/f_prime_y;
            if abs(f) > tolerance
                if counter > max_iterations % check for max iterations
                    new_C = JC_check(mu, x, y);
                    C_error = calc_error(new_C, C);
                    if C_error < 1e-12
                        zvc_result(end+1,:) = [x y -y new_C error_C counter];
                    else
                        break
                    end
                elseif d == 0 || r == 0 % check for dividing by 0
                    break
                elseif abs(f_prime_y) <= tolerance % check for zero derivative
                    break
                else
                    y = y - delta;
                    continue
                end
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y -y new_C error_C counter];
                break
            end
        end
    end
end
% for y = linspace(-1.5,1.5,1e3) % Find the curve for 0 < y < 1.21
for y = linspace(-1.5,1.5,1e3) % Find the curve for 0 < y < 1.21
    fprintf("y sweep: %f\n", y)
    if y == 0.752505
        fprintf("y = 0.752505!!\n")
    end
    for x = -1.5:0.1:1.5 % These are my guesses for x
        counter = 0;
        while 1
        counter = counter + 1;
        d = sqrt((x+mu)^2 + y^2);
        r = sqrt((x-1+mu)^2 + y^2);
        f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
        f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d^3 - 2*mu*(x-1+mu)/r^3;
        delta = f/f_prime_x;
        if abs(f) > tolerance
                if counter > max_iterations % check for max iterations
                    new_C = JC_check(mu, x, y);
                    C_error = calc_error(new_C, C);
                    if C_error < 1e-12
                        zvc_result(end+1,:) = [x y -y new_C error_C counter];
                    else
                        break
                    end
                elseif d == 0 || r == 0 % check for dividing by 0
                    break
                elseif abs(f_prime_x) <= tolerance % check for zero derivative
                    break
                else
                    x = x - delta;
                    continue
                end
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y -y new_C error_C counter];
                break
            end
        end
    end
end
zvc_table = array2table(zvc_result, 'VariableNames', {'x','y', '-y', 'JC', '% Error of C','Iterations'});

% Find the max error of JC
max_error = max(abs(zvc_result(:,4)));
fprintf("max error: %d\n", max_error)

% Locations of points when x=0 and y=0
x_1 = 0;
y_1 = 1.2722;

x_2 = 1.1159;
y_2 = 0;

%===================Plot the ZVC!===============
fig1 = figure('Name','ZVC');
earth = scatter(x_Earth, 0, 'blue', 'filled', 'SizeData', 200);
hold on
moon = scatter(x_Moon, 0, 'black', 'filled', 'SizeData', 100);
L1_plot = scatter(x_L1, y_L1, 'red', 'filled', 'SizeData', 10);
scatter(x_L2, y_L2, 'red', 'filled', 'SizeData', 10)
scatter(x_L3, y_L3, 'red', 'filled', 'SizeData', 10)
scatter(x_L4, y_L4, 'red', 'filled', 'SizeData', 10)
scatter(x_L5, y_L5, 'red', 'filled', 'SizeData', 10)
zvc_plot = scatter(zvc_table, 'x', 'y', 'SizeData', 2, 'MarkerFaceColor', '#53A1C9', 'MarkerEdgeColor', '#53A1C9');
scatter(zvc_table, 'x', '-y', 'SizeData', 2, 'MarkerFaceColor', '#53A1C9', 'MarkerEdgeColor', '#53A1C9');
hold off
xlim([-1.5 1.5])
ylim([-1.5 1.5])
axis square
xlabel("Non-dimensional X")
ylabel("Non-dimensional Y")
legend([earth, moon, L1_plot, zvc_plot], {'Earth', 'Moon', 'Eq. Points', 'ZVC'})
title({'The Zero Velocity Curves at z=0';['in the Earth-Moon System for Jacobi Constant = ', num2str(C)]})
box on
grid on
fontsize(14, 'points')