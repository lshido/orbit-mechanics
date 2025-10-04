% Script for modeling the orbit of a spacecraft and the corresponding planar ZVCs
% Problem Set C1 (c)
% Author: Lillian Shido
% Date: October 3, 2025

clear all

% Calculate C for the given position and velocity.
% Gravitational Parameters [km^3/s^2]
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

mu = mu_Moon/(mu_Earth + mu_Moon);
fprintf("mu %d\n", mu)

% Characteristic Length [km]
a_Moon = 384400; % around Earth
l_char = a_Moon;

% Calculate characteristic time
t_char = sqrt(l_char^3/(mu_Earth+mu_Moon));
fprintf("characteristic time: %d sec\n", t_char)

% position and velocity in NON-DIMENSIONAL units
r_vector = [-0.270 -0.420];
v_vector = [0.300 -1.018072];

% Calculate the Jacobi Constant for these ICs. 
% It should result in C = 3.2!!
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

% Locations of Libration Points
x_L1 = 8.3692e-01; % From Problem Set B
y_L1 = 0;

x_L2 = 1.1557e+00; % From Problem Set B
y_L2 = 0;

x_L3 = -1.0051; % From Problem Set B
y_L3 = 0;

x_L4 = cosd(60) - mu;
y_L4 = sind(60);

x_L5 = cosd(60) - mu;
y_L5 = -sind(60);

% To plot the ZVC for z = 0:
% Step 1: For a given x
% Step 2: Guess a y
% Step 3: Calculate Jacobi constant result
% Step 4: If Jacobi constant result is less than tolerance, break
% Step 5: If not, update y with Newton Raphson
% Step 6: Repeat

zvc_result = zeros(0,6);
tolerance = 1e-12;

%===================Calculate the ZVC around the Earth=====================
% Give an initial guess for y that's around the Earth)
for y = 0.1
    for x = [linspace(0,-1.21,1e3), linspace(0,1.21, 1e3)] %  Find the curve for -1.21 < x < 0
        counter = 0;
        while 1
            counter = counter + 1;
            d = sqrt((x+mu)^2 + y^2);
            r = sqrt((x-1+mu)^2 + y^2);
            f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
            f_prime_y = 2*y*( 1 - (1-mu)/d^3 - mu/r^3);
            delta = f/f_prime_y;
            if abs(f) > tolerance
                y = y - delta;
                continue
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y -y new_C error_C counter];
                counter = 0;
                break
            end
        end
    end
end

% Give an initial guess for x that's near the Earth)
for x = [-0.1, 0.1] % Try starting from the "outside" of the outer boundary
    for y = [linspace(0,-1.21,1e3), linspace(0,1.21,1e3)] % Find the curve for 0 < y < 1.21
        while 1
            counter = counter + 1;
            d = sqrt((x+mu)^2 + y^2);
            r = sqrt((x-1+mu)^2 + y^2);
            f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
            f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d^3 - 2*mu*(x-1+mu)/r^3;
            delta = f/f_prime_x;
            % Step 4:
            if abs(f) > tolerance
                x = x - delta;
                continue
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y -y new_C error_C counter];
                counter = 0;
                break
            end
        end
    end
end
%===================Calculate the ZVC around the Moon==================
% Give an initial guess for y that's around the Moon
for y = 0.01
    for x = [linspace(x_Moon-(1e-3), x_L1, 1e3), linspace(x_Moon, x_L2, 1e3)] %  Find the curve for -1.21 < x < 0
        counter = 0;
        while 1
            counter = counter + 1;
            d = sqrt((x+mu)^2 + y^2);
            r = sqrt((x-1+mu)^2 + y^2);
            f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
            f_prime_y = 2*y*( 1 - (1-mu)/d^3 - mu/r^3);
            delta = f/f_prime_y;
            if abs(f) > tolerance
                y = y - delta;
                continue
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y -y new_C error_C counter];
                fprintf("x: %d, y: %d\n", x, y)
                counter = 0;
                break
            end
        end
    end
end
% Give an initial guess for x that's around the moon
for x = x_Moon+1e-3;
    for y = [linspace(0,0.0999,1e3), linspace(0,-0.0999,1e3)] % Find the curve for 0 < y < 1.21
        counter = 0;
        while 1
            counter = counter + 1;
            d = sqrt((x+mu)^2 + y^2);
            r = sqrt((x-1+mu)^2 + y^2);
            f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
            f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d^3 - 2*mu*(x-1+mu)/r^3;
            delta = f/f_prime_x;
            if abs(f) > tolerance
                x = x - delta;
                continue
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y -y new_C error_C counter];
                counter = 0;
                break
            end
        end
    end
end
%===================Calculate the "outside" curve=====================
% Give an initial guess for y that starts "outside"
for y = 1.5 % Try starting from the "outside" of the outer boundary
    for x = [linspace(-1.21,0,1e3), linspace(0,1.21,1e3)] %  Find the curve for -1.21 < x < 0
        counter = 0;
        while 1
            counter = counter + 1;
            d = sqrt((x+mu)^2 + y^2);
            r = sqrt((x-1+mu)^2 + y^2);
            f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
            f_prime_y = 2*y*( 1 - (1-mu)/d^3 - mu/r^3);
            delta = f/f_prime_y;
            if abs(f) > tolerance
                y = y - delta;
                continue
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y -y new_C error_C counter];
                counter = 0;
                break
            end
        end
    end
end

% % Give an initial guess for x that starts "outside"
for x = [-1.5, 1.5]; % Try starting from the "outside" of the outer boundary
    for y = [linspace(-1.21,0,1e3), linspace(0,1.21,1e3)] % Find the curve for 0 < y < 1.21
        while 1
            counter = counter + 1;
            d = sqrt((x+mu)^2 + y^2);
            r = sqrt((x-1+mu)^2 + y^2);
            f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
            f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d^3 - 2*mu*(x-1+mu)/r^3;
            delta = f/f_prime_x;
            if abs(f) > tolerance
                x = x - delta;
                continue
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y -y new_C error_C counter];
                counter = 0;
                break
            end
        end
    end
end

zvc_table = array2table(zvc_result, 'VariableNames', {'x','y', '-y', 'JC', '% Error of C','Iterations'});
format short
disp(zvc_table)

% Find the max error of JC
max_error = max(abs(zvc_result(:,5)));
fprintf("max error: %d\n", max_error)

%==================================Propagate the Orbit=====================================
% Initial conditions
w0 = [r_vector(1);r_vector(2);v_vector(1);v_vector(2)];

% Set the final time in rads (non-dimensional time)
t_final = 100*pi;
fprintf("Non-dimensional Time: %d rad\n", t_final)
% Calculate the dimensional time
time_dim = t_final*t_char;
fprintf("Dimensional Time: %d sec\n", time_dim)
fprintf("Dimensional Time: %d years\n", time_dim/3600/24/365)

% Set the span of the integrator
tspan = [0 t_final];

% Set up the ODE
ode = @(t,w) [...
w(3);...
w(4);...
2*w(4) + w(1) - (1 - mu) * (w(1) + mu) / ((w(1) + mu)^2 + w(2)^2)^(3/2)...
- mu * (w(1) - 1 + mu) / ((w(1) - 1 + mu)^2 + w(2)^2)^(3/2);...
-2*w(3) + w(2) - (1 - mu) * w(2) / ((w(1) + mu)^2 + w(2)^2)^(3/2) - mu * w(2)/((w(1) - 1 + mu)^2 + w(2)^2)^(3/2)...
];
% Configure the tolerances
options = odeset('RelTol',1e-12,'AbsTol', 1e-14);
% Run the integrator!
[t,w] = ode45(ode, tspan, w0, options); 
% Pull out the values of the solution
x = w(:,1);
y = w(:,2);
v_x = w(:,3);
v_y = w(:,4);

%===================Plot!===============
fig1 = figure('Name','ZVC');
% Plot the Earth, Moon, and Libration Points in the system
earth = scatter(x_Earth, 0, 'blue', 'filled', 'SizeData', 200);
hold on
moon = scatter(x_Moon, 0, 'black', 'filled', 'SizeData', 100);
L1_plot = scatter(x_L1, y_L1, 'red', 'filled', 'SizeData', 50);
scatter(x_L2, y_L2, 'red', 'filled', 'SizeData', 50)
scatter(x_L3, y_L3, 'red', 'filled', 'SizeData', 50)
scatter(x_L4, y_L4, 'red', 'filled', 'SizeData', 50)
scatter(x_L5, y_L5, 'red', 'filled', 'SizeData', 50)
% Plot the ZVC
zvc_plot = scatter(zvc_table, 'x', 'y', 'SizeData', 10, 'MarkerFaceColor','#53A1C9', 'MarkerEdgeColor', '#53A1C9');
scatter(zvc_table, 'x', '-y', 'SizeData', 10, 'MarkerFaceColor','#53A1C9', 'MarkerEdgeColor', '#53A1C9');
% Plot the spacecraft trajectory
sc = plot(x,y);
sc.LineWidth = 1;
hold off
xlim([-1.5 1.5])
ylim([-1.5 1.5])
axis square
xlabel("Non-dimensional X")
ylabel("Non-dimensional Y")
legend([sc, earth, moon, L1_plot, zvc_plot], {'SC Trajectory', 'Earth', 'Moon', 'Eq. Points', 'ZVC'})
title({'The trajectory of a spacecraft in the Earth-Moon System';['for Jacobi Constant = ', num2str(C), ' orbiting for 3.73 years']})
box on
grid on
fontsize(14, 'points')