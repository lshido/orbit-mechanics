clear all

% Gravitational Parameters [km^3/s^2]
mu_Sun = 132712440017.99;
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

% Earth-Moon System
mu = mu_Moon/(mu_Earth + mu_Moon);
fprintf("Mu of Earth-Moon: %f\n", mu)
% Position of primary bodies
x_Earth = -mu;
x_Moon = 1-mu;
% Location of Earth-Moon L1 point
x_L1 = 8.3692e-01; % From Problem Set B
y_L1 = 0;
% Characteristic Length [km]
a_Moon = 384400; % around Earth
l_char = a_Moon;
fprintf("Characteristic Length: %f km \n", l_char)
% Calculate characteristic time
t_char = sqrt(a_Moon^3/(mu_Earth+mu_Moon));
fprintf("characteristic time: %d sec\n", t_char)

% Location of the other Earth-Moon Libration points
x_L2 = 1.1557e+00;
y_L2 = 0;

x_L3 = -1.0051;
y_L3 = 0;

x_L4 = cosd(60) - mu;
y_L4 = sind(60);

x_L5 = cosd(60) - mu;
y_L5 = -sind(60);

% Perturbations
xi_0 = 0.01;
xi_0_dim = xi_0*l_char;
fprintf("Dimensional xi_0: %f km\n", xi_0_dim)
eta_0 = 0;
x_0 = x_L1 + xi_0;
y_0 = 0 + eta_0;

% The solution to the perturbed orbit:
t_0 = 0; % non-dim time
t_end = 8*pi; % non-dim time
orbit_results = zeros(0,3);
% d, r, and partials all evaluated @ L1
d = sqrt((x_L1+mu)^2 + y_L1^2);
r = sqrt((x_L1-1+mu)^2 + y_L1^2);
U_xx = 1 - (1-mu)/d^3 - mu/r^3 + 3*(1-mu)*(x_L1+mu)^2/d^5 + 3*mu*(x_L1-1+mu)^2/r^5;
fprintf("Uxx: %f\n", U_xx)
U_yy = 1 - (1-mu)/d^3 - mu/r^3;
B_1 = 2 - (U_xx + U_yy)/2;
B_2_squared = -U_xx*U_yy;
s = (B_1 + (B_1^2 + B_2_squared)^(1/2))^(1/2);
B_3 = (s^2 + U_xx)/(2*s);
big_Lambda_1 = -B_1 + (B_1^2 + B_2_squared)^(1/2);
big_Lambda_2 = -B_1 - (B_1^2 + B_2_squared)^(1/2);
lambda_1 = sqrt(big_Lambda_1); % real
lambda_2 = -sqrt(big_Lambda_1); % real
lambda_3 = sqrt(big_Lambda_2); % imaginary
lambda_4 = -sqrt(big_Lambda_2); % imaginary
for t = t_0:0.1:t_end
    xi = xi_0*cos(s*(t-t_0)) + (eta_0/B_3)*sin(s*(t-t_0));
    x = xi + x_L1;
    eta = eta_0*cos(s*(t-t_0)) - B_3*xi_0*sin(s*(t-t_0)); 
    y = eta + y_L1; 
    orbit_results(end+1,:) = [t x y];
end
orbit_table = array2table(orbit_results, 'VariableNames', {'time', 'x', 'y'});
disp(t)

% Calculate the roots of the L1 point
fprintf("x_L: %f\n", x)
fprintf("lambda_1: %f + %fi\n", real(lambda_1), imag(lambda_1))
fprintf("lambda_2: %f + %fi\n", real(lambda_2), imag(lambda_2))
fprintf("lambda_3: %f + %fi\n", real(lambda_3), imag(lambda_3))
fprintf("lambda_4: %f + %fi\n", real(lambda_4), imag(lambda_4))

% Max distance from L1
% Where y is maximum
y_max = eta_0*cosd(90) - B_3*xi_0*sind(90) + y_L1;
fprintf('Max distance from L1 (non-dim): %f\n', y_max)
fprintf('Max distance from L1 (dim): %f\n', y_max*l_char)

% Period of the orbit around L1
P_around_L1 = 2*pi/s;
fprintf('Period around L1 (non-dim relative to the larger system): %f\n', P_around_L1)
fprintf('Period around L1 (dim): %f sec\n', P_around_L1*t_char)
fprintf('Period around L1 (dim): %f days\n', P_around_L1*t_char/3600/24)

% Initial velocities
v_x0 = eta_0*s/B_3;
v_y0 = -B_3*xi_0*s;
fprintf("initial velocity x-comp: %f\n", v_x0)
fprintf("initial velocity y-comp: %f\n", v_y0)
v_x0_dim = v_x0 * l_char / t_char * 1000;
v_y0_dim = v_y0 * l_char / t_char * 1000;
fprintf("initial velocity x-comp (dim): %f m/s\n", v_x0_dim)
fprintf("initial velocity y-comp (dim): %f m/s\n", v_y0_dim)

% Calculate the Jacobi Constant for these ICs. 
x_0 = x_L1 + xi_0;
y_0 = 0;
d_0 = sqrt((x_0+mu)^2 + y_0^2);
r_0 = sqrt((x_0-1+mu)^2 + y_0^2);
pseudo_U = (1-mu)/d_0 + mu/r_0 + (x_0^2+y_0^2)/2;
v_squared = v_x0^2 + v_y0^2;
C = 2*pseudo_U - v_squared;
fprintf("d_0: %f\n", d_0)
fprintf("r_0: %f\n", r_0)
fprintf("pseudo_U: %f\n", pseudo_U)
fprintf("v_squared: %f\n", v_squared)
fprintf("Jacobi Constant C: %d\n", C)

% To plot the ZVC for z = 0:
% Step 1: For a given x
% Step 2: Guess a y
% Step 3: Calculate Jacobi constant result
% Step 4: If Jacobi constant result is less than tolerance, break
% Step 5: If not, update y with Newton Raphson
% Step 6: Repeat

zvc_result = zeros(0,5);
tolerance = 1e-12;

%===================Calculate ZVC====================
% Give an initial guess for y that is "inside" %

for y = [-0.1, 0.1];
    for x = [linspace(0,-1.21,5e3), linspace(0,1.21,5e3)] %  Find the curve for -1.21 < x < 0
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
                zvc_result(end+1,:) = [x y new_C error_C counter];
                counter = 0;
                break
            end
        end
    end
end

% Give an initial guess for x that starts "inside"
for x = [-0.1, 0.1]; % Try starting from the "outside" of the outer boundary
    for y = [linspace(0,-1.21,5e3), linspace(0,1.21,5e3)] % Find the curve for 0 < y < 1.21
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
            % elseif counter > 1e3
            %     break
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y new_C error_C counter];
                counter = 0;
                break
            end
        end
    end
end

%===================Calculate the "outside" curve=====================%
% Give an initial guess for y that starts "outside"
for y = [-1.5, 1.5]; % Try starting from the "outside" of the outer boundary
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
                zvc_result(end+1,:) = [x y new_C error_C counter];
                counter = 0;
                break
            end
        end
    end
end

% Give an initial guess for x that starts "outside"
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
            % elseif counter > 1e3
            %     break
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y new_C error_C counter];
                counter = 0;
                break
            end
        end
    end
end

zvc_table = array2table(zvc_result, 'VariableNames', {'x','y', 'JC', '% Error of C','Iterations'});
format short
disp(zvc_table)

% Find the max error of JC
max_error = max(abs(zvc_result(:,4)));
fprintf("max error: %d\n", max_error)

%===================Plot the orbit===============
fig1 = figure("Name", "orbit");
earth = scatter(x_Earth, 0, 'blue', 'filled', 'SizeData', 200);
hold on
moon = scatter(x_Moon, 0, 'black', 'filled', 'SizeData', 100);
L1_plot = scatter(x_L1, y_L1, 'red', 'filled', 'SizeData', 10);
orbit = plot(orbit_table, 'x', 'y');
hold off
xlim([-0.5 1.5])
ylim([-1 1])
axis square
xlabel("Non-dimensional X")
ylabel("Non-dimensional Y")
legend([earth, moon, L1_plot, orbit], {'Earth', 'Moon', 'L1', 'Orbit'})
title({'The orbit around the L1 point in the Earth-Moon System'})
box on
grid on
fontsize(14, 'points')

%===================Plot the ZVC!===============%
fig1 = figure('Name','ZVC');
earth = scatter(x_Earth, 0, 'blue', 'filled', 'SizeData', 200);
hold on
moon = scatter(x_Moon, 0, 'black', 'filled', 'SizeData', 100);
L1_plot = scatter(x_L1, y_L1, 'red', 'filled', 'SizeData', 20);
scatter(x_L2, y_L2, 'red', 'filled', 'SizeData', 20)
scatter(x_L3, y_L3, 'red', 'filled', 'SizeData', 20)
scatter(x_L4, y_L4, 'red', 'filled', 'SizeData', 20)
scatter(x_L5, y_L5, 'red', 'filled', 'SizeData', 20)
zvc_plot = scatter(zvc_table, 'x', 'y', 'filled', 'SizeData', 5);
orbit = plot(orbit_table, 'x', 'y');
% plot(zvc_table, 'x', 'y')
hold off
xlim([-1.5 1.5])
ylim([-1.5 1.5])
axis square
xlabel("X [non-dim]")
ylabel("Y [non-dim]")
legend([earth, moon, L1_plot, zvc_plot, orbit], {'Earth', 'Moon', 'Eq. Points', 'ZVC', 'Orbit'})
title({'The Zero Velocity Curves at z=0';['in the Earth-Moon System for Jacobi Constant = ', num2str(C)]})
box on
grid on
fontsize(14, 'points')
