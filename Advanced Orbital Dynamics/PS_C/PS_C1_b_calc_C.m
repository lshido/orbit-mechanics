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

% Try calculating the C for (x,y) pairs
% Step 1: Given a set of x = [-1.5, 1.5] and y = [-1.5, 1.5]
% Step 2: Calculate the C value.
% Step 3: If it's within 1e-12 of the desired C value, save it. 
% Step 4: Else, skip the pair and move on. 

zvc_result = zeros(0,3);
tolerance = 1e-3;


% Step 1
for x = linspace(-1.5, 1.5, 1e4)
    fprintf("x: %d\n", x)
    for y = linspace(-1.5, 1.5, 1e4)
        % Step 3
        d = sqrt((x+mu)^2 + y^2);
        r = sqrt((x-1+mu)^2 + y^2);
        f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
        if abs(f) > tolerance
            continue
        else
            new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
            error_C = abs(new_C - C)/C*100;
            zvc_result(end+1,:) = [x y error_C];
        end
    end
end

zvc_table = array2table(zvc_result, 'VariableNames', {'x', 'y', '% Error of C'});
disp(zvc_table)

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
L1_plot = scatter(x_L1, y_L1, 'red', 'filled', 'SizeData', 50);
L2_plot = scatter(x_L2, y_L2, 'green', 'filled', 'SizeData', 50);
L3_plot = scatter(x_L3, y_L3, 'magenta', 'filled', 'SizeData', 50);
L4_plot = scatter(x_L4, y_L4, 'magenta', 'filled', 'SizeData', 50);
L5_plot = scatter(x_L5, y_L5, 'magenta', 'filled', 'SizeData', 50);
p1_plot = scatter(x_1, y_1, 'magenta', 'filled', 'SizeData', 50);
p2_plot = scatter(x_2, y_2, 'magenta', 'filled', 'SizeData', 50);
scatter(zvc_table, 'x', 'y')
hold off
xlim([-1.5 1.5])
ylim([-1.5 1.5])
axis square
xlabel("Non-dimensional X")
ylabel("Non-dimensional Y")
box on
grid on

