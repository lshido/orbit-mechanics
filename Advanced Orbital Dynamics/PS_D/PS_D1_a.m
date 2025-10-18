clear all

%================================DEFINE FUNCTIONS=====================

% Event to detect with trajectory crosses x-axis
function [value, isterminal, direction] = crossxEvent(t, w)
    value = w(2); % y-value
    isterminal = 1; % stop the integration
    direction = 0; % from either direction
    % you can add code here to do something when this event occurs like updating the ICs
end

% Calc error
function error = calc_error(actual, ideal)
    error =  abs(actual - ideal)/ideal;
end

%============================END DEFINE FUNCTIONS=====================

%====================DEFINE SYSTEM====================================
% Gravitational Parameters [km^3/s^2]
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

mu = mu_Moon/(mu_Earth + mu_Moon);

% Position of primary bodies
x_Earth = -mu;
x_Moon = 1-mu;

% Characteristic Length [km]
a_Moon = 384400; % around Earth
l_char = a_Moon;

% Calculate characteristic time
t_char = sqrt(l_char^3/(mu_Earth+mu_Moon));

% position and velocity in NON-DIMENSIONAL units
r_vector = [0.488 0.200];
v_vector = [-0.880 0.200];

%====================DEFINE SYSTEM====================================

%====================NUMERICAL INTEGRATOR=============================
% Initial conditions
w0 = [r_vector(1);r_vector(2);v_vector(1);v_vector(2)];

% Set the final time in rads (non-dimensional time)
t_final = 100*pi;
% Calculate the dimensional time
time_dim = t_final*t_char;

% Set the span of the integrator
tspan = [0 t_final];

% Set up the x ODE
ode = @(t,w) [...
w(3);...
w(4);...
2*w(4) + w(1) - (1 - mu) * (w(1) + mu) / ((w(1) + mu)^2 + w(2)^2)^(3/2)...
- mu * (w(1) - 1 + mu) / ((w(1) - 1 + mu)^2 + w(2)^2)^(3/2);...
-2*w(3) + w(2) - (1 - mu) * w(2) / ((w(1) + mu)^2 + w(2)^2)^(3/2) - mu * w(2)/((w(1) - 1 + mu)^2 + w(2)^2)^(3/2)...
];

% Configure the tolerances
options = odeset('Events', @(t, w) crossxEvent(t, w), 'RelTol',1e-12,'AbsTol', 1e-14);
% Run the integrator!
[t,w, te, we, ie] = ode45(ode, tspan, w0, options); 
% Pull out the values of the solution
x = w(:,1);
y = w(:,2);
v_x = w(:,3);
v_y = w(:,4);

% Calc error
error = calc_error(we(2), 1e-12);
%=================END NUMERICAL INTEGRATOR============================

%====================PRINT IMPORTANT NUMBERS==========================
fprintf("mu %d\n", mu)
fprintf("characteristic time: %f sec\n", t_char)
fprintf("characteristic length: %f sec\n", l_char)
fprintf("Initial x: %f km\n", r_vector(1)*l_char)
fprintf("Initial y: %f km\n", r_vector(2)*l_char)
fprintf("Initial v_x: %f m/s\n", v_vector(1)*l_char/t_char*1000)
fprintf("Initial v_y: %f m/s\n", v_vector(2)*l_char/t_char*1000)
fprintf("Non-dimensional event time: %d\n", te)
fprintf("Dimensional event time: %d sec\n", te*t_char)
fprintf("Dimensional event time: %d days\n", te*t_char/3600/24)
fprintf("Error of x-axis cross: %f\n", error)
%================END PRINT IMPORTANT NUMBERS==========================

%=====================Configure Plot==================================
fig1 = figure('Name','Traj');
earth = scatter(x_Earth, 0, 'blue', 'filled', 'SizeData', 200);
hold on
moon = scatter(x_Moon, 0, 'black', 'filled', 'SizeData', 100);
nonlinear_orbit = plot(x, y, 'Color', '#008000');
hold off
axis square
xlabel("x [non-dim]")
ylabel("y [non-dim]")
% xlim([-10e-3 5e-3])
% ylim([-4e-3 11e-3])
% xticks(-limit:.01:limit)
legend([earth, moon, nonlinear_orbit], {'Earth', 'Moon', 'Nonlinear Orbit'})
title({'Trajectory in Earth-Moon (x-y) (Lillian Shido)'})
box on
grid on
fontsize(14, 'points')
%=======================End Configure Plot============================