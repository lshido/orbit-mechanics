clear all 

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

%=====================DEFINE FUNCTIONS================================

% Event to detect with trajectory crosses x-axis
function [value, isterminal, direction] = crossxEvent(t, sv)
    value = sv(2); % y-value
    isterminal = 1; % stop the integration
    direction = 0; % from either direction
    % you can add code here to do something when this event occurs like updating the ICs
end

% Calc error
function error = calc_error(actual, ideal)
    error =  abs(actual - ideal)/ideal;
end


%==================END DEFINE FUNCTIONS===============================

% Set the span of the integrator
t_final = 1.5*pi;
tspan = [0 t_final];

% position and velocity in NON-DIMENSIONAL units
r_vector = [0.488 0.200];
v_vector = [-0.880 0.200];

%====================DEFINE SYSTEM====================================



%============================INTEGRATE THE EOM + STM=================================

% Set up the ODEs
function d_sv = odefun(t,sv,mu)

    d_sv = zeros(20,1);

    % EOM ODEs
    d_sv(1) = sv(3);
    d_sv(2) = sv(4);
    d_sv(3) = 2*sv(4) + sv(1) - (1 - mu) * (sv(1) + mu) / ((sv(1) + mu)^2 + sv(2)^2)^(3/2)...
- mu * (sv(1) - 1 + mu) / ((sv(1) - 1 + mu)^2 + sv(2)^2)^(3/2);
    d_sv(4) = -2*sv(3) + sv(2) - (1 - mu) * sv(2) / ((sv(1) + mu)^2 + sv(2)^2)^(3/2) - mu * sv(2)/((sv(1) - 1 + mu)^2 + sv(2)^2)^(3/2);
    
    % Calc the partials using the current x and y values
    d = sqrt((sv(1)+mu)^2 + sv(2)^2);
    r = sqrt((sv(1)-1+mu)^2 + sv(2)^2);
    U_xx = 1 - (1-mu)/d^3 - mu/r^3 + 3*(1-mu)*(sv(1)+mu)^2/d^5 + 3*mu*(sv(1)-1+mu)^2/r^5;
    U_yy = 1 - (1-mu)/d^3 - mu/r^3 + 3*(1-mu)*sv(2)^2/d^5 + 3*mu*sv(2)^2/r^5;
    U_xy = 3*(1-mu)*(sv(1)+mu)*sv(2)/d^5 + 3*mu*(sv(1)-1+mu)*sv(2)/r^5;
    
    % STM ODEs
    d_sv(5) = sv(13);
    d_sv(6) = sv(14);
    d_sv(7) = sv(15);
    d_sv(8) = sv(16);
    d_sv(9) = sv(17);
    d_sv(10) = sv(18);
    d_sv(11) = sv(19);
    d_sv(12) = sv(20);
    d_sv(13) = U_xx*sv(5) + U_xy*sv(9) + 2*sv(17);
    d_sv(14) = U_xx*sv(6) + U_xy*sv(10) + 2*sv(18);
    d_sv(15) = U_xx*sv(7) + U_xy*sv(11) + 2*sv(19);
    d_sv(16) = U_xx*sv(8) + U_xy*sv(12) + 2*sv(20);
    d_sv(17) = U_xy*sv(5) + U_yy*sv(9) - 2*sv(13);
    d_sv(18) = U_xy*sv(6) + U_yy*sv(10) - 2*sv(14);
    d_sv(19) = U_xy*sv(7) + U_yy*sv(11) - 2*sv(15);
    d_sv(20) = U_xy*sv(8) + U_yy*sv(12) - 2*sv(16);
end

% sv0 = [r_vector(1);r_vector(2);v_vector(1);v_vector(2);1;0;0;0;0;1;0;0;0;0;1;0;0;0;0;1];
sv0 = transpose([1:20]);
% set up the STM ODE
options = odeset('Events', @(t,sv) crossxEvent(t,sv), 'RelTol',1e-12,'AbsTol', 1e-14);
[t,sv, te, sve, ie] = ode45(@(t,sv) odefun(t,sv,mu), tspan, sv0, options);

% Pull out the values of the solution
x = sv(:,1);
y = sv(:,2);
v_x = sv(:,3);
v_y = sv(:,4);

error = calc_error(sv(2), 1e-12);

d_sv = odefun(t, sv0, mu);

%==========================END INTEGRATE THE EOM + STM===============================

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
xlim([-1.5 1.5])
ylim([-1.5 1.5])
% xticks(-limit:.01:limit)
legend([earth, moon, nonlinear_orbit], {'Earth', 'Moon', 'Nonlinear Orbit'})
title({'Trajectory in Earth-Moon (x-y) (Lillian Shido)'})
box on
grid on
fontsize(14, 'points')
%=======================End Configure Plot============================