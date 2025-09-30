% Use an integrator to solve the rotating frame EOMs
clear all
% Initial conditions
r0_vector = [-0.270 -0.420];
v0_vector = [0.300 -1.000];
w0 = [r0_vector(1);r0_vector(2);v0_vector(1);v0_vector(2)];

% Gravitational Parameters [km^3/s^2]
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

% mu = mu_Moon/(mu_Earth + mu_Moon);
mu = 1.2151e-02;

% Characteristic Length [km]
a_Moon = 384400; % around Earth
l_char = a_Moon;

% Calculate characteristic time
t_char = sqrt(l_char^3/(mu_Earth+mu_Moon));
fprintf("characteristic time: %d sec\n", t_char)

% Calculate dimensional time
mean_motion = sqrt((mu_Earth+mu_Moon)/(a_Moon)^3);
fprintf("Mean Motion: %d rad/sec\n", mean_motion)
Period = 2*pi/mean_motion;
fprintf("Period: %d sec\n", Period)
fprintf("Period: %d days\n", Period/3600/24)


% Locations of L1 and L2
x_L1 = 8.3692e-01;
y_L1 = 0;

x_L2 = 1.1557e+00;
y_L2 = 0;

x_L3 = -1.0051;
y_L3 = 0;

% tspan should be calculated with non-dim time
t_final = 1000*pi;

fprintf("Non-dimensional Time: %d rad\n", t_final)
% For non-dimensional time divide time_dim by characteristic time
time_dim = t_final*t_char;
fprintf("Dimensional Time: %d sec\n", time_dim)
fprintf("Dimensional Time: %d years\n", time_dim/3600/24/365)

tspan = [0 t_final];

ode = @(t,w) [...
w(3);...
w(4);...
2*w(4) + w(1) - (1 - mu) * (w(1) + mu) / ((w(1) + mu)^2 + w(2)^2)^(3/2)...
- mu * (w(1) - 1 + mu) / ((w(1) - 1 + mu)^2 + w(2)^2)^(3/2);...
-2*w(3) + w(2) - (1 - mu) * w(2) / ((w(1) + mu)^2 + w(2)^2)^(3/2) - mu * w(2)/((w(1) - 1 + mu)^2 + w(2)^2)^(3/2)...
];

options = odeset('RelTol',1e-8,'AbsTol', 1e-10);
[t,w] = ode45(ode, tspan, w0, options); 

x = w(:,1);
y = w(:,2);
v_x = w(:,3);
v_y = w(:,4);

d = sqrt((x+mu).^2 + y.^2);
r = sqrt((x-1+mu).^2 + y.^2);
x_y_sq = (x.^2+y.^2)./2;
term_1 = (1-mu)./d;
term_2 = mu./r;
pseudo_U = term_1 + term_2 + x_y_sq;
v_squared = w(:,3).^2 + w(:,4).^2;
C = 2*pseudo_U - v_squared;
C_correct = C(1,1);
C_diff = abs(C - C_correct)/C_correct*100;
% C_table = array2table(C)
% format longE
% disp(C_table)

fig2 = figure('Name', 'Jacobi C');
c_plot = scatter(t, C_diff);
hold on

% ode = @(t,w) [...
% w(3);...
% w(4);...
% 2*w(4) + w(1) - (1 - mu) * (w(1) + mu) / ((w(1) + mu)^2 + w(2)^2)^(3/2)...
% - mu * (w(1) - 1 + mu) / ((w(1) - 1 + mu)^2 + w(2)^2)^(3/2);...
% -2*w(3) + w(2) - (1 - mu) * w(2) / ((w(1) + mu)^2 + w(2)^2)^(3/2) - mu * w(2)/((w(1) - 1 + mu)^2 + w(2)^2)^(3/2)...
% ];

% options = odeset('RelTol',1e-12,'AbsTol', 1e-14);
% [t,w] = ode45(ode, tspan, w0, options); 

% x = w(:,1);
% y = w(:,2);
% v_x = w(:,3);
% v_y = w(:,4);

% d = sqrt((x+mu).^2 + y.^2);
% r = sqrt((x-1+mu).^2 + y.^2);
% x_y_sq = (x.^2+y.^2)./2;
% term_1 = (1-mu)./d;
% term_2 = mu./r;
% pseudo_U = term_1 + term_2 + x_y_sq;
% v_squared = w(:,3).^2 + w(:,4).^2;
% C = 2*pseudo_U - v_squared;
% C_correct = C(1,1);
% C_diff_2 = abs(C - C_correct)/C_correct*100;


% c2_plot = scatter(t, C_diff_2);


hold off
xlabel("Non-dimensional Time")
ylabel("Error Percentage")
% legend([c_plot, c2_plot], {"Abs. Tol: 1e-10, Rel. Tol: 1e-12", "Abs. Tol: 1e-14, Rel. Tol: 1e-16"})
legend([c_plot], {"Abs. Tol: 1e-12, Rel. Tol: 1e-14"})
title({"Jacobi Constant Error over time";"Problem B4(d.v), Lillian Shido"})
fontsize(14, "points")
%=============ANIMATE THE POSITIONS!==========================
% for k = 1:numel(t)
%     sc_marker.XData = w(k,1);
%     sc_marker.YData = w(k,2);
%     drawnow
%     % exportgraphics(gca,"parabola.gif",Append=true)
% end