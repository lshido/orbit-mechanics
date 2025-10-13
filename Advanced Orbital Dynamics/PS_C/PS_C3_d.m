clear all

%============================DEFINE FUNCTIONS=============================
% Calculate the location of L4
function [x_L4, y_L4] = calc_L4(mu)
    x_L4 = cosd(60) - mu;
    y_L4 = sind(60);
end

% Calculate eigenvalues for equilateral points:
function [lambda_1, lambda_2, lambda_3, lambda_4] = calc_eigenvalues_equilateral(mu)
    big_lambda_1 = 1/2*(-1 + (1 - 27*mu*(1-mu))^(1/2));
    big_lambda_2 = 1/2*(-1 - (1 - 27*mu*(1-mu))^(1/2));
    lambda_1 = sqrt(big_lambda_1);
    lambda_2 = -lambda_1;
    lambda_3 = sqrt(big_lambda_2);
    lambda_4 = -lambda_3;
end

% Calculate Uxx, Uyy, Uxy(=Uyx)
function [U_xx, U_yy, U_xy] = calc_U(mu, x_L, y_L)
    d = sqrt((x_L+mu)^2 + y_L^2);
    r = sqrt((x_L-1+mu)^2 + y_L^2);
    U_xx = 1 - (1-mu)/d^3 - mu/r^3 + 3*(1-mu)*(x_L+mu)^2/d^5 + 3*mu*(x_L-1+mu)^2/r^5;
    U_yy = 1 - (1-mu)/d^3 - mu/r^3 + 3*(1-mu)*y_L^2/d^5 + 3*mu*y_L^2/r^5;
    U_xy = 3*(1-mu)*(x_L+mu)*y_L/d^5 + 3*mu*(x_L-1+mu)*y_L/r^5;
end

%==========================END DEFINE FUNCTIONS===========================

%============================DEFINE SYSTEM================================
% Pluto-Charon System
G = 6.67430e-20; % Constant of Gravitation kg^-1 km^3 s^-2

% Gravitational Parameters [km^3/s^2]
mu_Pluto = 13024.6e18*G;
mu_Charon = 106.1;

% Characteristic Length [km]
a = 19600; % around Pluto
l_char = a;
% Calculate characteristic time
t_char = sqrt(a^3/(mu_Pluto+mu_Charon));


% Earth-Moon System
mu = mu_Pluto/(mu_Pluto + mu_Charon);
% Position of primary bodies
x_Pluto = -mu;
x_Charon = 1-mu;

% Calc L4 Libration point
[x_L4, y_L4] = calc_L4(mu);
%==========================END DEFINE SYSTEM==============================

%============================DEFINE PROBLEM===============================
xi_0 = 0.11;
eta_0 = 0;
xi_dot_0 = 0.005;
eta_dot_0 = 0;
t_end = 8*pi;

% xi_0: 0.010000, eta_0: 0.000000 [non-dim]
% xi_0: 196.000000, eta_0: 0.000000 [km]
% xi_dot_0: -0.005058, eta_dot_0: -0.002529 [non-dim]
% xi_dot_0: -1.128278, eta_dot_0: -0.564139 [m/s]

[lambda_1, lambda_2, lambda_3, lambda_4] = calc_eigenvalues_equilateral(mu);
g = 1-27*mu*(1-mu);

[U_xx, U_yy, U_xy] = calc_U(mu, x_L4, y_L4);

% Construct A Matrix
A = [0 0 1 0; 0 0 0 1; U_xx U_xy 0 2; U_xy U_yy -2 0];
[V,D] = eig(A)
 
a = real(lambda_2);
b = imag(lambda_2);

syms xi(C1,C2,t) eta(C3,C4,t)
xi(t) = C1*exp((a*t)*(cos(b*t)*real(V(1,4))-sin(b*t)*imag(V(1,4)))) + C2*exp((a*t)*(cos(b*t)*imag(V(1,4))+sin(b*t)*real(V(1,4))));
eta(t) = C1*exp((a*t)*(cos(b*t)*real(V(2,4))-sin(b*t)*imag(V(2,4)))) + C2*exp((a*t)*(cos(b*t)*imag(V(2,4))+sin(b*t)*real(V(2,4))));
% xi = C1*exp((a*t)*(cos(b*t)*real(V(1,4))-sin(b*t)*imag(V(1,4)))) + C2*exp((a*t)*(cos(b*t)*imag(V(1,4))+sin(b*t)*real(V(1,4))));
% eta = C1*exp((a*t)*(cos(b*t)*real(V(2,4))-sin(b*t)*imag(V(2,4)))) + C2*exp((a*t)*(cos(b*t)*imag(V(2,4))+sin(b*t)*real(V(2,4))));
xi_dot = diff(xi,t);
eta_dot = diff(eta,t);
xi_ddot = diff(xi,t,2);
eta_ddot = diff(eta,t,2);

xi_0 = 0.001;
eta_0 = 0;
xi_dot_0 = subs(xi_dot,t,0);
eta_dot_0 = subs(eta_dot,t,0);
xi_ddot_0 = subs(xi_ddot,t,0);
eta_ddot_0 = subs(eta_ddot,t,0);

syms C1 C2 C3 C4
eqns = [xi_ddot_0-2*eta_dot_0-U_xx*xi_0-U_xy*eta_0==0, eta_ddot_0+2*xi_dot_0-U_xy*xi_0-U_yy*eta_0==0];
S = vpasolve(eqns,[C1 C2 C3 C4])

% C1_new = 0.020;
% C2_new = 0.015;
% C3_new = 0.010;
% C4_new = 0.005;
C1_new = S.C1;
C2_new = S.C2;
C3_new = S.C3;
C4_new = S.C4;

% Update C values
xi_new = subs(xi, [C1,C2,C3,C4], [C1_new,C2_new,C3_new,C4_new]);
eta_new = subs(eta, [C1,C2,C3,C4], [C1_new,C2_new,C3_new,C4_new]);
xi_dot_new = subs(xi_dot, [C1,C2,C3,C4], [C1_new,C2_new,C3_new,C4_new]);
eta_dot_new = subs(eta_dot, [C1,C2,C3,C4], [C1_new,C2_new,C3_new,C4_new]);

% Get initial velocities
xi_dot_0_new = subs(xi_dot_new, t, 0);
eta_dot_0_new = subs(eta_dot_new, t, 0);

% Calculate the linear orbit for the equilateral points:
orbit_results = zeros(0,5);
for tau = 0:0.01:t_end
    xi_result = subs(xi_new, t, tau);
    eta_result = subs(eta_new, t, tau);
    x = xi_result + x_L4;
    y = eta_result + y_L4; 
    orbit_results(end+1,:) = [tau double(x) double(y) double(xi_result) double(eta_result)];
end
linear_table = array2table(orbit_results, 'VariableNames', {'time', 'x', 'y', 'xi', 'eta'});


%=========================END DEFINE PROBLEM==============================

%============================CALC NONLINEAR RESULT=========================
% Initial Conditions
x_0 = x_L4 + xi_0;
y_0 = y_L4 + eta_0;
w0 = [x_0;y_0;double(xi_dot_0_new);double(eta_dot_0_new)];

% tspan should be calculated with non-dim time
tspan = [0 t_end];

% Setting up the ODE
ode = @(t,w) [...
w(3);...
w(4);...
2*w(4) + w(1) - (1 - mu) * (w(1) + mu) / ((w(1) + mu)^2 + w(2)^2)^(3/2)...
- mu * (w(1) - 1 + mu) / ((w(1) - 1 + mu)^2 + w(2)^2)^(3/2);...
-2*w(3) + w(2) - (1 - mu) * w(2) / ((w(1) + mu)^2 + w(2)^2)^(3/2) - mu * w(2)/((w(1) - 1 + mu)^2 + w(2)^2)^(3/2)...
];
options = odeset('RelTol',1e-12,'AbsTol', 1e-14,'MaxStep', 1e-3);
[t,w] = ode113(ode, tspan, w0, options); 

% Retrieving the results from the integrator
% Position (dimensional)
x_nl = w(:,1);
y_nl = w(:,2);
% Position (non-dimensional)
xi_nl = w(:,1)-x_L4;
eta_nl = w(:,2)-y_L4;
% Velocity
v_x_nl = w(:,3);
v_y_nl = w(:,4);
%===========================END CALC NONLINEAR RESULT=========================

%==========================PRINT IMPORTANT NUMBERS========================

fprintf("========Pluto-Charon System==========\n")
fprintf("Characteristic Length: %f km \n", l_char)
fprintf("Characteristic time: %d sec\n", t_char)
fprintf("Characteristic time: %d days\n", t_char/3600/24)
fprintf("Mu of Pluto-Charon: %f\n", mu)
fprintf("g of Pluto-Charon: %f\n", g)
fprintf("x_L4: %f, y_L4: %f\n", x_L4, y_L4)
fprintf("========Eigenvalues==================\n")
fprintf("lambda_1: %f + %fi\n", real(lambda_1), imag(lambda_1))
fprintf("lambda_2: %f + %fi\n", real(lambda_2), imag(lambda_2))
fprintf("lambda_3: %f + %fi\n", real(lambda_3), imag(lambda_3))
fprintf("lambda_4: %f + %fi\n", real(lambda_4), imag(lambda_4))
fprintf("=========ICs=========================\n")
fprintf("xi_0: %f, eta_0: %f [non-dim]\n", xi_0, eta_0)
fprintf("xi_0: %f, eta_0: %f [km]\n", xi_0*l_char, eta_0*l_char)
fprintf("xi_dot_0: %f, eta_dot_0: %f [non-dim]\n", double(xi_dot_0_new), double(eta_dot_0_new))
fprintf("xi_dot_0: %f, eta_dot_0: %f [m/s]\n", double(xi_dot_0_new)*l_char/t_char*1000, double(eta_dot_0_new)*l_char/t_char*1000)
%======================END PRINT IMPORTANT NUMBERS========================

%======================CONFIGURE PLOTS====================================
fig1 = figure('Name', '1 :non-dim');
L4_plot = scatter(0, 0, 'red', 'filled', 'SizeData', 10);
hold on
linear_orbit = plot(linear_table, 'xi', 'eta', 'Color', '#800080');
nonlinear_orbit = plot(xi_nl, eta_nl, 'Color', '#008000');
hold off
% limit = 5*xi_0;
% ylim([-limit limit])
% xlim([-limit limit])
axis square
xlabel("\xi [non-dim]")
ylabel("\eta [non-dim]")
% xticks(-limit:.01:limit)
legend([L4_plot, linear_orbit, nonlinear_orbit], {'L4', 'Linear Orbit', 'Nonlinear Orbit'})
title({'L4 in Pluto-Charon (\xi-\eta) in'; ['\xi_0=', num2str(xi_0),', t=',num2str(t_end),'rad (Lillian Shido)']})
box on
grid on
fontsize(14, 'points')
%======================END CONFIGURE PLOTS====================================
