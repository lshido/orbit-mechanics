% Use an integrator to solve the rotating frame EOMs

% Initial conditions
r0_vector = [-0.270 -0.420];
v0_vector = [0.300 1.000];
w0 = [r0_vector(1);r0_vector(2);v0_vector(1);v0_vector(2)];

% Gravitational Parameters [km^3/s^2]
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

mu = 1.2151e-02;

% Characteristic Length [km]
a_Moon = 384400; % around Earth
l_char = a_Moon;

% Calculate characteristic time
t_char = sqrt(a_Moon^3/(mu_Earth+mu_Moon));
fprintf("characteristic time: %d sec\n", t_char)

% tspan should be calculated with non-dim time
tspan = [0 8*pi];

ode = @(t,w) [...
w(3);...
w(4);...
2*w(4) + w(1) - ((1-mu)*(w(1)+mu)/((w(1)+mu)^2 + w(2)^2)^(3/2)) - mu*(w(1)-1+mu)/(((w(1)-1+mu)^2 + w(2)^2)^(3/2));...
-2*w(3) + w(2) - ((1-mu)*w(2)/((w(1)+mu)^2 + w(2)^2)^(3/2)) - (mu*w(2)/(((w(1)-1+mu)^2 + w(2)^2)^(3/2)))...
];
options = odeset('RelTol',1e-12,'AbsTol', 1e-12, 'MaxStep',0.05);
% [t,w] = ode45(ode, tspan, w0, options); 
[t,w] = ode113(ode, tspan, w0, options); 

x = w(:,1);
y = w(:,2);

x_Earth = -mu;
x_Moon = 1-mu;

fig1 = figure('Name','static');
p = plot(x,y);
hold on
scatter(x_Earth, 0, 'black', 'filled', 'SizeData', 200)
scatter(x_Moon, 0, 'red', 'filled', 'SizeData', 100)
% sc_marker = scatter(NaN,NaN,[], "blue", "filled");
hold off
xlim([-1 1])
ylim([-1 1])
axis square
fontsize(14, "points")

%=============ANIMATE THE POSITIONS!==========================
% for k = 1:numel(t)
%     sc_marker.XData = w(k,1);
%     sc_marker.YData = w(k,2);
%     drawnow
%     % exportgraphics(gca,"parabola.gif",Append=true)
% end