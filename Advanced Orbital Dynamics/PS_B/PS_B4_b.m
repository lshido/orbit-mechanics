% spacecraft in the Earth-Moon system
% CR3BP

% Gravitational Parameters [km^3/s^2]
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

% Characteristic Length [km]
a_Moon = 384400; % around Earth
l_char = a_Moon;

% Calculate characteristic time
t_char = sqrt(a_Moon^3/(mu_Earth+mu_Moon));
fprintf("characteristic time: %d sec\n", t_char)

% position and velocity in NON-DIMENSIONAL units
r_vector = [-0.270 -0.420];
v_vector = [0.300 1.000];

r_dim_vector = r_vector * l_char *1000;
fprintf("dimensional position: %7e m\n", r_dim_vector)

v_dim_vector = v_vector * l_char/t_char*1000;
fprintf("dimensional vector: %7e m/s\n", v_dim_vector)

% Calculate Jacobi's constant for ZVC at each equilibrium point
% 0 = 2*pseudo_U - C
x_1 = 8.3692e-01;
y_1 = 0;

x_2 = 1.1557e+00
y_2 = 0;

x_3 = -1.0051;
y_3 = 0;

% mu of Earth-Moon system
mu_Earth_Moon = 1.2151e-02;
x_4 = cosd(60) - mu_Earth_Moon;
y_4 = sind(60);

x_5 = cosd(60) - mu_Earth_Moon;
y_5 = -sind(60);

% Create a table
System = ['L1';'L2';'L3';'L4';'L5'];
x_eq = [x_1;x_2;x_3;x_4;x_5];
y_eq = [y_1;y_2;y_3;y_4;y_5];
system_data = table(System,x_eq,y_eq)

C_results = zeros(0,7);
rows = height(system_data);
for row = 1:rows
    name = system_data{row,1};
    x = system_data{row,2};
    y = system_data{row,3};
    d = sqrt((x+mu_Earth_Moon)^2 + y^2);
    r = sqrt((x-1+mu_Earth_Moon)^2 + y^2);
    pseudo_U = (1-mu_Earth_Moon)/d + mu_Earth_Moon/r + (x^2+y^2)/2;
    C = 2*pseudo_U;
    C_results(end+1,:) = [convertCharsToStrings(name) x y d r pseudo_U C];
end

C_table = array2table(C_results, 'VariableNames', {'Eq. Point', 'x', 'y', 'd', 'r', 'pseudo U', 'C'});
format short
disp(C_table)