% Title: Newton-Raphson Algorithm to calculate values of gamma_2 and L2
% Author: Lillian Shido
% Date: September 27, 2025

clear all

% Characteristic length [km]
a_Earth = 149597898;
a_Moon = 384400; % around Earth
a_Titan = 1221865; % about Saturn
a_Phobos = 9376; % about Mars

% Gravitational Parameters [km^3/s^2]
mu_Sun = 132712440017.99;
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;
mu_Mars = 42828.314258067;
mu_Saturn = 37940626.061137;
mu_Phobos = 0.0007112;
mu_Titan = 8978.1382;

% Mass ratio of minor primary
mu_Earth_Moon = mu_Moon/(mu_Earth+mu_Moon);
mu_Sun_Earth = mu_Earth/(mu_Sun+mu_Earth);
mu_Saturn_Titan = mu_Titan/(mu_Saturn+mu_Titan);
mu_Mars_Phobos = mu_Phobos/(mu_Mars+mu_Phobos);

% Create a table
System = ["Earth-Moon";"Sun-Earth";"Saturn-Titan";"Mars-Phobos"];
mu_system = [mu_Earth_Moon;mu_Sun_Earth;mu_Saturn_Titan;mu_Mars_Phobos];
a_system = [a_Moon;a_Earth;a_Titan;a_Phobos];
system_data = table(System,mu_system,a_system)

% Newton-Raphson Algorithm

% Step 1: Pick a gamma
% Step 2: Plug into function
% Step 3: Check if result = abs(10^-12)
% Step 4: If not, update gamma using update equation and repeat
% Step 5: If it is, break and print gamma result.

L1_results = zeros(0,9);
tolerance = 1e-12;
counter = 0;
rows = height(system_data);
for row = 1:rows
    name = System(row,1);
    mu = system_data{row,2};
    a = system_data{row,3};
    % Step 1
    gamma = 0.2; % initial guess
    while 1
        counter = counter + 1;
        % Step 2
        f = 1-mu-gamma-((1-mu)/(1-gamma)^2)+(mu/gamma^2);
        f_prime = -1-(2*(1-mu)/(1-gamma)^3)-(2*mu/gamma^3);
        % result = 1-mu-gamma-((1-mu)/(1-gamma)^2)+(mu/gamma^2);
        % Step 3
        if abs(f) > tolerance
            % Step 4
            gamma = gamma - f/f_prime;
            % gamma = gamma - ((1-mu-gamma-((1-mu)/(1-gamma)^2)+(mu/gamma^2))/((-2*(1-mu)/(1-gamma)^3)-(2*mu/gamma^3)-1));
            continue
        else
            % Step 5
            x = 1 - mu - gamma;
            % Add check for accleration
            d_x = x*a+mu;
            r_x = x*a-1+mu;
            d = ((x*a+mu)^2)^(1/2);
            r = ((x*a-1+mu)^2)^(1/2);
            accel = -(1-mu)/d^3*d_x - mu/r^3*r_x;
            % Add check for partial wrt gamma
            partial = (1-mu)/(1-gamma)^2 - mu/(gamma)^2 - (1-mu-gamma);
            L1_results(end+1,:) = [string(name) mu gamma gamma*a gamma*100 x x*a accel partial];
            counter = 0;
            break
        end
    end
end
L1_table = array2table(L1_results,'VariableNames',{'System', 'mu of system', 'L1 gamma', 'L1 gamma_dim', 'percentage of a', 'x', 'x_dim', 'Acceleration', 'Partial'});
format shortE
disp(L1_table)


L3_results = zeros(0,9);
tolerance = 1e-12;
counter = 0;
rows = height(system_data);
for row = 1:rows
    name = System(row,:);
    mu = system_data{row,2};
    a = system_data{row,3};
    % Step 1
    gamma = 0.000001; % initial guess
    while 1
        counter = counter + 1;
        % Step 2
        f = ((1-mu)/gamma^2)+(mu/(gamma+1)^2)-mu-gamma;
        f_prime = -(2*(1-mu)/gamma^3)-(2*mu/(gamma+1)^3)-1;
        % Step 3
        if abs(f) > tolerance
            % Step 4
            gamma = gamma - (f/f_prime);
            continue
        else
            % Step 5
            x = - mu - gamma;
            % Add check for accleration
            d_x = x*a+mu;
            r_x = x*a-1+mu;
            d = ((x*a+mu)^2)^(1/2);
            r = ((x*a-1+mu)^2)^(1/2);
            accel = -(1-mu)/d^3*d_x - mu/r^3*r_x;
            % Add check for partial wrt gamma
            partial = -((1-mu)/gamma^2) - (mu/(gamma+1)^2) + (mu+gamma);
            L3_results(end+1,:) = [string(name) mu gamma gamma*a gamma*100 x x*a accel partial];
            counter = 0;
            break
        end
    end
end
L3_table = array2table(L3_results,'VariableNames',{'System', 'mu of system', 'L3 gamma', 'L3 gamma_dim', 'percentage of a', 'x', 'x_dim', 'Acceleration', 'Partial'});
format shortE
disp(L3_table)