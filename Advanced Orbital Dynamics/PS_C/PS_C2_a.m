% Gravitational Parameters [km^3/s^2]
mu_Sun = 132712440017.99;
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

% Earth-Moon System
mu = mu_Moon/(mu_Earth + mu_Moon);
% Locations of Earth-Moon Collinear Libration Points
x_L1 = 8.3692e-01; % From Problem Set B
y_L1 = 0;
x_L2 = 1.1557e+00; % From Problem Set B
y_L2 = 0;
x_L3 = -1.0051; % From Problem Set B
y_L3 = 0;

% Sun-Earth System
% mu = mu_Earth/(mu_Sun + mu_Earth);
% % Location of Sun-Earth Collinear Libration Points
% x_L1 = 0.99003; % From Problem Set B
% y_L1 = 0;
% x_L2 = 1.01; % From Problem Set B
% y_L2 = 0;
% x_L3 = -1; % From Problem Set B
% y_L3 = 0;

% Calculate the roots
y = 0; % collinear points
for x = [x_L1, x_L2, x_L3]
    d = sqrt((x+mu)^2 + y^2);
    r = sqrt((x-1+mu)^2 + y^2);
    U_xx = 1 - (1-mu)/d^3 - mu/r^3 + 3*(1-mu)*(x+mu)^2/d^5 + 3*mu*(x-1+mu)^2/r^5;
    U_yy = 1 - (1-mu)/d^3 - mu/r^3;
    B_1 = 2 - (U_xx + U_yy)/2;
    B_2_squared = -U_xx*U_yy;
    big_Lambda_1 = -B_1 + (B_1^2 + B_2_squared)^(1/2);
    big_Lambda_2 = -B_1 - (B_1^2 + B_2_squared)^(1/2);
    lambda_1 = sqrt(big_Lambda_1); % real
    lambda_2 = -sqrt(big_Lambda_1); % real
    lambda_3 = sqrt(big_Lambda_2); % imaginary
    lambda_4 = -sqrt(big_Lambda_2); % imaginary
    fprintf("x_L: %f\n", x)
    fprintf("lambda_1: %f + %fi\n", real(lambda_1), imag(lambda_1))
    fprintf("lambda_2: %f + %fi\n", real(lambda_2), imag(lambda_2))
    fprintf("lambda_3: %f + %fi\n", real(lambda_3), imag(lambda_3))
    fprintf("lambda_4: %f + %fi\n", real(lambda_4), imag(lambda_4))
end