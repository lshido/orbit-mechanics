% Givens
AU = 149597870.7; % km
a_Jupiter = 5.2*AU; % km
mu_Sun = 132712440017.99; % km^3/s^2
mu_Jupiter = 126712767.8578; % km^3/s^2
% mu_Jupiter = 0; % ignore since massless compared to Sun

% Determine the mean motion of Jupiter
mean_motion_Jupiter = sqrt((mu_Jupiter+mu_Sun)/a_Jupiter^3);
% mean_motion_Jupiter = sqrt((mu_Sun)/a_Jupiter^3); % 1/s
fprintf("mean motion of Jupiter: %.9e\n", mean_motion_Jupiter);
% Determine the period of Jupiter
Period_Jupiter = 2*pi/mean_motion_Jupiter/3600/24/365;
fprintf("Period of Jupiter: %.9e [years]\n", Period_Jupiter);

% Determine the mean motion of asteroid A.
mean_motion_Asteroid = (3/2)*mean_motion_Jupiter;
fprintf("mean motion of Asteroid A: %.9e\n", mean_motion_Asteroid);
% Determine the period of asteroid A.
Period_Asteroid = 2*pi/mean_motion_Asteroid/3600/24/365;
fprintf("Period of Asteroid A: %.9e [years]\n", Period_Asteroid);

% Determine the mean motion of comet C.
mean_motion_Comet = (3/4)*mean_motion_Jupiter;
fprintf("mean motion of Comet C: %.9e\n", mean_motion_Comet);
% Determine the period of comet C.
Period_Comet = 2*pi/mean_motion_Comet/3600/24/365;
fprintf("Period of Comet C: %.9e [years]\n", Period_Comet);

% Evaluate the ratios of mean motions
n_ratio_c_J = mean_motion_Comet/mean_motion_Jupiter;
n_ratio_a_J = mean_motion_Asteroid/mean_motion_Jupiter;
n_ratio_c_a = mean_motion_Comet/mean_motion_Asteroid;
fprintf("Ratio of Mean Motion Comet and Jupiter: %.9e\n", n_ratio_c_J);
fprintf("Ratio of Mean Motion Asteroid and Jupiter: %.9e\n", n_ratio_a_J);
fprintf("Ratio of Mean Motion Comet and Asteroid: %.9e\n", n_ratio_c_a);

% Evaluate the ratios of periods
P_ratio_c_J = Period_Comet/Period_Jupiter;
P_ratio_a_J = Period_Asteroid/Period_Jupiter;
P_ratio_c_a = Period_Comet/Period_Asteroid;
fprintf("Ratio of Period Comet and Jupiter: %.9e\n", P_ratio_c_J);
fprintf("Ratio of Period Asteroid and Jupiter: %.9e\n", P_ratio_a_J);
fprintf("Ratio of Period Comet and Asteroid: %.9e\n", P_ratio_c_a);