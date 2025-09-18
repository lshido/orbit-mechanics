% We don't know M or E, we just need to guess.
% So we start with E = M
% M = E-esin(E)
% we need to calculate E for each M.
% Step 1: For a given M, guess that E_guess is M
% STep 2: Plug in E_guess into E-esinE
% Step 3: Get M_guess
% step 4: calculate delta_M = M - M_guess
% step 5: calculate delta_E = delta_M / (1-ecos(E_guess))
% step 6: if delta_E > 10^-12,
% step 7: E_guess = E_guess + delta_E
% step 8: try new E_guess, repeat.
% step 9: if delta_E < 10^-12, break loop and print E.
% NOTE: perform calculation in rads, then convert final answer to degrees!

R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;
e1 = 0.999;

% Find mean motion n
mean_motion = sqrt(Gm_earth/(a^3));
fprintf('Mean motion (n): %.7e[rad/s]\n', mean_motion);

counter = 0;
tolerance = 10e-12;
results = zeros(0,5);

M = 342; %[deg]
E_guess = deg2rad(M);
while 1
    counter = counter + 1;
    M_guess = E_guess - e1*sin(E_guess);
    delta_M = deg2rad(M) - M_guess;
    delta_E = delta_M / (1 - e1*cos(E_guess));
    E_guess = E_guess + delta_E;
    results(end+1,:) = [counter rad2deg(M_guess) rad2deg(delta_M) rad2deg(delta_E) rad2deg(E_guess)];
    if abs(delta_E) > tolerance
        continue
    else
        counter = 0;
        break
    end
end

t = array2table(results,'VariableNames',{'Iteration' 'M_guess' 'delta_M', 'delta_E', 'E_guess'})

fig1 = figure('Name', 'A4b-2');
delta_E_plot = scatter(t, 'Iteration', 'delta_E', 'filled');
hold on;
E_guess_plot = scatter(t, 'Iteration', 'E_guess');
delta_E_plot.LineWidth = 2;
E_guess_plot.LineWidth = 2;
fontsize(14, "points")
legend('delta E', 'E guess')
hold off;
title("Values for delta E and Guessed E for M=342 degrees (Lillian Shido, PSA\_4b\_M342.m)")
ylabel("Degrees")
xlabel("Iteration Number")