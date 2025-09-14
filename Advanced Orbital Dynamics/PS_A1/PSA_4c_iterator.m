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
e = 0.999;

% Find mean motion n
mean_motion = sqrt(Gm_earth/(a^3));
fprintf('Mean motion (n): %.7e[rad/s]\n', mean_motion);

counter = 0;
tolerance = 10e-12;
results_i = zeros(0,4);
results_ii = zeros(0,4);
results_iii = zeros(0,4);

% M = 4 %[deg]
for M = 0:0.75:360
    E_guess = deg2rad(M);
    while 1
        counter = counter + 1;
        M_guess = E_guess - e*sin(E_guess);
        delta_M = deg2rad(M) - M_guess;
        delta_E = delta_M / (1 - e*cos(E_guess));
        E_guess = E_guess + delta_E;
        if abs(delta_E) > tolerance
            continue
        else
            results_i(end+1,:) = [M rad2deg(M_guess) rad2deg(E_guess) counter];
            counter = 0;
            break
        end
    end
end

t1 = array2table(results_i,'VariableNames',{'Target M [deg]', 'M_result [deg]', 'E [deg]', 'Iterations'})

fig1 = figure('Name', 'PS_A4_c');
i_plot = plot(t1, "Target M [deg]", "E [deg]");
i_plot.LineWidth = 2;
fontsize(14, "points")
hold on;
hold off;
title("E as a function of M (Lillian Shido, PSA\_4c\_iterator.m)")
ylabel("Eccentric Anomaly [degrees]")
xlim([0 360])
ylim([0 360])
yticks(0:30:360)
xticks(0:30:360)
xlabel("Mean Anomaly M [degrees]")