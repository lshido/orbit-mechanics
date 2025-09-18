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
e1 = 0.14;

% Find mean motion n
mean_motion = sqrt(Gm_earth/(a^3));
fprintf('Mean motion (n): %.7e[rad/s]\n', mean_motion);

counter = 0;
tolerance = 10e-12;
results_e1 = zeros(0,4);
results_e2 = zeros(0,4);

% M = 4 %[deg]
for M = 0:0.5:360
    E_guess = deg2rad(M);
    fprintf("E_initial: %d\n", rad2deg(E_guess));
    while 1
        counter = counter + 1;
        M_guess = E_guess - e1*sin(E_guess);
        delta_M = deg2rad(M) - M_guess;
        delta_E = delta_M / (1 - e1*cos(E_guess));
        E_guess = E_guess + delta_E;
        if abs(delta_E) > tolerance
            continue
        else
            results_e1(end+1,:) = [M rad2deg(M_guess) rad2deg(E_guess) counter];
            counter = 0;
            break
        end
    end
end

e2 = 0.999;
for M = 0:0.5:360
    E_guess = deg2rad(M);
    fprintf("E_initial: %d\n", rad2deg(E_guess));
    while 1
        counter = counter + 1;
        M_guess = E_guess - e2*sin(E_guess);
        delta_M = deg2rad(M) - M_guess;
        delta_E = delta_M / (1 - e2*cos(E_guess));
        E_guess = E_guess + delta_E;
        if abs(delta_E) > tolerance
            continue
        else
            results_e2(end+1,:) = [M rad2deg(M_guess) rad2deg(E_guess) counter];
            counter = 0;
            break
        end
    end
end

t1 = array2table(results_e1,'VariableNames',{'Target M [deg]', 'M_result [deg]', 'E [deg]', 'Iterations'})
t2 = array2table(results_e2,'VariableNames',{'Target M [deg]', 'M_result [deg]', 'E [deg]', 'Iterations'})


fig1 = figure('Name', 'PS_A4_b');
e1_plot = scatter(t1, "Target M [deg]", "Iterations");
hold on;
e2_plot = scatter(t2, "Target M [deg]", "Iterations", 'filled');
% r_m_plot = plot(t, "Target M [deg]", "radius in R_earth [R_earth]");
% r_e_plot.LineStyle = ":";
% e2_plot.LineStyle = "--";
% legend('e=0.14', 'e=0.999')
e1_plot.LineWidth = 2;
e2_plot.LineWidth = 2;
fontsize(14, "points")
hold off;
title("Iterations as a function of M, stepsize=0.5 (Lillian Shido, PSA\_4b\_iterator\_e.m)")
ylabel("Number of Iterations")
legend('e = 0.14', 'e = 0.999')
xlim([0 360])
xticks(0:30:360)
xlabel("Mean Anomaly M [degrees]")