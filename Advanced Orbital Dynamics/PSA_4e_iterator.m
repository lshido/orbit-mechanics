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
e = 0.99675;

step = 0.5;
counter = 0;
tolerance = 10e-12;
results_i = zeros(0,5);
results_ii = zeros(0,5);
results_iii = zeros(0,5);

% M = 4 %[deg]
for M = 0:step:360
    E_initial = deg2rad(M);
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
            results_i(end+1,:) = [M rad2deg(M_guess) rad2deg(E_guess) rad2deg(E_initial) counter];
            counter = 0;
            break
        end
    end
end

% for M = 0:step:360
%     E_initial = deg2rad(M) + e*sin(deg2rad(M));
%     E_guess = deg2rad(M) + e*sin(deg2rad(M));
%     while 1
%         counter = counter + 1;
%         M_guess = E_guess - e*sin(E_guess);
%         delta_M = deg2rad(M) - M_guess;
%         delta_E = delta_M / (1 - e*cos(E_guess));
%         E_guess = E_guess + delta_E;
%         if abs(delta_E) > tolerance
%             continue
%         else
%             results_ii(end+1,:) = [M rad2deg(M_guess) rad2deg(E_guess) rad2deg(E_initial) counter];
%             counter = 0;
%             break
%         end
%     end
% end

% for M = 0:step:360
%     E_initial = deg2rad(M) + e*sin(deg2rad(M)) + ((e^2)/2)*sin(2*deg2rad(M)) + (e^3)*((3/8)*sin(3*deg2rad(M))-(1/8)*sin(deg2rad(M)));
%     E_guess = deg2rad(M) + e*sin(deg2rad(M)) + ((e^2)/2)*sin(2*deg2rad(M)) + (e^3)*((3/8)*sin(3*deg2rad(M))-(1/8)*sin(deg2rad(M)));
%     while 1
%         counter = counter + 1;
%         M_guess = E_guess - e*sin(E_guess);
%         delta_M = deg2rad(M) - M_guess;
%         delta_E = delta_M / (1 - e*cos(E_guess));
%         E_guess = E_guess + delta_E;
%         if abs(delta_E) > tolerance
%             continue
%         else
%             results_iii(end+1,:) = [M rad2deg(M_guess) rad2deg(E_guess) rad2deg(E_initial) counter];
%             counter = 0;
%             break
%         end
%     end
% end

t1 = array2table(results_i,'VariableNames',{'Target M [deg]', 'M_result [deg]', 'E [deg]', 'E_initial [deg]', 'Iterations'})
t2 = array2table(results_ii,'VariableNames',{'Target M [deg]', 'M_result [deg]', 'E [deg]', 'E_initial [deg]','Iterations'})
t3 = array2table(results_iii,'VariableNames',{'Target M [deg]', 'M_result [deg]', 'E [deg]', 'E_initial [deg]','Iterations'})

fig1 = figure('Name', 'PS_A4_b');
E_plot = plot(t1, "Target M [deg]", "E [deg]");
hold on;
i_plot = plot(t1, "Target M [deg]", "E_initial [deg]");
ii_plot = plot(t2, "Target M [deg]", "E_initial [deg]");
iii_plot = plot(t3, "Target M [deg]", "E_initial [deg]");
i_plot.LineStyle = "-.";
ii_plot.LineStyle = ":";
iii_plot.LineStyle = "--";
i_plot.LineWidth = 2;
ii_plot.LineWidth = 2;
iii_plot.LineWidth = 2;
fontsize(14, "points")
legend('Converged E', 'Initial Guess for E (i): E_0 = M', 'Initial Guess for E (ii): 1st two terms', 'Initial Guess for E (iii): terms through order e^3')
hold off;
title("Converged E and Initial Guess for E as a function of M (Lillian Shido, PSA\_4d\_iterator.m)")
ylabel("Eccentric Anomaly E [degrees]")
xlim([0 360])
ylim([0 360])
xticks(0:30:360)
xlabel("Mean Anomaly M [degrees]")