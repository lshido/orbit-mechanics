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
e = 0.93;
a = 20 * R_earth;

% Find mean motion n
mean_motion = sqrt(Gm_earth/(a^3));
fprintf('Mean motion (n): %.7e[rad/s]\n', mean_motion);

counter = 0;
tolerance = 10e-12;
all_results = zeros(0,8);

for M = 0:1:360
    E_guess = deg2rad(M);
    while 1
        counter = counter + 1;
        time_since_periapsis = deg2rad(M) / mean_motion;
        M_guess = E_guess - e*sin(E_guess);
        delta_M = deg2rad(M) - M_guess;
        delta_E = delta_M / (1 - e*cos(E_guess));
        E_guess = E_guess + delta_E;
        if abs(delta_E) > tolerance
            continue
        else
            % calc True Anomaly
            TA = 2*(atan(sqrt((1+e)/(1-e))*tan(E_guess/2)));
            r = a*(1 - e*cos(E_guess));
            all_results(end+1,:) = [time_since_periapsis time_since_periapsis/3600 M rad2deg(M_guess) rad2deg(E_guess) counter mod(rad2deg(TA), 360) r/R_earth];
            counter = 0;
            break
        end
    end
end
t = array2table(all_results,'VariableNames',{'Time [sec]' 'Time [hours]' 'Target M [deg]', 'M_result [deg]', 'E [deg]', 'Iterations', 'True Anomaly [deg]', 'radius in R_earth [R_earth]'});
disp(t)

fig1 = figure('Name', 'PS_A3_anomalies');
r_ta_plot = plot(t, "True Anomaly [deg]", "radius in R_earth [R_earth]");
hold on;
r_e_plot = plot(t, "E [deg]", "radius in R_earth [R_earth]");
r_m_plot = plot(t, "Target M [deg]", "radius in R_earth [R_earth]");
r_e_plot.LineStyle = ":";
r_m_plot.LineStyle = "--";
legend
hold off;
title("Radius as a function of True Anomaly, E, and M")
ylabel("radius in R_earth [R_earth]")
xlim([0 360])
xticks(0:30:360)
xlabel("Angles [deg]")