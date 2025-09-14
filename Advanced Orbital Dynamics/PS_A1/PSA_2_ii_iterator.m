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
e = 0.85;
a = 7.5 * R_earth;

% Find mean motion n
mean_motion = sqrt(Gm_earth/(a^3));
fprintf('Mean motion (n): %.7e[rad/s]\n', mean_motion);

counter = 0;
tolerance = 10e-12;
all_results = zeros(0,5);

for M = 0:1:360
    E_guess = deg2rad(M);
    while 1
        counter = counter + 1;
        time_since_periapsis = deg2rad(M) / mean_motion;
        M_guess = E_guess - e*sin(E_guess);
        delta_M = deg2rad(M) - M_guess;
        delta_E = delta_M / (1 - e*cos(E_guess));
        if abs(delta_E) > tolerance
            E_guess = E_guess + delta_E;
            continue
        else
            % calc True Anomaly
            TA = 2*(atan(sqrt((1+e)/(1-e))*tan(E_guess/2)));
            all_results(end+1,:) = [time_since_periapsis time_since_periapsis/3600 M rad2deg(E_guess) mod(rad2deg(TA), 360)];
            counter = 0;
            break
        end
    end
end

t = array2table(all_results,'VariableNames',{'Time [sec]' 'Time [hours]' 'Mean Anomaly [deg]', 'Eccentric Anomaly [deg]', 'True Anomaly [deg]'});
head(t, 10)
tail(t, 10)

fig2 = figure('Name', 'PS_A2_anomalies');
ta_plot = plot(t, "Time [hours]", "True Anomaly [deg]");
hold on;
e_plot = plot(t, "Time [hours]", "Eccentric Anomaly [deg]");
m_plot = plot(t, "Time [hours]", "Mean Anomaly [deg]");
e_plot.LineStyle = ":";
m_plot.LineStyle = "--";
ta_plot.LineWidth = 2;
e_plot.LineWidth = 2;
m_plot.LineWidth = 2;
fontsize(14, "points")
legend ("True Anomaly \theta^*", "Eccentric Anomaly E", "Mean Anomaly M")
hold off;
title("Time history for True Anomaly, Eccentric Anomaly, and Mean Anomaly (Lillian Shido, PSA\_2\_ii\_iterator.m)")
ylabel("Angle [degrees]")
yticks(0:45:360)
ylim([0 360])