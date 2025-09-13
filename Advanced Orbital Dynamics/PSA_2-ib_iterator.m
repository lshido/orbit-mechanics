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
e = 0.85;
a = 7.5 * R_earth;

counter = 0;
tolerance = 10e-12 % [rad];
all_results = zeros(0,4);

for M = 0:10:360 % [deg]
    % Step 1
    E_guess = deg2rad(M);
    while 1
        counter = counter + 1;
        % Step 2
        M_guess = E_guess - e*sin(E_guess);
        % Step 3
        delta_M = deg2rad(M) - M_guess;
        % Step 4
        delta_E = delta_M / (1 - e*cos(E_guess));
        % Step 5
        if abs(delta_E) > tolerance
            % Step 6
            E_guess = E_guess + delta_E;
            continue
        else
            all_results(end+1,:) = [ M rad2deg(M_guess) rad2deg(E_guess) counter];
            counter = 0;
            break
        end
    end
end

t = array2table(all_results,'VariableNames',{'Target M [deg]', 'M_result [deg]', 'E [deg]', 'Iterations'});
disp(t)