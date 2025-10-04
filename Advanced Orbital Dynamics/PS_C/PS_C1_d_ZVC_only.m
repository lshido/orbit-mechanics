clear all

% Calculate C for the given position and velocity.
% Gravitational Parameters [km^3/s^2]
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

mu = mu_Moon/(mu_Earth + mu_Moon);
fprintf("mu %d\n", mu)

% mu_Earth_Moon = 1.2151e-02;

% position and velocity in NON-DIMENSIONAL units
% r_vector = [-0.270 -0.420];
% v_vector = [0.300 -1.000];

% x_0 = r_vector(1);
% y_0 = r_vector(2);
% d_0 = sqrt((x_0+mu)^2 + y_0^2);
% r_0 = sqrt((x_0-1+mu)^2 + y_0^2);
% pseudo_U = (1-mu)/d_0 + mu/r_0 + (x_0^2+y_0^2)/2;
% v_squared = norm(v_vector)^2;
% C = 2*pseudo_U - v_squared;
C = 3.150;
fprintf("Jacobi Constant C: %d\n", C)

% Position of primary bodies
x_Earth = -mu;
x_Moon = 1-mu;

% Locations of L1 and L2
x_L1 = 8.3692e-01;
y_L1 = 0;

x_L2 = 1.1557e+00;
y_L2 = 0;

x_L3 = -1.0051;
y_L3 = 0;

x_L4 = cosd(60) - mu;
y_L4 = sind(60);

x_L5 = cosd(60) - mu;
y_L5 = -sind(60);

% To plot the ZVC for z = 0:
% Step 1: For a given x
% Step 2: Guess a y
% Step 3: Calculate Jacobi constant result
% Step 4: If Jacobi constant result is less than tolerance, break
% Step 5: If not, update y with Newton Raphson
% Step 6: Repeat

zvc_result = zeros(0,6);
tolerance = 1e-12;

%==================Try to calculate one point=========================
% First things first, let's try and get a few points
x_vector = [0, 0, 0.1, 0, x_Moon];
y_vector = [0, 1, 0, 0.1, 0.1];
pairs = [x_vector; y_vector]'; % Transpose to get rows of x,y pairs
for k = 1:height(pairs)
    x = pairs(k,1);
    y = pairs(k,2);
    counter = 0;
    while 1
        counter = counter + 1;
        d = sqrt((x+mu)^2 + y^2);
        r = sqrt((x-1+mu)^2 + y^2);
        f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
        f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d^3 - 2*mu*(x-1+mu)/r^3;
        delta = f/f_prime_x;
        % Step 4:
        if abs(f) > tolerance
            x = x - delta;
            continue
        % elseif counter > 1e3
        %     break
        else
            new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
            error_C = abs(new_C - C)/C*100;
            zvc_result(end+1,:) = [x y -y new_C error_C counter];
            counter = 0;
            break
        end
    end
end
for k = 1:height(pairs)
    x = pairs(k,1);
    y = pairs(k,2);
    counter = 0;
    while 1
        counter = counter + 1;
        d = sqrt((x+mu)^2 + y^2);
        r = sqrt((x-1+mu)^2 + y^2);
        f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
        f_prime_y = 2*y*( 1 - (1-mu)/d^3 - mu/r^3);
        delta = f/f_prime_y;
        if abs(f) > tolerance
            y = y - delta;
            continue
        else
            new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
            error_C = abs(new_C - C)/C*100;
            zvc_result(end+1,:) = [x y -y new_C error_C counter];
            counter = 0;
            break
        end
    end
end

%===================Calculate the "inside" curve=====================
% Give an initial guess for y that's "inside"
% for y = [-0.1, 0.1];
%     for x = [linspace(0,-1.21,5e3), linspace(0,1.21,5e3)] %  Find the curve for -1.21 < x < 0
%         counter = 0;
%         while 1
%             counter = counter + 1;
%             d = sqrt((x+mu)^2 + y^2);
%             r = sqrt((x-1+mu)^2 + y^2);
%             f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
%             f_prime_y = 2*y*( 1 - (1-mu)/d^3 - mu/r^3);
%             delta = f/f_prime_y;
%             if abs(f) > tolerance
%                 y = y - delta;
%                 continue
%             else
%                 new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
%                 error_C = abs(new_C - C)/C*100;
%                 zvc_result(end+1,:) = [x y new_C error_C counter];
%                 counter = 0;
%                 break
%             end
%         end
%     end
% end

% Give an initial guess for x that starts "inside"
% for x = [-0.1, 0.1]; % Try starting from the "outside" of the outer boundary
%     for y = [linspace(0,-1.21,5e3), linspace(0,1.21,5e3)] % Find the curve for 0 < y < 1.21
%         while 1
%             counter = counter + 1;
%             d = sqrt((x+mu)^2 + y^2);
%             r = sqrt((x-1+mu)^2 + y^2);
%             f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
%             f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d^3 - 2*mu*(x-1+mu)/r^3;
%             delta = f/f_prime_x;
%             % Step 4:
%             if abs(f) > tolerance
%                 x = x - delta;
%                 continue
%             % elseif counter > 1e3
%             %     break
%             else
%                 new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
%                 error_C = abs(new_C - C)/C*100;
%                 zvc_result(end+1,:) = [x y new_C error_C counter];
%                 counter = 0;
%                 break
%             end
%         end
%     end
% end

%===================Calculate the "outside" curve=====================
% Give an initial guess for y that starts "outside"
% for y = [-1.5, 1.5]; % Try starting from the "outside" of the outer boundary
%     for x = [linspace(-1.21,0,1e3), linspace(0,1.21,1e3)] %  Find the curve for -1.21 < x < 0
%         counter = 0;
%         while 1
%             counter = counter + 1;
%             d = sqrt((x+mu)^2 + y^2);
%             r = sqrt((x-1+mu)^2 + y^2);
%             f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
%             f_prime_y = 2*y*( 1 - (1-mu)/d^3 - mu/r^3);
%             delta = f/f_prime_y;
%             if abs(f) > tolerance
%                 y = y - delta;
%                 continue
%             else
%                 new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
%                 error_C = abs(new_C - C)/C*100;
%                 zvc_result(end+1,:) = [x y new_C error_C counter];
%                 counter = 0;
%                 break
%             end
%         end
%     end
% end

% Give an initial guess for x that starts "outside"
% for x = [-1.5, 1.5]; % Try starting from the "outside" of the outer boundary
%     for y = [linspace(-1.21,0,1e3), linspace(0,1.21,1e3)] % Find the curve for 0 < y < 1.21
%         while 1
%             counter = counter + 1;
%             d = sqrt((x+mu)^2 + y^2);
%             r = sqrt((x-1+mu)^2 + y^2);
%             f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
%             f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d^3 - 2*mu*(x-1+mu)/r^3;
%             delta = f/f_prime_x;
%             if abs(f) > tolerance
%                 x = x - delta;
%                 continue
%             % elseif counter > 1e3
%             %     break
%             else
%                 new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
%                 error_C = abs(new_C - C)/C*100;
%                 zvc_result(end+1,:) = [x y new_C error_C counter];
%                 counter = 0;
%                 break
%             end
%         end
%     end
% end

zvc_table = array2table(zvc_result, 'VariableNames', {'x','y', '-y', 'JC', '% Error of C','Iterations'});
format short
disp(zvc_table)

% Find the max error of JC
max_error = max(abs(zvc_result(:,4)));
fprintf("max error: %d\n", max_error)

% Locations of points when x=0 and y=0
x_1 = 0;
y_1 = 1.2722;

x_2 = 1.1159;
y_2 = 0;

%===================Plot the ZVC!===============
fig1 = figure('Name','ZVC');
earth = scatter(x_Earth, 0, 'blue', 'filled', 'SizeData', 200);
hold on
moon = scatter(x_Moon, 0, 'black', 'filled', 'SizeData', 100);
L1_plot = scatter(x_L1, y_L1, 'red', 'filled', 'SizeData', 50);
scatter(x_L2, y_L2, 'red', 'filled', 'SizeData', 50)
scatter(x_L3, y_L3, 'red', 'filled', 'SizeData', 50)
scatter(x_L4, y_L4, 'red', 'filled', 'SizeData', 50)
scatter(x_L5, y_L5, 'red', 'filled', 'SizeData', 50)
% p1_plot = scatter(x_1, y_1, 'magenta', 'filled', 'SizeData', 50);
% p2_plot = scatter(x_2, y_2, 'magenta', 'filled', 'SizeData', 50);
zvc_plot = scatter(zvc_table, 'x', 'y', 'filled', 'SizeData', 10);
scatter(zvc_table, 'x', '-y')
hold off
xlim([-1.5 1.5])
ylim([-1.5 1.5])
axis square
xlabel("Non-dimensional X")
ylabel("Non-dimensional Y")
legend([earth, moon, L1_plot, zvc_plot], {'Earth', 'Moon', 'Eq. Points', 'ZVC'})
title({'The Zero Velocity Curves at z=0';['in the Earth-Moon System for Jacobi Constant = ', num2str(C)]})
box on
grid on
fontsize(14, 'points')