% EME 451 COMPUTATIONAL FLUID DYNAMICS
% ASSIGNMENT 1 
% HOMEWORK 1 PROBLEM 2
%%
% VINOD RAO A/L JAYAPRASADH (SCHOOL OF CHEMICAL ENGINEERING)
% MATRIC NUMBER: 158635
% NURUL AIN FAZWIN BINTI MOHAMAD SAKMAH (SCHOOL OF MECHANICAL ENGINEERING)
% MATRIC NUMBER: 153477
%%
% Code 1: Parameter Study Implementation

clc;
clear all;
close all;

% Given conditions
t0 = 0;
u0 = 1;
tEnd = 2;
delta_t = 0.01;
a = 2;
p = 1;

% Parameters to manipulate
changing_a = [2, 1/2, -1/2];
changing_delta_t = [1, 0.1, 0.01];
changing_p = [0.9, 1, 1.1];

% Time vector for default delta_t
T = t0:delta_t:tEnd;

% Analytical solution
U_True = u0 * exp(a * T);

% Plot styles and labels
plotstyles = {'-b', '-g', '-r', '--k'};
legends_p = {'p=0.9', 'p=1', 'p=1.1', 'Analytical'};
legends_a = {'a=2', 'a=1/2', 'a=-1/2', 'Analytical'};
legends_dt = {'delta t = 1', 'delta t = 0.1', 'delta t = 0.01', 'Analytical'};

% Manipulating variable 'p'
figure('Name', 'Manipulating variable p');
hold on
for j = 1:3
    U_p = solve_euler(u0, T, a, changing_p(j), delta_t);
    plot(T, U_p, plotstyles{j},'LineWidth', 1.5)
end
plot(T, U_True, plotstyles{4},'LineWidth', 1.0)
legend(legends_p, 'Location', 'northwest')
title('Manipulating variable p')
xlabel('time, t'); ylabel('u')
hold off

% Manipulating variable 'a'
figure('Name', 'Manipulating variable a');
hold on
for j = 1:3
    U_a = solve_euler(u0, T, changing_a(j), p, delta_t);
    plot(T, U_a, plotstyles{j},'LineWidth', 1.5)
end
plot(T, U_True, plotstyles{4},'LineWidth', 1.0)
legend(legends_a, 'Location', 'northwest')
title('Manipulating variable a')
xlabel('time, t'); ylabel('u')
hold off

% Manipulating variable 'delta t'
figure('Name', 'Manipulating variable delta t');
hold on
for j = 1:3
    changing_T = t0:changing_delta_t(j):tEnd;
    U_t = solve_euler(u0, changing_T, a, p, changing_delta_t(j));
    plot(changing_T, U_t, plotstyles{j},'LineWidth', 1.5)
end
plot(T, U_True, plotstyles{4},'LineWidth', 1.0)
legend(legends_dt, 'Location', 'northwest')
title('Manipulating variable delta t')
xlabel('time, t'); ylabel('u')
hold off

% Function to solve ODE using Euler's method
function U = solve_euler(u0, t, a, p, delta_t)
    U = zeros(size(t));
    U(1) = u0;
    for i = 1:length(t)-1
        U(i+1) = U(i) + delta_t * a * U(i)^p;
    end
end
