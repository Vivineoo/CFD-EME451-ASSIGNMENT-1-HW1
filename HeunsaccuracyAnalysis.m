% EME 451 COMPUTATIONAL FLUID DYNAMICS
% ASSIGNMENT 1 
% HOMEWORK 1 PROBLEM 2 (EXTRA CREDIT 10%)
%%
% VINOD RAO A/L JAYAPRASADH (SCHOOL OF CHEMICAL ENGINEERING)
% MATRIC NUMBER: 158635
% NURUL AIN FAZWIN BINTI MOHAMAD SAKMAH (SCHOOL OF MECHANICAL ENGINEERING)
% MATRIC NUMBER: 153477
%%
% Extra Credit (10%): Heun's Method Accuracy Analysis Implementation

% Clear workspace and figures
clc;
clear all;
close all;

% Given conditions
t0 = 0;
u0 = 1;
tEnd = 2;
a = 2;
p = 1;

% Analytical solution
U_True = u0*exp(a*tEnd);

% First Order Method
N1 = 1;
U_1(1) = u0;
error1 = 1;

% 1% and 0.01% accuracy for First Order
accuracies = [0.01, 0.0001];

for acc = 1:2
    while error1 > accuracies(acc)
        N1 = N1 + 1;
        delta_t1 = (tEnd-t0)/N1;
        for i = 1:N1
            U_1(i+1) = U_1(i)*(1+delta_t1^p*a);
        end
        error1 = abs((U_True-U_1(end))/U_True);
    end
    N1_values(acc) = N1;
    error1_values(acc) = error1;
    delta_t1_values(acc) = delta_t1;
end

% Second Order Method
N2 = 1;
U_2(1) = u0;
error2 = 1;

% 1% and 0.01% accuracy for Second Order
for acc = 1:2
    while error2 > accuracies(acc)
        N2 = N2 + 1;
        delta_t2 = (tEnd-t0)/N2;
        for i = 1:N2
            U_2(i+1) = U_2(i)*(1+delta_t2*a*(1+delta_t2*a/2));
        end
        error2 = abs((U_True-U_2(end))/U_True);
    end
    N2_values(acc) = N2;
    error2_values(acc) = error2;
    delta_t2_values(acc) = delta_t2;
end

% Heun's Method
N3 = 1;
U_3(1) = u0;
error3 = 1;

% 1% and 0.01% accuracy for Heun's Method
for acc = 1:2
    while error3 > accuracies(acc)
        N3 = N3 + 1;
        delta_t3 = (tEnd-t0)/N3;
        for i = 1:N3
            % Predictor (Euler step)
            k1 = a * U_3(i);
            u_pred = U_3(i) + delta_t3 * k1;
            % Corrector (Heun's formula)
            k2 = a * u_pred;
            U_3(i+1) = U_3(i) + (delta_t3/2) * (k1 + k2);
        end
        error3 = abs((U_True-U_3(end))/U_True);
    end
    N3_values(acc) = N3;
    error3_values(acc) = error3;
    delta_t3_values(acc) = delta_t3;

end

% Calculate Order of Accuracy for Heun's Method
OOA3a = log10(error3_values(1))/log10(delta_t3_values(1));
OOA3b = log10(error3_values(2))/log10(delta_t3_values(2));

% Displaying the results for Heun's Method
fprintf('\nHeuns Method Results:\n');
fprintf('---------------------------------------------------------------------------\n');
fprintf('Accuracy  Method    Interval Number    Time Step        Order of Accuracy\n');
fprintf('---------------------------------------------------------------------------\n');
fprintf('1%%        Heun      %d                 %.2e         %.4f\n', N3_values(1), delta_t3_values(1), OOA3a);
fprintf('0.01%%     Heun      %d                %.2e         %.4f\n', N3_values(2), delta_t3_values(2), OOA3b);
fprintf('---------------------------------------------------------------------------\n');

% Error Convergence plotting for all three methods
figure('Name', 'Error Convergence Comparison');
loglog(delta_t1_values, error1_values, 'bo-', ...
       delta_t2_values, error2_values, 'rs-', ...
       delta_t3_values, error3_values, 'gd--', 'LineWidth', 1.5);
grid on;
xlabel('log \Delta t');
ylabel('log Error');
title('Error Convergence Analysis - Comparison of Methods');
legend('First Order', 'Second Order', 'Heuns Method', 'Location', 'northwest');