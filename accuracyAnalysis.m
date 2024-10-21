% EME 451 COMPUTATIONAL FLUID DYNAMICS
% ASSIGNMENT 1 
% HOMEWORK 1 PROBLEM 2
%%
% VINOD RAO A/L JAYAPRASADH (SCHOOL OF CHEMICAL ENGINEERING)
% MATRIC NUMBER: 158635
% NURUL AIN FAZWIN BINTI MOHAMAD SAKMAH (SCHOOL OF MECHANICAL ENGINEERING)
% MATRIC NUMBER: 153477
%%
% Code 2: Accuracy Analysis Implementation

% Clear workspace and figures
clc;
close all;

% Given conditions
t0 = 0;
u0 = 1;
tEnd = 2;
a = 2;
p = 1;

% Analytical solution
U_True = u0*exp(a*tEnd);

% Initializing first order solutions
N1 = 1;
U_1(1) = u0;
error1 = 1;

% 1% accuracy
while error1>0.01
    N1 = N1+1;
    delta_t1 = (tEnd-t0)/N1;
    for i = 1:N1
        U_1(i+1) = U_1(i)*(1+delta_t1^p*a);
    end
    error1 = abs((U_True-U_1(end))/U_True);
end

N1a = N1;
error1a = error1;
delta_t1a = delta_t1;
OOA1a = log10(error1)/log10(delta_t1);

% 0.01% accuracy
while error1>0.0001
    N1 = N1+1;
    delta_t1 = (tEnd-t0)/N1;
    for i = 1:N1
        U_1(i+1) = U_1(i)*(1+delta_t1^p*a);
    end
    error1 = abs((U_True-U_1(end))/U_True);
end

N1b = N1;
error1b = error1;
delta_t1b = delta_t1;
OOA1b = log10(error1)/log10(delta_t1);

% Initializing second order solutions
N2 = 1;
U_2(1) = u0;
error2 = 1;

% 1% accuracy
while error2>0.01
    N2 = N2+1;
    delta_t2 = (tEnd-t0)/N2;
    for i = 1:N2
        U_2(i+1) = U_2(i)*(1+delta_t2*a*(1+delta_t2*a/2));
    end
    error2 = abs((U_True-U_2(end))/U_True);
end

N2a = N2;
error2a = error2;
delta_t2a = delta_t2;
OOA2a = log10(error2)/log10(delta_t2);

% 0.01% accuracy
while error2>0.0001
    N2 = N2+1;
    delta_t2 = (tEnd-t0)/N2;
    for i = 1:N2
        U_2(i+1) = U_2(i)*(1+delta_t2*a*(1+delta_t2*a/2));
    end
    error2 = abs((U_True-U_2(end))/U_True);
end

N2b = N2;
error2b = error2;
delta_t2b = delta_t2;
OOA2b = log10(error2)/log10(delta_t2);

% Displaying the results
fprintf('\nNumerical Methods Comparison:\n');
fprintf('---------------------------------------------------------------------------\n');
fprintf('Accuracy  Method    Interval Number    Time Step        Order of Accuracy\n');
fprintf('---------------------------------------------------------------------------\n');
fprintf('1%%        First     %d                %.2e           %.4f\n', N1a, delta_t1a, OOA1a);
fprintf('          Second    %d                 %.2e           %.4f\n', N2a, delta_t2a, OOA2a);
fprintf('---------------------------------------------------------------------------\n');
fprintf('0.01%%     First     %d              %.2e           %.4f\n', N1b, delta_t1b, OOA1b);
fprintf('          Second    %d                %.2e           %.4f\n', N2b, delta_t2b, OOA2b);
fprintf('---------------------------------------------------------------------------\n');

% Error Convergence plotting
figure('Name', 'Error Convergence');
dt_values = [delta_t1a delta_t1b; delta_t2a delta_t2b];
error_values = [error1a error1b; error2a error2b];
loglog(dt_values(1,:), error_values(1,:), 'bo-', ...
    dt_values(2,:), error_values(2,:), 'rs-', 'LineWidth', 1.5);
grid on;
xlabel('log \Delta t');
ylabel('log Error');
title('Error Convergence Analysis');
legend('First Order', 'Second Order', 'Location', 'northwest');