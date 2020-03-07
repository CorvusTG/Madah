function [] = Q01(T, N)
%solve the flame equation for T>0, y(0) = delta, with N (number of points
%as N), delta converges to 0
solveEquation(T, N, false);
title("results (Q01S01)");
%plot them both
figure

%solve the flame equation for T=2/delta, y(0) = delta, delta converges to 0
solveEquation(T, N, true);
title("results (Q01S02)");
end

function [] = solveEquation(T,N, flag)
%solve the flame equation
results = [];
deltas = [0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1];

if flag
    for i=1:length(deltas)
        results = [results, trapezoidal(deltas(i), 2/deltas(i), N)];
    end
else
    for i=1:length(deltas)
        results = [results, trapezoidal(deltas(i), T, N)];
    end
end
plot(deltas, results);
xlabel(strcat('delta = ', mat2str(deltas)));
ylabel('y(t)');
end

function [y_k_1] = trapezoidal(delta, T, N)
% This function receives a param delta,
%   and returns the vector of solutions to the ODE: y' = y^2-y^3 with
%   y(0)=0

% trapezoidal rule: y_k+1 = y_k + (h/2)*f(y_k)+(h/2)*f(y_k+1)
y0 = delta; % set initial value
h = (T/N); % discrete space
lambda = 0.000001; % desired accuracy
y_k = y0; % initial value
y_k_1 = do_iteration(y0,h); % do the first iteration of the rule

% run the trapezoidal rule loop
while abs(y_k_1 - y_k) >= lambda
    y_k = y_k_1;
    y_k_1 = do_iteration(y_k,h);
end
end

function [y] = do_iteration(y_last, h)
% does one iteration of the trapezoidal rule

% y' = y^2-Y^3, y(k+1)=y(k)+(h/2)(y(k)^2-y(k)^3)+(h/2)(y(k+1)^2-y(k+1)^3)
% => (h/2)y(k+1)^3-(h/2)y(k+1)^2+y(k)=y(k)+(h/2)(y(k)^2-y(k)^3)

d = y_last + (h/2)*(y_last^2-y_last^3); % d is the right side of the equation
y = [h/2 -h/2 1 -d]; % solve the third degree equation
roots_y = roots(y); % get the roots of the equation
y = roots_y(imag(roots_y) == 0); %take real answer
y = y(end);
end