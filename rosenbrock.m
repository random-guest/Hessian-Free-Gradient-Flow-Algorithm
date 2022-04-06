function [ f,gradf ] = rosenbrock(x)

% Description:
% Calculates Rosenbrock function and its gradient at given point x
%
% Input:
%  x: given point
%
% Output:
%  f: Rosenbrock function at point x
%  gradx: gradient of Rosenbrock fuction at point x
%
% Usage:
% [f,gradf] = rosenbrock(x)
% reference: https://www.sfu.ca/~ssurjano/rosen.html#:~:text=The%20Rosenbrock%20function%2C%20also%20referred,in%20a%20narrow%2C%20parabolic%20valley.

f     = 100*(x(1)^2 - x(2))^2 + (x(1)-1)^2;
gradf = [100*(2*(x(1)^2-x(2))*2*x(1)) + 2*(x(1)-1)...
         ;100*(-2*(x(1)^2-x(2)))];

end