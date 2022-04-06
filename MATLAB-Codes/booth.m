function [f,gradf] = booth(x)
% Description:
% Calculates booth function and its gradient at given point x
%
% Input:
%  x: given point
%
% Output:
%  f: booth function at point x
%  gradx: gradient of booth fuction at point x
%
% Usage:
% [f,gradf] = booth(x)
% Reference: https://www.sfu.ca/~ssurjano/booth.html
x_1 = x(1);
x_2 = x(2);
term1 = (x_1 + 2*x_2 - 7)^2;
term2 = (2*x_1 + x_2 - 5)^2;
f = term1 + term2;
g_x_1 = 10*x_1 + 8*x_2 -34; 
g_x_2 = 8*x_1 + 10*x_2 - 38;
gradf = [g_x_1,g_x_2]';
end