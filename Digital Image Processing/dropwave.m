function [f,gradf] = dropwave(x)
% Description:
% Calculates dropwave function and its gradient at given point x
%
% Input:
%  x: given point
%
% Output:
%  f: dropwave function at point x
%  gradx: gradient of dropwave fuction at point x
%
% Usage:
% [f,gradf] = dropwave(x)
% Reference: https://www.sfu.ca/~ssurjano/drop.html 
x_1 = x(1);
x_2 = x(2);
frac1 = 1 + cos(12*sqrt(x_1^2+x_2^2));
frac2 = 0.5*(x_1^2+x_2^2) + 2;
f = -frac1/frac2;
g_x1 = (12*x_1*sin(12*sqrt(x_1^2+x_2^2)))/(sqrt(x_1^2+x_2^2)*...
    ((x_1^2+x_2^2)/2+2))-(x_1*(-cos(12*sqrt(x_1^2+x_2^2))-1))/((x_1^2+x_2^2)/2+2)^2;
g_x2 = (12*x_2*sin(12*sqrt(x_2^2+x_1^2)))/(sqrt(x_2^2+x_1^2)*...
    ((x_2^2+x_1^2)/2+2))-(x_2*(-cos(12*sqrt(x_2^2+x_1^2))-1))/((x_2^2+x_1^2)/2+2)^2;
gradf = [g_x1,g_x2]';
end

