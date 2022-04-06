function [f,gradx] = shubert(x)
% Description:
% Calculates shubert function and its gradient at given point x
%
% Input:
%  x: given point
%
% Output:
%  f: shubert function at point x
%  gradx: gradient of shubert fuction at point x
%
% Usage:
% [f,gradf] = shubert(x)
% Reference: https://www.sfu.ca/~ssurjano/shubert.html
x_1 = x(1);
x_2 = x(2);
term1 = (cos(2*x_1+1)+2*cos(3*x_1+2)+3*cos(4*x_1+3)+4*cos(5*x_1+4)+5*cos(6*x_1+5));
term2 = (cos(2*x_2+1)+2*cos(3*x_2+2)+3*cos(4*x_2+3)+4*cos(5*x_2+4)+5*cos(6*x_2+5));
f = (term1)*(term2);
g1 = -2*(5*cos(6*x_2+5)+4*cos(5*x_2+4)+3*cos(4*x_2+3)+2*cos(3*x_2+2)+...
    cos(2*x_2+1))*(15*sin(6*x_1+5)+10*sin(5*x_1+4)+6*sin(4*x_1+3)+...
    3*sin(3*x_1+2)+sin(2*x_1+1));
g2 = -2*(5*cos(6*x_1+5)+4*cos(5*x_1+4)+3*cos(4*x_1+3)+2*cos(3*x_1+2)+...
    cos(2*x_1+1))*(15*sin(6*x_2+5)+10*sin(5*x_2+4)+6*sin(4*x_2+3)+...
    3*sin(3*x_2+2)+sin(2*x_2+1));
gradx = [g1,g2]';
end