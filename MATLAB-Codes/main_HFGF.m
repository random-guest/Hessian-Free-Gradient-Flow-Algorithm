function [f_x_final,n_iter,n_fev,n_gev] = main_HFGF(f,x0,tol,imax,eps,c,beta,amax)
%
% Author: Name:        Abdulkarim Atrash 
%         E-mail:      atrash.abdulkarim@metu.edu.tr
%         Address:     Middle East Technical University, Ankara, Turkey
%         Department:  Institute of Applied Mathematics, Scientific Computing Program
%
% Desctiption: 
% 
% An implementation of the Hessian Free Gradient Flow 
% algorithm presented in the following research paper:
% https://www.sciencedirect.com/science/article/pii/S0098135420303562  
% 
% Inputs: 
%     f    : The testing function under studies.
%     x0   :  Initial guess
%     tol  :  Tolerance
%     imax :  Max. # of iterations for the inner loop (solves Ax=b)
%     eps  :  Used to approx. the Hessian Matrix.
%     c    :  Parameter in finding the alpha that satisfies Armijo First Order Conditions
%     beta :  Parameter in finding the alpha that satisfies Armijo First Order Conditions
%     amax :  Maximum number of iterations in finding the alpha that satisfies Armijo First Order Conditions
%
% Outputs:
%
%     f_x_final     : Function Evaluated at the final xk found
%     n_iter        : The total number of iterations in the major loop
%     n_fev         : The number offunction evaluations
%     n_gev         : The number of gradient evaluations
%
% Usage: 
%   [f_x_final,n_iter,n_fev,n_gev] = main_HFGF(f,x0,tol,imax,eps,c,beta,amax)
%
%% First Initialize
p0 = zeros(2,1);   % Zero Vector, from truncated Newton Method, whose size is the dimension of the testing function (starting by 2D functions)
x0 = x0';          % Convert the initial guess to column vector
h0 = 10^(-3);      % step-size, A value between 0 and 1, chosen as specified in the paper
n_iter = 0;        % To get the number of iterations needed in the major loop
n_fev  = 0;        % To get the number of function evaluation
n_gev  = 0;        % To get the number of Gradient evaluation
%% Check the condition for the while loop
[~,g0] = f(x0);    % As we only need the gradient, where f & g  are the function & gradient evaluated at x0
N = norm(g0);      % Take the norm of g to check the condition of the while loop
n_gev = n_gev+1;   % Update the number of gradient evaluation
%% Start the loop
while N > tol 
    % Evalue the gradient at xk and at xk+eps*pk
    [fk,gk] = f(x0);
    n_fev = n_fev+1;                         % update the number of function evaluation
    n_gev = n_gev+1;                         % Update the number of gradient evaluation
    [~,gk_eps] = f(x0+eps*p0);
    n_fev = n_fev+1;                         % update the number of function evaluation
    n_gev = n_gev+1;                         % Update the number of gradient evaluation
    QkPk = (gk_eps - gk )/eps + (1/h0)*p0;   % Instead of multiplying with I, I transposed p0
    % Second Initializ
    r0 = QkPk + gk; d0 = -r0; pi = p0;
    % Start the inner for loop
    for i=1:imax
        % evaluate the gradient of the function at the new point, xk+eps*di
        [~,gk] = f(x0);                      % Note that this should be somewhere at the end!
        n_gev = n_gev+1;                     % Update the number of gradient evaluation
        n_fev = n_fev+1;                     % update the number of function evaluation
        [~,gk_eps] = f(x0+eps*d0);
        n_gev = n_gev+1;                     % Update the number of gradient evaluation
        n_fev = n_fev+1;                     % update the number of function evaluation
        Qkdi = (gk_eps - gk)/eps + (1/h0)*d0;
        alpha_i = (r0'*r0)/(d0'*Qkdi);
        p_i_plus_1 = pi + alpha_i*d0;
        r_i_plus_1 = r0 + alpha_i*d0 ;
        beta_i_plus_1  = (r_i_plus_1'*r_i_plus_1)/(r0'*r0);
        d_i_plus_1 = -r_i_plus_1 + beta_i_plus_1 * d0; % NOT USED !
        x_try = x0 + p_i_plus_1;
        % Evaluate the function at x_try
        [f_try,~] = f(x_try);
        n_fev = n_fev+1;                      % update the number of function evaluation
        n_gev = n_gev+1;                      % Update the number of gradient evaluation
        if f_try < fk
            % Check armijo first order condition, 
            % first find new alpha, named as alpha_Armijo
            [ alpha_Armijo ] = armijo(@dropwave, x_try, p_i_plus_1, alpha_i, c, beta, amax);
            if f_try < fk + alpha_Armijo*(gk'*p_i_plus_1)
                x0 = x_try;
                p0 = p_i_plus_1;
                h0 = 2*h0;  
                break
            elseif i == imax
                x0 = x_try;
                p0 = p_i_plus_1;
                h0 = 0.5*h0;
                break
            end
        elseif i == imax
            h0 = 0.5*h0;
            break       
        end
    end
    % Update the norm of the gradient 
    [~,gk] = f(x0);
    n_fev = n_fev+1;                    % update the number of function evaluation
    n_gev = n_gev+1;                    % Update the number of gradient evaluation
    N = norm(gk);
    n_iter = n_iter + 1;                % Update the number of iterations 
end
%% Evaluate the function understudies at x final
[f_x_final,~] = f(x0);
end