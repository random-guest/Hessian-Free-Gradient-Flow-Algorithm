%
% Author: Name:        Abdulkarim Atrash 
%         E-mail:      atrash.abdulkarim@metu.edu.tr
%         Address:     Middle East Technical University, Ankara, Turkey
%         Department:  Institute of Applied Mathematics, Scientific Computing Program
%
% Desctiption: 
% This Script will be responsible of creating the outputs in form of table
% for comparision with the reseacrh paper's ones.
% 
% Inputs: None , as they are defined bellow, if the user wishes to change
% any parameter, they can manually modify them down below
%
% Outputs:
% Table of results (to be compared against table 2 in the main reseacrh paper)
%
% Usage: Just run this script
%
%% Clear the working environment
clear; close all; clc;
%% Define the inputs for the main function, which is main_HFGF
x0           = [2,2];
tol          = 10^(-6);
imax         = 20;
eps          = 10^(-4);
c            = 10^(-4);
beta         = 0.5;
amax         = 100;
%% Call the three testing functions to get the results needed to construct the table
[f_dropwave,n_iter_dropwave,n_fev_dropwave,n_gev_dropwave]     = main_HFGF(@dropwave,x0,tol,imax,eps,c,beta,amax);
[f_shubert,n_iter_shubert,n_fev_shubert,n_gev_shubert]         = main_HFGF(@shubert,x0,tol,imax,eps,c,beta,amax);
[f_booth,n_iter_booth,n_fev_booth,n_gev_booth]                 = main_HFGF(@booth,x0,tol,imax,eps,c,beta,amax);
[f_rosenbrock,n_ite_rosenbrock,n_fev_rosenbrock,n_gev_rosenbrock]  = main_HFGF(@rosenbrock,x0,tol,imax,eps,c,beta,amax);
%% Define the components of the table
testing_Fns              = {'dropwave','shubert','booth','Rosenbrock'}; 
method                   = {'HFGF','HFGF','HFGF','HFGF'};
f_x                      = [f_dropwave,f_shubert,f_booth,f_rosenbrock];
n_it                     = [n_iter_dropwave,n_iter_shubert,n_iter_booth,n_ite_rosenbrock];
n_fv                     = [n_fev_dropwave,n_fev_shubert,n_fev_booth,n_fev_rosenbrock];
n_gv                     = [n_gev_dropwave,n_gev_shubert,n_gev_booth,n_gev_rosenbrock];
Function                 = testing_Fns';
Optimization_method      = method';
f_x_final                = f_x';
n_iter                   =  n_it';
n_fev                    = n_fv';
n_gev                    = n_gv';
%% Construct the Table
table(Function,Optimization_method,f_x_final,n_iter,n_fev,n_gev)
%%
