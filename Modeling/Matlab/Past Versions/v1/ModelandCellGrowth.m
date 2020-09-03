close all; clear all; clc
global tx N T R I pOpt

%TO DO - setup global parameters and grab the ones I need from here to be
%accessed in the growth model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   CHECK LIST                             %
% 1 - Change parameters in growth
% 2 - Lower bound and Upper Bound
% 3 - Data Used in growth objective
% 4 - Initial conditions  in Simulate Growth
% 5 - Get Parameter values
%Change the parameters that will be fitted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Make all changes here for optimization %%%%%%%%%%%%%%%%%
% Order of Parameters                                           %  
% 1_alpha 2_Thy 3_Thy_max 4_epsilon 5_a 6_c 7_b_R 8_mu 9_beta   %
% 10_g 11_b_T 12_d 13_e_T 14_e_R 15_f 16_kA 17_n 18_EntryNumber %
% 19_Notes 20_Naive 21_Activated 22_Treg 23_IL2 24_PreviousPset %
%                                                               %
p = GetParameters(10);

% Parameters to be optimized
mu = p(8); % Naive T cell production rate
beta = p(9); % Activation rate 
g = p(10); % Naive T cell death rate
b_T = p(11); %Activated T cell death rate
p0 = [mu, beta, g, b_T];

% optimize with fmincon
lb = [0,0,0,0]; % lower bound
ub = [50000,10,10, 30]; % upper bound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initial Conditions
N = p(20);
T = p(21);
R = p(22);
I = p(23);


tx = 0:432; %Maximum amount of time - 18 days

%Pick the variables that will be fitted to
%Change the parameter choice in the Growth function
%Change the parameters to be optimized in the Growth function

% optimize parameters

% no linear constraints
A = [];
b = [];
Aeq = [];
beq = [];
nlcon = [];

disp('Beginning Optimization...')
% options = optimoptions(@fmincon,'Algorithm','interior-point');
pOpt = fmincon(@GrowthObjective,p0,A,b,Aeq,beq,lb,ub,nlcon);


%Optimized Parameters

mu = pOpt(1);
beta = pOpt(2);
disp('The optimized parameters are')
disp(['mu: ' num2str(mu)])
disp(['beta: ' num2str(beta)])











