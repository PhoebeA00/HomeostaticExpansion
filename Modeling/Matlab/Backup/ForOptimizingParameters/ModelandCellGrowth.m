close all; clear all; clc
global tx N T R I pOpt p0 m
%global alpha epsilon a c b_R mu beta g b_T d e_T e_R f kA n
global epsilon c b_R mu beta g d f n


%TO DO - setup global parameters and grab the ones I need from here to be
%accessed in the growth model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   CHECK LIST                             %
% 1 - Change global vectors in: ModelandCellGrowth.m(Line 4),
% Growth.m(Line3)
%   1a - Add parameters that have been removed in previous experiment.
%        Copy and paste the commented parameter list
%   1b - Remove parameter for this experiment
% 2 - Change parameters to be optimized here
%   2a - Add parameters to p0 and set limits 
%   2b - Comment out the parameters to be explored in p0 below
%   2c - Remove % from parameters not being explored
% 3 - Choose the right set of parameters to be explored in the 
%       'GetParameters()' function (Line 27)   
% 4 - Change parameters to be optimized in Growth.m
% 5 - Depending on what you are optimizing for, choose the right
%       data sets for optimization in the GrowthObjective.m file
% 6 - Change the disp at the bottom to get the optimized values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                            %
p = GetParameters(22);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order of Parameters                                           %  
% 1_alpha 2_Thy 3_Thy_max 4_epsilon 5_a 6_c 7_b_R 8_mu 9_beta   %
% 10_g 11_b_T 12_d 13_e_T 14_e_R 15_f 16_kA 17_n 18_EntryNumber %
% 19_Notes 20_Naive 21_Activated 22_Treg 23_IL2 24_PreviousPset %
%     

% Parameters to be optimized
a = p(5); %2+rand(1,1)*(6 -2)
b_T = p(11);
e_T = p(13);
e_R = p(14);
kA = p(16);
alpha = p(1);

p0 = [a, b_T, e_T, e_R, kA, alpha];

% optimize with fmincon
lb = [0, 0, 0, 0, 0, 0] ; %[] lower bound
ub = [0.0000001, 1, 10, 10, 1000000000, 4]; %[] upper bound


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%alpha =    p(1);
%Thy =      Thymus ODE controls this parameters, look in Growth.m
%Thy_max =  ThyEq(432);Deprecated - K is now the Thy_max
epsilon =  p(4);
%a =        p(5);
c =        p(6);
b_R =      p(7);
mu =       p(8);
beta =     p(9);
g =        p(10);
%b_T =      p(11);
d =        p(12);
%e_T =      p(13);
%e_R =      p(14);
f =        p(15);
%kA =       p(16);
n =        p(17);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initial Conditions
N = p(20); %Naive T cells
T = p(21); %Activated T Cells
R = p(22); %T Regulatory Cells
I = p(23); %IL-2 Cytokine
m = 0.0023; %Average of the Thymus weight at day 0

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
a_opt = pOpt(1);
b_T_opt = pOpt(2);
e_T_opt = pOpt(3);
e_R_opt = pOpt(4);
kA_opt = pOpt(5);
alpha_opt = pOpt(6);


disp('The optimized parameters are')
disp(['a: ' num2str(a_opt)])
disp(['b_T: ' num2str(b_T_opt)])
disp(['e_T: ' num2str(e_T_opt)])
disp(['e_R: ' num2str(e_R_opt)])
disp(['kA: ' num2str(kA_opt)])
disp(['alpha: ' num2str(alpha_opt)])

