close all; clear all; clc
global tx N T R I pOpt p0
global alpha Thy Thy_max epsilon a c b_R mu beta g b_T d e_T e_R f kA n
global ThyEq 

%TO DO - setup global parameters and grab the ones I need from here to be
%accessed in the growth model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   CHECK LIST                             %
% 1 - Change parameters to be optimized here
% 2 - Change parameters to be optimized in Growth.m
% 2 - Lower bound and Upper Bound
% 3 - Data Used in growth objective
% 4 - Initial conditions  in Simulate Growth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order of Parameters                                           %  
% 1_alpha 2_Thy 3_Thy_max 4_epsilon 5_a 6_c 7_b_R 8_mu 9_beta   %
% 10_g 11_b_T 12_d 13_e_T 14_e_R 15_f 16_kA 17_n 18_EntryNumber %
% 19_Notes 20_Naive 21_Activated 22_Treg 23_IL2 24_PreviousPset %
%                                                               %
p = GetParameters(14);

% Parameters to be optimized
alpha = p(1); % Activation rate
p0 = alpha;

% optimize with fmincon
lb = 0; %[] lower bound
ub = 1000; %[] upper bound

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

%Values for the Thymus polynomial
p1 = -3.725e-09;
p2 = 2.419e-06;
p3 = -0.0001933;
p4 = 0.00398;

ThyEq = @(Mdlt) p1*Mdlt.^3 + p2*Mdlt.^2 + p3*Mdlt + p4;


%alpha =    p(1);
Thy =      p(2);
Thy_max =  ThyEq(432); %Maximum number of hours
epsilon =  p(4);
a =        p(5);
c =        p(6);
b_R =      p(7);
mu =       p(8);
beta =     p(9);
g =        p(10);
b_T =      p(11);
d =        p(12);
e_T =      p(13);
e_R =      p(14);
f =        p(15);
kA =       p(16);
n =        p(17);
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

alpha_1 = pOpt(1);
%b_T_1 = pOpt(2);
disp('The optimized parameters are')
disp(['alpha: ' num2str(alpha_1)])
%disp(['b_T: ' num2str(b_T_1)])










