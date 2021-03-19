close all; clear all; clc
global tx N T R I pOpt p0 m
%global alpha epsilon a c b_R mu beta g b_T d e_T e_R f kA n
global n d f

%----------------------------------------------------------------%
%------------------------CHECK LIST-------------------------%
%----------------------------------------------------------------%

%{
Where my code is located:
/home/jon/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab
1 - Change global vectors in:                                  
                      ModelandCellGrowth.m(Line 4),                              
                      Growth.m(Line3)                                            
  1a - Add parameters that have been removed.                  
       Copy and paste the commented parameter list
  1b - Remove parameter for this experiment
2 - Change parameters to be optimized here
  2a - Add parameters to p0 and set limits 
  2b - Comment out the parameters to be explored in p0 below
  2c - Remove % from parameters not being explored
3 - Change parameters to be optimized in Growth.m
4 - Depending on what you are optimizing for, choose the right
       data sets for optimization in the GrowthObjective.m file
5 - Change Paramter csv file to save in the right one
    5a: There are two here
    5b: When plotting there is one spot in GetParameters script

%}


% Parameter Ranges
alpha_min = 33914; %Thymic derived Tregs
alpha_max = 33914; 

a_min = 0.0001; %Self Replication rate for activated T cells
a_max = 0.2;

epsilon_min = 0;%Self Replication rate of Tregs
epsilon_max = 0.1895;

kA_min = 0; %Half suppression rate by Tregs
kA_max = 1000000; 

e_T_min = 100; %Consumption rate of T cells
e_T_max = 100; %activated T cells consume at max 1/10 of what they make

e_R_min = 1.3744E-09; %Consumption rate of Tregs
e_R_max = 2000;

g_min = 0;%Death rate of Naive
g_max = 0;

b_T_min = 0.2;%Death Rate of T cells
b_T_max = 0.5;

b_R_min = 0.2927;%Death Rate of Tregs
b_R_max = 0.2927;

mu_min =   100000; %Naive T production rate
mu_max = 100000; 

beta_min = 0.3; %activation rate
beta_max = 0.3;

c_min = 0.0206; %Naive differentiation to Tregs
c_max = 0.0206;

%j_min = 0; %Rate of desctruction of activated T cells
%j_max = 100; 

kB_min = 0; %half suppression rate of Treg death rate
kB_max = 100000;

%Randomizing the initial parameter choices
alpha = alpha_min + rand(1,1) * (alpha_max - alpha_min);
a =   a_min + rand(1,1) * (a_max - a_min);
kA =  kA_min + rand(1,1) * (kA_max - kA_min);
g = g_min + rand(1,1) * (g_max - g_min);
b_T = b_T_min + rand(1,1) * (b_T_max - b_T_min);
b_R = b_R_min + rand(1,1) * (b_T_max - b_T_min);
epsilon = epsilon_min + rand(1,1) * (epsilon_max - epsilon_min);
mu = mu_min + rand(1,1) * (mu_max - mu_min);
beta = beta_min + rand(1,1) * (beta_max - beta_min);
c = c_min + rand(1,1) * (c_max-c_min);
kB = kB_min + rand(1,1)*(kB_max-kB_min);
%Making sure that the consumption of Tregs is greater than that of
%activated T cell
e_T = e_T_min + rand(1,1) * (e_T_max - e_T_min);
e_R = e_R_min + rand(1,1) * (e_R_max - e_R_min);

while e_T > e_R
    e_T = e_T_min + rand(1,1) * (e_T_max - e_T_min);
    e_R = e_R_min + rand(1,1) * (e_R_max - e_R_min);
end

%------fmincon function arguments definitions-------
p0 = [alpha, a, kA, e_T, e_R, g, b_T, b_R, epsilon, mu, beta, c, kB];
lb = [alpha_min, a_min, kA_min, e_T_min, e_R_min, g_min, ...
    b_T_min, b_R_min, epsilon_min, mu_min, beta_min, c_min, kB_min]; %[] lower bound
ub = [alpha_max, a_max, kA_max, e_T_max, e_R_max, g_max, ...
    b_T_max, b_R_max, epsilon_max, mu_max, beta_max, c_max, kB_max]; %[] upper bound

% no linear constraints
A = [0 0 0 1 -1 0 0 0 0 0 0 0 0];
b = 0;
Aeq = [];
beq = [];
nlcon = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------Fixed Parameters-------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1; %Hill coefficient
d = 1000; %IL-2 production Rate
f = 1.38629; %IL-2 degradation Rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------Initial Conditions-----%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 10000; %Naive T cells
T = 0; %Activated T Cells
R = 300; %T Regulatory Cells
I = 0; %IL-2 Cytokine
m = 0.0023; %Average of the Thymus weight at day 0

tx = 0:432; %Maximum amount of time - 18 days

% optimize parameters
disp('Beginning Optimization...')
% options = optimoptions(@fmincon,'Algorithm','interior-point');
[pOpt, error] = fmincon(@GrowthObjective,p0,A,b,Aeq,beq,lb,ub,nlcon);
disp('...Ending Optimization')
%Optimized Parameters
alpha = pOpt(1);
a = pOpt(2);
kA = pOpt(3);
e_T= pOpt(4);
e_R= pOpt(5);
g= pOpt(6);
b_T= pOpt(7);
b_R= pOpt(8);
epsilon = pOpt(9);
mu = pOpt(10);
beta = pOpt(11);
c = pOpt(12);
kB = pOpt(13);

%-----Change this for saving files in a different location-----%
FileLocation = '../Data/ParameterRanges27.csv';
%---------------------------------------------------------------------------%
ParameterData = readtable(FileLocation);

if isempty(ParameterData.EntryNumber)
    EntryNumber = 1;
else
    EntryNumber = max(ParameterData.EntryNumber)+1;
end

parameters = [mu, beta, c, epsilon, n, d, f, alpha, a,... 
    kA, e_T, e_R, g, b_T, b_R, kB, error, EntryNumber];

dlmwrite(FileLocation ,parameters,'delimiter',',','-append');








