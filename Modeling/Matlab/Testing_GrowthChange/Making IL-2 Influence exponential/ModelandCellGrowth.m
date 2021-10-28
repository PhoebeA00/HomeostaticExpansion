close all; clc
%global tx N T R I ThyR RplR DiffR pOpt p0 m
global tx WTerror KOerror

% Parameter Ranges
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%-------------------- Naive -----------------------%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%

mu_min  = 0.415; %0.28 %Naive T production rate
mu_max = 0.415;

nK_min = 8890464;
nK_max = 8890464; % 85% of the maxmimum cd4CT of D56 mice, assuming 15% activation

z_min = 0.025; %Naive Self replication rate
z_max = 0.025;

g_min = 0.069944693966199;%Death rate of Naive
g_max = 0.069944693966199;

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%-------------------- Tregs -----------------------%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%

alpha_min = 0.0001; %Thymic derived Tregs
alpha_max = 0.03;

rK_min = 10459370; %mean of D56 Tregs
rK_max = 10459370; 

c_min = 0.002072865287658; %Naive Derived Tregs
c_max = 0.002072865287658;

epsilon_min = 0.018;%Self Replication rate of Tregs
epsilon_max = 0.018;

b_R_min = 0.057339990947437;%Death Rate of Tregs
b_R_max = 0.057339990947437;%0.79857;

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-%
%--------------------Activated T -----------------------%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-%

beta_min =  0.325;%0.102; %activation rate
beta_max = 0.325;

a_min = 0.006600026244555;%0.0076; %Self Replication rate for activated T cells
a_max = 0.006600026244555;

b_T_min = 0.034087191107056;%Death Rate of T cells
b_T_max = 0.034087191107056;

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-%
%-------------------- Consumption Rates -----------------------%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-%

e_T_min = 100; %Consumption rate of T cells
e_T_max = 100; %activated T cells consume at max 1/10 of what they make

e_R_min = 200;%393; %Consumption rate of Tregs
e_R_max = 200;

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-%
%--------------------  Suppression -----------------------%
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-%

kA_min = 314120; %Half suppression rate by Tregs
kA_max = 314120;

Ki_min = 0; %Half rate for activation suppression boost
Ki_max = 10;

j_min = 2.215568109662751e-06; %Rate of deactivation of activated T cells
j_max = 2.215568109662751e-06;

Kj_min = 0; % Half ratde for deactivation boost
Kj_max = 10;

kB_min = 0;%half suppression rate of Treg death rate
kB_max = 10;

n_min = 1; %Controls the sigmoidicity - Hill coefficient
n_max = 1;

d_min = 1000; %IL-2 Production
d_max = 1000; 

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
c = c_min + rand(1,1) * (c_max - c_min);
kB = kB_min + rand(1,1) * (kB_max - kB_min);
j = j_min + rand(1,1) * (j_max - j_min);
z = z_min + rand(1,1) * (z_max - z_min);
n = n_min + rand(1,1) * (n_max - n_min);
d = d_min + rand(1,1) * (d_max - d_min);
nK = nK_min + rand(1,1) * (nK_max - nK_min);
rK = rK_min + rand(1,1) * (rK_max - rK_min);
Ki = Ki_min + rand(1,1) * (Ki_max - Ki_min);
Kj = Kj_min + rand(1,1) * (Kj_max - Kj_min);

%Making sure that the consumption of Tregs is greater than that of
%activated T cell
e_T = e_T_min + rand(1,1) * (e_T_max - e_T_min);
e_R = e_R_min + rand(1,1) * (e_R_max - e_R_min);

while e_T > e_R
    e_T = e_T_min + rand(1,1) * (e_T_max - e_T_min);
    e_R = e_R_min + rand(1,1) * (e_R_max - e_R_min);
end

%------fmincon function arguments definitions-------
p0 = [alpha, a, kA, e_T, e_R, g, b_T, b_R, epsilon, mu, beta, c, kB, j, z, n, d, nK, rK, Ki, Kj];
lb = [alpha_min, a_min, kA_min, e_T_min, e_R_min, g_min, ...
    b_T_min, b_R_min, epsilon_min, mu_min, beta_min, c_min, kB_min, j_min, z_min, n_min, d_min,...
    nK_min, rK_min, Ki_min, Kj_min]; %[] lower bound
ub = [alpha_max, a_max, kA_max, e_T_max, e_R_max, g_max, ...
    b_T_max, b_R_max, epsilon_max, mu_max, beta_max, c_max, kB_max, j_max, z_max, n_max, d_max,...
    nK_max, rK_max, Ki_max, Kj_max]; %[] upper bound

% no linear constraints
A = [0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
b = 0;
Aeq = [];
beq = [];
nlcon = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------Initial Conditions-----%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tx = 0:432; %Maximum amount of time - 18 days


%-------------------------------------------------------------------------------------------%
%                 Preparing the data for the growth objective
%-------------------------------------------------------------------------------------------%

parameterStruct = {};
parameterStruct.WildType.CellData = readtable('../../RawData/ActivatedWTSpleen.csv');
parameterStruct.KnockOut.CellData = readtable('../../RawData/ActivatedKOSpleen.csv');
parameterStruct.WildType.ProlData = readtable('../../RawData/WTProl.csv');
parameterStruct.KnockOut.ProlData = readtable('../../RawData/KOProl.csv');

%Reducing the number of columns and organizing the populations to be
%iterated over
parameterStruct.WildType.CellData = parameterStruct.WildType.CellData(:,{'NaiveCT', 'ActivatedCD4CT', 'X4TregCT', ...
    'ThymicNaive', 'ActivatedNaiveCT', ...
    'ThymicDerivedTregsCT', 'NaiveDerivedTregsCT' ... 
    'hours'});
parameterStruct.KnockOut.CellData = parameterStruct.KnockOut.CellData(:,{'NaiveCT', 'ActivatedCD4CT', 'X4TregCT', ...
    'ThymicNaive', 'ActivatedNaiveCT', ...
    'ThymicDerivedTregsCT', 'NaiveDerivedTregsCT' ... 
    'hours'});
parameterStruct.WildType.ProlData = parameterStruct.WildType.ProlData(:,{ 'NaiveProlCT', 'ActivatedProlCT', 'X4TregProlCT', ...
    'hours'});
parameterStruct.KnockOut.ProlData = parameterStruct.KnockOut.ProlData(:,{ 'NaiveProlCT', 'ActivatedProlCT', 'X4TregProlCT', ...
    'hours'});

% optimize parameters
disp('Beginning Optimization...')
[pOptimized, error] = fmincon(@(p) GrowthObjective(p, parameterStruct), p0, A,b,Aeq,beq,lb,ub,nlcon);
disp('...Ending Optimization')
%Optimized Parameters

Genotype = [1, 2];
for i = Genotype
    PlottingResults(pOptimized, i)
end

PlottingEverything(pOptimized)

alpha = pOptimized(1);
a = pOptimized(2);
kA = pOptimized(3);
e_T = pOptimized(4);
e_R = pOptimized(5);
g = pOptimized(6);
b_T = pOptimized(7);
b_R = pOptimized(8);
epsilon = pOptimized(9);
mu = pOptimized(10);
beta = pOptimized(11);
c = pOptimized(12);
kB = pOptimized(13);
j = pOptimized(14);
z = pOptimized(15);
n = pOptimized(16);
d = pOptimized(17);
nK = pOptimized(18);
rK = pOptimized(19);
Ki = pOptimized(20);
Kj = pOptimized(21);
%%
%-----Change this for saving files in a different location-----%
FileLocation = './Data/ParameterSets.csv';
%---------------------------------------------------------------------------%
ParameterData = readtable(FileLocation);

if isempty(ParameterData.EntryNumber)
    EntryNumber = 1;
else
    EntryNumber = max(ParameterData.EntryNumber)+1;
end
%d = 1000; %IL-2 production Rate
f = 1.38629; %IL-2 degradation Rate
parameters = [mu, z, g, alpha, c, epsilon, b_R, beta, a, b_T,...
    e_T, e_R, kA, j, kB, n, d, f, nK, rK, Ki, Kj, error, WTerror, KOerror, EntryNumber];


dlmwrite(FileLocation ,parameters,'delimiter', ',',  'precision', 16,'-append');

disp(['Entry Number - ' num2str(EntryNumber)])












