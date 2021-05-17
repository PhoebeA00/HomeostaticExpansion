close all; clear all; clc
global tx N T R I ThyR RplR DiffR pOpt p0 m
%global alpha epsilon a c b_R mu beta g b_T d e_T e_R f kA n
global n d f

%----------------------------------------------------------------%
%------------------------CHECK LIST-------------------------%
%----------------------------------------------------------------%

% Parameter Ranges
alpha_min = 1533.3;%33914; %Thymic derived Tregs
alpha_max = 1533.3;%33914; 

a_min = 0.001; %Self Replication rate for activated T cells
a_max = 0.2;

epsilon_min = 0.16;%Self Replication rate of Tregs
epsilon_max = 0.2063;

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

c_min = 0.1302; %Naive differentiation to Tregs
c_max = 0.1302;

j_min = 0; %Rate of desctruction of activated T cells
j_max = 0.00001; 

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
c = c_min + rand(1,1) * (c_max - c_min);
kB = kB_min + rand(1,1)*(kB_max - kB_min);
j = j_min + rand(1,1)*(j_max - j_min);

%Making sure that the consumption of Tregs is greater than that of
%activated T cell
e_T = e_T_min + rand(1,1) * (e_T_max - e_T_min);
e_R = e_R_min + rand(1,1) * (e_R_max - e_R_min);

while e_T > e_R
    e_T = e_T_min + rand(1,1) * (e_T_max - e_T_min);
    e_R = e_R_min + rand(1,1) * (e_R_max - e_R_min);
end

%------fmincon function arguments definitions-------
p0 = [alpha, a, kA, e_T, e_R, g, b_T, b_R, epsilon, mu, beta, c, kB, j];
lb = [alpha_min, a_min, kA_min, e_T_min, e_R_min, g_min, ...
    b_T_min, b_R_min, epsilon_min, mu_min, beta_min, c_min, kB_min, j_min]; %[] lower bound
ub = [alpha_max, a_max, kA_max, e_T_max, e_R_max, g_max, ...
    b_T_max, b_R_max, epsilon_max, mu_max, beta_max, c_max, kB_max, j_max]; %[] upper bound

% no linear constraints
A = [0 0 0 1 -1 0 0 0 0 0 0 0 0 0];
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

ThyR = 300; % Thymic Derivied Tregs
RplR  = 0; % Self Replicating Tregs
DiffR = 0; % Naive Derived Tregs

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
j = pOpt(14);

Data = readtable('../RawData/ActivatedWTSpleen.csv');
%CellData = Data(:,{'NaiveCT', 'ActivatedCD4CT', 'X4TregFromThymusCT', ...                  
%     'hours'}); 
CellData = Data(:,{'NaiveCT', 'ActivatedCD4CT', 'X4TregCT', ...                  
     'ThymicDerivedTregsCT', 'X4TregProlCT', 'NaiveDerivedTregsCT' ... 
     'hours'}); 
ModelData = SimulateGrowth(pOpt);

PLT = figure(1);

%Activation of naive
ModelData(:,9) =beta.* ModelData(:,1).*(1./(1+(ModelData(:,3)./kA).^n)); 
%Activated T Cell Self Replication
ModelData(:,10) = a.*ModelData(:,2); 

subplot(3,4,1)
scatter(CellData.hours, CellData.NaiveCT)
hold on 
plot(tx, ModelData(:,1))
title('Naive T Cells')
hold off

subplot(3,4,2)
scatter(CellData.hours, CellData.ActivatedCD4CT)
hold on 
plot(tx, ModelData(:,2))
title('Activated T Cells')
hold off

subplot(3,4,3)
plot(tx, ModelData(:,9), 'DisplayName', 'From Naive')
hold on
plot(tx, ModelData(:,10), 'DisplayName', 'Self Replication')
title('Other Activated T Cells')
legend('Location','northwest')
hold off     

subplot(3,4,4)
plot(tx, ModelData(:,7))
title('IL-2')

subplot(3,4,5)
scatter(CellData.hours, CellData.X4TregCT)
hold on 
plot(tx, ModelData(:,3))
title('T Regulatory Cells')
hold off

subplot(3,4,6)
scatter(CellData.hours, CellData.ThymicDerivedTregsCT)
hold on
plot(tx, ModelData(:,4))
title('Thymic Derived')

subplot(3,4,7)
scatter(CellData.hours, CellData.X4TregProlCT)
hold on
plot(tx, ModelData(:,5))
title('Proliferating Tregs')

subplot(3,4,8)
scatter(CellData.hours, CellData.NaiveDerivedTregsCT)
hold on
plot(tx, ModelData(:,6))
title('Naive Derived')


Changing = {'*alpha',     alpha,     '   cells*hr−1';...
                    'a',           a,             '   hr−1';...
                    'kA',         kA,           '   cells';...
                    '*e_T',       e_T,         '   cells-1*hr−1';...
                    'e_R',       e_R,          '   cells-1*hr−1';...
                    '*g',          g,             '   hr−1';...
                    'b_T',       b_T,           '   hr−1';...
                    '*b_R',      b_R,          '   hr−1';...
                    'j',                j,             '    cell-1 *hour-1';};

columnname =   {'Parameters', '     Value          ', '           Units           '};
columnformat = {'char', 'numeric', 'char'}; 
uitable('Units','normalized',...
                 'Position', [0.12 0.04 0.3466 0.3],... % [ Horizontal Location, Verticle location,Right Line, Bottom Line]
                 'Data', Changing,...
                 'ColumnName', columnname,...
                 'ColumnFormat', columnformat,...
                 'RowName',[],...
                 'FontSize', 15,...
                 'ColumnWidth', {150 200 270});
             
             
Fixed =  {'*mu',       mu,         '   cells*hr−1';...
                '*beta',      beta,      '   hr−1';...   
                '*c',           c,           '   hr−1';...
                'epsilon',  epsilon,   '  hr−1';...
                '*n',           n,           '              -        ';...
                '*d',           d,           '   Molecules*cells-1*hr−1';...
                '*f',            f,           '   hr−1';...
                'kB'          kB           '   cells'};
            
columnname =   {'Parameters', '     Value          ', '           Units           '};
columnformat = {'char', 'numeric', 'char'}; 
uitable('Units','normalized',...
                 'Position', [0.57 0.04 0.394 0.3],... % [ Horizontal Location, Verticle location,Right Line, Bottom Line]
                 'Data', Fixed,...
                 'ColumnName', columnname,...
                 'ColumnFormat', columnformat,...
                 'RowName',[],...
                  'FontSize', 15,...
                 'ColumnWidth', {150 200 360});

PLT2 = figure(2);

%Hill suppression naive
ModelData(:,11) = (1./(1+(ModelData(:,3)./kA).^n));
%Hill suppression Treg death rate
ModelData(:,12) = (1./(1+(ModelData(:,7)./kB).^n));
%How many Activated T's are being destroyed
ModelData(:,13) = (j.*ModelData(:,3).*ModelData(:,2));

subplot(3,1,1)
plot(tx, ModelData(:,11))
title('Hill Value')
ylabel('Hill Value')

subplot(3,1,2)
plot(tx, ModelData(:,12))
title('Treg Death Suppression')
ylabel('Death Rate Suppression')

subplot(3,1,3)
plot(tx, ModelData(:,13))
title('Destroyed T Cells')


%%
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








