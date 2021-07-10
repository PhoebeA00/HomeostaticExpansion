close all; clc
%global tx N T R I ThyR RplR DiffR pOpt p0 m
global N T R ThyN ActN ThyR DiffR Nprol Tprol Rprol I m tx
global WhichGenotype

%1 = Wild Type, 2 = Knockout
WhichGenotype = 1;

% Parameter Ranges

%%%%%%%% Naive %%%%%%%%%%%%
mu_min =   339560; %Naive T production rate
mu_max = 339560; 

z_min = 0.0590; %Naive Self replication rate
z_max = 0.0590;

g_min = 0.03;%Death rate of Naive
g_max = 0.03;

%%%%%%%% Tregs %%%%%%%%%%%%
alpha_min = 3069.6;%33914; %Thymic derived Tregs
alpha_max = 3069.6; %3000;%33914; 

c_min = 0.0426; %Naive differentiation to Tregs
c_max = 0.0426;%0.0292; %0.05;

epsilon_min = 0.25;%Self Replication rate of Tregs
epsilon_max = 0.25;%0.2;

b_R_min = 0.7986;%Death Rate of Tregs
b_R_max = 0.7986;

%%%%%%%% Activated T %%%%%%%%%%%%
beta_min = 0.5122; %activation rate
beta_max = 0.5122;

a_min = 0.08; %Self Replication rate for activated T cells
a_max = 0.08;

b_T_min = 0.2698;%Death Rate of T cells
b_T_max = 0.2698;

%%%%%%%% Consumption Rates %%%%%%%%%%%%
e_T_min = 100; %Consumption rate of T cells
e_T_max = 100; %activated T cells consume at max 1/10 of what they make

e_R_min = 393;%1.3744E-09; %Consumption rate of Tregs
e_R_max = 393;

%%%%%%%% Suppression %%%%%%%%%%%%
kA_min = 448710; %Half suppression rate by Tregs
kA_max = 448710; 

j_min = 2.9788e-06; %Rate of desctruction of activated T cells
j_max = 2.9788e-06; 

kB_min = 3.9328; %half suppression rate of Treg death rate
kB_max = 3.9328;

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
kB = kB_min + rand(1,1)*(kB_max - kB_min);
j = j_min + rand(1,1)*(j_max - j_min);
z = z_min + rand(1,1)*(z_max - z_min);
n = n_min + rand(1,1)*(n_max - n_min);
d = d_min + rand(1,1)*(d_max - d_min);

%Making sure that the consumption of Tregs is greater than that of
%activated T cell
e_T = e_T_min + rand(1,1) * (e_T_max - e_T_min);
e_R = e_R_min + rand(1,1) * (e_R_max - e_R_min);

while e_T > e_R
    e_T = e_T_min + rand(1,1) * (e_T_max - e_T_min);
    e_R = e_R_min + rand(1,1) * (e_R_max - e_R_min);
end

%------fmincon function arguments definitions-------
p0 = [alpha, a, kA, e_T, e_R, g, b_T, b_R, epsilon, mu, beta, c, kB, j, z, n, d];
lb = [alpha_min, a_min, kA_min, e_T_min, e_R_min, g_min, ...
    b_T_min, b_R_min, epsilon_min, mu_min, beta_min, c_min, kB_min, j_min, z_min, n_min, d_min]; %[] lower bound
ub = [alpha_max, a_max, kA_max, e_T_max, e_R_max, g_max, ...
    b_T_max, b_R_max, epsilon_max, mu_max, beta_max, c_max, kB_max, j_max, z_max, n_max, d_max]; %[] upper bound

% no linear constraints
A = [0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0];
b = 0;
Aeq = [];
beq = [];
nlcon = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------Initial Conditions-----%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if WhichGenotype == 1
    
    N = 1027; %Naive T cells
    T = 252; %Activated T Cells
    R = 127; %T Regulatory Cells
    
    ThyN = 18; %Thymic Derived Naive Cells
    ActN = 252; % Activated Naive T Cells
    ThyR = 1418; % Thymic Derived Tregs
    DiffR = 1; %Naive Derived Tregs
    
    Nprol = 1008; %Self replicating naive T cells
    Tprol = 0; %Self replicating activated T cells
    Rprol = 127; %Self replicating Tregs
    
    I = 0; %IL-2 Cytokine
    m = 0.0023; %Average of the Thymus weight at day 0
elseif WhichGenotype == 2
    
    N = 1811; %Naive T cells
    T = 173; %Activated T Cells
    R = 30576; %T Regulatory Cells
    
    ThyN = 459; %Thymic Derived Naive Cells
    ActN = 141; % Activated Naive T Cells
    ThyR = 947; % Thymic Derived Tregs
    DiffR = 1; %Naive Derived Tregs
    
    Nprol = 1352; %Self replicating naive T cells
    Tprol = 32; %Self replicating activated T cells
    Rprol = 30576; %Self replicating Tregs
    
    I = 0; %IL-2 Cytokine
    m = 0.0023; %Average of the Thymus weight at day 0
end


tx = 0:432; %Maximum amount of time - 18 days

% optimize parameters
disp('Beginning Optimization...')
% options = optimoptions(@fmincon,'Algorithm','interior-point');
[pOpt, error] = fmincon(@GrowthObjective,p0,A,b,Aeq,beq,lb,ub,nlcon);
disp('...Ending Optimization')
%Optimized Parameters

PlottingResults(pOpt)
%{
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
z = pOpt(15);
n = pOpt(16);


ModelData = SimulateGrowth(pOpt);
ModelData = array2table(ModelData);

ModelData.Properties.VariableNames = {'NaiveCT' 'ActivCT' 'TregCT' ...
    'ThyDerivedNaive' 'ActivNaive' 'ThyDerivedTregs' 'NaiveDerivedTregs' ...
    'ProlNaive' 'ProlActiv' 'ProlTregs' ...
    'Il2' 'ThyWeight'};

%1 = WildType, 2 = Genotype
if WhichGenotype == 1
    CellData = readtable('../RawData/ActivatedWTSpleen.csv');
    ProlData = readtable('../RawData/WTProl.csv');
elseif WhichGenotype == 2
    CellData = readtable('../RawData/ActivatedKOSpleen.csv');
    ProlData = readtable('../RawData/KOProl.csv');
end

CellData = CellData(:,{'NaiveCT', 'ActivatedCD4CT', 'X4TregCT', ...
    'ThymicNaive', 'ActivatedNaiveCT', ...
    'ThymicDerivedTregsCT', 'NaiveDerivedTregsCT' ... 
    'hours'});

ProlData = ProlData(:,{ 'NaiveProlCT', 'ActivatedProlCT', 'X4TregProlCT', ...
    'hours'});

PLT = figure(1);
%Naive CT
subplot(3,4,1)
scatter(CellData.hours, CellData.NaiveCT)
hold on 
plot(tx, ModelData.NaiveCT)
title('Naive T Cells')
hold off
%Naive Prol
subplot(3,4,2)
scatter(ProlData.hours, ProlData.NaiveProlCT)
hold on 
plot(tx, ModelData.ProlNaive)
title('Proliferating Naive')
hold off
%Naive Thymic
subplot(3,4,3)
scatter(CellData.hours, CellData.ThymicNaive)
hold on 
plot(tx, ModelData.ThyDerivedNaive)
title('Thymic Naive')
hold off
%Activated CT
subplot(3,4,5)
scatter(CellData.hours, CellData.ActivatedCD4CT)
hold on 
plot(tx, ModelData.ActivCT)
title('Activated Counts')
hold off
%Activated Prol
subplot(3,4,6)
scatter(ProlData.hours, ProlData.ActivatedProlCT)
hold on 
plot(tx, ModelData.ProlActiv)
title('Activated Prol')
hold off
%Activated Naive derived
subplot(3,4,7)
scatter(CellData.hours, CellData.ActivatedNaiveCT)
hold on 
plot(tx, ModelData.ActivCT)
title('Activated Naive CT')
hold off
%IL-2
subplot(3,4,8)
plot(tx, ModelData.Il2)
title('IL-2')
hold off
%Treg CT
subplot(3,4,9)
scatter(CellData.hours, CellData.X4TregCT)
hold on 
plot(tx, ModelData.TregCT)
title('Treg CT')
hold off
%Treg Prol
subplot(3,4,10)
scatter(ProlData.hours, ProlData.X4TregProlCT)
hold on 
plot(tx, ModelData.ProlTregs)
title('Treg Prol')
hold off
%Thymic Tregs
subplot(3,4,11)
scatter(CellData.hours, CellData.ThymicDerivedTregsCT)
hold on 
plot(tx, ModelData.ThyDerivedTregs)
title('Thymic Tregs')
hold off
%Naive Derived Tregs
subplot(3,4,12)
scatter(CellData.hours, CellData.NaiveDerivedTregsCT)
hold on 
plot(tx, ModelData.NaiveDerivedTregs)
title('Naive Derived Tregs')
hold off

%Hill suppression naive
ModelData.HillNaive = (1./(1+(ModelData.NaiveCT./kA).^n));
%Hill suppression Treg death rate
ModelData.HillTregDeath = (1./(1+(ModelData.TregCT./kB).^n));
%How many Activated T's are being destroyed
ModelData.ActiveDestruction = (j.*ModelData.TregCT.*ModelData.ActivCT);


PLT2 = figure(2);

%Hill suppression naive
subplot(3,3,1)
plot(tx, ModelData.HillNaive)
title('Hill Value')
ylabel('Hill Value')

subplot(3,3,2)
plot(tx, ModelData.HillTregDeath)
title('Treg Death Suppression')
ylabel('Death Rate Suppression')

subplot(3,3,3)
plot(tx, ModelData.ActiveDestruction)
title('Destroyed T Cells')


Changing = {'*mu',       mu,         '   cells*hr−1';...
                    'z',             z,             'cells-1*hr-1';...
                    '*g',          g,             '   hr−1';...
                    
                    '*alpha',     alpha,     '   cells*hr−1';...
                    'c',           c,           '   hr−1';...
                    '*epsilon',  epsilon,   '  hr−1';...
                    '*b_R',      b_R,          '   hr−1';...
                    
                    'beta',      beta,      '   hr−1';...
                    '*a',           a,             '   hr−1';...
                    '*b_T',       b_T,           '   hr−1'};

columnname =   {'Parameters', '     Value          ', '           Units           '};
columnformat = {'char', 'numeric', 'char'}; 
uitable('Units','normalized',...
                 'Position', [0.12 0.2 0.3466 0.4],... % [ Horizontal Location, Verticle location,Right Line, Bottom Line]
                 'Data', Changing,...
                 'ColumnName', columnname,...
                 'ColumnFormat', columnformat,...
                 'RowName',[],...
                 'FontSize', 15,...
                 'ColumnWidth', {150 200 270});
             
             
Fixed =  {'*e_T',       e_T,         '   cells-1*hr−1';...
                '*e_R',       e_R,          '   cells-1*hr−1';...
                
                'kA',         kA,           '   cells';...
                'j',                j,             '    cell-2*hour-1';...
                '*kB',          kB,           '   cells';...
                
                'n',           n,           '              -        ';...
                '*d',           d,           '   Molecules*cells-1*hr−1';...
                '*f',            f,           '   hr−1'};
            
columnname =   {'Parameters', '     Value          ', '           Units           '};
columnformat = {'char', 'numeric', 'char'}; 
uitable('Units','normalized',...
                 'Position', [0.57 0.2 0.394 0.4],... % [ Horizontal Location, Verticle location, Right Line, Bottom Line]
                 'Data', Fixed,...
                 'ColumnName', columnname,...
                 'ColumnFormat', columnformat,...
                 'RowName',[],...
                  'FontSize', 15,...
                 'ColumnWidth', {150 200 360});
%}

%%

alpha = pOpt(1);
a = pOpt(2);
kA = pOpt(3);
e_T = pOpt(4);
e_R = pOpt(5);
g = pOpt(6);
b_T = pOpt(7);
b_R = pOpt(8);
epsilon = pOpt(9);
mu = pOpt(10);
beta = pOpt(11);
c = pOpt(12);
kB = pOpt(13);
j = pOpt(14);
z = pOpt(15);
n = pOpt(16);
d = pOpt(17);


%-----Change this for saving files in a different location-----%
FileLocation = '../Data/ParameterRanges28_WT.csv';
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
    e_T, e_R, kA, j, kB, n, d, f, error, EntryNumber];


dlmwrite(FileLocation ,parameters,'delimiter', ',', '-append');

disp(['Entry Number - ' num2str(EntryNumber)])












