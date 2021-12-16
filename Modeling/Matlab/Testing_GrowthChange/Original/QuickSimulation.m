close all; clc
global tx

FileName = 'Data/ParameterSets.csv';
EntryNumber = 27;
p = GetParameters(EntryNumber, FileName);

mu= 0.407;%Thymic Naive
z = p(2); %Prol Naive
g = p(3); %Naive Death
alpha = p(4); %Thymic Tregs
c = p(5); %Naive Derived Tregs
epsilon = p(6); %Treg Prol
b_R = p(7); %Treg Death
beta =p(8); %Activation Rate
a = p(9); %Activated Prol
b_T = p(10); %ActT Death
e_T = p(11); %ActT Consumption
e_R = p(12); %Treg Consumption
kA = p(13); %Beta Suppression
j = 6.2156e-07; %Deactivation
kB = p(15); %Treg Death Suppression
n = p(16);
d = p(17); %IL-2 production Rate
nK = 1835988; %Naive Carrying Capacity
rK = 10459000; %Treg Carrying Capacity
Ki = p(21);%Half rate for activation suppression boost
Kj = p(22);% Half rate for deactivation boost
dKO = p(23);

%{
Keeping this here if I want to replace the top parameter set with normal
values
mu= p(1);%Thymic Naive
z = p(2); %Prol Naive
g = p(3); %Naive Death
alpha = p(4); %Thymic Tregs
c = p(5); %Naive Derived Tregs
epsilon = p(6); %Treg Prol
b_R = p(7); %Treg Death
beta =p(8); %Activation Rate
a = p(9); %Activated Prol
b_T = p(10); %ActT Death
e_T = p(11); %ActT Consumption
e_R = p(12); %Treg Consumption
kA = p(13); %Beta Suppression
j = p(14); %Deactivation
kB = p(15); %Treg Death Suppression
n = p(16);
d = p(17); %IL-2 production Rate
nK = p(19); %Naive Carrying Capacity
rK = p(20); %Treg Carrying Capacity
Ki = p(21);%Half rate for activation suppression boost
Kj = p(22);% Half rate for deactivation boost
dKO = p(23);
%}

%Do not change this order
p0 = [alpha, a, kA, e_T, e_R, g, b_T, b_R, epsilon, mu, beta, c, kB, j, z, n, d, nK, rK, Ki, Kj, dKO];

tx = 0:1000; %Maximum amount of time - 18 days

Genotype = [1, 2];

for i = Genotype
    PlottingResults(p0, i)
end


% PlottingEverything(p0)

%%
global  WTerror KOerror

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


GrowthObjective(p0, parameterStruct);

%%

% Saving parameters to parameterset.csv


error = WTerror + KOerror;
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
    e_T, e_R, kA, j, kB, n, d, f, nK, rK, Ki, Kj, dKO, error, WTerror, KOerror, EntryNumber];


dlmwrite(FileLocation ,parameters,'delimiter', ',',  'precision', 16, '-append');

disp(['Entry Number - ' num2str(EntryNumber)])

%%

% Saving the simulation results to compare to data in R

ModelDataWT = SimulateGrowth(p0, 1);
csvwrite('~/my.work/PhD/HomestaticExpansionProject/Code/Stats plots and data management/ModelOutputWT.csv', ModelDataWT)


ModelDataKO = SimulateGrowth(p0, 2);
csvwrite('~/my.work/PhD/HomestaticExpansionProject/Code/Stats plots and data management/ModelOutputKO.csv', ModelDataKO)


