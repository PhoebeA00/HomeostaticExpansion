close all; clc
global tx

FileName = 'Data/ParameterSets.csv';
EntryNumber = 21;
p = GetParameters(EntryNumber, FileName);

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
Kj = 0.001;% Half rate for deactivation boost

% mu= 0.175;%Thymic Naive
% z = 0.019; %Prol Naive
% g = p(3); %Naive Death
% alpha = p(4); %Thymic Tregs
% c = p(5); %Naive Derived Tregs
% epsilon = p(6); %Treg Prol
% b_R = p(7); %Treg Death
% beta =0.102; %Activation Rate
% a = p(9); %Activated Prol
% b_T = p(10); %ActT Death
% e_T = p(11); %ActT Consumption
% e_R = p(12); %Treg Consumption
% kA = 314120; %Beta Suppression
% j = 3.975834604883678e-07; %Deactivation
% kB = p(15); %Treg Death Suppression
% n = p(16);
% d = p(17); %IL-2 production Rate
% nK = p(19); %Naive Carrying Capacity
% rK = p(20); %Treg Carrying Capacity
% Ki = 3;%Half rate for activation suppression boost
% Kj = p(22);% Half rate for deactivation boost


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
%}

%Do not change this order
p0 = [alpha, a, kA, e_T, e_R, g, b_T, b_R, epsilon, mu, beta, c, kB, j, z, n, d, nK, rK, Ki, Kj];

tx = 0:432; %Maximum amount of time - 18 days

Genotype = [1, 2];

for i = Genotype
    PlottingResults(p0, i)
end


% PlottingEverything(p0)

%%

GrowthObjective(p0);
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
    e_T, e_R, kA, j, kB, n, d, f, nK, rK, Ki, Kj, error, WTerror, KOerror, EntryNumber];


dlmwrite(FileLocation ,parameters,'delimiter', ',', '-append');

disp(['Entry Number - ' num2str(EntryNumber)])
