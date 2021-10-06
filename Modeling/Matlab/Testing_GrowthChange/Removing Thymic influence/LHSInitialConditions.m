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



%---------------------------------------------------------------------%
%                       Initial Conditions
%---------------------------------------------------------------------%

N = 1027; %Naive T cells
T = 252; %Activated T Cells
R = 128; %T Regulatory Cells (10% of all T cells)

ThyN = 18; %Thymic Derived Naive Cells
ActN = 252; % Activated Naive T Cells
ThyR = 24; % Thymic Derived Tregs
DiffR = 1; %Naive Derived Tregs

Nprol = 1008; %Self replicating naive T cells
Tprol = 0; %Self replicating activated T cells
Rprol = 510; %Self replicating Tregs

I = 0.0001; %IL-2 Cytokine
m = 0.0023; %Average of the Thymus weight at day 0



%%
% Choose a percentage here
% choose witch parameters are going to be changed
% Have the algorithm take care of only the values that you have chosen
% 

%{
    Defining the Initial Conditions
    N - Naive T cells
    T - Activated T Cells
    R - T Regulatory Cells
    
    ThyN - Thymic Derived Naive Cells
    ActN - Activated Naive T Cells
    ThyR - Thymic Derived Tregs
    DiffR - Naive Derived Tregs
    
    Nprol - Self replicating naive T cells
    Tprol - Self replicating activated T cells
    Rprol - Self replicating Tregs
    
    I - IL-2 Cytokine
%}

%Initial Conditions
KeySets = {'N', 'T', 'R', 'ThyN', 'ActN', 'ThyR', 'DiffR', 'Nprol', 'Tprol', 'Rprol', 'I'};
Values = {1027, 252, 128, 18, 252, 24, 1, 1008, 10, 510, 0.0001};
InitCond = containers.Map(KeySets, Values);
tx = 0:432; %Maximum amount of time - 18 days

%-------------------------------------------------------------------------------------------------------%
%                                       Only make changes here                                              % 
CondKeys = {'T', 'ActN', 'Tprol'};
SampleSize = 10;
PctChange = 0.25; %What percentage should the initial conditions vary?

%--------------------------------------------------------------------------------------------------------%


%Generating Samples
par = length(CondKeys); %number of parameters being sampled
LHSVector = lhsdesign(SampleSize,par, 'criterion','maximin','smooth','off');
%%
%Generating Samples
%Defining the sampling Map
LHSConditions = containers.Map('KeyType', 'char', 'ValueType','any');
%Adding all of the randomized values to the map
for sample = 1:length(CondKeys)
    
    CurrentKey = CondKeys(sample); 
    CurrentKey = CurrentKey{1}; %Extracting the string key
    value = values(InitCond, CondKeys(sample)); %Extracting the key's value
    
    % calculate the min and max
    minI = value{1}-r*value{1};
    maxI = value{1}+r*value{1};   
    
    % use the inverse CDF to calculate the initial conditions spread
    LHSCDFinv = unifinv(LHSVector(:,sample),minI,maxI);
    
    % save the new array to the map associated with the right key
    LHSConditions(CurrentKey) = LHSCDFinv;
end 




%%
%Renaming
 S_LHSSamples = X(:,1);
 E_LHSSamples = X(:,2);
 I_LHSSamples = X(:,3);

% using inverse CDF to get the random value.
r = 0.05;
minS = S_0-r*S_0;
maxS = S_0+r*S_0;
S_vec = unifinv(S_LHSSamples,minS,maxS);

minI = I_0-r*I_0;
maxI = I_0+r*I_0;
I_vec = unifinv(I_LHSSamples,minI,maxI);



%Storage Vector for the cellular populations we want to test


LHSOutput = zeros(length(tx), SampleSize); %Each column will hold one iteration




%

%%
values(InitCond, CondKeys)
%%
Keyset = {'N', 'T', 'R'};
values = [1027 252 128];

InitC = containers.Map(Keyset, values);
%%
keyset = {'N', 'T'};
a = values(InitC, keyset);


