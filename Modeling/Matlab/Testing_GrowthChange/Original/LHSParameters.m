clc; clear; close all;
%-------------------------------------------------------------------------------------------------------%
%                                       Only make changes here
% Here are the choices {'N', 'T', 'R', 'ThyN', 'ActN', 'ThyR', 'DiffR', 'Nprol', 'Tprol', 'Rprol', 'I'};
CondKeys = {'Kj'};
SampleSize = 100;
PctChange = 1; %What percentage should the initial conditions vary?
EntryNumber = 27;
%--------------------------------------------------------------------------------------------------------%

%Preparing 

FileName = 'Data/ParameterSets.csv';
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
Kj = p(22);% Half rate for deactivation boost
dKO = p(23); %Production rate of IL-2 KO 


%Do not change this order


%Generating map of Initial Conditions
KeySets = {'mu', 'z', 'g', 'alpha', 'c', 'epsilon', 'b_R', 'beta', 'a', 'b_T', 'e_T', 'e_R', 'kA',...
    'j', 'kB', 'n', 'd', 'nK', 'rK', 'Ki', 'Kj', 'dKO'};
Values = {mu, z, g, alpha, c, epsilon, b_R, beta, a, b_T, e_T, e_R, kA, ...
    j, kB, n, d, nK, rK, Ki, Kj, dKO};
% Values = {p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9), p(10), p(11), p(12), p(13),...
%     p(14), p(15), p(16), p(17), p(19), p(20), p(21), p(22), p(23)};
ParamValues = containers.Map(KeySets, Values);
tx = 0:1000; %Maximum amount of time - 18 days

%Generating LHS Samples
par = length(CondKeys); %number of parameters being sampled
LHSVector = lhsdesign(SampleSize,par, 'criterion','maximin','smooth','off');

%Defining the sampling Map
LHSParaSamples = containers.Map('KeyType', 'char', 'ValueType','any');

%Generating Samples
for sample = 1:par
    
    CurrentKey = CondKeys(sample); 
    CurrentKey = CurrentKey{1}; %Extracting the string key
    value = values(ParamValues, CondKeys(sample)); %Extracting the key's value
    
    % calculate the min and max
    minI = value{1}-PctChange*value{1};
    maxI = value{1}+PctChange*value{1};   
    
    % use the inverse CDF to calculate the initial conditions spread
    LHSCDFinv = unifinv(LHSVector(:,sample),minI,maxI);
    
    % save the new array to the map associated with the right key
    LHSParaSamples(CurrentKey) = LHSCDFinv;
end 

% Generating duplicate values for parameters that have not been
% through the LHS sampling. This is to make solving of future iterations
% easier

remove(ParamValues , CondKeys); %removes keys that were used for LHS

for static = 1:length(keys(ParamValues))
    ParamKeys = keys(ParamValues);
    
    %Grab keys
    keystring = ParamKeys(static);
    keystring = keystring{1}; %Grab just the key string
    value = ParamValues (keystring);  %Grabbing the keys value

    a = repelem(value, SampleSize); %duplicating the values
    LHSParaSamples(keystring) = a';
    
end


%-----------------------------------------------------------------------------------------------%
%                                     Preparing for Simulations
%-----------------------------------------------------------------------------------------------%


% Storage Vectors for all the cellular growths per sample size. Each column
% in these will represent one iteration of the sameple size
NaiveCTWT = zeros(length(tx),SampleSize);
ActTCTWT = zeros(length(tx),SampleSize);
TregCTWT = zeros(length(tx),SampleSize);
ThyNWT = zeros(length(tx),SampleSize);
ActNWT = zeros(length(tx),SampleSize);
ThyRWT = zeros(length(tx),SampleSize);
DiffRWT = zeros(length(tx),SampleSize);
NprolWT = zeros(length(tx),SampleSize);
TprolWT = zeros(length(tx),SampleSize);
RprolWT = zeros(length(tx),SampleSize);
IWT = zeros(length(tx),SampleSize);

%IL-2 KO results
NaiveCTKO = zeros(length(tx),SampleSize);
ActTCTKO = zeros(length(tx),SampleSize);
TregCTKO = zeros(length(tx),SampleSize);
ThyNKO = zeros(length(tx),SampleSize);
ActNKO = zeros(length(tx),SampleSize);
ThyRKO = zeros(length(tx),SampleSize);
DiffRKO = zeros(length(tx),SampleSize);
NprolKO = zeros(length(tx),SampleSize);
TprolKO = zeros(length(tx),SampleSize);
RprolKO = zeros(length(tx),SampleSize);
IKO = zeros(length(tx),SampleSize);

Genotype = [1, 2]; %Will run simulation for both genotypes

%Prepping the Parameter values for interation
mu_LHS = LHSParaSamples('mu'); %Thymic Naive
z_LHS = LHSParaSamples('z'); %Prol Naive
g_LHS = LHSParaSamples('g'); %Naive Death
alpha_LHS = LHSParaSamples('alpha'); %Thymic Tregs
c_LHS = LHSParaSamples('c'); %Naive Derived Tregs
epsilon_LHS = LHSParaSamples('epsilon'); %Treg Prol
b_R_LHS = LHSParaSamples('b_R'); %Treg Death
beta_LHS = LHSParaSamples('beta'); %Activation Rate
a_LHS = LHSParaSamples('a'); %Activated Prol
b_T_LHS = LHSParaSamples('b_T'); %ActT Death
e_T_LHS = LHSParaSamples('e_T'); %ActT Consumption
e_R_LHS = LHSParaSamples('e_R'); %Treg Consumption
kA_LHS = LHSParaSamples('kA'); %Beta Suppression
j_LHS = LHSParaSamples('j'); %Deactivation
kB_LHS = LHSParaSamples('kB'); %Treg Death Suppression
n_LHS = LHSParaSamples('n');
d_LHS = LHSParaSamples('d'); %IL-2 production Rate
nK_LHS = LHSParaSamples('nK'); %Naive Carrying Capacity
rK_LHS = LHSParaSamples('rK'); %Treg Carrying Capacity
Ki_LHS = LHSParaSamples('Ki');%Half rate for activation suppression boost
Kj_LHS = LHSParaSamples('Kj');% Half rate for deactivation boost
dKO_LHS = LHSParaSamples('dKO'); %Production rate of IL-2 KO 

IterationNumber = 0;

%Running Simulation
for iter = 1:SampleSize
    % Setting up the Init Conditions
    N = 1027;
    T = 252;
    R = 128;
    ThyN = 18;
    ActN = 252;
    ThyR = 24;
    DiffR = 1;
    Nprol = 1008;
    Tprol = 0;
    Rprol = 510;
    I = 0.0001;
    m = 0.0023; %Too lazy to remove this from every where rn
        
    T0 = [N T R ...
    ThyN ActN ThyR DiffR ... 
    Nprol Tprol Rprol ...
    I m ];
    
    %Setting up Parameters
    mu = mu_LHS(iter);
    z = z_LHS(iter);
    g = g_LHS(iter);
    alpha = alpha_LHS(iter);
    c = c_LHS(iter);
    epsilon = epsilon_LHS(iter);
    b_R = b_R_LHS(iter);
    beta = beta_LHS(iter);
    a = a_LHS(iter);
    b_T = b_T_LHS(iter);
    e_T = e_T_LHS(iter);
    e_R = e_R_LHS(iter);
    kA = kA_LHS(iter);
    j = j_LHS(iter);
    kB = kB_LHS(iter);
    n = n_LHS(iter);
    d = d_LHS(iter);
    nK = nK_LHS(iter);
    rK = rK_LHS(iter);
    Ki = Ki_LHS(iter);
    Kj = Kj_LHS(iter);
    dKO = dKO_LHS(iter);
    
    
    p0 = [alpha, a, kA, e_T, e_R, g, b_T, b_R, epsilon, mu, beta, c, kB, j, z, n, d, nK, rK, Ki, Kj, dKO];
    
    IterationNumber = IterationNumber + 1;
    disp(['Iteration Number: ', num2str(IterationNumber)])
    for Gene = Genotype
        
        %------------------This is where the solving happens----------------%
        ModelData = SimulateGrowthLHS(T0, tx,p0, Gene );
        %------------------------------------------------------------------------------%
        
        if Gene == 1        
            NaiveCTWT(:,iter) = ModelData(:,1);
            ActTCTWT(:,iter) = ModelData(:,2);
            TregCTWT(:,iter)= ModelData(:,3);
            %Non Proliferating Trackers
            ThyNWT(:,iter) = ModelData(:,4);
            ActNWT(:,iter) = ModelData(:,5);
            ThyRWT(:,iter) = ModelData(:,6);
            DiffRWT(:,iter) = ModelData(:,7);
            % Proliferating Cell Trackers
            NprolWT(:,iter) = ModelData(:,8);
            TprolWT(:,iter) = ModelData(:,9);
            RprolWT(:,iter) = ModelData(:,10);
            % Other
            IWT(:,iter) = ModelData(:,11);
            
        elseif Gene ==2
            NaiveCTKO(:,iter) = ModelData(:,1);
            ActTCTKO(:,iter) = ModelData(:,2);
            TregCTKO(:,iter)= ModelData(:,3);
            %Non Proliferating Trackers
            ThyNKO(:,iter) = ModelData(:,4);
            ActNKO(:,iter) = ModelData(:,5);
            ThyRKO(:,iter) = ModelData(:,6);
            DiffRKO(:,iter) = ModelData(:,7);
            % Proliferating Cell Trackers
            NprolKO(:,iter) = ModelData(:,8);
            TprolKO(:,iter) = ModelData(:,9);
            RprolKO(:,iter) = ModelData(:,10);
            % Other
            IKO(:,iter) = ModelData(:,11);
        end
        
    end
    
end

%-----------------------------------------------------------------------------------------------%
%                               Preparing Data for Plotting
%-----------------------------------------------------------------------------------------------%

%Saving WT Data
CellularData(:,:,1,1) = NaiveCTWT;
CellularData(:,:,2,1) = ActTCTWT;
CellularData(:,:,3,1) = TregCTWT;
CellularData(:,:,4,1) = ThyNWT;
CellularData(:,:,5,1) = ActNWT;
CellularData(:,:,6,1) = ThyRWT;
CellularData(:,:,7,1) = DiffRWT;
CellularData(:,:,8,1) = NprolWT;
CellularData(:,:,9,1) = TprolWT;
CellularData(:,:,10,1) = RprolWT;
CellularData(:,:,11,1) = IWT;
%Saving KO Data
CellularData(:,:,1,2) = NaiveCTKO;
CellularData(:,:,2,2) = ActTCTKO;
CellularData(:,:,3,2) = TregCTKO;
CellularData(:,:,4,2) = ThyNKO;
CellularData(:,:,5,2) = ActNKO;
CellularData(:,:,6,2) = ThyRKO;
CellularData(:,:,7,2) = DiffRKO;
CellularData(:,:,8,2) = NprolKO;
CellularData(:,:,9,2) = TprolKO;
CellularData(:,:,10,2) = RprolKO;
CellularData(:,:,11,2) = IKO;

%------------ All Statistical calculations are done here------------------%
StatsOfCells = CalculateTheFillRanges(CellularData);

%----------------------Plotting the LHS Results-----------------------------%
PlottingLHSResults(StatsOfCells, tx)



