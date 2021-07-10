global N T R ThyN ActN ThyR DiffR Nprol Tprol Rprol I m tx
global WhichGenotype

FileName = '../Data/ParameterRanges28_WT.csv';
EntryNumber = 11;
p = GetParameters(EntryNumber, FileName);

WhichGenotype = 1;

mu= p(1);
z = p(2);
g = p(3);
alpha = p(4);
c = p(5);
epsilon = p(6); 
b_R = p(7);
beta =p(8);
a = p(9);
b_T = p(10);
e_T = p(11);
e_R = p(12);
kA = p(13);
j = p(14);
kB = p(15);
n = p(16);
%d = 0;
d = p(17);

%{
mu= p(1);
z = p(2);
g = p(3);
alpha = p(4);
c = p(5);
epsilon = p(6); 
b_R = p(7);
beta =p(8);
a = p(9);
b_T = p(10);
e_T = p(11);
e_R = p(12);
kA = p(13);
j = p(14);
kB = p(15);
n = p(16);
d = p(17);
%}

p0 = [alpha, a, kA, e_T, e_R, g, b_T, b_R, epsilon, mu, beta, c, kB, j, z, n, d];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------Initial Conditions-----%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if WhichGenotype == 1
    
    N = 1027; %Naive T cells
    T = 252; %Activated T Cells
    %R = 29792; %T Regulatory Cells
    R = 0; %T Regulatory Cells
    
    ThyN = 18; %Thymic Derived Naive Cells
    ActN = 252; % Activated Naive T Cells
    ThyR = 1418; % Thymic Derived Tregs
    DiffR = 1; %Naive Derived Tregs
    
    Nprol = 1008; %Self replicating naive T cells
    Tprol = 0; %Self replicating activated T cells
    Rprol = 29732; %Self replicating Tregs
    
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

PlottingResults(p0)


