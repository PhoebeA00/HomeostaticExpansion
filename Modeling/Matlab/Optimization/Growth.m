function dTdt = Growth(t, x, p, i)
%global alpha epsilon a c b_R mu beta g b_T d e_T e_R f kA n


alpha   = p(1);
a       = p(2);
kA      = p(3);
e_T     = p(4);
e_R     = p(5);
g        = p(6);
b_T     = p(7);
b_R     = p(8);
epsilon = p(9);
mu = p(10);
beta = p(11);
c = p(12);
kB = p(13);
j = p(14);
z = p(15);
n = p(16);
d = p(17);
n1 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------Fixed Parameters-------%
% These should never be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%d = 1000; %IL-2 production Rate
f = 1.38629; %IL-2 degradation Rate

%Initial conditions
N = x(1); %Naive T Cells
T = x(2); %Activated T Cells
R = x(3); %Regulatory T Cells

ThyN = x(4); %Thymic Dervied Naive
ActN =  x(5);  % Activated Naive T Cells
ThyR = x(6); % Thymic Derived Tregs
DiffR = x(7); % Naive Derived Tregs

Nprol = x(8); % Self replicating naive T cells
Tprol = x(9); % Self replicating activated T cells
Rprol = x(10); % Self replicating Tregs

I = x(11); %IL-2 Cytokine
m = x(12); %Mass of Thymus


%Thymus Parameters
K = 0.074896;
lambda = 0.016932;

% Naive and sub populations  
dNdt = mu*(m/K)-beta*N*(1/(1+(R/kA)^n)) + z*N - c*N - g*N;
dNproldt = z*N - beta*Nprol*(1/(1+(R/kA)^n)) - c*Nprol - g*Nprol ; %Self replicating naive T cells
dThyNdt = mu*(m/K) - beta*ThyN*(1/(1+(R/kA)^n)) - c*ThyN- g*ThyN; %Thymic derived naive T cells

% Activated T Cells and sub populations
dTdt = beta*N*(1/(1+(R/kA)^n)) + a*T - j*R*T - b_T*T;
dActNdt = beta*N*(1/(1+(R/kA)^n)) - b_T*ActN; %Activated Naive T cells
dTproldt = a*T - b_T*Tprol; %Self replicating activated T cells

% Treg and Treg sub-populations
dRdt = alpha*(m/K) + epsilon*R + c*N - b_R*R*(1/(1+(I/kB)^n1));
dThyRdt = alpha*(m/K) - b_R*ThyR*(1/(1+(I/kB)^n1)); %Thymic Dervied Tregs
dRproldt = epsilon*R - b_R*Rprol*(1/(1+(I/kB)^n1)); %Self Replicating Tregs
dDiffRdt = c*N - b_R*DiffR*(1/(1+(I/kB)^n1)); %Naive Derived Tregs

% IL-2
dIdt = d*T - e_T*I*T - e_R*I*R - f*I;
% Thymic Mass
dmdt = lambda * m * (1 - (m/K));

%{
disp('Before Returning the Data')
disp('Could it be the data?')
disp('naive')
disp(dNdt)
disp('Activated')
disp(dTdt)
disp('AllTregs')
disp(dRdt)
disp('thyR')
disp(dThyRdt)
disp('RplR')
disp(dRplRdt)
disp('dDiffRdt')
disp(dDiffRdt)
disp('Thymic mass')
disp(dmdt)
disp('IL-2')
disp(dIdt)
%}

%The apostrophe returns a transposed matrix of [4x1]
%dTdt = [dNdt, dTdt, dRdt, dThyRdt, dRplRdt, dDiffRdt, dIdt, dmdt]';
dTdt = [dNdt, dTdt, dRdt, ... %1, 2, 3,
    dThyNdt, dActNdt, dThyRdt, dDiffRdt, ... %4, 5, 6, 7
    dNproldt, dTproldt, dRproldt, ... %8, 9, 10
    dIdt, dmdt]'; %11, 12



