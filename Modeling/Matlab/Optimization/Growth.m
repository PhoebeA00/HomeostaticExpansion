function dTdt = Growth(t, x, p, i)
%global alpha epsilon a c b_R mu beta g b_T d e_T e_R f kA n
global d f n

alpha   = p(1);
a       = p(2);
kA      = p(3);
e_T     = p(4);
e_R     = p(5);
g       = p(6);
b_T     = p(7);
b_R     = p(8);
epsilon = p(9);
mu = p(10);
beta = p(11);
c = p(12);
kB = p(13);



%Order of Parameters
% 1_alpha 2_Thy 3_Thy_max 4_epsilon 5_a 6_c 7_b_R 8_mu 9_beta 10_g 11_b_T 
% 12_d 13_e_T 14_e_R 15_f 16_kA 17_n 18_EntryNumber 19_Notes 20_Naive
% 21_Activated 22_Treg 23_IL2 24_PreviousPset

%Initial conditions
N = x(1);
T = x(2);
R = x(3);
I = x(4);
m = x(5); %Mass of Thymus

%Thymus Parameters
K = 0.074896;
lambda = 0.016932;

%Cell  
dNdt = mu*(m/K)-beta*N*(1/(1+(R/kA)^n)) - c*N - g*N;
    
dTdt = beta*N*(1/(1+(R/kA)^n)) + a*T - b_T*T;

dRdt = alpha*(m/K) + epsilon*R + c*N - b_R*R*(1/(1+(I/kB)^n));

dIdt = d*T - e_T*I*T - e_R*I*R - f*I;

dmdt = lambda * m * (1 - (m/K));

%The apostrophe returns a transposed matrix of [4x1]
dTdt = [dNdt, dTdt, dRdt, dIdt, dmdt]';



