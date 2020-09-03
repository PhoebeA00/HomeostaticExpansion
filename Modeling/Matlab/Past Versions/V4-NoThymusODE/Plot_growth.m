function dTdt = Plot_growth(t, x, p, i)

%Order of Parameters
% 1_alpha 2_Thy 3_Thy_max 4_epsilon 5_a 6_c 7_b_R 8_mu 9_beta 10_g 11_b_T 
% 12_d 13_e_T 14_e_R 15_f 16_kA 17_n 18_EntryNumber 19_Notes 20_Naive
% 21_Activated 22_Treg 23_IL2 24_PreviousPset

%Thymus Equation

%Values for the Thymus polynomial
p1 = -3.725e-09;
p2 = 2.419e-06;
p3 = -0.0001933;
p4 = 0.00398;

ThyEq = @(Mdlt) p1*Mdlt.^3 + p2*Mdlt.^2 + p3*Mdlt + p4;



%Change parameters that will be fitted accordingly using the p parameter
alpha = p(1);
Thy = ThyEq(i); %Current Thymus value
Thy_max = ThyEq(432); %The maximum value
epsilon = p(4);
a = p(5);
c = p(6);
b_R = p(7);
mu = p(8);
beta = p(9);
g = p(10);
b_T = p(11);
d = p(12);
e_T = p(13);
e_R = p(14);
f = p(15);
kA = p(16);
n = p(17);
 


%Initial conditions
N = x(1);
T = x(2);
R = x(3);
I = x(4);

%Cell  
dNdt = mu*(Thy/Thy_max)-beta*N*(1/(1+(R/kA)^n))- c*N - g*N;
    
dTdt = beta*N*(1/(1+(R/kA)^n)) + a*I*T - b_T*T;

dRdt = alpha*(Thy/Thy_max) + epsilon*a*I*R + c*N - b_R*R;

dIdt = d*T - e_T*I*T - e_R*I*R - f*I;
%The apostrophe returns a transposed matrix of [4x1]
dTdt = [dNdt, dTdt, dRdt, dIdt]';



