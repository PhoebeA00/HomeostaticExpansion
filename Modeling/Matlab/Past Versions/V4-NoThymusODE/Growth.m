function dTdt = Growth(t, x, p, i)
global alpha Thy Thy_max epsilon a c b_R mu beta g b_T d e_T e_R f kA n
global ThyEq
%Getting all parameter values

%%%%%%%% Make all changes here for optimization %%%%%%%%%%%%%%%%%
%                                                               %
%Prmt = GetParameters(10); %Adjust parameter choice here %%%%%%%%%
alpha = p(1);
mu = p(2);                                                      %
%
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Order of Parameters
% 1_alpha 2_Thy 3_Thy_max 4_epsilon 5_a 6_c 7_b_R 8_mu 9_beta 10_g 11_b_T 
% 12_d 13_e_T 14_e_R 15_f 16_kA 17_n 18_EntryNumber 19_Notes 20_Naive
% 21_Activated 22_Treg 23_IL2 24_PreviousPset

 
%Change parameters that will be fitted accordingly using the p parameter
 %alpha = 0; %Prmt(1);
 %Thy = 1;%Prmt(2);
 %Thy_max = 1;%Prmt(3);
 %epsilon = 0;%Prmt(4);
 %a = 0.09; %Prmt(5);
 %c = 0;%Prmt(6);
 %b_R = 0; %Prmt(7);
 %mu = Prmt(8);
 %beta = Prmt(9);
 %g = 0.01;%Prmt(10);
 %b_T = 20;%Prmt(11);
 %d = 500;%Prmt(12);
 %e_T = 300;%Prmt(13);
 %e_R = 800;%Prmt(14);
 %f = 0.1;%Prmt(15);
 %kA = 1;%Prmt(16);
 %n = 0;%Prmt(17);

%Initial conditions
N = x(1);
T = x(2);
R = x(3);
I = x(4);


Thy = ThyEq(i);
%disp(['Thymus: ' num2str(Thy)])

%Cell  
dNdt = mu*(Thy/Thy_max)-beta*N*(1/(1+(R/kA)^n))- c*N - g*N;
    
dTdt = beta*N*(1/(1+(R/kA)^n)) + a*I*T - b_T*T;

dRdt = alpha*(Thy/Thy_max) + epsilon*a*I*R + c*N - b_R*R;

dIdt = d*T - e_T*I*T - e_R*I*R - f*I;
%The apostrophe returns a transposed matrix of [4x1]
dTdt = [dNdt, dTdt, dRdt, dIdt]';



