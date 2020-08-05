function dTdt = Growth(t, x, p)
%Getting all parameter values
%Prmt = GetParameters
Prmt = p;

%Order of Parameters
% 1_alpha 2_Thy 3_Thy_max 4_epsilon 5_a 6_c 7_b_R 8_mu 9_beta 10_g 11_b_T 
% 12_d 13_e_T 14_e_R 15_f 16_kA 17_n 18_EntryNumber 19_Notes 20_Naive
% 21_Activated 22_Treg 23_IL2 24_PreviousPset
%disp(['alpha' num2str(Prmt(1))])
%disp(['Thy' num2str(Prmt(2))])
%disp(['Thy_max' num2str(Prmt(3))])
%disp(['epsilon' num2str(Prmt(4))])
%disp(['a' num2str(Prmt(5))])

%Change parameters that will be fitted accordingly using the p parameter
 alpha = Prmt(1);
 Thy = Prmt(2);
 Thy_max = Prmt(3);
 epsilon = Prmt(4);
 a = Prmt(5);
 c = Prmt(6);
 b_R = Prmt(7);
 mu = Prmt(8);
 beta = Prmt(9);
 g = Prmt(10);
 b_T = Prmt(11);
 d = Prmt(12);
 e_T = Prmt(13);
 e_R = Prmt(14);
 f = Prmt(15);
 kA = Prmt(16);
 n = Prmt(17);

%v = [N T R I]
%dd = [dN/dt; dT/dt; dR/dt; dI/dt]
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



