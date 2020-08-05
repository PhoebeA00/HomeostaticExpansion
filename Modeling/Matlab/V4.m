%This is an example of a function.
function dd=V4(t,v)
global alpha  Thy  Thy_max epsilon a c b_R mu beta g b_T d e_T e_R f kA n

%v = [N T R I]
%dd = [dN/dt; dT/dt; dR/dt; dI/dt]

N = v(1);
T = v(2);
R = v(3);
I = v(4);

dd=[mu*(Thy/Thy_max)-beta*N*(1/(1+(R/kA)^n))- c*N - g*N;
    beta*N*(1/(1+(R/kA)^n)) + a*I*T - b_T*T;
    alpha*(Thy/Thy_max) + epsilon*a*I*R + c*N - b_R*R;
    d*T - e_T*I*T - e_R*I*R - f*I];