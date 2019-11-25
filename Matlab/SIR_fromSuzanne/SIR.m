%This is an example of a function.
function dd=SIR(t,d)
global beta gamma

%d = [S I R]
%dd = [dS/dt; dI/dt; dR/dt]

dd=[-beta*d(2)*d(1);
		beta*d(2)*d(1)-gamma*d(2);
        gamma*d(2)];