function ModelData = Plot_Simulation(ParameterNumber)

tx = 1:432; %Max hours in the simulation
ModelData = zeros(length(tx),0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%     Choose parameters here      %%%%%%%%%%%%%%%%%
%                                                               %
p = GetParameters(ParameterNumber);                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Change parameters that will be fitted accordingly using the p parameter
%alpha = p(1);
%Thy = p(2);
%Thy_max = p(3);
%epsilon = p(4);
%a = p(5);
%c = p(6);
%b_R = p(7);
%mu = p(8);
%beta = p(9);
%g = p(10);
%b_T = p(11);
%d = p(12);
%e_T = p(13);
%e_R = p(14);
%f = p(15);
%kA = p(16);
%n = p(17);


%Initial Conditions
N = p(20);
T = p(21);
R = p(22);
I = p(23);
m = 0.0023; %Average of the Thymus weight starting at day 0

T0 = [N,T,R,I,m];

%First Entry of the data frame
ModelData(1,1) = T0(1);
ModelData(1,2) = T0(2);
ModelData(1,3) = T0(3);
ModelData(1,4) = T0(4);
ModelData(1,5) = T0(5);

disp('Simulating...')
for i = 1:length(tx)-1
    ts = [tx(i),tx(i+1)];%Time step by hour
    sol = ode15s(@(t,x)Plot_growth(t,x,p,i),ts,T0);
    T0 = [sol.y(1,end),sol.y(2,end),sol.y(3,end),sol.y(4,end),sol.y(5,end)];
    ModelData(i+1,1) = T0(1); %Naive
    ModelData(i+1,2) = T0(2); %Activated T Cells
    ModelData(i+1,3) = T0(3); %Tregs
    ModelData(i+1,4) = T0(4); %IL-2
    ModelData(i+1,5) = T0(5); %Thymus Weight
end