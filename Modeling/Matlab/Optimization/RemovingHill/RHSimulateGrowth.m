function ModelData = RHSimulateGrowth(p)
global tx N T R I m

ModelData = zeros(length(tx),0);

T0 = [N,T,R,I,m];

%First Entry of the data frame
ModelData(1,1) = T0(1);
ModelData(1,2) = T0(2);
ModelData(1,3) = T0(3);
ModelData(1,4) = T0(4);
ModelData(1,5) = T0(5); %%The average of the Thymus weight at age 0

for i = 1:length(tx)-1
    
    ts = [tx(i),tx(i+1)];
    sol = ode15s(@(t,x)RHGrowth(t,x,p,i),ts,T0); %i is passed for the Thymus equation
    T0 = [sol.y(1,end),sol.y(2,end),sol.y(3,end),sol.y(4,end),sol.y(5,end)];
    ModelData(i+1,1) = T0(1); % Naive
    ModelData(i+1,2) = T0(2); % Activated T cells
    ModelData(i+1,3) = T0(3); % Tregs
    ModelData(i+1,4) = T0(4); % IL-2
    ModelData(i+1,5) = T0(5); % Thymus Weight
end