function ModelData = SimulateGrowth(p)
global tx N T R I m ThyR RplR DiffR

ModelData = zeros(length(tx),0);

T0 = [N,T,R, ThyR, RplR, DiffR, I, m];

%First Entry of the data frame
ModelData(1,1) = T0(1); % Naive
ModelData(1,2) = T0(2); % Activated T cells
ModelData(1,3) = T0(3); % Tregs
ModelData(1,4) = T0(4); % Thymic Derivied Tregs
ModelData(1,5) = T0(5); % Self Replicating Tregs
ModelData(1,6) = T0(6); % Naive Derived Tregs
ModelData(1,7) = T0(7); % IL-2
ModelData(1,8) = T0(8); % Thymus mass


for i = 1:length(tx)-1
    ts = [tx(i),tx(i+1)];
    %options1 = odeset('NonNegative', 1);
    sol = ode15s(@(t,x)Growth(t,x,p,i),ts,T0); %i is passed for the Thymus equation
    T0 = [sol.y(1,end),sol.y(2,end),sol.y(3,end),sol.y(4,end),sol.y(5,end),...
        sol.y(6,end), sol.y(7,end), sol.y(8,end)];
    ModelData(i+1,1) = T0(1); % Naive
    ModelData(i+1,2) = T0(2); % Activated T cells
    ModelData(i+1,3) = T0(3); % Tregs
    ModelData(i+1,4) = T0(4); % Thymic Derivied Tregs
    ModelData(i+1,5) = T0(5); % Self Replicating Tregs
    ModelData(i+1,6) = T0(6); % Naive Derived Tregs
    ModelData(i+1,7) = T0(7); % IL-2
    ModelData(i+1,8) = T0(8); % Thymus Weight
end
%save('SimmyTest', 'ModelData') %Used this to trouble shoot a problem