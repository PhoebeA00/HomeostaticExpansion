function ModelData = SimulateGrowth(p)
global tx N T R I

ModelData = zeros(length(tx),0);

T0 = [N,T,R,I];

%First Entry of the data frame
ModelData(1,1) = T0(1);
ModelData(1,2) = T0(2);
ModelData(1,3) = T0(3);
ModelData(1,4) = T0(4);

%disp('Simulating...')

for i = 1:length(tx)-1
    
    %{
    Instance = rem(i, 25);
    if Instance == 0
        disp(['Instance = ', num2str(i)])
    end
    %}
        
    ts = [tx(i),tx(i+1)];
    sol = ode15s(@(t,x)Growth(t,x,p),ts,T0);
    T0 = [sol.y(1,end),sol.y(2,end),sol.y(3,end),sol.y(4,end)];
    ModelData(i+1,1) = T0(1);
    ModelData(i+1,2) = T0(2);
    ModelData(i+1,3) = T0(3);
    ModelData(i+1,4) = T0(4);
end