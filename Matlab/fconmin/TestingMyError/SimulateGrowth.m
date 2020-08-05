%function dTdt = SimulateGrowth(p)

close all; clear all; clc
 
tx = 1:500;
ModelData = zeros(length(tx),0);

p = GetParameters(10); %Adjust parameter choice here %%%%%%%%%%%%%%%%

%Initial Conditions
N = p(20);
T = p(21);
R = p(22);
I = p(23);

T0 = [N,T,R,I];
T0 = T0./3; % I think I did this in the python code because it makes it match the data better

%First Entry of the data frame
ModelData(1,1) = T0(1);
ModelData(1,2) = T0(2);
ModelData(1,3) = T0(3);
ModelData(1,4) = T0(4);


for i = 1:length(tx)-1
    %disp(['HEY HERE IS THE I', num2str(i)])
    ts = [tx(i),tx(i+1)];
    sol = ode15s(@(t,x)Growth(t,x,p),ts,T0);
    T0 = [sol.y(1,end),sol.y(2,end),sol.y(3,end),sol.y(4,end)];
    ModelData(i+1,1) = T0(1);
    ModelData(i+1,2) = T0(2);
    ModelData(i+1,3) = T0(3);
    ModelData(i+1,4) = T0(4);
end

figure(1)

subplot(2,2,1)
plot(tx,ModelData(:,1))
ylabel('Naive T Cells')

subplot(2,2,2)
plot(tx,ModelData(:,2))
ylabel('Activated T Cells')

subplot(2,2,3)
plot(tx,ModelData(:,3))
ylabel('Tregs')

subplot(2,2,4)
plot(tx,ModelData(:,4))
ylabel('Il-2')



