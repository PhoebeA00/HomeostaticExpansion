close all; clear all; clc
global MaxTime

MaxTime = 0:432; %Maximum amount of hours in data

%Parameters to be optimized
K = 0.04;
lambda = 0.016932;
p0 = [K, lambda];%p0 will constantly change until the algorithm finds the
%least rsquare

%Lower and Upper Bounds
lb = [0.00, 0.016932];
ub = [0.09, 0.016932];

% no linear constraints
A = [];
b = [];
Aeq = [];
beq = [];
nlcon = [];

disp('Beginning Optimization...')
% options = optimoptions(@fmincon,'Algorithm','interior-point');
pOpt = fmincon(@ThymObjective,p0,A,b,Aeq,beq,lb,ub,nlcon);


%Optimized Parameters

k = pOpt(1);
lambda = pOpt(2);
disp('The optimized parameters are')
disp(['K: ' num2str(k)])
disp(['lambda: ' num2str(lambda)])


%%
ThymusData = readtable('../Data/ThymusData.csv');
ModelData = SimulateThymus(pOpt);

scatter(ThymusData.Hours, ThymusData.Weight)
hold on
plot(MaxTime, ModelData)
title('Thymus weight in grams')
ylabel('Weight in grams')
xlabel('Time in hours')
legend('Data', 'Logistical Growth')







