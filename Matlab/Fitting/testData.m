
close all
clear all;

data = importdata('JonathanDataToRead.csv')

hours = data(:,1);
cells = data(:,2)*10^6;

%Revised Model:
% dNdt = \alpha N
alpha = 0.012;
Nt    = (20000)*exp(alpha*hours);
semilogy(hours,cells,'o')
hold
semilogy(hours,Nt,'-r');

figure
plot(hours,cells,'o')
hold
plot(hours,Nt,'-r')


%Johnathan's Model
mu = 1244;
g  = 0.1;
c0 = 20000 - mu/g;

% dT/dt = \alpha N
% N(t) = (mu/g)*(Thy(t)/ThyFull) + c0 e^(-gt) 
Nt = (mu/g)*ones(size(cells)) + c0*(exp(-g*hours))

%plot(hours,Nt,'-r')

%dN/dt        = mu*Size(thymus) - g*N
% dN/dt + g*N = mu*

%Integrating factor:
% Multiply both sides by exp(g*t)
% dN/dt e^(gt) + g*N*e^(gt) = mu e^(gt)
% ( N(t) e^(gt) )' = mu e^(gt)
% N(t) e^(gt) = (mu/g)e^(gt) + c0
% N(t) = (mu/g) + c0 e^(-gt) 

%N(0) = 20000
%N(0) = (mu/g) + c0 --> c0 = (20000 - mu/g)