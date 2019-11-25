

%Clear is always good to run at the beginning (removes all variables from
%memory
clear;
close all;

%to = initial time
to = 0;

%tf = final time
tf = 10;

%Initial Conditions of your cellular populations
yo = [33*2 0.25*2]; %Initial Conditions [T(0), I(0)]

%Steady State Values: 
%I(ss) = 0.2500
%T(ss) = 33.3333


%SOLVE the differential equation
[t y] = ode45('yprf',[to tf],yo);
plot(y(:,2),y(:,1))
hold
plot(yo(2),yo(1),'ro')
title('Equation 1')
ylabel('Activated T Cells')
xlabel('IL-2')

figure
plot(t,y(:,1))
title('Activated T Cells')

figure
plot(t,y(:,2))
title('IL-2')

