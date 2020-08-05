%SIR disease Model

%%%%%
% This Matlab Code solves an SIR model and fits the parameters to the
% data
% This is matlab SCRIPT. (it does not say function)
%%%%%

%%%%%
%Dependencies: 
%function SIR.m 
%%%%%%

%The first step in using
clear all; close all;

global beta gamma

%initial parameters and initial conditions
%initial parameters
day=[1:1:14];
infected=[3 6 25 73 222 294 258 237 191 125 69 27 11 4];

%This is setting choices for the initial conditions
beta=.0025; %What are these units? 
gamma=0.5;  %What are these units? 
Init=[762; 1; 0];	%initial conditions

%This is where you SOLVE the differential equation
options = odeset('RelTol',1e-6,'AbsTol',1e-6); %These are options for the ode solver
[t d]=ode45('SIR',[0 15],Init,options); %This solves the Differential Equation
                                        %ode45 = Runga Kutta 4,5 Method
                                        %(You should know a bit about
                                        %this.)
%From ode45, you get back:
%t <- vector of all the times. 
%d <- the states of the system at these times 
%      Our system is 3 variables: d(times,#numStates)
%      d(i,1) = Susceptibles at time t(i)
%      d(10,2) = The number of infected at time t(2)

TotalIntervals=length(d);
for i=1:TotalIntervals
    popln(i)=d(i,1)+d(i,2)+d(i,3);
end
%total population=susceptible at time t+ infected +recovered
figure(1)
set(gcf,'Color','white')
set(0,'defaultaxesfontsize',22);
plot(t,d(:,1),'-b','LineWidth',3);  %Blue, Susceptible
hold on
plot(t,d(:,2),'k--','LineWidth',3); %Black , INfected -- dashed
hold on
plot(t,d(:,3),'g--','LineWidth',3); %Green, Recovered -- dashed
hold on
plot(t,popln,'m--','LineWidth',3) %Magenta Total Population
xlabel('Time','FontWeight','Bold')
ylabel('States','FontWeight','Bold')
legend('S','I','R');
title('SIR Model, \beta=0.0025, \gamma=0.3','FontWeight','Bold')
grid on

%This creates figure 2; Figure 2 shows the "data"
figure(2)
set(gcf,'Color','white')
set(0,'defaultaxesfontsize',22);
plot(day,infected,'-b','LineWidth',3);
xlabel('Time (days)','FontWeight','Bold')
ylabel('Infected','FontWeight','Bold')
title('Data','FontWeight','Bold')
grid on

%Figure 3 plots the "data" with the model output"
figure(3)
set(gcf,'Color','white')
set(0,'defaultaxesfontsize',22);
plot(day,infected,'-b','LineWidth',3);
hold on
plot(t,d(:,2),'k--','LineWidth',3);
xlabel('Time (days)','FontWeight','Bold')
ylabel('Infected','FontWeight','Bold')
title('Comparison','FontWeight','Bold')
legend('Data','Infected, model')
grid on

P3=polyfit(day,infected,3);
P4=polyfit(day,infected,4);
P5=polyfit(day,infected,5);

P3t=polyval(P3,t);
P4t=polyval(P4,t);
P5t=polyval(P5,t);
options = odeset('RelTol',1e-6,'AbsTol',1e-6); %These are options for the ode solver
[t d]=ode45('SIR',[0 15],Init,options); 
Ipch=pchip(day,infected,t);

figure(4)
set(gcf,'Color','white')
set(0,'defaultaxesfontsize',22);
plot(day,infected,'*b','MarkerSize',15);
hold on
plot(t,P3t,'k--','LineWidth',3);
hold on
plot(t,P4t,'g--','LineWidth',3);
hold on
plot(t,P5t,'r--','LineWidth',3);
hold on
plot(t,Ipch,'--m','LineWidth',3);
xlabel('Time (days)','FontWeight','Bold')
ylabel('Infected','FontWeight','Bold')
title('Comparison','FontWeight','Bold')
legend('Data','P3','P4','P5','pchip')
grid on

p=[beta; gamma];
% err=0;
% while(err>1e-4)
%     err = odefit(t,Ipch,p)
% end