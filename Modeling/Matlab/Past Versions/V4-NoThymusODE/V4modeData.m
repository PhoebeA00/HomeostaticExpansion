%System of ODE's that simulate the adaptive Immune system

%%%%%
% This Matlab Code solves an Immun model
% This is matlab SCRIPT. (it does not say function)
%%%%%

%%%%%
%Dependencies: 
%function V4.m 
%%%%%%


%The first step in using
clear all; close all;

global alpha  Thy  Thy_max epsilon a c b_R mu beta g b_T d e_T e_R f kA n

%initial parameters and initial conditions
%initial parameters
day=[1:1:100];


%%%%%%%%%%%%
%  Thymus  %
%%%%%%%%%%%%
alpha = 66;  %------------ T Regulatory Cells
mu = 12427;   %------------ Naive T cells
Thy = 1;        %------------ Size of the thymus
Thy_max = 1;    %------------ Max size of the thymus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Naive T cell Differentiation Rates   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 598;       %-------------To T regulatory Cells
beta = 289;      %------------ To activated T cells

%%%%%%%%%%%
%  Tregs  %
%%%%%%%%%%%
epsilon = 0.0328;    %------------T regulatory cell Self replication
n       = 1;    %------------hill coefficient
kA      = 700;   %------------halfSaturationRate 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IL-2 Cytokine Expression and Consumption  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 1000;       %-----------T Cell Expression
a = 1.102;        %-----------Activated T cells
e_T = 3000;     %---------- T Cell Consumption Rate
e_R = 8000;     %---------- T Reg Consumption Rate

%%%%%%%%%%%%%%%%%%
%  Death Rates   %
%%%%%%%%%%%%%%%%%%
g = 126;       %-----------Naive T cells
b_T = 126;      %-----------Activated T cells
b_R = 55;      %-----------Regulatory T Cells
f = 0.1;          %-------------IL-2 Cytokine


%Initial Conditions Init = [N T R I]
Init=[50; 0; 0; 0];	%initial conditions

%This is where you SOLVE the differential equation
options = odeset('RelTol',1e-6,'AbsTol',1e-6); %These are options for the ode solver
[t v]=ode45('V4',[0 100],Init,options); %This solves the Differential Equation
                                        %ode45 = Runga Kutta 4,5 Method
                                        %(You should know a bit about
                                        %this.)
                                       
                                        
                                        
%From ode45, you get back:
%t <- vector of all the times. 
%v <- the states of the system at these times 
%      Our system is 3 variables: d(times,#numStates)
%      d(i,1) = Susceptibles at time t(i)
%      d(10,2) = The number of infected at time t(2)

%TotalIntervals=length(v);
%for i=1:TotalIntervals
 %   popln(i)=v(i,1)+v(i,2)+v(i,3)+v(i,4);
%end

figure(1)
set(gcf,'Color','white')
set(0,'defaultaxesfontsize',22);
plot(t,v(:,1), '-b', 'LineWidth', 3); %Blue, Naive T cells
xlabel('Time','FontWeight','Bold')
ylabel('Cells','FontWeight','Bold')
title('Naive Kinetics','FontWeight','Bold')
grid on

figure(2)
set(gcf,'Color','white')
set(0,'defaultaxesfontsize',22);
plot(t,v(:,2), '-k', 'LineWidth', 3); %Blue, Activated T cells
xlabel('Time','FontWeight','Bold')
ylabel('Cells','FontWeight','Bold')
title('Activated T Kinetics','FontWeight','Bold')
grid on

figure(3)
set(gcf,'Color','white')
set(0,'defaultaxesfontsize',22);
plot(t,v(:,3), '-g', 'LineWidth', 3); %Green, Tregs
xlabel('Time','FontWeight','Bold')
ylabel('Cells','FontWeight','Bold')
title('Treg Kinetics','FontWeight','Bold')
grid on

figure(4)
set(gcf,'Color','white')
set(0,'defaultaxesfontsize',22);
plot(t,v(:,4), '-p', 'LineWidth', 3); %, IL-2 Cytokine
xlabel('Time','FontWeight','Bold')
ylabel('Cells','FontWeight','Bold')
title('IL-2 Cytokine','FontWeight','Bold')
grid on


%total population=susceptible at time t+ infected +recovered
figure(1)
set(gcf,'Color','white')
set(0,'defaultaxesfontsize',22);
plot(t,v(:,1),'-b','LineWidth',3);  %Blue, Naive T cells
hold on
plot(t,v(:,2),'k--','LineWidth',3); %Black , Activated T Cells -- dashed
hold on
plot(t,v(:,3),'g--','LineWidth',3); %Green, T Regulatory Cells-- dashed
hold on
plot(t,v(:,4),'p--', 'LineWidth',3); %Purple, IL-2 Cytokines-- dashed
%plot(t,popln,'m--','LineWidth',3) %Magenta Total Population
xlabel('Time','FontWeight','Bold')
ylabel('Cells','FontWeight','Bold')
legend('N','T','R', 'I');
title('Immune Kinetics','FontWeight','Bold')
grid on


