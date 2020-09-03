%%
clear all
ThymusData = readtable('ThymusData.csv');
ThymusDataFinal = ThymusData(:,{'Age', 'Weight'});
%scatter(24*ThymusDataFinal.Age, ThymusDataFinal.Weight)

time = 24*(ThymusDataFinal{:,'Age'}); %Changing age to hourly
weight = ThymusDataFinal{:,'Weight'};


[fit, stats] = fit(time,weight,'poly3'); 

%Best fit Poly4


%% Prepping for plot

p1 = -3.725e-09;
p2 = 2.419e-06;
p3 = -0.0001933;
p4 = 0.00398;


Mdlt = 0:1:432;
FitEq = @(Mdlt) p1*Mdlt.^3 + p2*Mdlt.^2 + p3*Mdlt + p4;
%FitEq = inline('p1*Mdlt^4 + p2*Mdlt^3 + p3*Mdlt^2 + p4*Mdlt + p5');

%%
figure
plot(time, weight, 'o');
hold on

Mdlt = 0:1:432;
plot(Mdlt, FitEq(Mdlt), 'r')


