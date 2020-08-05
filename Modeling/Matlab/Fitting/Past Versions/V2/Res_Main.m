%% This is the main script for calculating the Residuals of all the data
%Gets residuals for all data types
close all; clear all; clc

Data = readtable('ActivatedDataForModel_WTSPLEEN.csv');
CellData = Data(:,{'NaiveT_Cells', 'Activated_Cells', ...                  
     'X4Treg_Cells', 'hours'});
tx = 1:432;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Make all changes here for optimization %%%%%%%%%%%%%%%%%
%                                                               %
DataUsed = [1, 2, 3]; % Naive = 1, Activ = 2, Treg = 3, hours = 4   %
%The data location is equivalent in the ModelData df            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ModelData = Plot_Simulation(15);

%Setting up the hours
DataHours = unique(CellData.hours);

Residuals = table;
%ModelData = SimulateGrowth(p0);

%Resetting the variables
N_Residuals = 0;
T_Residuals = 0;
R_Residuals = 0; 

for i = 1:length(DataHours)%Cycles through each unique hour

    hour = DataHours(i);
    CellDataIndex = CellData.hours == hour; %all the data for the hour
    CellDataForResid = CellData{CellDataIndex,:};%all the data for the hour
    SimulationValueForResid = ModelData(hour,:); %Row data from that hour   
    
    for DataNumber = DataUsed
        if DataNumber == 1
            N_Residuals = Res_Calculate(SimulationValueForResid, ...
                                        CellDataForResid, DataNumber);
        elseif DataNumber == 2
            T_Residuals = Res_Calculate(SimulationValueForResid, ...
                                        CellDataForResid, DataNumber);
        else 
            R_Residuals = Res_Calculate(SimulationValueForResid, ...
                                        CellDataForResid, DataNumber);
        end        
    end
    %If is here to make sure that there no mistakes in the lengths
    if length(N_Residuals) ~= length(R_Residuals) || ...
                length(N_Residuals) ~= length(T_Residuals)
        error('Lengths Do not match ')
    end
               
    for j = 1:length(N_Residuals)
        row = table(N_Residuals(j), T_Residuals(j),...
            R_Residuals(j), hour);
        Residuals = [Residuals; row];
    end
end

Residuals.Properties.VariableNames = {'N', 'T', 'R', 'hour'};

%%

figure(1)

%Simulation and Data

subplot(2,3,1)
scatter(CellData.hours, CellData.NaiveT_Cells)
hold on 
plot(tx, ModelData(:,1))
title('Naive T Cells')

subplot(2,3,2)
scatter(CellData.hours, CellData.Activated_Cells)
hold on
plot(tx, ModelData(:,2))
title('Activated T cells')

subplot(2,3,3)
scatter(CellData.hours, CellData.X4Treg_Cells)
hold on
plot(tx,ModelData(:,3))
title('T Regulatory Cells')

%Residuals

subplot(2,3,4)
scatter(Residuals.hour, Residuals.N)
hold on
yline(0)
title('Naive T cells')

subplot(2,3,5)
scatter(Residuals.hour, Residuals.T)
hold on
yline(0)
title('Activated T Cells')

subplot(2,3,6)
scatter(Residuals.hour, Residuals.R)
hold on
yline(0)
title('T Regulator Cells')

%LEARN TO USE FIT

