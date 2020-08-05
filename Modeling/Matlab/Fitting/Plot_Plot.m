%% - Naive T cells and Activated T cells
%close all; clear all; clc
%Choose Parameter number
ModelD = Plot_Simulation(13);

Data = readtable('ActivatedDataForModel_WTSPLEEN.csv');
CellData = Data(:,{'NaiveT_Cells', 'Activated_Cells', ...                  
     'X4Treg_Cells', 'hours'});
 
tx = 1:432;

 %Initial results
figure(2)

subplot(1,2,1)
scatter(CellData.hours, CellData.NaiveT_Cells)
hold on
plot(tx,ModelD(:,1))
title('Naive T Cells')
xlabel('Hours')

subplot(1,2,2)
scatter(CellData.hours, CellData.Activated_Cells)
hold on
plot(tx,ModelD(:,2))
title('Activated T Cells')
xlabel('Hours')

%% Parameter list and comparing it to the data

ModelData = Plot_Simulation(15);

Data = readtable('ActivatedDataForModel_WTSPLEEN.csv');
CellData = Data(:,{'NaiveT_Cells', 'Activated_Cells', ...                  
     'X4Treg_Cells', 'hours'});
 
tx = 1:432;


figure(1)

subplot(1,3,1)
scatter(CellData.hours, CellData.NaiveT_Cells)
hold on 
plot(tx, ModelData(:,1))
title('Naive T Cells')

subplot(1,3,2)
scatter(CellData.hours, CellData.Activated_Cells)
hold on
plot(tx, ModelData(:,2))
title('Activated T cells')

subplot(1,3,3)
scatter(CellData.hours, CellData.X4Treg_Cells)
hold on
plot(tx,ModelData(:,3))
title('T Regulatory Cells')














