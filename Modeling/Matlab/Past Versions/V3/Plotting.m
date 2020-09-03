%% - Plotting Naive and Activated normal and optimized
global p0 pOpt

Data = readtable('ActivatedDataForModel_WTSPLEEN.csv');
CellData = Data(:,{'NaiveT_Cells', 'Activated_Cells', ...                  
     'X4Treg_Cells', 'hours'});
 
ModelData = SimulateGrowth(p0);
OptModel = SimulateGrowth(pOpt);

%Initial results
figure(1)

subplot(2,2,1)
scatter(CellData.hours, CellData.NaiveT_Cells)
hold on
plot(tx,ModelData(:,1))
title('Naive T Cells and NonOpt Model')
xlabel('Hours')

subplot(2,2,2)
scatter(CellData.hours, CellData.Activated_Cells)
hold on
plot(tx,ModelData(:,2))
title('Activated T Cells and NonOpt Model')
xlabel('Hours')

subplot(2,2,3)
scatter(CellData.hours, CellData.NaiveT_Cells)
hold on
plot(tx, OptModel(:,1))
title('NaiveT and Optimized Model')

subplot(2,2,4)
scatter(CellData.hours, CellData.Activated_Cells)
hold on
plot(tx,OptModel(:,2))
title('Activated T and Optimized Model')

%% - Plotting without optimization

%Initial results
figure(1)

subplot(2,2,1)
scatter(CellData.hours, CellData.NaiveT_Cells)
hold on
plot(tx,ModelData(:,1))
title('Naive T Cells and NonOpt Model')
xlabel('Hours')

subplot(2,2,2)
scatter(CellData.hours, CellData.Activated_Cells)
hold on
plot(tx,ModelData(:,2))
title('Activated T Cells and NonOpt Model')
xlabel('Hours')

subplot(2,2,3)
scatter(CellData.hours, CellData.NaiveT_Cells)
hold on
plot(tx, OptModel(:,1))
title('NaiveT and Optimized Model')

subplot(2,2,4)
scatter(CellData.hours, CellData.Activated_Cells)
hold on
plot(tx,OptModel(:,2))
title('Activated T and Optimized Model')

%% Plotting Treg data
global p0 pOpt

Data = readtable('ActivatedDataForModel_WTSPLEEN.csv');
CellData = Data(:,{'NaiveT_Cells', 'Activated_Cells', ...                  
     'X4Treg_Cells', 'hours'});
 
ModelData = SimulateGrowth(p0);
OptModel = SimulateGrowth(pOpt);

figure(2)

subplot(1,2,1)
scatter(CellData.hours, CellData.X4Treg_Cells)
hold on
plot(tx,ModelData(:,3))
title('Tregs and NonOpt Model')
xlabel('Hours')

subplot(1,2,2)
scatter(CellData.hours, CellData.X4Treg_Cells)
hold on
plot(tx,OptModel(:,3))
title('Tregs and Optimized Model')
xlabel('Hours')





