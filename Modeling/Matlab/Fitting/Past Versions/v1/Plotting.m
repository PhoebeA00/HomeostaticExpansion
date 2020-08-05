Data = readtable('ActivatedDataForModel_WTSPLEEN.csv');
CellData = Data(:,{'NaiveT_Cells', 'Activated_Cells', ...                  
     'X4Treg_Cells', 'hours'});

%%
ModelData = SimulateGrowth(p0);
OptModel = SimulateGrowth(pOpt);



%%

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