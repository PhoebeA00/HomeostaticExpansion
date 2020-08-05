function Residuals = Res_Calculate(SimulationValueForResid,... 
    CellDataForResid, DataNumber)

CellData = CellDataForResid(:,DataNumber); %1:N 2:T 3:R
SimData = SimulationValueForResid(DataNumber);
Res = []; %Setting it up 
disp('Cycle Begins')
count = 1;
for h = 1:length(CellData)
       disp(['Cycle number: ', num2str(h) ]) 
    CellValue = CellData(h);
    Residual = (SimData - CellValue);
    
    Res = [Res Residual];
end
Residuals = Res;