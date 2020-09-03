function obj= GrowthObjective(p)
%load data
Data = readtable('./RawData/ActivatedDataForModel_WTSPLEEN.csv');
CellData = Data(:,{'NaiveT_Cells', 'Activated_Cells', ...                  
     'X4Treg_Cells', 'hours'}); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Make all changes here for optimization %%%%%%%%%%%%%%%%%
%                                                               %
DataUsed = [1, 3]; % Naive = 1, Activ = 2, Treg = 3, hours = 4   %
%The data location is equivalent in the ModelData df            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ModelData = SimulateGrowth(p);

%Setting up the hours
DataHours = unique(CellData.hours);

Rsquare = 0;
%ModelData = SimulateGrowth(p0);
for i = DataUsed
    %Calculates each data used at a time
    %Cell data numbers match 
    Cells = CellData(:,[i,4]); %Grabs data and hours
    SimulationData = ModelData(:,i);
    
    
    for j = 1:length(DataHours)
        
        %{
        colnm = ThymusData.Properties.VariableNames;
        disp(' ')
        disp(' ')
        disp(['Population = ', colnm{i}])
        disp(['Working on hour: ', num2str(DataHours(j))])
        %}
        hour = DataHours(j);
        CellDataIndex = Cells.hours == hour;
        CellDataForRSqr = Cells{CellDataIndex,1}; %Do not need hours anymore
        SimulationValue = SimulationData(hour);
        
        
        for h = 1:length(CellDataForRSqr)
            %disp(['Simulation Data: ', num2str(SimulationValue)])
            %disp(['Hour: ', num2str(hour), ' Data: ', ...
            %    num2str(CellDataForRSqr(h))])
            %disp(' ')
            
            CellValue = CellDataForRSqr(h);
            RSquareValue = (SimulationValue - CellValue).^2;
            %RSquareValue = ((SimulationValue - CellValue)./CellValue).^2;
            Rsquare = Rsquare + RSquareValue;
           
        end

    end
    
end
disp(['Objective = ' num2str(Rsquare)])
obj = Rsquare;

    
