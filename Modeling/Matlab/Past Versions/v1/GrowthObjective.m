function obj= GrowthObjective(p)
%load data
Data = readtable('ActivatedDataForModel_WTSPLEEN.csv');
CellData = Data(:,{'NaiveT_Cells', 'Activated_Cells', ...                  
     'X4Treg_Cells', 'hours'}); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Make all changes here for optimization %%%%%%%%%%%%%%%%%
%                                                               %
DataUsed = [1, 2]; % Naive = 1, Activ = 2, Treg = 3, hours = 4   %
%The data location is equivalent in the ModelData df            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mu = p(1);
beta = p(2);
g = p(3);

disp(['mu: ' num2str(mu)])
disp(['beta: ' num2str(beta)])
disp(['g: ' num2str(g)])

ModelData = SimulateGrowth(p);

%Setting up the hours
DataHours = unique(CellData.hours);

Rsquare = 0;
%ModelData = SimulateGrowth(p0);
for i = DataUsed
    %disp(['Data number = ', num2str(i)])
    
    Cells = CellData(:,[i,4]);
    SimulationData = ModelData(:,i);
    
    
    for j = 1:length(DataHours)
        colnm = CellData.Properties.VariableNames;
        %{
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

    
