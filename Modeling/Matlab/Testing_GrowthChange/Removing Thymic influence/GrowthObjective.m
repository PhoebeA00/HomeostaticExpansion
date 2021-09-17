function obj= GrowthObjective(p)
global Objective WTerror KOerror
%load data
Gntype = [1, 2];

Rsquare = 0;
% nWT = 0;
% nKO = 0;

for Genotype = Gntype
    %1 = WildType, 2 = IL-2 KO
    if Genotype == 1
        CellData = readtable('../../RawData/ActivatedWTSpleen.csv');
        ProlData = readtable('../../RawData/WTProl.csv');
    elseif Genotype == 2
        CellData = readtable('../../RawData/ActivatedKOSpleen.csv');
        ProlData = readtable('../../RawData/KOProl.csv');
    end

    CellData = CellData(:,{'NaiveCT', 'ActivatedCD4CT', 'X4TregCT', ...
        'ThymicNaive', 'ActivatedNaiveCT', ...
        'ThymicDerivedTregsCT', 'NaiveDerivedTregsCT' ... 
        'hours'});

    ProlData = ProlData(:,{ 'NaiveProlCT', 'ActivatedProlCT', 'X4TregProlCT', ...
        'hours'});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Numbers for the loop will choose the right data
    % The order of which the data/model is set up is critically important
    %                                                               
    DataUsed = [1, 2, 3, 4, 5, 6, 7]; 
    % 1 = Naive T cells
    % 2 = Activated CD4 cells
    % 3 = All Tregs
    % 4 = Thymic Derivied Naive T Cells
    % 5 = Activated Naive T Cells
    % 6 = Thymic Derived Tregs
    % 7 = Naive Derived Tregs
    % 8 = Hours
    %The data location is equivalent in the ModelData df            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ModelData = SimulateGrowth(p, Genotype);

    %Setting up the hours
    DataHours = unique(CellData.hours);
    %Calculating R squared from the NON-proliferating populations
    
    
%     disp('**************HERES THE GENOTYPE LOOP****************')
%     disp(['Genotype = ', num2str(Genotype)])
%     disp(['Rsquare = ', num2str(Rsquare)])
%     disp('***********************************************************')
   
    
    for i = DataUsed
        %Calculates each data used at a time
        %Cell data numbers match 
        Cells = CellData(:,[i,8]); %Grabs data and hours
        SimulationData = ModelData(:,i);
        

        for j = 1:length(DataHours)
            %Pull out all of the data needed; condition is hour match
            hour = DataHours(j);      
            CellDataIndex = Cells.hours == hour;
            CellDataForRSqr = Cells{CellDataIndex,1}; %Do not need hours anymore
            SimulationValue = SimulationData(hour);


            for h = 1:length(CellDataForRSqr)

                CellValue = CellDataForRSqr(h);
                %RSquareValue = (SimulationValue - CellValue).^2;
                RSquareValue = ((SimulationValue - CellValue)/CellValue).^2;
                Rsquare = Rsquare + RSquareValue;
                
                if hour < 289
                    RSquareValue = RSquareValue/10000;
                    Rsquare = Rsquare + RSquareValue;
                else
                    RSquareValue = RSquareValue*10000;
                    Rsquare = Rsquare + RSquareValue;
                end   
              %{  
                if Genotype == 1
                    nWT = nWT +1;
                elseif Genotype == 2
                    nKO = nKO + 1;
                end
                %}
                
                
%                 disp(['Data Number = ' num2str(i)])
%                 disp(['Hour = ' num2str(hour)])
%                 disp(['Cell Value = ' num2str(CellValue)])
%                 disp(['SimulationValue = ' num2str(SimulationValue)])
%                 disp(['RSquareValue = ' num2str(RSquareValue)])
%                 disp(['Rsquare = ' num2str(Rsquare)])
%                 disp(' ')
              

            end

        end

    end
    
%     disp('**************START OF PROL CALCULATIONS****************')
%     disp(['Genotype = ', num2str(Genotype)])
%     disp(['Rsquare = ' num2str(Rsquare)])
%     disp(' ')
%     disp('***********************************************************')
    
    
    %Calculating R squared from the proliferating populations
    ProlDataUsed = [1, 2, 3];
    % 1 = Self replicating naive T cells
    % 2 = Self replicating activated T cells
    % 3 = Self replicating Tregs
    % 4 = Hours
    
    %{
    if Genotype == 2
        disp('KO Objective before Proliferation data')
        disp(['Total Objective = ' num2str(Rsquare)])
    end
    %}
    
    for i = ProlDataUsed
        %Calculates each data used at a time
        %Cell data numbers match 
        Cells = ProlData(:,[i,4]); %Grabs data and hours
        SimulationData = ModelData(:,i+7); %The +7 makes sure that the right data in the ModelData is chosen
%                 disp(['Data Number = ' num2str(i)])
%                 disp(['Hour = ' num2str(hour)])
%                 disp(['Cell Value = ' num2str(CellValue)])
%                 disp(['SimulationValue = ' num2str(SimulationValue)])
%                 disp(['RSquareValue = ' num2str(RSquareValue)])
%                 disp(['Rsquare = ' num2str(Rsquare)])
%                 disp(' ')
%                 

        for j = 1:length(DataHours)
            hour = DataHours(j);      
            CellDataIndex = Cells.hours == hour;
            CellDataForRSqr = Cells{CellDataIndex,1}; %Do not need hours anymore
            SimulationValue = SimulationData(hour);
            
            


            for h = 1:length(CellDataForRSqr)

                CellValue = CellDataForRSqr(h);
                if (CellValue ~= 1)
                    %RSquareValue = (SimulationValue - CellValue).^2;
                    RSquareValue = ((SimulationValue - CellValue)/CellValue).^2;
                    
                    if hour < 289
                        RSquareValue = RSquareValue/10000;
                        Rsquare = Rsquare + RSquareValue;
                    else
                        RSquareValue = RSquareValue*10000;
                        Rsquare = Rsquare + RSquareValue;
                    end                  
                    
                    
                    %{
                    if Genotype == 1
                        nWT = nWT +1;
                    elseif Genotype == 2
                        nKO = nKO + 1;
                    end
                    %}
                end
                
                
                
                
%                 disp(['Data Number = ' num2str(i)])
%                 disp(['Hour = ' num2str(hour)])
%                 disp(['Cell Value = ' num2str(CellValue)])
%                 disp(['SimulationValue = ' num2str(SimulationValue)])
%                 disp(['RSquareValue = ' num2str(RSquareValue)])
%                 disp(['Rsquare = ' num2str(Rsquare)])
%                 disp(' ')
%                 
                
            end

        end
        
    end
    
    if Genotype == 1
        WTerror = Rsquare;        
    elseif Genotype == 2
        KOerror = Rsquare - WTerror;
    end
    
end

disp(['Total Objective = ' num2str(Rsquare)])
disp(['WT Objective = ' num2str(WTerror)])
disp(['KO Objective = ' num2str(KOerror)])
%disp(['WT n = ' num2str(nWT)])
%disp(['KO n = ' num2str(nKO)])
disp(' ')

Objective = Rsquare;
obj = Rsquare;
    
