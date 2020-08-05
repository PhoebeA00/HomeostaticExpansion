%function obj= objective(p)
%load data
CellData = readtable('Data.csv');

%Setting up data files 
    %You can control the order in the columsn here.
    %I am going to make the first column the hours and second column the Cell
    %data. This will be the standard
NaiveCellData = CellData(:,{'hour','number'});

%Will save ModelData in a double variable
    %Each order of entry will be indicative of an hour. 
    %Model(3) - Will be hour 4
ModelNaive = [5,7,9,13,18,20,30];


%Setting up the hours
DataHours = unique(NaiveCellData.hour);

%Pulling by hour condition for both Cell Data and model data
%for i = DataHours
CellIndex = find(NaiveCellData.hour == 5);
CellDataForRSqr = NaiveCellData(CellIndex,:); %Pullind data based on index

RSquare = 0; %All of the Rsquare values will be added to this
%First Attempt of the for loop
for i = 1:length(DataHours)
    hour = DataHours(i);
    ModelHourData = ModelNaive(hour);
    CellDataIndex = find(NaiveCellData.hour == hour);
    CellDataForRSqr = NaiveCellData(CellDataIndex,:);
    %disp(['The hour is: ',num2str(hour)])
    %disp(CellDataForRSqr.number)
    %disp(['Length of this data: ',num2str(length(CellDataForRSqr.number))])
    for j = 1:length(CellDataForRSqr.number)
        CellValue = CellDataForRSqr.number(j);
        ModelValue = ModelNaive(hour);
        RSquareValue = (CellValue - ModelValue).^2;
        RSquare = RSquare + RSquareValue;
        %disp(['The hour is: ',num2str(hour)])
        %disp(['The data is:', num2str(CellRSqrData)])
    end
end
disp(['Final R Square Value: ', num2str(RSquare)])
