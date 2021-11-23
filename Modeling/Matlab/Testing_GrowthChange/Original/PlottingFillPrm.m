function PlottingFillPrm(ModelRates, tx, PlotType)

%{
The order of the plots 
subplot(3,4, n)
1 = NaiveCT -           StatsOfCells(:,:,1,1) = NaiveCTWT;
2 = Prol Naive -        StatsOfCells(:,:,8,1) = NprolWT;
3 = ThyNaive -          StatsOfCells(:,:,4,1) = ThyNWT;
4 = EMPTY
5 = ActivatedCT -     StatsOfCells(:,:,2,1) = ActTCTWT;
6 = ActTProl -           StatsOfCells(:,:,9,1) = TprolWT;
7 = ActN -                StatsOfCells(:,:,5,1) = ActNWT;
8 = IL-2 -                  StatsOfCells(:,:,11,1) = IWT;
9 = TregCT -             StatsOfCells(:,:,3,1) = TregCTWT;
10 =  TregProl -        StatsOfCells(:,:,9,1) = TprolWT;
11 = ThyTregs -        StatsOfCells(:,:,6,1) = ThyRWT;
12 = NaiveTregs -     StatsOfCells(:,:,7,1) = DiffRWT;

Location of the Stats StatsOfCells(:,x,-,-)
1 - Means
2 - Standard Deviation
3 - +1 std
4 - -1 std
5 - Lowest 10%ile
6 - Highes 90%ile
%}

%Choosing the plot type

if PlotType == "Percentile"
    UpperLimits = 6;
    LowerLimits = 5;
elseif PlotType == "Std"
    UpperLimits = 3;
    LowerLimits = 4;
end


Genotype = [1,2];
   
%Setting up variables for saving the figure    
loc = './Plots/LHS_Parameters/';
plt1 = append(loc, 'NaiveBalance_LHS_prmt.png');

%Figure Positioning
left = 0;
bottom = 400;
width = 1800;
height = 1050;


PLT1 = figure('Visible', 'off');
set(PLT1,'Position',[left bottom width height]); %This sents the position of the figure itself


strWT = '#9B13A2';
strKO = '#09E600';

colorWT = sscanf(strWT(2:end),'%2x%2x%2x',[1 3])/255;
colorKO = sscanf(strKO(2:end),'%2x%2x%2x',[1 3])/255;

%Proliferating Naive
subplot(1,1,1)

plot(tx, ModelRates(2).NaiveProl(:,1,1),'black-', 'LineWidth', 4)%Plotting the means
hold on
plot(tx, ModelRates(2).NaiveProl(:,1,2),'black--', 'LineWidth', 4)%Plotting the means
tx2 = [tx, fliplr(tx)];
%WT fill
inBetween = [ModelRates(2).NaiveProl(:,UpperLimits,1)', ...
    fliplr(ModelRates(2).NaiveProl(:,LowerLimits,1)')];
fill(tx2, inBetween, colorWT, 'FaceAlpha',0.25);
%KO fill
inBetween = [ModelRates(2).NaiveProl(:,UpperLimits,2)', ...
    fliplr(ModelRates(2).NaiveProl(:,LowerLimits,2)')];
fill(tx2, inBetween, colorKO, 'FaceAlpha',0.25);

title('Prol Naive T Cells')
ylabel('Cell Numbers')
xlabel('Age in Hours')
legend('WT', 'IL-2 KO', 'WT', 'IL-2 KO')
hold off




saveas(PLT1, sprintf(plt1))






























































end