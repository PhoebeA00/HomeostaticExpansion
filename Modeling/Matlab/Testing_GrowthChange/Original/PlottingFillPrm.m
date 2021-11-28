function PlottingFillPrm(ModelRates, StatsOfCells, tx, PlotType)

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

%-----------------Setting up variables for plotting-----------------------%

if PlotType == "Percentile"
    UpperLimits = 6;
    LowerLimits = 5;
elseif PlotType == "Std"
    UpperLimits = 3;
    LowerLimits = 4;
end

% Genotype = [1,2];

tx2 = [tx, fliplr(tx)]; %Used for the fill in ranges

%Saving the hex colors that I chose
strWT = '#9B13A2';
strKO = '#09E600';
%Converting the hex to the RGB that matlab uses
colorWT = sscanf(strWT(2:end),'%2x%2x%2x',[1 3])/255;
colorKO = sscanf(strKO(2:end),'%2x%2x%2x',[1 3])/255;

%The opacity of the fill
FaceAlpha = 0.25;

%LineWidth
LineWidth = 1.5;

%Setting up variables for saving the figure    
loc = './Plots/LHS_Parameters/';
plt1 = append(loc, 'Figure_4.png');

%Figure Positioning
left = 0;
bottom = 400;
width = 1800;
height = 1050;


PLT1 = figure('Visible', 'off');
set(PLT1,'Position',[left bottom width height]); %This sents the position of the figure itself

%------------------End of variables-------------------------------------------------%



%%%%%%%%%%%%%%
%Treg Death Rate
%%%%%%%%%%%%%%
subplot(2,4,2)

plot(tx, ModelRates(2).TregDeath(:,1,1),'b-', 'LineWidth', LineWidth)%Plotting the means
hold on
plot(tx, ModelRates(2).TregDeath(:,1,2),'b--', 'LineWidth', LineWidth)%Plotting the means

%WT fill
inBetweenWT = [ModelRates(2).TregDeath(:,UpperLimits,1)', ...
    fliplr(ModelRates(2).TregDeath(:,LowerLimits,1)')];
fill(tx2, inBetweenWT, colorWT, 'FaceAlpha',FaceAlpha);
%KO fill
inBetweenKO = [ModelRates(2).TregDeath(:,UpperLimits,2)', ...
    fliplr(ModelRates(2).TregDeath(:,LowerLimits,2)')];
fill(tx2, inBetweenKO, colorKO, 'FaceAlpha',FaceAlpha);

title('Treg Death Rates')
ylabel('Cell Numbers')
xlabel('Age in Hours')
legend('WT', 'IL-2 KO', 'WT', 'IL-2 KO')
hold off

%%%%%%%%%%%%%%%%%%%
%Zoomed in Treg Death Rate
%%%%%%%%%%%%%%%%%%%
subplot(2,4,6)

plot(tx, ModelRates(2).TregDeath(:,1,1),'b-', 'LineWidth', LineWidth)%Plotting the means
hold on
plot(tx, ModelRates(2).TregDeath(:,1,2),'b--', 'LineWidth', LineWidth)%Plotting the means
%WT fill
inBetweenWT = [ModelRates(2).TregDeath(:,UpperLimits,1)', ...
    fliplr(ModelRates(2).TregDeath(:,LowerLimits,1)')];
fill(tx2, inBetweenWT, colorWT, 'FaceAlpha',FaceAlpha);
%KO fill
inBetweenKO = [ModelRates(2).TregDeath(:,UpperLimits,2)', ...
    fliplr(ModelRates(2).TregDeath(:,LowerLimits,2)')];
fill(tx2, inBetweenKO, colorKO, 'FaceAlpha',FaceAlpha);

title('Zoomed in Treg Death Rates')
ylabel('Cell Numbers')
xlabel('Age in Hours')
legend('WT', 'IL-2 KO', 'WT', 'IL-2 KO')
xlim([0 120])
ylim([0 30])
hold off

%%%%%%%%
%Prol Tregs
%%%%%%%%
subplot(2,4,3)

plot(tx, ModelRates(2).ProlTregs(:,1,1),'b-', 'LineWidth', LineWidth)%Plotting the means
hold on
plot(tx, ModelRates(2).ProlTregs(:,1,2),'b--', 'LineWidth', LineWidth)%Plotting the means

%WT fill
inBetweenWT = [ModelRates(2).ProlTregs(:,UpperLimits,1)', ...
    fliplr(ModelRates(2).ProlTregs(:,LowerLimits,1)')];
fill(tx2, inBetweenWT, colorWT, 'FaceAlpha',FaceAlpha);
%KO fill
inBetweenKO = [ModelRates(2).ProlTregs(:,UpperLimits,2)', ...
    fliplr(ModelRates(2).ProlTregs(:,LowerLimits,2)')];
fill(tx2, inBetweenKO, colorKO, 'FaceAlpha',FaceAlpha);

title('Treg Death Rates')
ylabel('Cell Numbers')
xlabel('Age in Hours')
legend('WT', 'IL-2 KO', 'WT', 'IL-2 KO')
hold off

%%%%%%%%%%%%%
%Zoomed in Prol Tregs
%%%%%%%%%%%%%
subplot(2,4,7)

plot(tx, ModelRates(2).ProlTregs(:,1,1),'b-', 'LineWidth', LineWidth)%Plotting the means
hold on
plot(tx, ModelRates(2).ProlTregs(:,1,2),'b--', 'LineWidth', LineWidth)%Plotting the means

%WT fill
inBetweenWT = [ModelRates(2).ProlTregs(:,UpperLimits,1)', ...
    fliplr(ModelRates(2).ProlTregs(:,LowerLimits,1)')];
fill(tx2, inBetweenWT, colorWT, 'FaceAlpha',FaceAlpha);
%KO fill
inBetweenKO = [ModelRates(2).ProlTregs(:,UpperLimits,2)', ...
    fliplr(ModelRates(2).ProlTregs(:,LowerLimits,2)')];
fill(tx2, inBetweenKO, colorKO, 'FaceAlpha',FaceAlpha);

title('Zoomed in Proliferating Tregs')
ylabel('Cell Numbers')
xlabel('Age in Hours')
legend('WT', 'IL-2 KO', 'WT', 'IL-2 KO')
xlim([0 60])
ylim([1.5 7])
hold off


%%%%%%%%%%%%
%Treg Hill Results
%%%%%%%%%%%%
subplot(2,4,5)

plot(tx, ModelRates(2).TregDeathHill(:,1,1),'b-', 'LineWidth', LineWidth)%Plotting the means
hold on
plot(tx, ModelRates(2).TregDeathHill(:,1,2),'b--', 'LineWidth', LineWidth)%Plotting the means

%WT fill
inBetweenWT = [ModelRates(2).TregDeathHill(:,UpperLimits,1)', ...
    fliplr(ModelRates(2).TregDeathHill(:,LowerLimits,1)')];
fill(tx2, inBetweenWT, colorWT, 'FaceAlpha',FaceAlpha);
%KO fill
inBetweenKO = [ModelRates(2).TregDeathHill(:,UpperLimits,2)', ...
    fliplr(ModelRates(2).TregDeathHill(:,LowerLimits,2)')];
fill(tx2, inBetweenKO, colorKO, 'FaceAlpha',FaceAlpha);

title('Treg Death Capacity')
ylabel('Cell Numbers')
xlabel('Age in Hours')

hold off

%%%%%%%%%%%%%%
%Thymic Tregs Results
%%%%%%%%%%%%%%
subplot(2,4,5)

plot(tx, ModelRates(2).TregDeathHill(:,1,1),'b-', 'LineWidth', LineWidth)%Plotting the means
hold on
plot(tx, ModelRates(2).TregDeathHill(:,1,2),'b--', 'LineWidth', LineWidth)%Plotting the means

%WT fill
inBetweenWT = [ModelRates(2).TregDeathHill(:,UpperLimits,1)', ...
    fliplr(ModelRates(2).TregDeathHill(:,LowerLimits,1)')];
fill(tx2, inBetweenWT, colorWT, 'FaceAlpha',FaceAlpha);
%KO fill
inBetweenKO = [ModelRates(2).TregDeathHill(:,UpperLimits,2)', ...
    fliplr(ModelRates(2).TregDeathHill(:,LowerLimits,2)')];
fill(tx2, inBetweenKO, colorKO, 'FaceAlpha',FaceAlpha);

title('Treg Death Capacity')
ylabel('Cell Numbers')
xlabel('Age in Hours')

hold off



%-----------------Plotting from StatsOfCells DF------------------------%
% Here are what each of the multidimensional slots represents
% (A, B, C, D)
% A - Each row represents an hour
% B - The different stats are here
% C - The cellular populations
% D - Genotypes

% Note: I didn't learn how to properly use data structs at the time when I
% was coding the function that calculates StatsOfCells.If I have time, I'll change
% the structure.

%%%%%
% IL-2
%%%%%
subplot(2,4,1)

plot(tx, StatsOfCells(:,1,11,1),'b-', 'LineWidth', LineWidth)%Plotting the means WT
hold on
plot(tx, StatsOfCells(:,1,11,2),'b--', 'LineWidth', LineWidth)%Plotting the means KO
%WT Fill
inBetweenWT = [StatsOfCells(:,UpperLimits,11,1)', fliplr(StatsOfCells(:,LowerLimits,11,1)')];
fill(tx2, inBetweenWT , colorWT, 'FaceAlpha',FaceAlpha);
%KO Fill
inBetweenKO = [StatsOfCells(:,UpperLimits,11,2)', fliplr(StatsOfCells(:,LowerLimits,11,2)')]; % KO
fill(tx2, inBetweenKO , colorKO, 'FaceAlpha',FaceAlpha);

title('IL-2 Cytokine')
ylabel('Cell Numbers')
xlabel('Age in Hours')
hold off

%%%%%%%%
% Total Tregs
%%%%%%%%
subplot(2,4,4)

plot(tx, StatsOfCells(:,1,3,1),'b-', 'LineWidth', LineWidth)%Plotting the means
hold on
plot(tx, StatsOfCells(:,1,3,2),'b--', 'LineWidth', LineWidth)%Plotting the means

%WT Fill
inBetweenWT = [StatsOfCells(:,UpperLimits,3,1)', fliplr(StatsOfCells(:,LowerLimits,3,1)')];
fill(tx2, inBetweenWT , colorWT, 'FaceAlpha',FaceAlpha);

%KO FIll
inBetweenKO = [StatsOfCells(:,UpperLimits,3,2)', fliplr(StatsOfCells(:,LowerLimits,3,2)')];
fill(tx2, inBetweenKO , colorKO, 'FaceAlpha',FaceAlpha);

title('Total Tregs')
ylabel('Cell Numbers')
xlabel('Age in Hours')
hold off

%%%%%%%%%%%%%%
% zoomed in Total Tregs
%%%%%%%%%%%%%%
subplot(2,4,8)

plot(tx, StatsOfCells(:,1,3,1),'b-', 'LineWidth', LineWidth)%Plotting the means
hold on
plot(tx, StatsOfCells(:,1,3,2),'b--', 'LineWidth', LineWidth)%Plotting the means

%WT Fill
inBetweenWT = [StatsOfCells(:,UpperLimits,3,1)', fliplr(StatsOfCells(:,LowerLimits,3,1)')];
fill(tx2, inBetweenWT , colorWT, 'FaceAlpha',FaceAlpha);

%KO FIll
inBetweenKO = [StatsOfCells(:,UpperLimits,3,2)', fliplr(StatsOfCells(:,LowerLimits,3,2)')];
fill(tx2, inBetweenKO , colorKO, 'FaceAlpha',FaceAlpha);

title('Zoomed in Total Tregs')
ylabel('Cell Numbers')
xlabel('Age in Hours')
xlim([0 60])
ylim([130 400])
hold off




saveas(PLT1, sprintf(plt1))






























































end