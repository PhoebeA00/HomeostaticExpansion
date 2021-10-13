function PlottingLHSResults(StatsOfCells, tx)

%{
Figuring out the order of my Celullar pages for plotting
(3,4, n)
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
Genotype = [1,2];
for gene = Genotype
    %Selecting the proper files and file names to save figure
    %1 = WildType, 2 = Genotype
    
    if gene == 1
        CellData = readtable('../../RawData/ActivatedWTSpleen.csv');
        ProlData = readtable('../../RawData/WTProl.csv');
        Gntype = "WT";

    elseif gene == 2
        CellData = readtable('../../RawData/ActivatedKOSpleen.csv');
        ProlData = readtable('../../RawData/KOProl.csv');
        Gntype = "KO";
    end

    %Setting up data for plotting
    CellData = CellData(:,{'NaiveCT', 'ActivatedCD4CT', 'X4TregCT', ...
        'ThymicNaive', 'ActivatedNaiveCT', ...
        'ThymicDerivedTregsCT', 'NaiveDerivedTregsCT' ... 
        'hours'});

    ProlData = ProlData(:,{ 'NaiveProlCT', 'ActivatedProlCT', 'X4TregProlCT', ...
        'hours'});
    
    %Setting up variables for saving the figure    
    loc = './Plots/LHS_InitialConditions/';
    plt = append(loc, Gntype, '_LHSinit.png');
    
    %Figure Positioning
    left = 0;
    bottom = 400;
    width = 1800;
    height = 1050;
    
    PLT = figure('Visible', 'off');
    set(PLT,'Position',[left bottom width height]); %This sents the position of the figure itself
    
    %----------------------------------------------------------------------------------------------%
    %                                           Naive T Cells
    %----------------------------------------------------------------------------------------------%
    
    
    %Total Naive T Cells
    subplot(3,4,1)
    scatter(CellData.hours, CellData.NaiveCT, 'blue')
    ylim([0, 4000000])
    hold on
    plot(tx, StatsOfCells(:,1,1,gene),'b-', 'LineWidth', 1.5)%Plotting the means
    tx2 = [tx, fliplr(tx)];
    inBetween = [StatsOfCells(:,3,1,gene)', fliplr(StatsOfCells(:,4,1,gene)')];
    fill(tx2, inBetween, [0.6 0.1 1.0], 'FaceAlpha',0.2);
    title('Naive T Cells')
    ylabel('Cell Numbers')
    xlabel('Age in Hours')
    hold off
    
    %Proliferating Naive
    subplot(3,4,2)
    scatter(ProlData.hours, ProlData.NaiveProlCT, 'blue')
    ylim([0, 4000000])
    hold on
    plot(tx, StatsOfCells(:,1,8,gene),'b-', 'LineWidth', 1.5)%Plotting the means
    tx2 = [tx, fliplr(tx)];
    inBetween = [StatsOfCells(:,3,8,gene)', fliplr(StatsOfCells(:,4,8,gene)')];
    fill(tx2, inBetween, [0.6 0.1 1.0], 'FaceAlpha',0.2);
    title('Prol Naive T Cells')
    ylabel('Cell Numbers')
    xlabel('Age in Hours')
    hold off
    
    %Thymic Derived Naive T Cells
    subplot(3,4,3)
    scatter(CellData.hours, CellData.ThymicNaive, 'blue')
    ylim([0, 4000000])
    hold on
    plot(tx, StatsOfCells(:,1,4,gene),'b-', 'LineWidth', 1.5)%Plotting the means
    tx2 = [tx, fliplr(tx)];
    inBetween = [StatsOfCells(:,3,4,gene)', fliplr(StatsOfCells(:,4,4,gene)')];
    fill(tx2, inBetween, [0.6 0.1 1.0], 'FaceAlpha',0.2);
    title('Thymic Derived Naive')
    ylabel('Cell Numbers')
    xlabel('Age in Hours')
    hold off
    
    %----------------------------------------------------------------------------------------------%
    %                                       Activated T Cells
    %----------------------------------------------------------------------------------------------%
    
    %Total Activated T Cells
    subplot(3,4,5)
    scatter(CellData.hours, CellData.ActivatedCD4CT, 'blue')
    ylim([0, 8000000])
    hold on
    plot(tx, StatsOfCells(:,1,2,gene),'b-', 'LineWidth', 1.5)%Plotting the means
    tx2 = [tx, fliplr(tx)];
    inBetween = [StatsOfCells(:,3,2,gene)', fliplr(StatsOfCells(:,4,2,gene)')];
    fill(tx2, inBetween, [0.6 0.1 1.0], 'FaceAlpha',0.2);
    title('Total Activated T')
    ylabel('Cell Numbers')
    xlabel('Age in Hours')
    hold off
    
    %Prol ActT
    subplot(3,4,6)
    scatter(ProlData.hours, ProlData.ActivatedProlCT, 'blue')
    ylim([0, 8000000])
    hold on
    plot(tx, StatsOfCells(:,1,9,gene),'b-', 'LineWidth', 1.5)%Plotting the means
    tx2 = [tx, fliplr(tx)];
    inBetween = [StatsOfCells(:,3,9,gene)', fliplr(StatsOfCells(:,4,9,gene)')];
    fill(tx2, inBetween, [0.6 0.1 1.0], 'FaceAlpha',0.2);
    title('Prol ActT')
    ylabel('Cell Numbers')
    xlabel('Age in Hours')
    hold off
    
    % Naive Derived Activated T
    subplot(3,4,7)
    scatter(CellData.hours, CellData.ActivatedNaiveCT, 'blue')
    ylim([0, 8000000])
    hold on
    plot(tx, StatsOfCells(:,1,5,gene),'b-', 'LineWidth', 1.5)%Plotting the means
    tx2 = [tx, fliplr(tx)];
    inBetween = [StatsOfCells(:,3,5,gene)', fliplr(StatsOfCells(:,4,5,gene)')];
    fill(tx2, inBetween, [0.6 0.1 1.0], 'FaceAlpha',0.2);
    title('Naive Derived ActT')
    ylabel('Cell Numbers')
    xlabel('Age in Hours')
    hold off
    
    % IL-2
    subplot(3,4,8)
    plot(tx, StatsOfCells(:,1,11,gene),'b-', 'LineWidth', 1.5)%Plotting the means
    tx2 = [tx, fliplr(tx)];
    inBetween = [StatsOfCells(:,3,11,gene)', fliplr(StatsOfCells(:,4,11,gene)')];
    fill(tx2, inBetween, [0.6 0.1 1.0], 'FaceAlpha',0.2);
    title('IL-2 Cytokine')
    ylabel('Cell Numbers')
    xlabel('Age in Hours')
    hold off
    
    %----------------------------------------------------------------------------------------------%
    %                                       Regulatory T Cells
    %----------------------------------------------------------------------------------------------%
    
    % Total Tregs
    subplot(3,4,9)
    scatter(CellData.hours, CellData.X4TregCT, 'blue')
    ylim([0, 820000])
    hold on
    plot(tx, StatsOfCells(:,1,3,gene),'b-', 'LineWidth', 1.5)%Plotting the means
    tx2 = [tx, fliplr(tx)];
    inBetween = [StatsOfCells(:,3,3,gene)', fliplr(StatsOfCells(:,4,3,gene)')];
    fill(tx2, inBetween, [0.6 0.1 1.0], 'FaceAlpha',0.2);
    title('Total Tregs')
    ylabel('Cell Numbers')
    xlabel('Age in Hours')
    hold off
    
    %Prol Treg
    subplot(3,4,10)
    scatter(ProlData.hours, ProlData.X4TregProlCT, 'blue')
    ylim([0, 820000])
    hold on
    plot(tx, StatsOfCells(:,1,3,gene),'b-', 'LineWidth', 1.5)%Plotting the means
    tx2 = [tx, fliplr(tx)];
    inBetween = [StatsOfCells(:,3,3,gene)', fliplr(StatsOfCells(:,4,3,gene)')];
    fill(tx2, inBetween, [0.6 0.1 1.0], 'FaceAlpha',0.2);
    title('Prol Tregs')
    ylabel('Cell Numbers')
    xlabel('Age in Hours')
    hold off
    
    
    % ThymicDerivedTregsCT
    subplot(3,4,11)
    scatter(CellData.hours, CellData.ThymicDerivedTregsCT, 'blue')
    ylim([0, 820000])
    hold on
    plot(tx, StatsOfCells(:,1,6,gene),'b-', 'LineWidth', 1.5)%Plotting the means
    tx2 = [tx, fliplr(tx)];
    inBetween = [StatsOfCells(:,3,6,gene)', fliplr(StatsOfCells(:,4,6,gene)')];
    fill(tx2, inBetween, [0.6 0.1 1.0], 'FaceAlpha',0.2);
    title('Thymic Derived Tregs')
    ylabel('Cell Numbers')
    xlabel('Age in Hours')
    hold off
    
    % ThymicDerivedTregsCT
    subplot(3,4,12)
    scatter(CellData.hours, CellData.NaiveDerivedTregsCT, 'blue')
    ylim([0, 820000])
    hold on
    plot(tx, StatsOfCells(:,1,7,gene),'b-', 'LineWidth', 1.5)%Plotting the means
    tx2 = [tx, fliplr(tx)];
    inBetween = [StatsOfCells(:,3,7,gene)', fliplr(StatsOfCells(:,4,7,gene)')];
    fill(tx2, inBetween, [0.6 0.1 1.0], 'FaceAlpha',0.2);
    title('Naive Derived Tregs CT')
    ylabel('Cell Numbers')
    xlabel('Age in Hours')
    hold off
    
    saveas(PLT, sprintf(plt))
end