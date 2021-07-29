function PlottingResults(pOpt, Genotype)
global tx

%This is to get the parameter values on the figure
alpha = pOpt(1);
a = pOpt(2);
kA = pOpt(3);
e_T= pOpt(4);
e_R= pOpt(5);
g= pOpt(6);
b_T= pOpt(7);
b_R= pOpt(8);
epsilon = pOpt(9);
mu = pOpt(10);
beta = pOpt(11);
c = pOpt(12);
kB = pOpt(13);
j = pOpt(14);
z = pOpt(15);
n = pOpt(16);

if Genotype == 1
    d = pOpt(17);
elseif Genotype == 2
    d = 0;
end



ModelData = SimulateGrowth(pOpt, Genotype);
ModelData = array2table(ModelData);

ModelData.Properties.VariableNames = {'NaiveCT' 'ActivCT' 'TregCT' ...
    'ThyDerivedNaive' 'ActivNaive' 'ThyDerivedTregs' 'NaiveDerivedTregs' ...
    'ProlNaive' 'ProlActiv' 'ProlTregs' ...
    'Il2' 'ThyWeight'};

%Selecting the proper files and file names to save figure
%1 = WildType, 2 = Genotype
if Genotype == 1
    CellData = readtable('../RawData/ActivatedWTSpleen.csv');
    ProlData = readtable('../RawData/WTProl.csv');
    Gntype = "WT";
    Order = "1_";
    
elseif Genotype == 2
    CellData = readtable('../RawData/ActivatedKOSpleen.csv');
    ProlData = readtable('../RawData/KOProl.csv');
    Gntype = "KO";
    Order = "2_";
end

%Setting up data for plotting
CellData = CellData(:,{'NaiveCT', 'ActivatedCD4CT', 'X4TregCT', ...
    'ThymicNaive', 'ActivatedNaiveCT', ...
    'ThymicDerivedTregsCT', 'NaiveDerivedTregsCT' ... 
    'hours'});

ProlData = ProlData(:,{ 'NaiveProlCT', 'ActivatedProlCT', 'X4TregProlCT', ...
    'hours'});

%Setting file names
loc = '../Plots/';

plt =  append (loc, Order, Gntype, '_1_pops.png');
plt2 = append (loc, Order, Gntype, '_2_Hills.png');


%Figure Positioning
left = 0;
bottom = 400;
width = 1800;
height = 1050;

%Begin Plotting
PLT = figure('Visible', 'off');
set(PLT,'Position',[left bottom width height]); %This sents the position of the figure itself

%Naive CT
subplot(3,4,1)
scatter(CellData.hours, CellData.NaiveCT)
ylim([0, 4000000])
hold on 
plot(tx, ModelData.NaiveCT)
title('Naive T Cells')
hold off
%Naive Prol
subplot(3,4,2)
scatter(ProlData.hours, ProlData.NaiveProlCT)
ylim([0, 1000000])
hold on 
plot(tx, ModelData.ProlNaive)
title('Proliferating Naive')
hold off
%Naive Thymic
subplot(3,4,3)
scatter(CellData.hours, CellData.ThymicNaive)
ylim([0, 3500000])
hold on 
plot(tx, ModelData.ThyDerivedNaive)
title('Thymic Naive')
hold off

%Calculating Treg Frequency
TregFreq = ModelData.TregCT ./(ModelData.NaiveCT+ModelData.ActivCT);
%TregFreq (TregFreq > 0.1) = NaN;
subplot(3,4,4)
plot(tx, TregFreq)
title('Treg frequency')

%Activated CT
subplot(3,4,5)
scatter(CellData.hours, CellData.ActivatedCD4CT)
ylim([0, 8000000])
hold on 
plot(tx, ModelData.ActivCT)
title('Activated Counts')
hold off
%Activated Prol
subplot(3,4,6)
scatter(ProlData.hours, ProlData.ActivatedProlCT)
ylim([0, 2000000])
hold on 
plot(tx, ModelData.ProlActiv)
title('Activated Prol')
hold off
%Activated Naive derived
subplot(3,4,7)
scatter(CellData.hours, CellData.ActivatedNaiveCT)
ylim([0, 6000000])
hold on 
plot(tx, ModelData.ActivCT)
title('Activated Naive')
hold off
%IL-2
subplot(3,4,8)
plot(tx, ModelData.Il2)
title('IL-2')
hold off
%Treg CT
subplot(3,4,9)
scatter(CellData.hours, CellData.X4TregCT)
ylim([0, 700000])
hold on 
plot(tx, ModelData.TregCT)
title('Treg CT')
hold off
%Treg Prol
subplot(3,4,10)
scatter(ProlData.hours, ProlData.X4TregProlCT)
ylim([0, 600000])
hold on 
plot(tx, ModelData.ProlTregs)
title('Treg Prol')
hold off
%Thymic Tregs
subplot(3,4,11)
scatter(CellData.hours, CellData.ThymicDerivedTregsCT)
ylim([0, 35000])
hold on 
plot(tx, ModelData.ThyDerivedTregs)
title('Thymic Tregs')
hold off
%Naive Derived Tregs
subplot(3,4,12)
scatter(CellData.hours, CellData.NaiveDerivedTregsCT)
ylim([0, 400000])
hold on 
plot(tx, ModelData.NaiveDerivedTregs)
title('Naive Derived Tregs')
hold off

saveas(PLT, sprintf(plt))
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------Fixed Parameters---------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%d = 1000; %IL-2 production Rate
f = 1.38629; %IL-2 degradation Rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------Calculating---------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n1 = 1;

%Hill suppression naive
ModelData.HillNaive = (1./(1+(ModelData.TregCT./kA).^n));
%Hill suppression Treg death rate
ModelData.HillTregDeath = (1./(1+(ModelData.Il2./kB).^n1));
%How many Activated T's are being destroyed
ModelData.ActiveDestruction = (j.*ModelData.TregCT.*ModelData.ActivCT);


PLT2 = figure('Visible', 'off');
set(PLT2,'Position',[left bottom width height]); %This sents the position of the figure itself

%Hill suppression naive
subplot(3,3,1)
plot(tx, ModelData.HillNaive)
title('Hill Value')
ylabel('Hill Value')

subplot(3,3,2)
plot(tx, ModelData.HillTregDeath)
title('Treg Death Suppression')
ylabel('Death Rate Suppression')

subplot(3,3,3)
plot(tx, ModelData.ActiveDestruction)
title('Destroyed T Cells')


Changing = {'mu',       mu,         '   cells*hr−1';...
                    'z',             z,             'cells-1*hr-1';...
                    '*g',          g,             '   hr−1';...
                    
                    'alpha',     alpha,     '   cells*hr−1';...
                    'c',           c,           '   hr−1';...
                    'epsilon',  epsilon,   '  hr−1';...
                    'b_R',      b_R,          '   hr−1';...
                    
                    'beta',      beta,      '   hr−1';...
                    'a',           a,             '   hr−1';...
                    'b_T',       b_T,           '   hr−1'};

columnname =   {'Parameters', '     Value          ', '           Units           '};
columnformat = {'char', 'numeric', 'char'}; 
uitable('Units','normalized',...
                 'Position', [0.12 0.2 0.3466 0.4],... % [ Horizontal Location, Verticle location,Right Line, Bottom Line]
                 'Data', Changing,...
                 'ColumnName', columnname,...
                 'ColumnFormat', columnformat,...
                 'RowName',[],...
                 'FontSize', 15,...p
                 'ColumnWidth', {150 200 270});
             
             
Fixed =  {'*e_T',       e_T,         '   cells-1*hr−1';...
                '*e_R',       e_R,          '   cells-1*hr−1';...
                
                'kA',         kA,           '   cells';...
                '*j',                j,             '    cell-2*hour-1';...
                'kB',          kB,           '   cells';...
                
                '*n',           n,           '              -        ';...
                '*d',           d,           '   Molecules*cells-1*hr−1';...
                '*f',            f,           '   hr−1'};
            
columnname =   {'Parameters', '     Value          ', '           Units           '};
columnformat = {'char', 'numeric', 'char'}; 
uitable('Units','normalized',...
                 'Position', [0.57 0.2 0.394 0.4],... % [ Horizontal Location, Verticle location, Right Line, Bottom Line]
                 'Data', Fixed,...
                 'ColumnName', columnname,...
                 'ColumnFormat', columnformat,...
                 'RowName',[],...
                  'FontSize', 15,...
                 'ColumnWidth', {150 200 360});            

saveas(PLT2, sprintf(plt2))
%}
end