function PlottingEverything(pOpt)
global tx

%{
FileName = '../Data/ParameterSetAll.csv';
p = GetParameters(EntryNumber, FileName);

%This is to get the parameter values on the figure
mu= p(1);
z = p(2);
g = p(3);
alpha = p(4);
c = p(5);
epsilon = p(6); 
b_R = p(7);
beta =p(8);
a = p(9);
b_T = p(10);
e_T = p(11);
e_R = p(12);
kA = p(13);
j = p(14);
kB = p(15);
n = p(16);
n1 = 1;
f = 1.38629;
%}

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
n1 = 1;
f = 1.38629;

%Different Il-2 production Values
dWT = pOpt(17);
dKO = 500;

%Thymus Parameters
K = 0.074896;
lambda = 0.016932;

%Setting up the right order of things
%p0 = [alpha, a, kA, e_T, e_R, g, b_T, b_R, epsilon, mu, beta, c, kB, j, z, n, dWT];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------Organizing Data---------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ModelDataWT = SimulateGrowth(pOpt, 1);
ModelDataKO = SimulateGrowth(pOpt, 2);
ModelDataWT = array2table(ModelDataWT);
ModelDataKO = array2table(ModelDataKO);

ModelDataWT.Properties.VariableNames = {'NaiveCT' 'ActivCT' 'TregCT' ...
    'ThyDerivedNaive' 'ActivNaive' 'ThyDerivedTregs' 'NaiveDerivedTregs' ...
    'ProlNaive' 'ProlActiv' 'ProlTregs' ...
    'Il2' 'ThyWeight'};
ModelDataKO.Properties.VariableNames = {'NaiveCT' 'ActivCT' 'TregCT' ...
    'ThyDerivedNaive' 'ActivNaive' 'ThyDerivedTregs' 'NaiveDerivedTregs' ...
    'ProlNaive' 'ProlActiv' 'ProlTregs' ...
    'Il2' 'ThyWeight'};

%Selecting the proper files and file names to save figure
%1 = WildType, 2 = Genotype

CellDataWT = readtable('../RawData/ActivatedWTSpleen.csv');
ProlDataWT = readtable('../RawData/WTProl.csv');
  
CellDataKO = readtable('../RawData/ActivatedKOSpleen.csv');
ProlDataKO = readtable('../RawData/KOProl.csv');

%Setting up data for plotting
CellDataWT = CellDataWT(:,{'NaiveCT', 'ActivatedCD4CT', 'X4TregCT', ...
    'ThymicNaive', 'ActivatedNaiveCT', ...
    'ThymicDerivedTregsCT', 'NaiveDerivedTregsCT' ... 
    'hours'});
CellDataKO = CellDataKO(:,{'NaiveCT', 'ActivatedCD4CT', 'X4TregCT', ...
    'ThymicNaive', 'ActivatedNaiveCT', ...
    'ThymicDerivedTregsCT', 'NaiveDerivedTregsCT' ... 
    'hours'});

ProlDataWT = ProlDataWT(:,{ 'NaiveProlCT', 'ActivatedProlCT', 'X4TregProlCT', ...
    'hours'});
ProlDataKO = ProlDataKO(:,{ 'NaiveProlCT', 'ActivatedProlCT', 'X4TregProlCT', ...
    'hours'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------Calculating---------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-=-=-=-=-=-=-=-=-=-=-=%
%          Fig 2 Things
%-=-=-=-=-=-=-=-=-=-=-=%

%---Hill suppression naive---%
ModelDataWT.HillNaive = (1./(1+(ModelDataWT.TregCT./kA).^n));
ModelDataKO.HillNaive = (1./(1+(ModelDataKO.TregCT./kA).^n));

%Hill suppression Treg death rate
ModelDataWT.HillTregDeath = (1./(1+(ModelDataWT.Il2./kB).^n1));
ModelDataKO.HillTregDeath = (1./(1+(ModelDataKO.Il2./kB).^n1));

%How many Activated T's are being destroyed
ModelDataWT.ActiveDestruction = (j.*ModelDataWT.TregCT.*ModelDataWT.ActivCT);
ModelDataKO.ActiveDestruction = (j.*ModelDataKO.TregCT.*ModelDataKO.ActivCT);

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%            Converting Suppression
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%

%---Activation Suppression---%
ModelDataWT.ActivationSuppression = 1 - ModelDataWT.HillNaive;
ModelDataKO.ActivationSuppression = 1 - ModelDataKO.HillNaive;

%---Treg Death Rate Suppression---%
ModelDataWT.TregDeathSuppression = 1 - ModelDataWT.HillTregDeath;
ModelDataKO.TregDeathSuppression = 1 - ModelDataKO.HillTregDeath;

%-=-=-=-=-=-=-=-=-=-=-=%
%               Naive
%-=-=-=-=-=-=-=-=-=-=-=%

%---Gains---%
ModelDataWT.NaiveGain = mu.*(ModelDataWT.ThyWeight./K) + (z.*ModelDataWT.NaiveCT);
ModelDataKO.NaiveGain = mu.*(ModelDataKO.ThyWeight./K) + (z.*ModelDataKO.NaiveCT);

%---Loss---%
ModelDataWT.NaiveLoss =  ((beta.*ModelDataWT.NaiveCT).*(1./(1.+(ModelDataWT.TregCT./kA).^n))) + ...
    + (c.*ModelDataWT.NaiveCT) + (g.*ModelDataWT.NaiveCT); 
ModelDataKO.NaiveLoss =  ((beta.*ModelDataKO.NaiveCT).*(1./(1.+(ModelDataKO.TregCT./kA).^n))) + ...
    + (c.*ModelDataKO.NaiveCT) + (g.*ModelDataKO.NaiveCT); 

%---Thymic Production---%
ModelDataWT.ThymicProductionIR = mu.*(ModelDataWT.ThyWeight./K);
ModelDataKO.ThymicProductionIR = mu.*(ModelDataKO.ThyWeight./K);

%---Naive Self Replication---%
ModelDataWT.NaiveProlIR = (z.*ModelDataWT.NaiveCT);
ModelDataKO.NaiveProlIR = (z.*ModelDataKO.NaiveCT);

%---Activated Naive T Cells---%
ModelDataWT.ActivatedNaiveIR = ((beta.*ModelDataWT.NaiveCT).*(1./(1.+(ModelDataWT.TregCT./kA).^n)));
ModelDataKO.ActivatedNaiveIR = ((beta.*ModelDataKO.NaiveCT).*(1./(1.+(ModelDataKO.TregCT./kA).^n)));

%---Differentiating Naive to Tregs---%
ModelDataWT.NaiveToTregIR = c.*ModelDataWT.NaiveCT;
ModelDataKO.NaiveToTregIR = c.*ModelDataKO.NaiveCT;

%---Death Rate---%
ModelDataWT.NaiveDeathIR = g.*ModelDataWT.NaiveCT;
ModelDataKO.NaiveDeathIR = g.*ModelDataKO.NaiveCT;

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%            Regulatory T Cells
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%

%---Gains---%
ModelDataWT.TregGain = (alpha*(ModelDataWT.ThyWeight./K)) + (epsilon.*ModelDataWT.TregCT)...
    + c.*(ModelDataWT.NaiveCT);
ModelDataKO.TregGain = (alpha*(ModelDataKO.ThyWeight./K)) + (epsilon.*ModelDataKO.TregCT)...
    + c.*(ModelDataKO.NaiveCT);

%---Loss---%
ModelDataWT.TregLoss = b_R.*ModelDataWT.TregCT.*(1./(1+(ModelDataWT.Il2./kB).^n1));
ModelDataKO.TregLoss = b_R.*ModelDataKO.TregCT.*(1./(1+(ModelDataKO.Il2./kB).^n1));

%---Thymic Production---%
ModelDataWT.TregThymicProductionIR = (alpha*(ModelDataWT.ThyWeight./K));
ModelDataKO.TregThymicProductionIR = (alpha*(ModelDataKO.ThyWeight./K));

%---Treg Self Replication---%
ModelDataWT.TregProlIR = epsilon.*ModelDataWT.TregCT;
ModelDataKO.TregProlIR = epsilon.*ModelDataKO.TregCT;

%---Naive to Treg---%
ModelDataWT.NaiveToTregIR = c.*ModelDataWT.NaiveCT;
ModelDataKO.NaiveToTregIR = c.*ModelDataKO.NaiveCT;

%---Death Rate---%
ModelDataWT.TregLoss = b_R.*ModelDataWT.TregCT.*(1./(1+(ModelDataWT.Il2./kB).^n1));
ModelDataKO.TregLoss = b_R.*ModelDataKO.TregCT.*(1./(1+(ModelDataKO.Il2./kB).^n1));

%---Calculating Treg Frequency---%
TregFreqWT = ModelDataWT.TregCT ./(ModelDataWT.NaiveCT+ModelDataWT.ActivCT);
TregFreqKO = ModelDataKO.TregCT ./(ModelDataKO.NaiveCT+ModelDataKO.ActivCT);

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%            Activated T Cells
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%

%---Gains---%
ModelDataWT.TGains =  ((beta.*ModelDataWT.NaiveCT).*(1./(1.+(ModelDataWT.TregCT./kA).^n))) +...
    (a.*ModelDataWT.ActivCT);
ModelDataKO.TGains =  ((beta.*ModelDataKO.NaiveCT).*(1./(1.+(ModelDataKO.TregCT./kA).^n))) +...
    (a.*ModelDataKO.ActivCT);

%---Loss---%
ModelDataWT.TLoss = (j.*ModelDataWT.TregCT.*ModelDataWT.ActivCT)+(b_T.*ModelDataWT.ActivCT);
ModelDataKO.TLoss = (j.*ModelDataKO.TregCT.*ModelDataKO.ActivCT)+(b_T.*ModelDataKO.ActivCT);

%---Activation of Naive T Cells ---%
ModelDataWT.ActivatedNaiveIR = ((beta.*ModelDataWT.NaiveCT).*(1./(1.+(ModelDataWT.TregCT./kA).^n)));
ModelDataKO.ActivatedNaiveIR = ((beta.*ModelDataKO.NaiveCT).*(1./(1.+(ModelDataKO.TregCT./kA).^n)));

%---Self Replicating---%
ModelDataWT.TProlIR = a.*ModelDataWT.ActivCT;
ModelDataKO.TProlIR = a.*ModelDataKO.ActivCT;

%---Death---%
ModelDataWT.TDeathIR = b_T.*ModelDataWT.ActivCT;
ModelDataKO.TDeathIR = b_T.*ModelDataKO.ActivCT;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------Saving Plots---------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


loc = '../Plots/';

plt1 = append(loc, '3_', 'HillsComarison.png');
plt2 = append(loc, '4_', 'Comparison.png');
plt3 = append(loc, '5_', "NaiveBalance.png");
plt4 = append(loc, '6_', 'TregBalance.png');
plt5 = append(loc,'7_', 'ActTBalance.png');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------Figures--------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%               Figure Parameters
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%Figure Positioning
left = 0;
bottom = 400;
width = 1800;
height = 1050;

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%                 Figure 1 - Hills Comparison
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%

PLT1 = figure('Visible', 'off');
set(PLT1,'Position',[left bottom width height]); %This sents the position of the figure itself

%---Activation Suppression---%
subplot(3,3,1)
plot(tx, ModelDataWT.ActivationSuppression)
hold on
plot(tx, ModelDataKO.ActivationSuppression)
hold off
title('Percent of Activation Suppression')
ylabel('Hill Value')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Treg Death Suppression---%
subplot(3,3,2)
plot(tx, ModelDataWT.TregDeathSuppression)
hold on
plot(tx, ModelDataKO.TregDeathSuppression)
hold off
title('Percent of Treg Death Suppression')
ylabel('Death Rate Suppression')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Removed Activated T Cells---%
subplot(3,3,3)
plot(tx, ModelDataWT.ActiveDestruction)
hold on
plot(tx, ModelDataKO.ActiveDestruction)
hold off
title('Removed Activated T Cells')
legend ({'WT', 'KO'}, 'Location','northwest')


Changing = {'*mu',       mu,         '   cells*hr−1';...
                    '*z',             z,             'cells-1*hr-1';...
                    '*g',          g,             '   hr−1';...
                    
                    '*alpha',     alpha,     '   cells*hr−1';...
                    '*c',           c,           '   hr−1';...
                    '*epsilon',  epsilon,   '  hr−1';...
                    '*b_R',      b_R,          '   hr−1';...
                    
                    '*beta',      beta,      '   hr−1';...
                    '*a',           a,             '   hr−1';...
                    '*b_T',       b_T,           '   hr−1'};

columnname =   {'Parameters', '     Value          ', '           Units           '};
columnformat = {'char', 'numeric', 'char'}; 
uitable('Units','normalized',...
                 'Position', [0.12 0.2 0.3466 0.4],... % [ Horizontal Location, Verticle location,Right Line, Bottom Line]
                 'Data', Changing,...
                 'ColumnName', columnname,...
                 'ColumnFormat', columnformat,...
                 'RowName',[],...
                 'FontSize', 15,...
                 'ColumnWidth', {150 200 270});
             
             
Fixed =  {'*e_T',       e_T,         '   cells-1*hr−1';...
                '*e_R',       e_R,          '   cells-1*hr−1';...
                
                'kA',         kA,           '   cells';...
                'j',                j,             '    cell-2*hour-1';...
                'kB',          kB,           '   cells';...
                
                '*n',           n,           '              -        ';...
                '*d',           '1000/500',           '   Molecules*cells-1*hr−1';...
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
             
 saveas(PLT1, sprintf(plt1))

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%        Figure 2 - Comparing Models
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%

%Begin Plotting
PLT2 = figure('Visible', 'off');
set(PLT2,'Position',[left bottom width height]); %This sents the position of the figure itself

%---Naive CT---%
subplot(3,4,1)
%ylim([0, 4000000])
plot(tx, ModelDataWT.NaiveCT)
hold on 
plot(tx, ModelDataKO.NaiveCT)
title('Naive T Cells')
hold off
legend ({'WT', 'KO'}, 'Location','northwest')

%---Naive Prol---%
subplot(3,4,2)
plot(tx, ModelDataWT.ProlNaive)
hold on 
plot(tx, ModelDataKO.ProlNaive)
title('Proliferating Naive')
hold off
legend ({'WT', 'KO'}, 'Location','northwest')

%---Naive Thymic---%
subplot(3,4,3)
plot(tx, ModelDataWT.ThyDerivedNaive)
hold on 
plot(tx, ModelDataKO.ThyDerivedNaive)
title('Thymic Naive')
hold off
legend ({'WT', 'KO'}, 'Location','northwest')

%---Treg Frequency---%
subplot(3,4,4)
plot(tx, TregFreqWT)
hold on
plot(tx, TregFreqKO)
hold off
title('Treg frequency')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Activated CT---%
subplot(3,4,5)
plot(tx, ModelDataWT.ActivCT)
hold on 
plot(tx, ModelDataKO.ActivCT)
title('Activated Counts')
hold off
legend ({'WT', 'KO'}, 'Location','northwest')

%---Activated Prol---%
subplot(3,4,6)
plot(tx, ModelDataWT.ProlActiv)
hold on 
plot(tx, ModelDataKO.ProlActiv)
title('Activated Prol')
hold off
legend ({'WT', 'KO'}, 'Location','northwest')

%---Activated Naive derived---%
subplot(3,4,7)
plot(tx, ModelDataWT.ActivCT)
hold on 
plot(tx, ModelDataKO.ActivCT)
title('Activated Naive')
hold off
legend ({'WT', 'KO'}, 'Location','northwest')

%---IL-2---%
subplot(3,4,8)
plot(tx, ModelDataWT.Il2)
hold on
plot(tx, ModelDataKO.Il2)
title('IL-2')
hold off
legend ({'WT', 'KO'}, 'Location','northwest')

%---Treg CT---%
subplot(3,4,9)
plot(tx, ModelDataWT.TregCT)
hold on 
plot(tx, ModelDataKO.TregCT)
title('Treg CT')
hold off
legend ({'WT', 'KO'}, 'Location','northwest')

%---Treg Prol---%
subplot(3,4,10)
plot(tx, ModelDataWT.ProlTregs)
hold on 
plot(tx, ModelDataKO.ProlTregs)
title('Treg Prol')
hold off
legend ({'WT', 'KO'}, 'Location','northwest')

%---Thymic Tregs---%
subplot(3,4,11)
plot(tx, ModelDataWT.ThyDerivedTregs)
hold on 
plot(tx, ModelDataKO.ThyDerivedTregs)
title('Thymic Tregs')
hold off
legend ({'WT', 'KO'}, 'Location','northwest')

%---Naive Derived Tregs---%
subplot(3,4,12)
plot(tx, ModelDataWT.NaiveDerivedTregs)
hold on 
plot(tx, ModelDataKO.NaiveDerivedTregs)
title('Naive Derived Tregs')
hold off

saveas(PLT2, sprintf(plt2))

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%                Figure 3 - Naive
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%

PLT3 = figure('Visible', 'on');
set(PLT3,'Position',[left bottom width height]); %This sents the position of the figure itself

%---Naive Gain---%
subplot(3,4,2)
plot(tx, ModelDataWT.NaiveGain)
hold on
plot(tx, ModelDataKO.NaiveGain)
hold off
title('Naive Gain')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Naive Loss---%
subplot(3,4,3)
plot(tx, ModelDataWT.NaiveLoss)
hold on
plot(tx, ModelDataKO.NaiveLoss)
hold off
title('Naive Loss')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Thymic Production Rate---%
subplot(3,4,5)
plot(tx, ModelDataWT.ThymicProductionIR)
hold on
plot(tx, ModelDataKO.ThymicProductionIR)
hold off
title('Thymic Production Rate')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Naive Prolifeartion---%
subplot(3,4,6)
plot(tx, ModelDataWT.NaiveProlIR)
hold on
plot(tx, ModelDataKO.NaiveProlIR)
hold off
title('Naive Prolifeartion')
legend ({'WT', 'KO'}, 'Location','northwest')


%---Acitvation of Naive T Cells---%
subplot(3,4,7)
plot(tx, ModelDataWT.ActivatedNaiveIR)
hold on
plot(tx, ModelDataKO.ActivatedNaiveIR)
hold off
title('Acitvation of Naive T Cells')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Naive Differentiation to Tregs---%
subplot(3,4,8)
plot(tx, ModelDataWT.NaiveToTregIR)
hold on
plot(tx, ModelDataKO.NaiveToTregIR)
hold off
title('Differentiation to Tregs')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Naive Death Rate---%
subplot(3,4,9)
plot(tx, ModelDataWT.NaiveDeathIR)
hold on
plot(tx, ModelDataKO.NaiveDeathIR)
hold off
title('Dying Naive T Cells')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Percent of Suppression---%
subplot(3,4,12)
plot(tx, ModelDataWT.ActivationSuppression)
hold on
plot(tx, ModelDataKO.ActivationSuppression)
hold off
title('Percent of Suppression')
legend ({'WT', 'KO'}, 'Location','northwest')

saveas(PLT3, sprintf(plt3))


%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%             Figure 4 - Treg Balance
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
PLT4 = figure('Visible', 'on');
set(PLT4,'Position',[left bottom width height]); %This sents the position of the figure itself

%---Treg Gain---%
subplot(3,4,2)
plot(tx, ModelDataWT.TregGain)
hold on
plot(tx, ModelDataKO.TregGain)
hold off
title('Treg Gain')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Treg Loss---%
subplot(3,4,3)
plot(tx, ModelDataWT.TregLoss)
hold on
plot(tx, ModelDataKO.TregLoss)
hold off
title('Treg Loss')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Thymic Treg---%
subplot(3,4,5)
plot(tx, ModelDataWT.TregThymicProductionIR)
hold on
plot(tx, ModelDataKO.TregThymicProductionIR)
hold off
title('Thymic Treg')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Replicating Tregs---%
subplot(3,4,6)
plot(tx, ModelDataWT.TregProlIR)
hold on
plot(tx, ModelDataKO.TregProlIR)
hold off
title('Replicating Tregs')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Naive Derived Tregs---%
subplot(3,4,7)
plot(tx, ModelDataWT.NaiveToTregIR)
hold on
plot(tx, ModelDataKO.NaiveToTregIR)
hold off
title('Peripherally Derived Tregs')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Dying Tregs---%
subplot(3,4,8)
plot(tx, ModelDataWT.TregLoss)
hold on
plot(tx, ModelDataKO.TregLoss)
hold off
title('Dying Tregs')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Treg Death rate suppression---%
subplot(3,4,10)
plot(tx, ModelDataWT.TregDeathSuppression)
hold on
plot(tx, ModelDataKO.TregDeathSuppression)
hold off
title('Treg Death Suppression')
legend ({'WT', 'KO'}, 'Location','northwest')

%---Treg Frequency---%
subplot(3,4,11)
plot(tx, TregFreqWT)
hold on
plot(tx, TregFreqKO)
hold off
title('Treg Frequency')
legend ({'WT', 'KO'}, 'Location','northwest')

saveas(PLT4, sprintf(plt4))

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%             Figure 5 - Activation Balance
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
PLT5 = figure('Visible', 'off');
set(PLT5,'Position',[left bottom width height]); %This sents the position of the figure itself

%--AcT Gains---%
subplot(3,3,1)
plot(tx, ModelDataWT.TGains)
hold on
plot(tx, ModelDataKO.TGains)
hold off
title('AcT Gains')
legend ({'WT', 'KO'}, 'Location','northwest')

%--AcT Gains---%
subplot(3,3,2)
plot(tx, ModelDataWT.TLoss)
hold on
plot(tx, ModelDataKO.TLoss)
hold off
title('AcT Loss')
legend ({'WT', 'KO'}, 'Location','northwest')

%--Naive Derived AcT---%
subplot(3,3,4)
plot(tx, ModelDataWT.ActivatedNaiveIR)
hold on
plot(tx, ModelDataKO.ActivatedNaiveIR)
hold off
title('Naive Derived AcT')
legend ({'WT', 'KO'}, 'Location','northwest')

%--Proliferating AcT---%
subplot(3,3,5)
plot(tx, ModelDataWT.TProlIR)
hold on
plot(tx, ModelDataKO.TProlIR)
hold off
title('Proliferating AcT')
legend ({'WT', 'KO'}, 'Location','northwest')

%--Dying AcT---%
subplot(3,3,7)
plot(tx, ModelDataWT.TDeathIR)
hold on
plot(tx, ModelDataKO.TDeathIR)
hold off
title('Dying AcT')
legend ({'WT', 'KO'}, 'Location','northwest')

%--Removed AcT---%
subplot(3,3,8)
plot(tx, ModelDataWT.ActiveDestruction)
hold on
plot(tx, ModelDataKO.ActiveDestruction)
hold off
title('Removed AcT')
legend ({'WT', 'KO'}, 'Location','northwest')

saveas(PLT5, sprintf(plt5))

PLT6  = figure();
set(PLT6,'Position',[left bottom width height]); %This sents the position of the figure itself

%---Dying Tregs---%
plot(tx, ModelDataWT.TregLoss)
hold on
plot(tx, ModelDataKO.TregLoss)
hold off
title('Dying Tregs')
legend ({'WT', 'KO'}, 'Location','northwest')

end














