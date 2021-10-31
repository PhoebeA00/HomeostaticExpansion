function Ratesresults = CalculatingRatesFromLHSSampling(ModelData)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------Calculating---------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-=-=-=-=-=-=-=-=-=-=-=%
%               Naive
%-=-=-=-=-=-=-=-=-=-=-=%

%---Gains---%
ModelDataWT.NaiveGain = (mu.*ModelDataWT.NaiveCT.*(1-(ModelDataWT.NaiveCT./nK))) + (z.*ModelDataWT.NaiveCT);
ModelDataKO.NaiveGain = (mu.*ModelDataKO.NaiveCT.*(1-(ModelDataKO.NaiveCT./nK))) + (z.*ModelDataKO.NaiveCT);

%---Loss---%
ModelDataWT.NaiveLoss =  ((beta.*ModelDataWT.NaiveCT).*(1./(1+((ModelDataWT.TregCT.*(ModelDataWT.Il2./(Ki + ModelDataWT.Il2)))./kA).^n))) + ...
    + (c.*ModelDataWT.NaiveCT) + (g.*ModelDataWT.NaiveCT); 
ModelDataKO.NaiveLoss =  ((beta.*ModelDataKO.NaiveCT).*(1./(1+((ModelDataKO.TregCT.*(ModelDataKO.Il2./(Ki + ModelDataKO.Il2)))./kA).^n))) + ...
    + (c.*ModelDataKO.NaiveCT) + (g.*ModelDataKO.NaiveCT); 

%---Thymic Production---%
ModelDataWT.ThymicProductionIR = (mu.*ModelDataWT.NaiveCT.*(1-(ModelDataWT.NaiveCT./nK)));
ModelDataKO.ThymicProductionIR = (mu.*ModelDataKO.NaiveCT.*(1-(ModelDataKO.NaiveCT./nK)));

%---Naive Self Replication---%
ModelDataWT.NaiveProlIR = (z.*ModelDataWT.NaiveCT);
ModelDataKO.NaiveProlIR = (z.*ModelDataKO.NaiveCT);

%---Activated Naive T Cells---%
ModelDataWT.ActivatedNaiveIR = ((beta.*ModelDataWT.NaiveCT).*(1./(1+((ModelDataWT.TregCT.*(ModelDataWT.Il2./(Ki + ModelDataWT.Il2)))./kA).^n)));
ModelDataKO.ActivatedNaiveIR = ((beta.*ModelDataKO.NaiveCT).*(1./(1+((ModelDataKO.TregCT.*(ModelDataKO.Il2./(Ki + ModelDataKO.Il2)))./kA).^n)));

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
ModelDataWT.TregGain = (alpha.*ModelDataWT.TregCT.*(1-(ModelDataWT.TregCT./rK))) ...
    + (epsilon.*ModelDataWT.TregCT) + c.*(ModelDataWT.NaiveCT);
ModelDataKO.TregGain = (alpha.*ModelDataKO.TregCT.*(1-(ModelDataKO.TregCT./rK)))...
    + (epsilon.*ModelDataKO.TregCT) + c.*(ModelDataKO.NaiveCT);

%---Loss---%
ModelDataWT.TregLoss = b_R.*ModelDataWT.TregCT.*(1./(1+(ModelDataWT.Il2./kB).^n1));
ModelDataKO.TregLoss = b_R.*ModelDataKO.TregCT.*(1./(1+(ModelDataKO.Il2./kB).^n1));

%---Thymic Production---%
ModelDataWT.TregThymicProductionIR = (alpha.*ModelDataWT.TregCT.*(1-(ModelDataWT.TregCT./rK)));
ModelDataKO.TregThymicProductionIR = (alpha.*ModelDataWT.TregCT.*(1-(ModelDataWT.TregCT./rK)));

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
ModelDataWT.TGains =  ((beta.*ModelDataWT.NaiveCT).*(1./(1+((ModelDataWT.TregCT.*(ModelDataWT.Il2./(Ki + ModelDataWT.Il2)))./kA).^n))) +...
    (a.*ModelDataWT.ActivCT);
ModelDataKO.TGains =  ((beta.*ModelDataKO.NaiveCT).*(1./(1+((ModelDataKO.TregCT.*(ModelDataKO.Il2./(Ki + ModelDataKO.Il2)))./kA).^n))) +...
    (a.*ModelDataKO.ActivCT);

%---Loss---%
ModelDataWT.TLoss = (j.*ModelDataWT.TregCT.*ModelDataWT.ActivCT)+(b_T.*ModelDataWT.ActivCT);
ModelDataKO.TLoss = (j.*ModelDataKO.TregCT.*ModelDataKO.ActivCT)+(b_T.*ModelDataKO.ActivCT);

%---Activation of Naive T Cells ---%
ModelDataWT.ActivatedNaiveIR = ((beta.*ModelDataWT.NaiveCT).*(1./(1+((ModelDataWT.TregCT.*(ModelDataWT.Il2./(Ki + ModelDataWT.Il2)))./kA).^n)));
ModelDataKO.ActivatedNaiveIR = ((beta.*ModelDataKO.NaiveCT).*(1./(1+((ModelDataKO.TregCT.*(ModelDataKO.Il2./(Ki + ModelDataKO.Il2)))./kA).^n)));

%---Self Replicating---%
ModelDataWT.TProlIR = a.*ModelDataWT.ActivCT;
ModelDataKO.TProlIR = a.*ModelDataKO.ActivCT;

%---Death---%
ModelDataWT.TDeathIR = b_T.*ModelDataWT.ActivCT;
ModelDataKO.TDeathIR = b_T.*ModelDataKO.ActivCT;

%-=-=-=-=-=-=-=-=-=-=-=%
%          Fig 2 Things
%-=-=-=-=-=-=-=-=-=-=-=%

%---Hill suppression naive---%
ModelDataWT.ActivationSuppression = (1./(1+((ModelDataWT.TregCT.*(ModelDataWT.Il2./(Ki + ModelDataWT.Il2)))./kA).^n));
ModelDataKO.ActivationSuppression = (1./(1+((ModelDataKO.TregCT.*(ModelDataKO.Il2./(Ki + ModelDataKO.Il2)))./kA).^n));

%Hill suppression Treg death rate
ModelDataWT.TregDeathSuppression = (1./(1+(ModelDataWT.Il2./kB).^n1));
ModelDataKO.TregDeathSuppression = (1./(1+(ModelDataKO.Il2./kB).^n1));

%How many Activated T's are being destroyed
ModelDataWT.ActiveDestruction = (j.*ModelDataWT.TregCT.*ModelDataWT.ActivCT) ...
    .* (ModelDataWT.Il2./(Ki + ModelDataWT.Il2));
ModelDataKO.ActiveDestruction = (j.*ModelDataKO.TregCT.*ModelDataKO.ActivCT) ...
    .* (ModelDataKO.Il2./(Ki + ModelDataKO.Il2));

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
%                   Interleukin - 2
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%

%---Gain---%
ModelDataWT.ILGain = d.*ModelDataWT.ActivCT;
ModelDataKO.ILGain = d.*ModelDataKO.ActivCT;

%---Loss---%
ModelDataWT.ILLoss = e_T.*ModelDataWT.Il2.*ModelDataWT.ActivCT...
    + e_R.*ModelDataWT.Il2.*ModelDataWT.TregCT + f.*ModelDataWT.Il2;
ModelDataKO.ILLoss = e_T.*ModelDataKO.Il2.*ModelDataKO.ActivCT...
    + e_R.*ModelDataKO.Il2.*ModelDataKO.TregCT + f.*ModelDataKO.Il2;

%---ActT Consumption---%
ModelDataWT.AcTConsumption = e_T.*ModelDataWT.Il2.*ModelDataWT.ActivCT;
ModelDataKO.AcTConsumption = e_T.*ModelDataKO.Il2.*ModelDataKO.ActivCT;

%---Treg Consumption
ModelDataWT.TregConsumption = e_R.*ModelDataWT.Il2.*ModelDataWT.TregCT;
ModelDataKO.TregConsumption = e_R.*ModelDataKO.Il2.*ModelDataKO.TregCT;

%---Death of IL-2---%
ModelDataWT.IL2Death = f.*ModelDataWT.Il2;
ModelDataKO.IL2Death = f.*ModelDataKO.Il2;

%---Activation Suppression Bosst---%
ModelDataWT.ActivationSuprBoost = (ModelDataWT.Il2./(Ki + ModelDataWT.Il2));
ModelDataKO.ActivationSuprBoost = (ModelDataKO.Il2./(Ki + ModelDataKO.Il2));

%---Deactivation Boost---%
ModelDataWT.DeactivationBoost = (ModelDataWT.Il2./(Kj + ModelDataWT.Il2));
ModelDataKO.DeactivationBoost = (ModelDataKO.Il2./(Kj + ModelDataKO.Il2));
end