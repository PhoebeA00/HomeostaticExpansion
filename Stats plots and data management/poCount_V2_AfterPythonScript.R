#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#Before continuing                                                                                             #
#Run this in bash:                                                                                             #
#~/my.work/PhD/HomestaticExpansionProject/Code/Stats\ plots\ and\ data\ management/CalculatingActivatedTCells.py
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

##############################
#                            #
# Prepping data for Modeling #
#                            #
##############################

ActivationData = read.csv("~/my.work/PhD/HomestaticExpansionProject/ModelData/ActivatedCD4pop.csv")

#Creating the naive T cells
ActivationData$NaiveCT = ActivationData$NoTregCD4CT - ActivationData$ActivatedCD4CT

#Deprecated, but keeping here just in case something comes up because I made this mistake
#ActivationData$AllTregs = ActivationData$NaiveDerivedTregsCT + ActivationData$X4TregFromThymusCT

#Saving the data to Modeldata for plotting
write.csv(ActivationData, "~/my.work/PhD/HomestaticExpansionProject/ModelData/ActivatedCD4pop.csv")

#Preparing data for fitting
ActivationData$NaiveDerivedTregsCT = ActivationData$NaiveDerivedTregsCT * 10**6
ActivationData$ThymicDerivedTregsCT = ActivationData$ThymicDerivedTregsCT * 10**6
ActivationData$NaiveCT = ActivationData$NaiveCT * 10**6
ActivationData$ActivatedCD4CT = ActivationData$ActivatedCD4CT * 10**6
ActivationData$X4TregCT = ActivationData$X4TregCT * 10**6
ActivationData$hours = ActivationData$intage * 24
ActivationData$hours[ActivationData$hours == 0] <- 1 #the 0 won't work with my matlab code

#Removing Thymus for the fitting process
ActivatedWTSpleen = subset(ActivationData, Genotype == 'WT' & Organ == 'Spleen')
ActivaedKOSpleen = subset(ActivationData, Genotype == 'KO' & Organ == 'Spleen')
#SavingData for modeling
write.csv(ActivatedWTSpleen, '~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedWTSpleen.csv')
write.csv(ActivaedKOSpleen, '~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedKOSpleen.csv')
