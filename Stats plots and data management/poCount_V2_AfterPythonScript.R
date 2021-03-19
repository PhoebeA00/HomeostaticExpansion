#Run this in bash:
#~/my.work/PhD/HomestaticExpansionProject/Code/Stats\ plots\ and\ data\ management/CalculatingActivatedTCells.py



#install.packages("reticulate")
# library(reticulate)
# 
# import("pandas")
# py_run_file('~/my.work/PhD/HomestaticExpansionProject/Code/Stats plots and data management/CalculatingActivatedTCells.py')


################################################
#
#       Here is a slot for sourcing the LookingForCD4Significance script that makes
#       ~/my.work/PhD/HomestaticExpansionProject/ModelData/ActivatedCD4pop.csv
# import("pandas")
# py_run_file('~/my.work/PhD/HomestaticExpansionProject/Code/Stats plots and data management/CalculatingActivatedTCells.py')
#
###################################################

#ActivationData = read.csv("~/my.work/PhD/HomestaticExpansionProject/ModelData/ActivatedCD4pop.csv")

#Making Naive 
#ActivationData$NaiveTCT = ActivationData$CD4CT - (ActivationData$X4TregCT + ActivationData$ActivatedCD4CT)
#write.csv(ActivationData, "~/my.work/PhD/HomestaticExpansionProject/ModelData/ActivatedCD4pop.csv")


##############################
#                            #
# Prepping data for Modeling #
#                            #
##############################

ActivationData = read.csv("~/my.work/PhD/HomestaticExpansionProject/ModelData/ActivatedCD4pop.csv")

#Creating the naive T cells all Tregs
ActivationData$NaiveCT = ActivationData$NoTregCD4CT - ActivationData$ActivatedCD4CT
ActivationData$AllTregs = ActivationData$NaiveDerivedTregsCT + ActivationData$X4TregFromThymusCT

#Saving the data to Modeldata for plotting
write.csv(ActivationData, "~/my.work/PhD/HomestaticExpansionProject/ModelData/ActivatedCD4pop.csv")

#Preparing data for fitting
ActivationData$NaiveDerivedTregsCT = ActivationData$NaiveDerivedTregsCT * 10**6
ActivationData$X4TregFromThymusCT = ActivationData$X4TregFromThymusCT * 10**6
ActivationData$NaiveCT = ActivationData$NaiveCT * 10**6
ActivationData$ActivatedCD4CT = ActivationData$ActivatedCD4CT * 10**6
ActivationData$AllTregs = ActivationData$AllTregs * 10**6
ActivationData$hours = ActivationData$intage * 24
ActivationData$hours[ActivationData$hours == 0] <- 1 #the 0 won't work with my matlab code

#Removing WT and Thymus for the fitting process
ActivatedWTSpleen = subset(ActivationData, Genotype == 'WT' & Organ == 'Spleen')
ActivaedKOSpleen = subset(ActivationData, Genotype == 'KO' & Organ == 'Spleen')
#SavingData for modeling
write.csv(ActivatedWTSpleen, '~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedWTSpleen.csv')
write.csv(ActivaedKOSpleen, '~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedKOSpleen.csv')
