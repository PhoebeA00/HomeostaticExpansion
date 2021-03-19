#This script 
library(ggplot2)
require(plyr)
#kristen hogquist
#Ellenbnn robread.csv2()ey
#pc = read.csv('/media/jon/Seagate Expansion Drive/ModelData/Version8_Batchprocsesing.csv', header = T)
#pc = read.csv('/ModelData/Version8_Batchprocsesing.csv', header = T)
#pc = read.csv('/ModelData/Version8_Batchprocsesing_0-9.csv', header = T)
pop = read.csv('~/my.work/PhD/HomestaticExpansionProject/ModelData/Version8_Batchprocsesing.csv', 
               header = T, blank.lines.skip = TRUE)
####For my laptop#####
#setwd('/media/jon/Seagate Expansion Drive/ModelData/PlotsAll/')
##Remove extra data
#pop=pop[1:125,]

####For the lab computer#####
#setwd('/ModelData/PlotsAll/')
#### Not having the numbers as numeric is really screwing everything up
# pop[37:56] <- lapply(pop[37:56], as.numeric)
#pop[15:41] <- lapply(pop[15:41], as.numeric)
pop$Genotype = as.character(pop$Genotype)
pop$Organ = as.character(pop$Organ)
pop$intage = pop$Age
pop$Age = as.factor(pop$Age)#This made plotting categories possible
#pop$Age = as.character(pop$Age) #This has to be treated as a character in order to group the right rows properly
# pop[37:56] <- lapply(pop[37:56], factor)

#Need the ratios
pop$Bratio = pop$Bcells_events / pop$Singlets2_events
pop$CD4Ratio = pop$pCD4_events / pop$Singlets2_events
pop$CD8Ratio = pop$pCD8_events / pop$Singlets2_events
pop$DPRatio = pop$pDP_events / pop$Singlets2_events
pop$DNRatio = pop$DN_events / pop$Singlets2_events
pop$X4TregRatio = pop$x4Tregs_events / pop$CD4_events
pop$X8TregRatio = pop$x8Tregs_events / pop$CD8_events
pop$TCRbRatio = pop$TCRb_events / pop$Singlets2_events
pop$CD4ProlRatio = pop$CD4prol_events / pop$CD4_events
pop$CD8ProlRatio = pop$CD8Prol_events / pop$CD8_events
pop$BProlRatio = pop$BcellProl_events / pop$Bcells_events
pop$X4TregProlRatio = pop$x4tregProl_events / pop$x4Tregs_events
pop$X8TregProlRatio = pop$x8tregsProl / pop$x8Tregs_events

#The huge numbers don't behave well in ggplot
#pop$TotalLiveCountInMillions = pop$TotalLiveCountInMillions * 10**6
pop$Bct = pop$Bratio * pop$TotalLiveCountInMillions
pop$CD4CT = pop$CD4Ratio * pop$TotalLiveCountInMillions
pop$CD8ct = pop$CD8Ratio * pop$TotalLiveCountInMillions
pop$DPct = pop$DPRatio *pop$TotalLiveCountInMillions
pop$DNct = pop$DNRatio * pop$TotalLiveCountInMillions
pop$X4TregCT = pop$X4TregRatio * pop$CD4CT
pop$X8TregCT = pop$X8TregRatio * pop$CD8ct
pop$TCRbCT = pop$TCRbRatio * pop$TotalLiveCountInMillions
#The rest after here are based on the cell numbers count
###CD4's
pop$CD4ProlCT = pop$CD4ProlRatio * pop$CD4CT
pop$CD8ProlCT = pop$CD8ProlRatio * pop$CD8ct
pop$BprolCT = pop$BProlRatio * pop$Bct
pop$X4TregProlCT = pop$X4TregProlRatio * pop$X4TregCT
pop$X8TregProlCT = pop$X8TregProlRatio * pop$X8TregCT


####################################################
# Calculating the Tregs derived from Naive T cells #
####################################################

library(dplyr)
library("Rmisc")
library("reshape2")
#Need to remove dates that are not symmetrical for subtractions - Dates are 2/25/2018 and 12/6/2017
#There is a missing thymus and only one thymus for one of the dates - Incomplete data

popNoInc <- pop[!(pop$intage == 18 & pop$expDate=="2/25/2018" | pop$expDate=="12/6/2017"),]

#Works fine with day 0, I don't remember what the problem was
#day 0 screwing with the data
#popNoInc = popNoInc[!(popNoInc$Age == "0"),]

#A d56 spleen doesn't have a thymus partner
popNoInc = popNoInc[!(popNoInc$FileID == "JA022518WK8M1WTS"),]

#Grouping and then subtracting Thymus treg freq from Spleen

#The numbers don't work out too well with this format.
popNoInc = popNoInc %>%
  group_by(Age, Genotype, expDate) %>% #expDate is to subtract only by the thymus from the same experiment
  mutate(NaiveDerivedTregsRatio = X4TregRatio - X4TregRatio[Organ == 'Thymus'])
#Tregs derived from the Thymus
popNoInc$ThymusDerivedTregRatio = popNoInc$X4TregRatio - popNoInc$NaiveDerivedTregsRatio
#Replacing the Negative Values with 0
popNoInc$NaiveDerivedTregsRatio[popNoInc$NaiveDerivedTregsRatio < 0] <- 0

#Calculating the Counts for modeling
popNoInc$NaiveDerivedTregsCT = popNoInc$CD4CT * popNoInc$NaiveDerivedTregsRatio
popNoInc$X4TregFromThymusCT = popNoInc$CD4CT * popNoInc$ThymusDerivedTregRatio
popNoInc$NoTregCD4CT = popNoInc$CD4CT - (popNoInc$NaiveDerivedTregsCT + popNoInc$X4TregFromThymusCT)


# popNoInc saved into a csv below

################################################
#                                              #
#Prepping Genevieves Data                      #
#                                              #
################################################


#Fixing the CD69 data that Genevieve gave me
#empty strings need to be filled with NA
CD69df = read.csv("~/my.work/PhD/HomestaticExpansionProject/T cell Activation Summary_Jon.csv", header = FALSE, na.strings=c("","NA"))
#has a weird first row with one entry saying "CD69.Data.Summary", so I'm removing it, and the row with names
CD69df = CD69df[-(1:2),]
#Removing an empty column
CD69df[[15]] = NULL
# Using actual Column names now
colnames(CD69df) <- c("PerformedBy", "Date", "Notebook", "Page", "Mouse Tag",
                      "TissuesUsed", "TubeNumbers", "Genotype", "Age", "Notes", 
                      "CD4_pct", "CD8_pct", "CD4CD69_pct", "CD8CD69_pct",
                      "CD4CD44CD62L_pct", "CD8CD44CD62L_pct")

#Now removing the rows that have empty data under the CD4_pct column
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
CD69df = completeFun(CD69df, "CD4_pct")
# Numbers aren't being read as numeric because the first lines are messed up
cl4 = which( colnames(CD69df)=="CD4_pct" )
CD69df[cl4:ncol(CD69df)] = lapply(CD69df[cl4:ncol(CD69df)], as.character)
CD69df[cl4:ncol(CD69df)] = lapply(CD69df[cl4:ncol(CD69df)], as.numeric)

#Age seems to have a problem, it read age 4 as the highest when I made a plot
CD69df["Age"] = lapply(CD69df["Age"], as.character)
CD69df["Age"] = lapply(CD69df["Age"], as.numeric)


# Replacing the IL-2 HET with WT
CD69df$Genotype[CD69df$Genotype == "IL-2-HET"] = "WT"



####################
#
# Saving my Results
#
####################


write.csv(pop, '~/my.work/PhD/HomestaticExpansionProject/ModelData/AfterCalculations.csv')
write.csv(CD69df, '~/my.work/PhD/HomestaticExpansionProject/ModelData/CD69DataFromGen.csv')
write.csv(popNoInc, "~/my.work/PhD/HomestaticExpansionProject/ModelData/NaiveTregDifferentiation.csv")

#Run the python script: CalculatingActivatedTCells.py, then the poCount_V2_AfterPythonScript.R  