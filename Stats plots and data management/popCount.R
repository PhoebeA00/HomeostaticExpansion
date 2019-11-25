library(ggplot2)
require(plyr)
#kristen hogquist
#Ellenbnn robread.csv2()ey
pc = read.csv('/ModelData/FinalVersion6Norm.csv', header = T)
####For my laptop#####
#setwd('/media/jon/Seagate Expansion Drive/ModelData/PlotsAll/')

####For the lab computer#####
setwd('/ModelData/PlotsAll/')
pop = pc 
#### Not having the numbers as numeric is really screwing everything up
# pop[37:56] <- lapply(pop[37:56], as.numeric)
pop[15:52] <- lapply(pop[15:51], as.numeric)
pop$Genotype = as.character(pop$Genotype)
pop$Organ = as.character(pop$Organ)
pop$Age = as.factor(pop$Age)#This made plotting categories possible
#pop$Age = as.character(pop$Age) #This has to be treated as a character in order to group the right rows properly
# pop[37:56] <- lapply(pop[37:56], factor)

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


pop$CD4PnACT = pop$CD4PnARatio * pop$CD4CT
pop$CD4Neither = pop$CD4NeitherRatio * pop$CD4CT
pop$CD4ProlCT = pop$CD4ProlRatio * pop$CD4CT + pop$CD4PnARatio
pop$CD4ActivCT = pop$CD4ActivRatio * pop$CD4CT + pop$CD4PnACT
####### CD8's

pop$CD8PnACT = pop$CD8PnARatio  * pop$CD8ct
pop$CD8Neither = pop$CD8NeitherRatio * pop$CD8ct
pop$CD8ProlCT = pop$CD8ProlRatio * pop$CD8ct + pop$CD8PnACT
pop$CD8ActivCT = pop$CD8ActivRatio * pop$CD8ct + pop$CD8PnACT

####### Bcells

pop$BPnACT = pop$BPnARatio * pop$Bct
pop$BNeither = pop$BNeitherRatio * pop$Bct
pop$BprolCT = pop$BProlRatio * pop$Bct + pop$BPnACT
pop$BactivCT = pop$BActivRatio * pop$Bct + pop$BPnACT

####### 4Tregs

pop$X4TregProlCT = pop$X4TregProlRatio * pop$X4TregCT
pop$X4TregActivCT = pop$X4TregActivRatio * pop$X4TregCT
pop$X4TregPnACT = pop$X4TregPnARatio * pop$X4TregCT
pop$X4TregNeitherCT = pop$X4TregNeither * pop$X4TregCT
####### 8tregs
pop$X8TregProlCT = pop$X8TregProlRatio * pop$X8TregCT
pop$X8TregActivCT = pop$X8TregActivRatio * pop$X8TregCT
pop$X8TregPnACT = pop$X8TregPnARatio * pop$X8TregCT
pop$X8TregNeitherCT = pop$X8TregNeither * pop$X8TregCT
####### TCRb ADD THIS LATER
# pop[14:56] <- lapply(pop[14:56], as.numeric)
#pop[, 14:56][is.na(pop[,14:56])] <- 0
#d21 = subset(dat, Age == '21')
# colnames(pop)
# "CD4ActivCT"              
# "CD4PnACT"
# d18 = subset(pop, Age == 18)

#################3
###### Plots Version 2
###################

#function(popl, TotalCell, TotalRatio, fileName, PageTitle)
PlotPops(pop, "CD4CT", "CD4Ratio", "CD4 T cells", "10-CD4.pdf")
PlotPops(pop, "CD8ct", "CD8Ratio", "CD8 T cells", "12-CD8.pdf")
PlotPops(pop, "Bct", "Bratio", "B cells", "13-BCells.pdf")
PlotPops(pop, "X4TregCT", "X4TregRatio", "CD4 Tregs", "14-Tregs.pdf")
PlotPops(pop, "X8TregCT", "X8TregRatio", "CD8 T regs", "15-TTvregs.pdf")
PlotPops(pop, "TCRbCT", "TCRbRatio", "TCRb", "16-zzTcrb.pdf")
PlotPops(pop, "CD4ProlCT", "CD4ProlRatio", "CD4 Proliferation", "17-CD4Prol.pdf")
PlotPops(pop, "CD4ActivCT", "CD4ActivRatio", "CD4 Activation", "18-CD4Act.pdf")
PlotPops(pop, "CD4PnACT", "CD4PnARatio", "CD4 PnA", "19-CD4PnA.pdf")
PlotPops(pop, "CD4Neither", "CD4Neither", "CD4 Neither PnA", "20-CD4Neither.pdf")
PlotPops(pop, "CD8ProlCT", "CD8ProlRatio", "CD8 Proliferation", "21-CD8Prol.pdf")
PlotPops(pop, "CD8ActivCT", "CD8ActivRatio", "CD8 Activation", "22-CD8Act.pdf")
PlotPops(pop, "CD8PnACT", "CD8PnARatio", "CD8 PnA", "23-CD8PnA.pdf")
PlotPops(pop, "CD8Neither", "CD8NeitherRatio", "CD8 Neither PnA", "24-CD8Neither.pdf")
PlotPops(pop, "BprolCT", "BProlRatio", "B Proliferation", "25-BProl.pdf")
PlotPops(pop, "BactivCT", "BActivRatio", "B Activation", "26-BActiv.pdf")
PlotPops(pop, "BPnACT", "BPnARatio", "B PnA", "27-BPnA.pdf")
PlotPops(pop, "BNeither", "BNeither", "B Neither PnA", "28-BNeither.pdf" )
PlotPops(pop, "X4TregProlCT", "X4TregProlRatio", "CD4 Treg Proliferation", "29-X4tregProli.pdf")
PlotPops(pop, "X4TregActivCT", "X4TregActivRatio", "CD4 Treg Activation", "30-X4Activation.pdf")
PlotPops(pop, "X4TregPnACT", "X4TregPnARatio", "CD4 Treg PnA", "31-cd4pna.pdf")
PlotPops(pop, "X4TregNeitherCT" , "X4TregNeither", "CD4 Treg Neither", "32-4tregneither.pdf" )
PlotPops(pop, "X8TregProlCT", "X8TregProlRatio", "CD8 Treg Proliferation", "33-8tregprol.pdf")
PlotPops(pop, "X8TregActivCT", "X8TregActivRatio", "CD8 Treg Activation", "34-8tregactiv.pdf")
PlotPops(pop, "X8TregPnACT", "X8TregPnARatio", "CD8 Treg PnA", "35-8tregPna.pdf")
PlotPops(pop, "X8TregNeitherCT", "X8TregNeither", "CD8 Treg Neither", "36-8tregNeither" )
setwd('/ModelData/')
write.csv(pop, "AfterCalulations.csv")
