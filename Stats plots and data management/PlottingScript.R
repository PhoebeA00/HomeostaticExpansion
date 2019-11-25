#This script is to plot all of the calculations that were made in "popCount_V2.R
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



PlotPops(pop, "CD8ProlCT", "CD8ProlRatio", "CD8 Proliferation", "21-CD8Prol.pdf")


PlotPops(pop, "BprolCT", "BProlRatio", "B Proliferation", "25-BProl.pdf")

PlotPops(pop, "X4TregProlCT", "X4TregProlRatio", "CD4 Treg Proliferation", "29-X4tregProli.pdf")

PlotPops(pop, "X8TregProlCT", "X8TregProlRatio", "CD8 Treg Proliferation", "33-8tregprol.pdf")

setwd('/ModelData/')
write.csv(pop, "AfterCalulations.csv")
