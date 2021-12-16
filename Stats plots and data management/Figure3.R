library(ggplot2)
library(ggpubr)

ModeldataWT = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Stats plots and data management/ModelOutputWT.csv")
ModeldataKO = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Stats plots and data management/ModelOutputKO.csv")
WTData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedWTSpleen.csv")
ProlWTData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/WTProl.csv")
KOData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedKOSpleen.csv")
ProlKOData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/KOProl.csv")

WTData$hours = WTData$hours / 24
KOData$hours = KOData$hours / 24
ProlWTData$hours = ProlWTData$hours / 24
ProlKOData$hours = ProlKOData$hours / 24

colnames(ModeldataWT) = c("NaiveCT", "ActTCT", "TregCT", "ThyNaive", "ActTNaive", "ThyTregs",
                          "TregNaive", "ProlNaive", "ProlActT", "ProlTreg", "IL-2", "ThymWeigth")
colnames(ModeldataKO) = c("NaiveCT", "ActTCT", "TregCT", "ThyNaive", "ActTNaive", "ThyTregs",
                          "TregNaive", "ProlNaive", "ProlActT", "ProlTreg", "IL-2", "ThymWeigth")

ModeldataWT$time = 0:431
ModeldataKO$time = 0:431
ModeldataWT$time = ModeldataWT$time / 24
ModeldataKO$time = ModeldataKO$time / 24

###########################################
#  WT - Naive
##############################################

# NAIVE T CELLS

NaiveCTWT = ggplot(WTData, aes(x=hours, y=NaiveCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=NaiveCT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Total Naive T Cell Count(WT)", x = "Age in days", y = "Cell Counts")

#PROLIFERATING NAIVE
ProlNaiveWT = ggplot(ProlWTData, aes(x=hours, y=NaiveProlCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=ProlNaive), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Proliferating Naive T Cells(WT)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,4000000))

#THYMIC NAIVE
ThymicNaiveWT = ggplot(WTData, aes(x=hours, y=ThymicNaive)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=ThyNaive), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Thymic Derived Naive T Cells(WT)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,4000000))

#############################################
# IL-2 KO Naive
#############################################


# NAIVE T CELLS

NaiveCTKO = ggplot(KOData, aes(x=hours, y=NaiveCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=NaiveCT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Total Naive T Cell Count(KO)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,4000000))

#PROLIFERATING NAIVE
ProlNaiveKO = ggplot(ProlKOData, aes(x=hours, y=NaiveProlCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=ProlNaive), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Proliferating Naive T Cells(KO)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,4000000))

#THYMIC NAIVE
ThymicNaiveKO = ggplot(KOData, aes(x=hours, y=ThymicNaive)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=ThyNaive), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Thymic Derived Naive T Cells(KO)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,4000000))



#################################################
#################################################
########### B - Activated T Cells ###############
#################################################
#################################################

# Total Activated

ActTCD4CTWT = ggplot(WTData, aes(x=hours, y=ActivatedCD4CT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=ActTCT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Total Activated T Cell Count(WT)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))

#PROLIFERATING ACTIVATED T
ProlActTWT = ggplot(ProlWTData, aes(x=hours, y=ActivatedProlCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=ProlActT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Proliferating Activated T Cells(WT)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))

#Naive Derived ActT
NaiveActTWT = ggplot(WTData, aes(x=hours, y=ActivatedNaiveCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=ActTNaive), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Naive Derived Activated T Cells (WT)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))

#############################################
# IL-2 KO Naive
#############################################

# Total Activated

ActTCD4CTKO = ggplot(KOData, aes(x=hours, y=ActivatedCD4CT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=ActTCT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Total Activated T Cell Count (KO)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))

#PROLIFERATING ACTIVATED T
ProlActTKO = ggplot(ProlKOData, aes(x=hours, y=ActivatedProlCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=ProlActT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Proliferating Activated T Cells (KO)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))

#Naive Derived ActT
NaiveActTKO = ggplot(KOData, aes(x=hours, y=ActivatedNaiveCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=ActTNaive), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Naive Derived Activated T Cells (KO)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))


###########################################
#  WT - Acitve
##############################################

# Total Activated

ActTCD4CTWT = ggplot(WTData, aes(x=hours, y=ActivatedCD4CT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=ActTCT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Total Activated T Cell Count(WT)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))

#PROLIFERATING ACTIVATED T
ProlActTWT = ggplot(ProlWTData, aes(x=hours, y=ActivatedProlCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=ProlActT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Proliferating Activated T Cells(WT)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))

#Naive Derived ActT
NaiveActTWT = ggplot(WTData, aes(x=hours, y=ActivatedNaiveCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=ActTNaive), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Naive Derived Activated T Cells (WT)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))

#############################################
# IL-2 KO Naive
#############################################

# Total Activated

ActTCD4CTKO = ggplot(KOData, aes(x=hours, y=ActivatedCD4CT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=ActTCT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Total Activated T Cell Count (KO)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))

#PROLIFERATING ACTIVATED T
ProlActTKO = ggplot(ProlKOData, aes(x=hours, y=ActivatedProlCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=ProlActT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Proliferating Activated T Cells (KO)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))

#Naive Derived ActT
NaiveActTKO = ggplot(KOData, aes(x=hours, y=ActivatedNaiveCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=ActTNaive), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Naive Derived Activated T Cells (KO)", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,7075000))

#####################################
#####################################
########### C - Tregs ###############
#####################################
#####################################

###########################################
#  WT - Tregs
##############################################

#Total Tregs
TregCTWT = ggplot(WTData, aes(x=hours, y=X4TregCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=TregCT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Total Treg Counts", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,665000))

# Thymic Tregs

ThymicTregWT = ggplot(WTData, aes(x=hours, y=ThymicDerivedTregsCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=ThyTregs), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Thymic Derived Tregs", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,665000))

#Naive Derived Tregs
NaiveTregWT = ggplot(WTData, aes(x=hours, y=NaiveDerivedTregsCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=TregNaive), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Peripherally Derived Tregs", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,665000))

#Proliferating Tregs
ProlTregWT = ggplot(ProlWTData, aes(x=hours, y=X4TregProlCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataWT, aes(x = time, y=ProlTreg), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Proliferating Tregs", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,665000))

#############################################
# IL-2 KO Tregs
#############################################


#Total Tregs
TregCTKO = ggplot(KOData, aes(x=hours, y=X4TregCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=TregCT), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Total Treg Counts", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,665000))

# Thymic Tregs

ThymicTregKO = ggplot(KOData, aes(x=hours, y=ThymicDerivedTregsCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=ThyTregs), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Thymic Derived Tregs", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,665000))

#Naive Derived Tregs
NaiveTregKO = ggplot(KOData, aes(x=hours, y=NaiveDerivedTregsCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=TregNaive), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Peripherally Derived Tregs", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,665000))

#Proliferating Tregs
ProlTregKO = ggplot(ProlKOData, aes(x=hours, y=X4TregProlCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line", linetype="dotted")+
  geom_line(data = ModeldataKO, aes(x = time, y=ProlTreg), colour = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Proliferating Tregs", x = "Age in days", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,820000))



########################
########################
##### The ggsave #######
########################
########################
#---------------#
#------3A-------#
#---------------#
NaivePlots = ggarrange(NaiveCTWT, ProlNaiveWT, ThymicNaiveWT, NaiveCTKO, ProlNaiveKO, ThymicNaiveKO,
                       labels = c("A"),
                       ncol = 3, nrow = 2)

ggsave(file = "~/my.work/PhD/HomestaticExpansionProject/Figures/ForPaper/Figure3/Figure3A.pdf", NaivePlots,
       height = 8,
       width = 18)

#---------------#
#------3B-------#
#---------------#
ActTPlots = ggarrange(ActTCD4CTWT, ProlActTWT, NaiveActTWT, ActTCD4CTKO, ProlActTKO, NaiveActTKO,
                      labels = c("B"),
                      ncol = 3, nrow = 2)

ggsave(file = "~/my.work/PhD/HomestaticExpansionProject/Figures/ForPaper/Figure3/Figure3B.pdf", ActTPlots,
       height = 8,
       width = 18)

#---------------#
#------3C-------#
#---------------#

TregPlots = ggarrange(TregCTWT, ThymicTregWT, NaiveTregWT, ProlTregWT, TregCTKO, ThymicTregKO, NaiveTregKO, ProlTregKO,
                      labels = c("C"),
                      ncol = 4, nrow = 2)

ggsave(file = "~/my.work/PhD/HomestaticExpansionProject/Figures/ForPaper/Figure3/Figure3C.pdf", TregPlots,
       height = 8,
       width = 18)


