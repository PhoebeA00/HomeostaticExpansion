library(ggplot2)
library(ggpubr)

ModeldataWT = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Stats plots and data management/ModelOutputWT.csv")
ModeldataKO = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Stats plots and data management/ModelOutputKO.csv")
WTData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedWTSpleen.csv")
ProlWTData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/WTProl.csv")
KOData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedKOSpleen.csv")
ProlKOData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/KOProl.csv")


colnames(ModeldataWT) = c("NaiveCT", "ActTCT", "TregCT", "ThyNaive", "ActTNaive", "ThyTregs",
                          "TregNaive", "ProlNaive", "ProlActT", "ProlTreg", "IL-2", "ThymWeigth")
colnames(ModeldataKO) = c("NaiveCT", "ActTCT", "TregCT", "ThyNaive", "ActTNaive", "ThyTregs",
                          "TregNaive", "ProlNaive", "ProlActT", "ProlTreg", "IL-2", "ThymWeigth")

ModeldataWT$time = 0:431
ModeldataKO$time = 0:431

ggplot(WTData, aes(x=Age, y=NaiveCT)) + geom_point() +
  stat_summary(fun=mean, colour="red", geom="line")

###########################################
#  WT - Naive
##############################################

# NAIVE T CELLS

NaiveCTWT = ggplot(WTData, aes(x=hours, y=NaiveCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line")+
  geom_line(data = ModeldataWT, aes(x = time, y=NaiveCT), colour = "purple")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Total Naive T Cell Count(WT)", x = "Age in hours", y = "Cell Counts")

#PROLIFERATING NAIVE
ProlNaiveWT = ggplot(ProlWTData, aes(x=hours, y=NaiveProlCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line")+
  geom_line(data = ModeldataWT, aes(x = time, y=ProlNaive), colour = "purple")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Proliferating Naive T Cells(WT)", x = "Age in hours", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,4000000))

#THYMIC NAIVE
ThymicNaiveWT = ggplot(WTData, aes(x=hours, y=ThymicNaive)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line")+
  geom_line(data = ModeldataWT, aes(x = time, y=ThyNaive), colour = "purple")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Thymic Derived Naive T Cells(WT)", x = "Age in hours", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,4000000))

#############################################
# IL-2 KO Naive
#############################################


# NAIVE T CELLS

NaiveCTKO = ggplot(KOData, aes(x=hours, y=NaiveCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line")+
  geom_line(data = ModeldataKO, aes(x = time, y=NaiveCT), colour = "purple")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Total Naive T Cell Count(KO)", x = "Age in hours", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,4000000))

#PROLIFERATING NAIVE
ProlNaiveKO = ggplot(ProlKOData, aes(x=hours, y=NaiveProlCT)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line")+
  geom_line(data = ModeldataKO, aes(x = time, y=ProlNaive), colour = "purple")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Proliferating Naive T Cells(KO)", x = "Age in hours", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,4000000))

#THYMIC NAIVE
ThymicNaiveKO = ggplot(KOData, aes(x=hours, y=ThymicNaive)) + geom_point() +
  stat_summary(fun=mean, colour="black", geom="line")+
  geom_line(data = ModeldataKO, aes(x = time, y=ThyNaive), colour = "purple")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  labs(titles = "Thymic Derived Naive T Cells(KO)", x = "Age in hours", y = "Cell Counts")+
  scale_y_continuous(limits = c(0,4000000))

NaivePlots = ggarrange(NaiveCTWT, ProlNaiveWT, ThymicNaiveWT, NaiveCTKO, ProlNaiveKO, ThymicNaiveKO,
                    labels = c("A"),
                    ncol = 3, nrow = 2)

ggsave(file = "~/my.work/PhD/HomestaticExpansionProject/Figures/ForPaper/Figure3A.png", NaivePlots,
       height = 8,
       width = 18)
