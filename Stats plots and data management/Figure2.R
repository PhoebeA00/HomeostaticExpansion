library(scales)
library(ggplot2)
library(scales)
library(tidyr)
library(ggpubr)

#Paramters that control the plots
dotSize = 2 
#Makign a function that change the decimal places ofy axis values
scaleFUN <- function(x) sprintf("%.1f", x)
#Creating Y axis label for D-F
YlabelDF = expression(CD44^"+"~CD62L^"-"~(10^6))
YlabelA =  expression("%"~CD69^"+")
YlabelB = expression("%"~CD44^"+"~CD62L^"-")
YlabelC = expression("%"~Foxp3^"+")
###########################################33
# Percentage Difference between WT and IL-KO
#############################################


#--------------#
# Figure 2A
#--------------#

CD69Data = read.csv('~/my.work/PhD/HomestaticExpansionProject/ModelData/CD69DataFromGen.csv')
#Removing day 0, because it is always weird
CD69Data = subset(CD69Data, Age > 0)
CD69Data$CD4CD69_pct = CD69Data$CD4CD69_pct
CD69 = ggplot(CD69Data, aes(Age, CD4CD69_pct, color = Genotype)) + 
  scale_color_manual(values = c("#DB3D2F", "#3027E6"))+
  geom_point(position = position_dodge(1), size = dotSize)+
  scale_x_continuous(breaks = c(0,5,10,15,18))+
  scale_y_continuous(breaks = seq(0,max(CD69Data$CD4CD69_pct, na.rm = TRUE),
                                  length.out =5),
                     labels = scaleFUN)+
  labs(titles = "Early Activation", x = "Age in Days", y =YlabelA)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.position = "none",
        axis.title.x = element_text(colour="black", size=20),
        axis.title.y = element_text(colour = "black", size = 20),
        plot.title = element_text(lineheight=.8, size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  stat_summary(aes(group=Genotype, color = Genotype), fun=mean, geom="line")

#--------------#
# Figure 2B
#--------------#

ActivData = read.csv('~/my.work/PhD/HomestaticExpansionProject/ModelData/TCellActivationSummary_EdittedinR.csv')

CD62L = ggplot(ActivData, aes(Age, pct_CD4_CD44_pos_CD62L_neg, color = Genotype)) + 
  geom_point(position = position_dodge(1), size = dotSize)+
  scale_color_manual(labels = c('IL-2 KO', "WT"), values = c("#DB3D2F", "#3027E6"))+
  scale_x_continuous(breaks = c(0,5,10,15,18))+
  scale_y_continuous(breaks = seq(0,max(ActivData$pct_CD4_CD44_pos_CD62L_neg, na.rm = TRUE),
                                  length.out =5),
                     labels = scaleFUN)+
  labs(titles = "Activated T Cells", x = "Age in Days", y = YlabelB)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.position = "none",
        axis.title.x = element_text(colour="black", size=20),
        axis.title.y = element_text(colour = "black", size = 20),
        plot.title = element_text(lineheight=.8, size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  stat_summary(aes(group=Genotype, color = Genotype), fun=mean, geom="line")

#--------------#
# Figure 2C
#--------------#

TregDF = read.csv('~/my.work/PhD/HomestaticExpansionProject/ModelData/AfterCalculations.csv')
#Removing day 0, because it is always weird
TregDF = subset(TregDF, Age > 0 & Age < 50)
#So percentages show up right
TregDF$X4TregRatio = TregDF$X4TregRatio * 100   

Treg = ggplot(TregDF, aes(Age, X4TregRatio, color = Genotype)) + 
  geom_point(position = position_dodge(1), size = dotSize)+
  scale_color_manual(values = c("#DB3D2F", "#3027E6"))+
  scale_x_continuous(breaks = c(0,5,10,15,18))+
  scale_y_continuous(breaks = seq(0,max(TregDF$X4TregRatio, na.rm = TRUE),
                                  length.out =5),
                     labels = scaleFUN)+
  labs(titles = "Regulatory T Cell Percentage", x = "Age in Days", y = YlabelC)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.position = "none",
        axis.title.x = element_text(colour="black", size=20),
        axis.title.y = element_text(colour = "black", size = 20),
        plot.title = element_text(lineheight=.8, size = 20),
        axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  stat_summary(aes(group=Genotype, color = Genotype), fun=mean, geom="line")


#######################################
# Activation and Proliferating Data
#######################################

WTProl = read.csv('~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/WTProl.csv')
KOProl = read.csv('~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/KOProl.csv')
ActivatedWTSpleen = read.csv('~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedWTSpleen.csv')
ActivatedKOSpleen = read.csv('~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedKOSpleen.csv')


#--------------#
# Figure 2D
#--------------#
ActivatedWTSpleen$ActivatedCD4CTWT = ActivatedWTSpleen$ActivatedCD4CT

#Setting up Data
ActivatedCD4CTWT = subset(ActivatedWTSpleen, select = c("ActivatedCD4CTWT", "Age"))
LongActivatedCD4CTWT = gather(ActivatedCD4CTWT, variable, value, -Age)

ActivatedCD4CTKO = subset(ActivatedKOSpleen, select = c("ActivatedCD4CT", "Age"))
LongActivatedCD4CTKO = gather(ActivatedCD4CTKO, variable, value, -Age)
LongActivatedCD4CT = rbind(LongActivatedCD4CTWT, LongActivatedCD4CTKO)

LongActivatedCD4CT$variable = factor(LongActivatedCD4CT$variable, levels = c("ActivatedCD4CTWT", "ActivatedCD4CT"), labels = c("WT", "KO"))



ActTCT = ggplot(LongActivatedCD4CT, aes(x = Age, y = value, color = variable)) + 
  geom_point(position = position_dodge(1), size = dotSize)+
  stat_summary(aes(group=variable, color = variable), fun=mean, geom="line")+
  labs(titles = "Total Activated T Cells", x = "Age in Days", y = YlabelDF)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.position = "none",
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        #axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  # scale_y_continuous(breaks = seq(0,3880000, length.out = 5))+
  scale_y_continuous(limits=c(0, 7072000), breaks = c(0, 1768000, 3536000, 5304000, 7072000), labels = c(0, 1.7, 3.5, 5.3, 7.0))+
  scale_x_continuous(breaks = c(0,5,10,15,18))+
  scale_color_manual(values = c("#3027E6", "#DB3D2F"))+
  guides(color = guide_legend("Legend"))


#######################################
#           Naive Derived Activated T Cells
#######################################
ActivatedWTSpleen$ActivatedNaiveCTWT = ActivatedWTSpleen$ActivatedNaiveCT
#Setting up Data
ActivatedNaiveWT = subset(ActivatedWTSpleen, select = c("ActivatedNaiveCTWT", "Age"))
LongActivatedNaiveWT = gather(ActivatedNaiveWT, variable, value, -Age)

ActivatedNaiveKO = subset(ActivatedKOSpleen, select = c("ActivatedNaiveCT", "Age"))
LongActivatedNaiveKO = gather(ActivatedNaiveKO, variable, value, -Age)
LongActivatedNaive = rbind(LongActivatedNaiveWT, LongActivatedNaiveKO)

LongActivatedNaive$variable = factor(LongActivatedNaive$variable, levels = c("ActivatedNaiveCTWT", "ActivatedNaiveCT"), labels = c("WT", "KO"))

NaiveActT = ggplot(LongActivatedNaive, aes(x = Age, y = value, color = variable)) + 
  geom_point(position = position_dodge(1), size = dotSize)+
  stat_summary(aes(group=variable, color = variable), fun=mean, geom="line")+
  labs(titles = "Naive Derived Activated T Cells", x = "Age in Days", y = YlabelDF)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.position = "none",
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        #axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  # scale_y_continuous(breaks = seq(0,3880000, length.out = 5))+
  scale_y_continuous(limits=c(0, 7072000), breaks = c(0, 1768000, 3536000, 5304000, 7072000), labels = c(0, 1.7, 3.5, 5.3, 7.0))+
  scale_x_continuous(breaks = c(0,5,10,15,18))+
  scale_color_manual(values = c("#3027E6", "#DB3D2F"))+
  guides(color = guide_legend("Legend"))


#####################################################
#         Proliferating Activated T Cells
#####################################################
WTProl$ActivatedProlCTWT = WTProl$ActivatedProlCT

#Setting up Data
ActivatedProlCTWT = subset(WTProl, select = c("ActivatedProlCTWT", "Age"))
LongActivatedProlCTWT = gather(ActivatedProlCTWT, variable, value, -Age)

ActivatedProlCTKO = subset(KOProl, select = c("ActivatedProlCT", "Age"))
LongActivatedProlCTKO = gather(ActivatedProlCTKO, variable, value, -Age)
LongActivatedProlCT = rbind(LongActivatedProlCTWT, LongActivatedProlCTKO)

LongActivatedProlCT$variable = factor(LongActivatedProlCT$variable, levels = c("ActivatedProlCTWT", "ActivatedProlCT"), labels = c("WT", "KO"))

ProlActT = ggplot(LongActivatedProlCT, aes(x = Age, y = value, color = variable)) + 
  geom_point(position = position_dodge(1), size = dotSize)+
  stat_summary(aes(group=variable, color = variable), fun=mean, geom="line")+
  labs(titles = "Proliferating Activated T Cells", x = "Age in Days", y = YlabelDF)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.position = "none",
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        #axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  # scale_y_continuous(breaks =seq(0,2138000, length.out = 5 ))+
  scale_y_continuous(limits=c(0, 7072000), breaks = c(0, 1768000, 3536000, 5304000, 7072000), labels = c(0, 1.7, 3.5, 5.3, 7.0))+
  scale_x_continuous(breaks = c(0,5,10,15,18))+
  scale_color_manual(values = c("#3027E6", "#DB3D2F"))+
  guides(color = guide_legend("Legend"))


a = ggarrange(CD69, CD62L, Treg, ActTCT, NaiveActT, ProlActT,
          labels = c("A", "B", "C", "D", "E", "F" ),
          ncol = 3, nrow = 2)
# height - 1559 width = 837
ggsave(file = "~/my.work/PhD/HomestaticExpansionProject/Figures/ForPaper/Figure2/Figure2.pdf", a,
       height = 8,
       width = 15)

#------------------#
# For the Legned
#------------------#

JustForLegend = ggplot(LongActivatedProlCT, aes(x = Age, y = value, color = variable)) + 
  geom_point(position = position_dodge(1), size = dotSize)+
  stat_summary(aes(group=variable, color = variable), fun=mean, geom="line")+
  labs(titles = "Proliferating Activated T Cells", x = "Age in Days", y = YlabelDF)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.background = (element_rect(colour= "black", fill = "white")),
        axis.title.x = element_text( colour="black", size=20),
        axis.title.y = element_text( colour = "black", size = 20),
        plot.title = element_text(lineheight=.8,  size = 20),
        #axis.ticks.length=unit(.25, "cm"),
        text = element_text(size=20))+
  # scale_y_continuous(breaks =seq(0,2138000, length.out = 5 ))+
  scale_y_continuous(limits=c(0, 7072000), breaks = c(0, 1768000, 3536000, 5304000, 7072000), labels = c(0, 1.7, 3.5, 5.3, 7.0))+
  scale_x_continuous(breaks = c(0,5,10,15,18))+
  scale_color_manual(labels = c("WT", "IL-2 KO"), values = c("#3027E6", "#DB3D2F"))+
  guides(color = guide_legend(guide_legend(title = "Genotype")))

ggsave(file = "~/my.work/PhD/HomestaticExpansionProject/Figures/ForPaper/Figure2/legend.pdf", JustForLegend,
       height = 8,
       width = 15)
