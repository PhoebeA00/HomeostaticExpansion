WTData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedWTSpleen.csv")
ProlWTData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/WTProl.csv")
KOData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/ActivatedKOSpleen.csv")
ProlKOData = read.csv("~/my.work/PhD/HomestaticExpansionProject/Code/Modeling/Matlab/RawData/KOProl.csv")
pop = read.csv("/home/jon/my.work/PhD/HomestaticExpansionProject/ModelData/AfterCalculations.csv")

WTData$hours



# Mu - Thymic naive
ThymicNaive = summary(lm(data=WTData, ThymicNaive~hours))$coefficients[2]

# z - Self replication
NaiveProlCT = summary(lm(data=ProlWTData, NaiveProlCT~hours))$coefficients[2]

# alpha - Thymic naive
ThymicDerivedTregsCT = summary(lm(data=WTData, ThymicDerivedTregsCT~hours))$coefficients[2]

# c - Naive Dervied Tregs
NaiveDerivedTregsCT = summary(lm(data=WTData, NaiveDerivedTregsCT~hours))$coefficients[2]

# epsilon - self replicating Tregs
X4TregProlCT = summary(lm(data=ProlWTData, X4TregProlCT~hours))$coefficients[2]

# beta - activation rate of naive t cells
ActivatedNaiveCT = summary(lm(data=WTData, ActivatedNaiveCT~hours))$coefficients[2]

# a - self replication rate of activated T cells
ActivatedProlCT = summary(lm(data=ProlWTData, ActivatedProlCT~hours))$coefficients[2]

parameters = c("\\mu$", "$z$", "$\\alpha$", "$c$", "$\\epsilon$", "$\\beta$", 
               "$a$")
rates = c(ThymicNaive, NaiveProlCT, ThymicDerivedTregsCT, NaiveDerivedTregsCT,X4TregProlCT,
          ActivatedNaiveCT, ActivatedProlCT)

Rates.data = data.frame(parameters, rates)

write.csv(Rates.data, "~/my.work/PhD/HomestaticExpansionProject/ModelData/LinearRegressionRates.csv")

# # rK - Carrying capacity for Tregs
# popday56 = subset(pop, intage > 40)
# Tregs_rK = mean(popday56$X4TregCT) * 10**6


pd <- position_dodge(1)
ggplot(WTData,aes(hours, ThymicNaive)) +
  # stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm')+
  geom_point(position=pd, size=3)+
  labs(title = "Naive T Cells - 5005 cells/hour", y = "Total counts", x = "Age in Days")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text( size = 20, face = "bold")
  )

