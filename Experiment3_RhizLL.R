rm(list=ls())
require(stringr)
require(dplyr)
require(tidyr)
library(agricolae)
require(nlme)
require(ggplot2)
require(grid)
require(broom)
require(glmm)
library(glmm)

setwd("~/Documents/Mitchell Lab/Manuscripts/Growth Chamber Work/Data/Experiment 3")
Key<-read.csv(file="Experiment3_Key.csv", sep=",", header=T)
df_original<-read.csv(file="Experiment3.csv", sep=",", header=T)
df_merged<-merge(x = df_original, y = Key, by = "Plant_ID", all.x = TRUE)

df_merged$Plant_ID<-as.factor(df_merged$Plant_ID)
df_merged$Plant_Leaf <- paste(df_merged$Plant_ID,"_",df_merged$Leaf_ID)
df_merged<-na.omit(df_merged) #this also removed the IDs of leaves that died
df_merged$RhizoctoniaLL<-as.numeric(levels(df_merged$RhizoctoniaLL)[df_merged$RhizoctoniaLL])
df_merged<-df_merged %>% na.omit()
df_merged$Endo_Cat<-as.factor(paste(df_merged$Endo_Treatment, df_merged$Endo_Confirmed, sep="_"))

#subset just Rhizoc treatments
df_rhiz <- df_merged[df_merged$Parasite_Treatment == "CR" | df_merged$Parasite_Treatment == "Rhizoctonia" ,]

#excluding the plants that were inoculated with endophyte but it wasn't detected
df_rhiz <- df_rhiz[df_rhiz$Endo_Cat == "Inoculated_Confirmed" | df_rhiz$Endo_Cat == "Free_Free" ,] 

#Exclude Day 0
df_rhiz2<-df_rhiz[df_rhiz$DAI>0,]

#Only include the oldest leaf (the one I inoculated)
df_rhiz2<-df_rhiz2[df_rhiz2$Leaf_ID<2,]


#0s are when Rhizoctonia inoculation failed, only a few so excluded
df_rhiz2<-df_rhiz2[df_rhiz2$RhizoctoniaLL>0,] 
df_rhiz2$logRhizLL<-log(df_rhiz2$RhizoctoniaLL)

df_rhiz2<-df_rhiz2[df_rhiz2$DAI<12,]
df_rhiz2$Plant_Leaf<-as.factor(df_rhiz2$Plant_Leaf)
#Center Time
library(dplyr)
detach(package::plyr)
df_mean<-df_rhiz2 %>% 
  group_by(Plant_Leaf)%>%
  summarize(meanDAI=mean(DAI))

df_rhiz3<-merge(x=df_rhiz2, y=df_mean, by="Plant_Leaf", all.x=TRUE)
df_rhiz3$cDAI<-df_rhiz3$DAI-df_rhiz3$meanDAI

#Modeling
B1 <- gls(RhizoctoniaLL ~ 1 + Chamber*Shelf+cDAI*Parasite_Treatment*Endo_Cat,
            method = "REML", data = df_rhiz3)
B2 <- lme(RhizoctoniaLL ~1 + Chamber*Shelf+cDAI*Parasite_Treatment*Endo_Cat, 
          data = df_rhiz3, method="REML",
          random = ~1 | Plant_Leaf)
B5 <- lme(RhizoctoniaLL ~ 1 + Chamber*Shelf+cDAI*Parasite_Treatment*Endo_Cat, 
          data = df_rhiz3,
          control=list(opt="optim",maxIter=5000000, msMaxIter=500000), #let it run longer than the defaults)
          random = ~1 + cDAI | Plant_Leaf, method="REML")
B6 <- lme(RhizoctoniaLL ~ 1 + Chamber*Shelf+cDAI*Parasite_Treatment*Endo_Cat, 
          data = df_rhiz3,
          control=list(opt="optim",maxIter=5000000, msMaxIter=500000), #let it run longer than the defaults)
          random = ~1 + cDAI | Plant_ID/Plant_Leaf, method="REML")


AIC(B1, B2, B5, B6)
anova(B5, B6) #Ahh only one leaf per plant so redundant, going with B5
#going forward with B5


B5ml <- lme(RhizoctoniaLL ~ 1 + Chamber*Shelf+cDAI*Parasite_Treatment*Endo_Cat, 
            data = df_rhiz3,
            control=list(opt="optim",maxIter=5000000, msMaxIter=500000), #let it run longer than the defaults)
            random = ~1 + cDAI | Plant_Leaf, method="ML")

resids.fig(B5ml, df_rhiz3) #not homoscedastic, and lots of 0s

df_rhiz3$sqLL<-sqrt(df_rhiz3$RhizoctoniaLL)
df_rhiz3$CS<-as.factor(paste(df_rhiz3$Chamber, df_rhiz3$Shelf))
B5mlsq <- lme(sqLL ~ CS+cDAI*Parasite_Treatment*Endo_Treatment, 
          data = df_rhiz3,
          control=list(opt="optim",maxIter=5000000, msMaxIter=500000), #let it run longer than the defaults)
          random = ~1 + cDAI | Plant_Leaf, method="ML")



levels()
y_title<-expression(paste(italic("Rhizoctonia"), " lesion length (cm)"))
x_title<-expression(paste("Days after ", italic("Rhizoctonia"), " inoculation"))
cbPalette <- c(  "#E69F00","#56B4E9","#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
B5mlsq %>% augment() %>% 
  ggplot(aes(x=DAI,y=RhizoctoniaLL, shape=Parasite_Treatment, lty=Endo_Cat, color=Parasite_Treatment))+
  theme_bw() +
  scale_color_manual(values=cbPalette,
                     name="Inoculation Treatment",
                     breaks=c("CR", "Rhizoctonia"),
                     labels=c("Coinfection", "Single Infection"))+
  #geom_line(data =  df_rhiz2, aes(x=DAI, y=logRhizLL, group = Plant_Leaf, color=Parasite_Treatment, alpha=0.00001))+
  geom_line(aes(group=Plant_Leaf), size=0.5, alpha=0.2)+
  geom_smooth(aes(y=(.fixed)^2),size=3, method = "lm", se = F)+
  #facet_wrap(~Endo_Cat)+
  #geom_line(aes(y=.fixed),size=3)+
  labs(x=x_title, y=y_title)+
  #guides(color=guide_legend(title="Parasite Treatment"))+
  guides(lty=FALSE)+
  scale_linetype_manual(values=c("solid", "dotted"))+
  theme(axis.text=element_text(size=12),
        axis.title= element_text(size=12),
        axis.ticks=element_line(size=1),
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        legend.background = element_blank(), legend.position = "bottom")
  
  
resids.fig <- function(mod, df) {
    residdf <- dplyr::mutate(df, resids = residuals(mod, type = 'normalized'),
                             fits = fitted(mod))
    fig2 <-ggplot(residdf, aes(x = fits, y = resids)) + geom_point() +
      labs(x = 'Fitted values', y = '')
    
    fig3 <- ggplot(residdf) + stat_qq(aes(sample = resids)) +
      labs(x = 'Theoretical Quantiles', y = 'Sample Quantiles')
    
    # qqline plot = FALSE, according to James should work
    
    fig4 <- ggplot(residdf, aes(x = resids)) + geom_histogram(aes(y=..density..), colour = 'grey50') +
      labs(x = 'Residuals', y = 'Frequency') + scale_y_continuous(expand = c(0, 0)) +
      stat_function(fun = dnorm, color = "red", args = list(mean = mean(residdf$resids),
                                                            sd = sd(residdf$resids)))
    grid.draw(rbind(ggplotGrob(fig2), ggplotGrob(fig3), ggplotGrob(fig4), size = 'first'))
    
    return(summary(mod))
}

  
variables.fig <- function(df, mod, variable){
    df %>% 
      mutate(resids = residuals(mod, type = 'normalized')) %>% 
      ggplot(aes_string(x=variable, y="resids", color=variable)) +
      xlab(variable) +
      ylab("residuals") +
      geom_boxplot() +
      geom_point(position=position_jitter(h=0, w=0.4))
}

ggplot(df_rhiz2, aes(x=DAI, y=exp(logRhizLL), color=Parasite_Treatment))+
  #geom_point(aes(color= timing), size = 0.9) + xlab("DAI") + ylab("Rs Lesion Length")+
  geom_line(aes(group=Plant_Leaf, alpha=0.01))+
  geom_smooth(aes(group = Parasite_Treatment), size=3, se = FALSE, method = lm)+
  facet_wrap(~Endo_Treatment)
#geom_line(aes(group = uniqueID, color=timing), alpha = 0.65)+theme_bw()+scale_x_continuous(breaks =seq(0, 12, 2))+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(size = 20),
      axis.text=element_text(size=20),
      axis.title= element_text(face="bold", size=16),
      axis.ticks=element_line(size=2),
      legend.title = element_text(size=16, face="bold"), 
      legend.text=element_text(size=12),
      legend.background = element_blank())


