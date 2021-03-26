require(stringr)
require(dplyr)
require(tidyr)
library(agricolae)
require(nlme)
require(ggplot2)
require(grid)
require(broom)

GCdat<-tbl_df(read.csv('~/Documents/Mitchell Lab/Manuscripts/Growth Chamber Work/Data/Experiment1/Experiment1.csv', na.strings = c("NA", ""))) 
GCdat$fEndophyte<-as.factor(GCdat$Endophyte)

summary(GCdat) 

GCdat.long <- gather(GCdat, key = day, value = LL, 
                    X2DAI.020117, X3DAI.020217, X4DAI.020317, X7DAI.020617, X8DAI.020717, X9DAI.020817 ) %>% 
  mutate(DAI = as.numeric(substr(day, 2, 2)))
summary(GCdat.long)


#remove LL rows with missing values (they mess up the AUDPC calculation)
GCdat.long <- GCdat.long[!is.na(GCdat.long$LL),]


#lets look at the residuals
#this function takes the name of the model and the dataframe that was used to generate it
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

ggplot(GCdat.long, aes(x=DAI, y=LL))+
  geom_point() +
  geom_smooth()

ggplot(GCdat.long, aes(x=DAI, y=LL))+
  geom_point() +
  geom_line(aes(color = Plant_ID))

GCdat.long2 <- mutate(GCdat.long,
                    fEndophyte=as.character(Endophyte),
                    fStrain = as.character(Rs_Strain))  
GCdat.long2<-GCdat.long2[GCdat.long2$DAI<11,]
GCdat.long2$fStrain<-as.factor(GCdat.long2$fStrain)
GCdat.long2$fDAI<-as.factor(GCdat.long2$DAI)
GCdat.long2$fEndophyte<-as.factor(GCdat.long2$fEndophyte)

#Center Time
library(dplyr)
detach(package::plyr)
df_mean<-GCdat.long2 %>% 
  group_by(Plant_ID)%>%
  summarize(meanDAI=mean(DAI))

GCdat.long2a<-merge(x=GCdat.long2, y=df_mean, by="Plant_ID", all.x=TRUE)
GCdat.long2a$cDAI<-GCdat.long2a$DAI-GCdat.long2a$meanDAI
Rmod.gls <- gls(LL ~ Chamber_Shelf + cDAI*fEndophyte*fStrain, data=GCdat.long2a, method="REML") 
# It may be worth putting an autoregressive correlation function in your model since you're repeatedly measuring the 
# same individuals over time. 

GCdat.long2a$fEndophyte <- factor(GCdat.long2a$fEndophyte, levels = c("EI", "EF"))

# Start with a full model, don't forget to nest leaf within plant, since you have multiple nested observations per leaf
Rmod.full<-lme(LL ~  Chamber_Shelf + cDAI*fStrain*fEndophyte, 
           random=~cDAI|Plant_ID, 
           method="REML",
           data=GCdat.long2a,
           control=list(opt="optim",maxIter=5000000, msMaxIter=500000)) #let it run longer than the defaults
resids.fig(Rmod.full, GCdat.long2a) 
anova(Rmod.full)
# I've modified the random intercepts model to include uniqueID nested in plants, and the autoregressive correlation function
Rmod.ri <- lme(LL ~ Chamber_Shelf + cDAI*fEndophyte*fStrain, 
               random=~1|Plant_ID,
               method="REML",
               data=GCdat.long2a,
               control=list(opt="optim",maxIter=5000000, msMaxIter=500000))


  
anova(Rmod.gls, Rmod.full, Rmod.ri)
# random slopes is a significant improvement to the model (p<0.05)
# This says that the growth rate of parasites varies among leaves.

anova(Rmod.full)

Rmod.rs2 <- update(Rmod.full, weights = varIdent(form =~ 1 | fEndophyte*fStrain)) #that's a pretty serious weights function.
# may be worth considering if you really need the endophyte part of it.
anova(Rmod.full, Rmod.rs2) #it improves the model substantially

resids.fig(Rmod.rs2, GCdat.long2)
# to me, it doesn't look like we've solved the problem of heteroskedasticity here.
Rmod.rs3 <- update(Rmod.rs2, weights = varIdent(form =~ 1 | cDAI*fEndophyte*fStrain)) #that's a pretty serious weights function.

variables.fig(GCdat.long2, Rmod.rs2, "fEndophyte")

variables.fig(GCdat.long2, Rmod.rs2, "fStrain") #this looks better to me

anova(Rmod.rs2) #still no differences

#evaluate fixed effects using ML
Rmod.ml<-update(Rmod.full, method="ML")
require(car)
anova(Rmod.ml)
#no three way interaction
Rmod2 <- update(Rmod.ml, .~. -cDAI:fEndophyte:fStrain)
anova(Rmod.ml, Rmod2)
anova(Rmod2)
Rmod3 <- update(Rmod2, .~. -fEndophyte:fStrain)
anova(Rmod2, Rmod3)
anova(Rmod3)

##Plots of model results
cbPalette <- c( "#56B4E9","grey69", "#0072B2", "#56B4E9", "#E69F00", "#E69F00","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(broom.mixed)
y_title<-expression(paste(italic("Rhizoctonia"), " lesion length (cm)"))
x_title<-expression(paste("Days after ", italic("Rhizoctonia"), " inoculation"))
augment(Rmod3) %>% ggplot(aes(x=DAI,y=LL, group=fEndophyte,color=fEndophyte))+
  #geom_point(size=2, alpha=0.5, position="jitter")+
  geom_line(aes(group=Plant_ID), alpha=0.2)+
  geom_smooth(aes(y=.fixed, color=fEndophyte),size=3, method = "lm", se = F)+#this is what you want to change.
  labs(x=x_title, y=y_title)+
  scale_color_manual(values=cbPalette,
                     breaks=c("EF", "EI"),
                     labels=c("Endophyte absent", "Endophyte present"))+
  theme_bw()+
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10))+
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title.x= element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.ticks=element_line(size=2),
        legend.position="bottom",
        legend.title = element_blank(), 
        legend.text=element_text(size=18),
        legend.background = element_blank())



