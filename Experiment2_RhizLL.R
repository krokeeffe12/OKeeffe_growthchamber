require(stringr)
require(dplyr)
require(tidyr)
library(agricolae)
require(nlme)
require(ggplot2)
require(grid)
require(broom)

#23 is last day for Rs symptoms in Rs2Cc1
#32 is last day for Col symptoms in Rs1Cc2

######read in data to use
setwd("~/Documents/Mitchell Lab/Manuscripts/Growth Chamber Work/Data/Experiment 2")
GCdat<-(tbl_df(read.csv('~/Downloads/GC_len_correction_24daysFINAL_edited.csv', na.strings = c("NA", "")))) 


summary(GCdat) 
GCdat$Plant_ID = as.factor(GCdat$plant_ID) 


######reshape data for easier analysis
GCdat.long <- gather(GCdat, key = day, value = Disease.Len, contains("DAI.len.")) %>%
  #Reshuffle columns so that columns for each day after infection become rows and the measures are under new column called LL for "lesion length"
  mutate(DAI = substr(day, 9,11)) #Makes DAI column from the days

GCdat.long$uniqueID = paste(GCdat.long$plant_ID, GCdat.long$leaf_ID, sep="")
GCdat.long <- GCdat.long[!is.na(GCdat.long$Disease.Len),] #remove rows with missing values (they mess up the AUDPC calculation)
GCdat.long$Disease.Len = as.numeric(GCdat.long$Disease.Len) #make disease lengths numbers
GCdat.long$DAI = as.numeric(GCdat.long$DAI)

head(GCdat.long)

####longitudinal model of rhizoctonia lesions: timing EFFECT
#creating dataframe and model 
GCdat.long.l = GCdat.long[GCdat.long$disease_type == "L",]
GCdat.long.l <- GCdat.long.l[!is.na(GCdat.long.l$Disease.Len),]
GCdat.long.l <- GCdat.long.l[GCdat.long.l$timing == "rs1cc1" | GCdat.long.l$timing == "rs1cc2" | GCdat.long.l$timing == "rs2cc1"| GCdat.long.l$timing == "rs",]
GCdat.long.l <- GCdat.long.l[GCdat.long.l$DAI <=12,]
GCdat.long.l$CS<-paste(GCdat.long.l$chamber_ID, GCdat.long.l$shelf_ID)
GCdat.long.l$CS<-as.factor(GCdat.long.l$CS)
#No random effects
mod<-gls(Disease.Len ~ CS + DAI*timing, data = GCdat.long.l, method="REML")
#Random intercepts
modL = lme(Disease.Len ~ CS + DAI*timing, random = ~1|plant_ID, method="REML", data = GCdat.long.l)
AIC(mod, modL)
anova(mod, modL) #modL is better

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
#Random slopes
modLa = lme(Disease.Len ~ CS + DAI*timing, random = ~DAI|plant_ID/leaf_ID, 
            method="REML", data = GCdat.long.l,
            control=list(opt="optim", maxIter=5000000, msMaxIter=5000000))
AIC(modL, modLa)
anova(modL, modLa) #modLa is better
resids.fig(modLa, GCdat.long.l) #heteroscedasticity issue

modLa.anova = anova(modLa)
modLaML = lme(Disease.Len ~ CS + DAI*timing, random = ~DAI|plant_ID/leaf_ID, 
            method="ML", data = GCdat.long.l,
            control=list(opt="optim", maxIter=5000000, msMaxIter=5000000))

modLaML.anova = anova(modLaML)
resids.fig(modLaML, GCdat.long.l) #normal, but not homoscedastic

#let's not force through 0
GCdat.long.l2 = subset(GCdat.long.l, DAI!=0)

#modL3=no random effects
modL3 = gls(Disease.Len ~ CS + DAI:timing, method="REML", data = GCdat.long.l2)
modL4 = lme(Disease.Len ~ CS + DAI:timing, random=~DAI|plant_ID/uniqueID, data=GCdat.long.l2,
            method="REML",
            control=list(opt="optim", maxIter=5000000, msMaxIter=5000000))
AIC(modL3, modL4) #modL4 is better
modL4.anova = anova(modL4)

resids.fig(modL4, GCdat.long.l2) #normal not homoscedastic
#log transformation

#Center Time
library(dplyr)
detach(package::plyr)
df_mean<-GCdat.long.l2 %>% 
  group_by(uniqueID)%>%
  summarize(meanDAI=mean(DAI))

GCdat.long.l2a<-merge(x=GCdat.long.l2, y=df_mean, by="uniqueID", all.x=TRUE)
GCdat.long.l2a$cDAI<-GCdat.long.l2a$DAI-GCdat.long.l2a$meanDAI

modL4a = lme(log1p(Disease.Len) ~ CS+cDAI*timing, random=~cDAI|plant_ID/uniqueID, data=GCdat.long.l2a,
            method="ML",
            control=list(opt="optim", maxIter=5000000, msMaxIter=5000000))
resids.fig(modL4a, GCdat.long.l2) #better
anova(modL4a)
summary(modL4a)
require(lsmeans)
library(lsmeans)

levels(GCdat.long.l2a$timing)
lsmeans(modL4a, pairwise~timing, adjust="tukey", transform="log")
emtrends(modL4a, pairwise ~ timing, var = "cDAI")
emmeans (modL4a, pairwise~DAI*timing, adjust="tukey")


GCdat.long.l2a$timing<-GCdat.long.l2a$timing %>% droplevels
levels(GCdat.long.l2a$timing)
GCdat.long.l2a$timing <- factor(GCdat.long.l2a$timing, levels = c("rs1cc1", "rs1cc2", "rs2cc1", "rs"))
modL4a_relevel = lme(log1p(Disease.Len) ~ CS+cDAI*timing, random=~cDAI|plant_ID/uniqueID, data=GCdat.long.l2a,
             method="ML",
             control=list(opt="optim", maxIter=5000000, msMaxIter=5000000))
summary(modL4a_relevel)
lsmeans(modL4a_relevel, pairwise~timing, adjust="tukey", transform="log")
emtrends(modL4a_relevel, pairwise ~ timing, var = "cDAI")

require(emmeans)
emmeans(modL4a, pairwise ~ timing )
#only simultaneous is significantly different


######LINEAR MODEL WITH RANDOM EFFECTS
cbPalette <- c("#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
y_title<-expression(paste(italic("Rhizoctonia"), " lesion length (cm)"))
x_title<-expression(paste("Days after ", italic("Rhizoctonia"), " inoculation"))
augment(modL4a_relevel) %>% ggplot(aes(x=DAI,y=Disease.Len, group=timing,color=timing))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(aes(group=uniqueID), alpha=0.3)+
  labs(x=x_title, y=y_title)+
  #geom_point(size=2, alpha=0.5, position="jitter")+
  geom_smooth(aes(y=exp(.fixed), color=timing),size=3, method = "lm", se = F)+ #this is what you want to change.
  scale_color_manual(values=cbPalette,
                     name="Treatment",
                     breaks=c("rs", "rs1cc1", "rs1cc2", "rs2cc1"),
                     labels=c("Single Inoculation", "Coinfection: Simultaneous", "Coinfection: Rhiz. First", "Coinfection:Rhiz. Second"))+
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title.x= element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.ticks=element_line(size=2),
        legend.position="bottom",
        legend.title = element_blank(), 
        legend.text=element_text(size=12),
        legend.background = element_blank())+
  guides(colour = guide_legend(nrow = 2))

