require(survival)

setwd("~/Documents/Mitchell Lab/Manuscripts/Growth Chamber Work/Data/Experiment 3")
Key<-read.csv(file="Experiment_Key_R.csv", sep=",", header=T)
df_original<-read.csv(file="~/Documents/Mitchell Lab/Manuscripts/Growth Chamber Work/Data/ERC_Inoc_Data_R4.csv", sep=",", header=T)
df_merged<-merge(x  = df_original, y = Key, by = "Plant_ID", all.x = TRUE)

df_merged$Plant_ID<-as.factor(df_merged$Plant_ID)
df_merged$Plant_Leaf <- paste(df_merged$Plant_ID,"_",df_merged$Leaf_ID)
df_merged$RhizoctoniaLL<-as.numeric(df_merged$RhizoctoniaLL)

df_merged$Censor<-0
df_merged$Censor[df_merged$Leaf_Length=="dead"] <- 1
df_merged$Exposed<-0
df_merged$test<-df_merged$DAI-df_merged$Exposed
df_merged<-df_merged %>% na.omit()

names(df_merged)
new<-df_merged[,1:3]
new<-cbind(new, df_merged$DAI)
new<-cbind(new, df_merged[,14:20])
new$DAI<-new$`df_merged$DAI`

library(tidyr)
data_wide <- spread(new, DAI, Censor)
data_wide

w <- reshape(data_wide, 
             timevar = "DAI",
             idvar = c("Plant_ID", "Chamber_Shelf.x", "Leaf_ID", "Endo_Treatment", "Parasite_Treatment", "Plant_Leaf", "Exposed"),
             direction = "wide")

write.csv(w,'w.csv')
w<-read.csv(file="~/Documents/Mitchell Lab/Surv_Onlyoldest.csv", sep=",", header=T)

surv0 <- survfit(Surv(DeathDay,Censor)~1,data=w)
plot(surv0)
str(surv0)

#Create data frame with the time, estimated survival function, CI, and censors
survplotframe <- data.frame(time=surv0$time,surv=surv0$surv,
                            upper=surv0$upper,lower=surv0$lower,
                            censor=surv0$n.censor)

require(ggplot2)
theme_set(theme_bw())
survplot <- ggplot(survplotframe,aes(x=time,y=surv))
survplot+geom_step()+geom_step(aes(y=upper),linetype=2
)+geom_step(aes(y=lower),linetype=2)+geom_point(
  data=survplotframe[survplotframe$censor>0,]
  ,shape=3,size=4)


w<-w[w$Endo_Cat=="FreeFree"|w$Endo_Cat=="InoculatedConfirmed",]
w$Endo_Cat<-w$Endo_Cat %>% droplevels()
levels(w$Endo_Cat)<-c("FreeFree", "InoculatedConfirmed")
w$Parasite_Endo<-as.factor(paste(w$Parasite_Treatment,w$Endo_Cat, sep="_"))
summary(w$Parasite_Endo)

surv1 <- survfit(Surv(DeathDay,Censor)~Parasite_Treatment,data=w)
survplotframe<-data.frame(time=surv1$time,surv=surv1$surv,upper=surv1$upper,lower=surv1$lower,Parasite_Treatment=unlist(sapply(factor(c("Colletotrichum", "CONTROL", "CR", "Rhizoctonia")),function(x)rep(x,surv1$strata[x]))),censor=surv1$n.censor)
survplot1 <- ggplot(survplotframe,aes(x=time,y=surv,color=Parasite_Treatment))
survplot1+geom_step()+
  geom_step(aes(y=upper), linetype=2)+
  geom_step(aes(y=lower), linetype=2)+
  theme(legend.position="bottom", legend.title=element_blank())+
  xlab("Days After Inoculation")+
  ylab("Proportion Surviving")

#Compute log rank test/Endo
deathdiff<-survdiff(Surv(DeathDay,Censor)~Parasite_Treatment, data=w)
deathdiff #significant
#Add rho argument, which changes weighting of the observation on the survival history
#rho 1 gets Gehan-Wilcoxon test, which weights survivorship by the number surviving which focuses the statistic on early survivorship
deathdiff1<-survdiff(Surv(DeathDay,Censor)~Parasite_Treatment,data=w, rho=1)
deathdiff1 #more significant

#Cox proportional hazards model/Parasite_Treatment
deathcox <- coxph(Surv(DeathDay,Censor)~as.factor(Parasite_Treatment),data=w)
summary(deathcox) 
#The exponentiated coefficients represent the hazard ratio between the Parasite Treatment compared to the reference Treatment, Colletotrichum. 
#So there is 2.48 x higher hazard in CR than Colletotrichum only, for example

#Let's take care of some random effects
deathcox2 <- coxph(Surv(DeathDay,Censor)~as.factor(Parasite_Treatment)+cluster(Chamber_Shelf.x),data=w)
summary(deathcox2)
#decreases significance, only CR is significantly different from the reference.
w$Parasite_Treatment = relevel(w$Parasite_Treatment, ref = "CONTROL")
w2<-w[w$Censor>0,]
#w2<-w2[w2$Parasite_Treatment!="CONTROL",]
#w3<-w2[w2$Parasite_Treatment!="Colletotrichum",]
deathcox3 <- coxph(Surv(DeathDay,Censor)~as.factor(Parasite_Treatment)+cluster(Chamber_Shelf.x),data=w)
summary(deathcox3)
#significant difference with CR holds

#Baseline hazard function not estimated by the Cox model but it's possible to get an estimate by applying 
#the survfit on the Cox proportional hazards object
#have to do separately for each factor
cox1<-survfit(deathcox3,newdata=list(Parasite_Treatment=factor("CONTROL")))
cox2 <- survfit(deathcox3,newdata=list(Parasite_Treatment=factor("Colletotrichum")))
cox3 <- survfit(deathcox3,newdata=list(Parasite_Treatment=factor("Rhizoctonia")))
cox4 <- survfit(deathcox3,newdata=list(Parasite_Treatment=factor("CR")))

#Plot as comparison to K-M estimates
#new variable that indicates if it's a Cox estimate or not
survplotframe$cox<-"Kaplan-Meier"
#plot against K-M estimates
survplotframe <- rbind(survplotframe,data.frame(time=cox1$time,surv=cox1$surv,upper=cox1$upper,lower=cox1$upper, Parasite_Treatment="CONTROL", cox="Cox",
                                                censor=0))
survplotframe <- rbind(survplotframe,data.frame(time=cox2$time,surv=
                                                  cox2$surv,upper=cox2$upper,lower=cox2$lower,
                                                Parasite_Treatment="Colletotrichum",censor=0,cox="Cox"))
survplotframe <- rbind(survplotframe,data.frame(time=cox3$time,surv=
                                                  cox3$surv,upper=cox3$upper,lower=cox3$lower,
                                                Parasite_Treatment="Rhizoctonia",censor=0,cox="Cox"))
survplotframe <- rbind(survplotframe,data.frame(time=cox4$time,surv=
                                                  cox4$surv,upper=cox4$upper,lower=cox4$lower,
                                                Parasite_Treatment="CR",censor=0,cox="Cox"))

survplotframe = survplotframe[survplotframe$cox!="Kaplan-Meier",]

cbPalette <- c("#D55E00","#999999","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#levels(survplotframe$Parasite_Treatment)<-c("CONTROL", "Colletotrichum", "Rhizoctonia", "CR")
survplot1 <- ggplot(survplotframe,aes(x=time,y=surv,color=Parasite_Treatment))
Parasite<-survplot1+geom_step(size=1.2)+geom_point(data=survplotframe[survplotframe$censor>0,],
                                                   shape=3,size=4) + scale_linetype_manual(values=c("solid", "dashed"))+
  xlab("Days after Inoculation")+ylab("Proportion of Leaves Surviving")+scale_y_continuous(breaks = seq(0, 1, 0.2))+
  theme_bw()+
  geom_rug(data=w2, aes(x=DeathDay, y=0, color=Parasite_Treatment), length= unit(0.05, "npc"), inherit.aes = F, sides="b", position="jitter")+
  scale_color_manual(values=cbPalette, 
                     name="Parasite Inoculation\nTreatment",
                     breaks=c("CONTROL", "Colletotrichum", "Rhizoctonia", "CR"),
                     labels=c("Control", "Colletotrichum\nonly", "Rhizoctonia\nonly", "Coinfection"))+
  ylim(0,1)+
  theme(axis.text=element_text(size=12),
        axis.title= element_text(face="bold", size=12),
        axis.ticks=element_line(size=1),
        legend.text=element_text(size=12),
        legend.background = element_blank(),
        legend.position = "bottom", legend.title=element_text(size=12, face="bold"))+
  guides(color=guide_legend(nrow=2,byrow=TRUE))

surv3 <- survfit(Surv(DeathDay,Censor)~Parasite_Endo,data=w)
survplotframe1<-data.frame(time=surv3$time,surv=surv3$surv,upper=surv3$upper,lower=surv3$lower,
                           Parasite_Endo=unlist(sapply(factor(c("Colletotrichum_FreeFree", "CONTROL_InoculatedConfirmed", "CR_FreeFree", "Rhizoctonia_FreeFree","Colletotrichum_InoculatedConfirmed", "CR_InoculatedConfirmed", "Rhizoctonia_InoculatedConfirmed", "CONTROL_FreeFree")),
                                                       function(x)rep(x,surv3$strata[x]))),censor=surv3$n.censor)
survplot <- ggplot(survplotframe1,aes(x=time,y=surv,color=Parasite_Endo))
survplot+geom_step()+
  geom_step(aes(y=upper), linetype=2)+
  geom_step(aes(y=lower), linetype=2)+
  theme(legend.position="bottom")+
  xlab("Days After Inoculation")+
  ylab("Proportion Surviving")

#Compute log rank test/Endo
deathdiff<-survdiff(Surv(DeathDay,Censor)~Parasite_Endo, data=w)
deathdiff #significant
#Add rho argument, which changes weighting of the observation on the survival history
#rho 1 gets Gehan-Wilcoxon test, which weights survivorship by the number surviving which focuses the statistic on early survivorship
deathdiff1<-survdiff(Surv(DeathDay,Censor)~Parasite_Endo,data=w, rho=1)
deathdiff1 #more significant

w1<-w[!(w$Parasite_Endo=="CONTROL_FreeFree" | w$Parasite_Endo=="CONTROL_InoculatedConfirmed"| w$Parasite_Endo=="Colletotrichum_FreeFree"| w$Parasite_Endo=="Colletotrichum_InoculatedConfirmed"),]
w1$Parasite_Endo<-w1$Parasite_Endo %>%droplevels()
levels(w1$Parasite_Endo)
#Cox proportional hazards model/Parasite_Treatment
deathcox2 <- coxph(Surv(DeathDay,Censor)~as.factor(Parasite_Endo),data=w1)
summary(deathcox2) 
#The exponentiated coefficients represent the hazard ratio between the Parasite Treatment compared to the reference Treatment, Colletotrichum. 
#So there is 2.48 x higher hazard in CR than Colletotrichum only, for example

#Let's take care of some random effects
deathcox2 <- coxph(Surv(DeathDay,Censor)~as.factor(Parasite_Endo)+cluster(Chamber_Shelf.x),data=w1)
summary(deathcox2)
#decreases significance, only CR is significantly different from the reference.
w1$Parasite_Endo = relevel(w1$Parasite_Endo, ref = "Rhizoctonia_InoculatedConfirmed")
deathcox3 <- coxph(Surv(DeathDay,Censor)~as.factor(Parasite_Endo)+cluster(Chamber_Shelf.x),data=w1)
summary(deathcox3)
#significant difference with CR holds

#Baseline hazard function not estimated by the Cox model but it's possible to get an estimate by applying 
#the survfit on the Cox proportional hazards object
#have to do separately for each factor
cox2 <- survfit(deathcox2,newdata=list(Parasite_Endo=factor("Rhizoctonia_FreeFree")))
cox3 <- survfit(deathcox2,newdata=list(Parasite_Endo=factor("CR_FreeFree")))
cox5 <- survfit(deathcox2,newdata=list(Parasite_Endo=factor("Rhizoctonia_InoculatedConfirmed")))
cox6 <- survfit(deathcox2,newdata=list(Parasite_Endo=factor("CR_InoculatedConfirmed")))

#Plot as comparison to K-M estimates
#new variable that indicates if it's a Cox estimate or not
survplotframe1$cox<-"Kaplan-Meier"
survplotframe1$linetype<-"blank"
survplotframe1$Parasite<-"blank"
#plot against K-M estimates
survplotframe1 <- rbind(survplotframe1,data.frame(time=cox2$time,surv=
                                                    cox2$surv,upper=cox2$upper,lower=cox2$lower,
                                                  Parasite_Endo="Rhizoctonia_FreeFree",censor=0,cox="Cox", linetype="solid", Parasite="Rhizoctonia"))
survplotframe1 <- rbind(survplotframe1,data.frame(time=cox3$time,surv=
                                                    cox3$surv,upper=cox3$upper,lower=cox3$lower,
                                                  Parasite_Endo="CR_FreeFree",censor=0,cox="Cox", linetype="solid", Parasite="CR"))

survplotframe1 <- rbind(survplotframe1,data.frame(time=cox5$time,surv=
                                                    cox5$surv,upper=cox5$upper,lower=cox5$lower,
                                                  Parasite_Endo="Rhizoctonia_InoculatedConfirmed",censor=0,cox="Cox", linetype="dotted", Parasite="Rhizoctonia"))
survplotframe1 <- rbind(survplotframe1,data.frame(time=cox6$time,surv=
                                                    cox6$surv,upper=cox6$upper,lower=cox6$lower,
                                                  Parasite_Endo="CR_InoculatedConfirmed",censor=0,cox="Cox", linetype="dotted", Parasite="CR"))

survplotframe1 = survplotframe1[survplotframe1$cox!="Kaplan-Meier",]
survplotframe1$linetype<-as.factor(survplotframe1$linetype)
levels(survplotframe1$linetype)<-c("Endophyte Present", "Endophyte Absent")
survplotframe1$Parasite<-as.factor(survplotframe1$Parasite)
levels(survplotframe1$Parasite)<-c("Coinfection", "Rhizoctonia\nonly")
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
survplotframe1$Endophyte<-survplotframe1$linetype
survplot <- ggplot(survplotframe1,aes(x=time,y=surv,lty=Endophyte,color=Parasite))
w2<-w1[w1$Censor>0,]
w2$Parasite<-w2$Parasite_Treatment
w2$Parasite<-w2$Parasite %>% droplevels
w2$Endophyte<-w2$Endo_Cat
library(plyr)
w2$Endophyte<-revalue(w2$Endophyte, c("FreeFree"="Endophyte Absent", "InoculatedConfirmed"="Endophyte Present"))
w2$Parasite<-revalue(w2$Parasite, c("Rhizoctonia"="Rhizoctonia\nonly", "CR"="Coinfection"))
P_Endo<-survplot+geom_step(size=1.2)+geom_point(data=survplotframe1[survplotframe1$censor>0,],
                                                shape=3,size=4) + scale_linetype_manual(values=c("solid", "dotdash"))+
  xlab("Days after Inoculation")+ylab("Proportion of Leaves Surviving")+scale_y_continuous(breaks = seq(0, 1, 0.2))+
  theme_bw()+
  ylim(0,1)+
  geom_rug(data=w2, aes(x=DeathDay, y=0, color=Parasite, lty=Endophyte), length= unit(0.05, "npc"), inherit.aes = F, sides="b", position="jitter")+
  scale_color_manual(values=cbPalette, labels=c("Coinfection", "Rhizoctonia\nonly"))+
  theme(axis.text=element_text(size=12),
        axis.title= element_text(face="bold", size=12),
        axis.ticks=element_line(size=1),
        legend.text=element_text(size=12),
        legend.background = element_blank(),
        legend.position = "bottom", legend.title=element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE), linetype=guide_legend(nrow=2,byrow=TRUE))

library(cowplot)
pdf(file = "~/Desktop/GC_Fig4_031221.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches
prow <- plot_grid( Parasite + theme(legend.position="bottom"),
                   P_Endo + theme(legend.position="bottom"),
                   align = 'vh',
                   labels = c("A", "B"),
                   hjust = -1,
                   axis="l",
                   nrow = 1
)
prow
dev.off()
legend_b <- get_legend(LONG + theme(legend.position="bottom", legend.title=element_blank()))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
p
plot_grid(long, audps, labels = c("A", "B"))






