#############################################################################
## This script was developed by Dr. Erin Sauer for Sauer et al. (submitted) #
# "Male pathology regardless of behaviour drives transmission in an         #
#  avian host-pathogen system"                                              #
# The script fits all statistical models described in the manuscript        #
# Data were collected in the canary - MG experiment described in the        #
#     manuscript, unless otherwise noted.                                   #
#############################################################################

library(lme4)
library(car)
library(tidyverse)
library(mgcv)

######## Eye scores #############################
es <- read.csv("Eye score data.csv")
str(es)
es$Day <- as.numeric(es$Day) 
es$Sex <- as.factor(es$Sex)
#only MG exposed birds are included in this analysis
es <- subset(es, Treatment == "MG")  

es.m <- gamm(ES ~ Sex + s(Day), niterPQL=100, family=poisson(),
            correlation=corCAR1(form=~Day|ID), data = es)
plot(es.m$gam, pages=1)
summary(es.m$gam)

esfig <- ggplot(es, aes(y=ES, x=Day, color=Sex)) +
   geom_smooth()+
   geom_point()+
   theme_bw()+
   coord_cartesian(ylim=c(0, 6))+
   geom_jitter(height = .01) +
   theme(axis.text=element_text(size=20, color="black"),
         panel.border = element_blank(),
         axis.title=element_text(size=25),
         legend.position="top",
         legend.text=element_text(size=20),
         legend.title=element_blank())+
   ylab("Total eye score")
ggsave("esfig.png", width = 2000, height = 2000, units = "px", dpi=300)

#############   House Finch data for comparison   ###########################
# This data is a subset from Adelman et al 2015 doi: 10.1098/rspb.2015.1429 #
#  and is publicly available. Here, only eye scores from the index birds    #
#  are used.                                                                #
#############################################################################

HOFIes <- read.csv("Adelman et al 2015.csv")
str(HOFIes)
HOFIes$sex <- as.factor(HOFIes$sex)

ggplot(HOFIes, aes(y=es, x=day, color=sex)) +
  geom_smooth()+
  geom_point()+
  theme_bw()+
  coord_cartesian(ylim=c(0, 6))+
  geom_jitter(height = .01) +
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=25),
        legend.position="right",
        legend.text=element_text(size=20),
        legend.title=element_blank())+
  ylab("Total eye score")+xlab("Day")

mHOFI <- gamm(es ~ sex + s(day , k=3), family=poisson(),
           correlation=corCAR1(form=~day|id), data = HOFIes)
plot(mHOFI$gam, pages=1) #this plot lets us visually assess the fit
summary(mHOFI$gam)

################## MG load #######################
mg <- read.csv("MG load.csv")
str(mg)
mg$Day <- as.numeric(mg$Day) #make day numeric

# this analysis only uses MG exposed birds and, 
#  for the purpose of useing a linear mode, does not include
#  days 0 (no infections) or day 35 when most birds had recovered
mg <- subset(mg, Treatment == "MG")
mg <- subset(mg, Day != 0)
mg <- subset(mg, Day != 35)

m3 <- lmer(logLoad ~ Sex * Day + (1|ID), mg)
summary(m3)
Anova(m3)

#this code mades a database of predicted data from the model and plots it
points5 <- predict(m3, mg)
unname(points5)
pps5 <- as.data.frame(cbind(points5, mg))

mgfig<-ggplot(pps5, aes(y=points5, x=Day, color=Sex)) +
  #geom_boxplot()+
  geom_smooth(method="lm")+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=25),
        legend.position="top",
        legend.text=element_text(size=20),
        legend.title=element_blank())+
  ylab("log10(MG load)")

ggsave("mgfig.png", width = 2000, height = 2000, units = "px", dpi=300)

######## Body condition #############################
# this section has both mass and fat score analyses #
#####################################################

bc <- read.csv("Body condition data.csv")
str(bc) 
bc$Day <- as.numeric(bc$Day) #make day numeric

######## fat score ##########
bc$ST <- revalue(bc$ST, c("FemaleControl"="Female:Control", 
                          "MaleControl"="Male:Control",
                          "MaleMG"="Male:MG",
                          "FemaleMG"="Female:MG"))

bc$STorder <- revalue(bc$ST, c("Female:Control"=2, 
                               "Male:Control"=4,
                               "Male:MG"=3,
                               "Female:MG"=1))
bc$STorder <- as.numeric(bc$STorder)
#raw data plot
ggplot(bc, aes(y=Fat, x=Day, color=ST)) +
  geom_smooth(method="lm")+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))+
  ylab("Fat Score")

fatm2 <- lmer(Fat ~ Sex * Day * Treatment + (1|ID),  bc)
summary(fatm2)
Anova(fatm2)

# make predictive plot which is used in supplement
points3 <- predict(fatm2, bc)
View(points3)
unname(points3)
pps3 <- as.data.frame(cbind(points3, bc))
pps3$ST <- as.factor(pps3$ST)
pps3 <- pps3 %>% mutate(ST = fct_reorder(ST, STorder))
#plot!
fatfig <- ggplot(pps3, aes(y=points3, x=Day, color=ST)) +
  geom_smooth(method="lm")+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="top",
        legend.text=element_text(size=13),
        legend.title=element_blank())+
  ylab("Fat Score")
fatfig
ggsave("fatfig.png", width = 2000, height = 2000, units = "px", dpi=300)

######## body mass ########
#raw data plot
ggplot(bc, aes(y=Mass, x=Day, color=ST)) +
  geom_smooth(method="lm")+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_blank())+
  ylab("Mass (g)")

massm2 <- lmer(Mass ~ Sex * Day * Treatment + (1|ID), bc)
summary(massm2)
Anova(massm2)

#make predictive plot for supplement
points4 <- predict(massm2, bc)
#View(points4)
unname(points4)
pps4 <- as.data.frame(cbind(points4, bc))
str(pps4)

pps4 <- pps4 %>% mutate(ST = fct_reorder(ST, STorder))

massfig <- ggplot(pps4, aes(y=points4, x=Day, color=ST)) +
  geom_smooth(method="lm")+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(size=15, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=20),
        legend.position="top",
        legend.text=element_text(size=13),
        legend.title=element_blank())+
  ylab("Mass (g)")
bodycfig <- ggarrange(fatfig,massfig,ncol=1, nrow=2, 
                  labels=c("A","B"), font.label=list(size=25))
ggsave("bodycfig.png", 
       width = 2000, height = 3000, units = "px", dpi=300)

######### Antibody data ######################################
ab <- read.csv("antibody data.csv")
str(ab)
ab$Day <- as.numeric(ab$Day) 
#this analysis only uses MG exposed birds
abMG <- subset(ab, Treatment == "MG") 

m1 <- lmer(antibody ~ Sex * Day + (1|ID), abMG)
summary(m1)
Anova(m1)
#make predictive plot for manuscript
points <- predict(m1, abMG)
View(points)
unname(points)
pps <- as.data.frame(cbind(points, abMG))

#plot!
abfig <- ggplot(ab, aes(y=antibody, x=Day, color=Sex)) +
  geom_smooth(method=lm)+
  geom_point()+
  theme_bw()+coord_cartesian(ylim=c(0, 0.115))+ 
  #geom_hline(yintercept=0.0229, linetype="dashed", color="black", size=.3)+
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title.y =element_text(size=20),
        axis.title.x=element_text(size=25),
        legend.position="top",
        legend.text=element_text(size=20),
        legend.title=element_blank())+
        #legend.text=element_text(size=15),
        #legend.title=element_text(size=20))+
  ylab("MG-specific antibody (OD)")
abfig
ggsave("abfig.png", width = 2000, height = 2000, units = "px", dpi=300)

###### fig2 panel ######################################
Fig2 <- ggarrange(esfig,mgfig,abfig, ncol=3, nrow=1, 
                  labels=c("A","B","C"), font.label=list(size=25))
ggsave("Fig2.png", 
       width = 4000, height = 1500, units = "px", dpi=300)

###### WBC counts ########################################
wbc <- read.csv("WBC counts.csv")
str(wbc)
#analysis of change in relative abundance of leukocytes 
# from pre-exposure to one week post exposure
wbc7 <- subset(wbc, Day < 8)
colnames(wbc7)
wbc7 <- wbc7[,c(2:6,8:11)]
wbc7 <- wbc7 %>% group_by(ST,ID) %>% 
  pivot_wider(names_from = Day, values_from = c(Eosinphil, Heterophil,Lymphocyte,Monocyte)) %>%
  mutate(Epc = (Eosinphil_7 - Eosinphil_0)/Eosinphil_0*100,
         Hpc = (Heterophil_7 - Heterophil_0)/Heterophil_0*100,
         Lpc = (Lymphocyte_7 - Lymphocyte_0)/Lymphocyte_0*100,
         Mpc = (Monocyte_7 - Monocyte_0)/Monocyte_0*100,
         Eosinophil = (Eosinphil_7 - Eosinphil_0),
         Heterophil = (Heterophil_7 - Heterophil_0),
         Lymphocyte = (Lymphocyte_7 - Lymphocyte_0),
         Monocyte = (Monocyte_7 - Monocyte_0),
         HLR_0=(Heterophil_0/Lymphocyte_0),
         HLR_7=(Heterophil_7/Lymphocyte_7),
         dHLR=(HLR_7-HLR_0)) 
View(wbc7) 
str(wbc7)
wbc7 <- do.call(data.frame, lapply(wbc7, function(x) {
  replace(x, is.infinite(x) | is.na(x), NA)
})
)
View(wbc7)
#fix some NAs that should be 0
wbc7[1,16] <- 0
wbc7[15,16] <- 0
wbc7[12,16] <- 0

wbcM <- manova(cbind(Eosinophil,Heterophil,Lymphocyte,Monocyte)~
                 Sex * Treatment, data=wbc7)
summary(wbcM)
summary.aov(wbcM)

colnames(wbc7)
wbcplot <- wbc7[,c(4,17:20)]
wbcplot <- wbcplot %>% pivot_longer(!ST, names_to = "WBC", values_to = "values")
wbcplot$ST <- revalue(wbcplot$ST, c("FemaleControl"="Female:Control", 
                              "MaleControl"="Male:Control",
                              "MaleMG"="Male:MG",
                              "FemaleMG"="Female:MG"))

wbcplot$STorder <- revalue(wbcplot$ST, c("Female:Control"=2, 
                                   "Male:Control"=4,
                                   "Male:MG"=3,
                                   "Female:MG"=1))
wbcplot$STorder <- as.numeric(wbcplot$STorder)
wbcdf <- ggplot(wbcplot, aes(y=values, x=ST, color=WBC)) +
 geom_boxplot()+
  theme_bw()+
  #coord_cartesian(ylim=c(-102, 250))+
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=25),
        legend.position="top",
        legend.text=element_text(size=20),
        legend.title=element_blank())+
  ylab("Change in count post exposure")+xlab("")

wbcdf
ggsave("wbcdf.png", 
       width = 2500, height = 2000, units = "px", dpi=300)

#pre exposure
wbc0 <- subset(wbc, Day ==0)
str(wbc0)
wbcM0 <- manova(cbind(Eosinphil,Heterophil,Lymphocyte,Monocyte)~
                     Sex , data=wbc0)
sum1 <- summary(wbcM0)
sum2 <- summary.aov(wbcM0)

ggplot(wbc0, aes(y=HLR, x=Sex)) +
  geom_boxplot()+
  theme_bw()+
  #coord_cartesian(ylim=c(-102, 250))+
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=25),
        legend.position="top",
        legend.text=element_text(size=20),
        legend.title=element_blank())+
  ylab("Eosinphil")+xlab("")

###### hematocrit ###########
hema <- read.csv("hematocrit.csv")
str(hema)
hema$Sex <- as.factor(hema$Sex)
#this analysis did not include values prior to exposure
hema <- subset(hema, Hem >0)
hema$ST <- revalue(hema$ST, c("ControlFemale"="Female:Control", 
                          "ControlMale"="Male:Control",
                          "MGMale"="Male:MG",
                          "MGFemale"="Female:MG"))

hema$STorder <- revalue(hema$ST, c("Female:Control"=2, 
                               "Male:Control"=4,
                               "Male:MG"=3,
                               "Female:MG"=1))
hema$STorder <- as.numeric(hema$STorder)

#MK 114 was excluded from this analysis because this bird had 
# abnormally high hematocrit both before and after exposure.
# (>80% compared to all other birds <70%, nearly all between 50-70)
hema <- subset(hema, ID!="MK 114")

hem1 <- lmer(Hem ~ Sex * Day * Treatment + (1|ID), hema)
summary(hem1)
Anova(hem1)
#making predictive plot fig. S3
hempred <- predict(hem1, hema)
unname(hempred)
hempred <- as.data.frame(cbind(hempred, hema))

hempred <- hempred %>% mutate(ST = fct_reorder(ST, STorder))
hemaplot <- ggplot(hempred, aes(y=hempred, x=Day, color = ST)) +
  #geom_boxplot()+
  #geom_line()+
  geom_smooth(method="lm")+
  geom_point()+
  #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5) +
  theme_bw()+
  theme(axis.text=element_text(size=20, color="black"),
        panel.border = element_blank(),
        axis.title=element_text(size=25),
        legend.position="top",
        legend.text=element_text(size=13),
        legend.title=element_blank())+
  ylab("% Hematocrit")
ggsave("hemaplot.png", 
       width = 2000, height = 1500, units = "px", dpi=300)
