# Bd and sperm quality - pre and post chytrid
# L. aurea
# Sept-24, Rose Upton

library(lme4)
library(emmeans)
library(glmmTMB)
library(tidyverse)
library(lmerTest)
library(readxl)
library(bbmle) ## for AICtab
library(gridExtra) # Multiple ggplots per page
library(ggeffects)
library(DHARMa)
library(ggplot2)


# From Ben Bolker's glmm faq page - for assessing overdispersion using Pearson residuals
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

setwd("C://Users//Rose Upton//OneDrive - The University Of Newcastle//Documents//UoN//Students//Anne//Forum_chytrid_bellfrogs")
data1=read_excel("2024_Bd_sperm_quality.xlsx",sheet="Sheet1")

## For models
# Calculate number of sperm/ml and sperm/sample
data1 = data1 %>% mutate(sperm_ml = (sperm*dil*10000/primary_sq)) %>%
  mutate(sperm_sample = sperm_ml * vol / 1000)

# # Scaling factor that converts sperm count to sperm per sample
# data1 = mutate(data1, scale_sample=(dil*10000/primary_sq) * vol / 1000)

# Scaling factor that converts sperm count to sperm per ml
data1 = mutate(data1,scale_ml=(dil*10000/primary_sq))

# Also reciprocal of the scaling factor per ml and log of reciprocal
data1 = data1 %>% mutate(scale_ml_inv = 1/scale_ml) %>%
  mutate(log_scale_ml_inv = log(scale_ml_inv))
names(data1)

# Data Filtering
# Aim 1: Pre-heat treatment only = data 2
data2a=filter(data1, heat_treat !="post")
data2=filter(data2a, heat_treat !="post2")

# Aim 2: First two time points only = data 3
data3=filter(data1, heat_treat !="post2")

# Aim 3: All time points, Chytrid pos only = data 4
data4=filter(data1, chytrid !="neg")

# Remove frogs from first two time points that don't have final time point
# For data exploration
data5a=filter(data4, chip !="120")
data5b=filter(data5a, chip !="31C")
data5c=filter(data5b, chip !="4F9")
data5d=filter(data5c, chip !="A6F")
data5=filter(data5d, chip !="B75")

## Data Exploration to inform model design
ggplot(data=data1, aes(x=heat_treat, y=spconc, color=chytrid))+geom_point()
ggplot(data=data1, aes(x=heat_treat, y=spconc, color=severity))+geom_point()
ggplot(data=data1, aes(x=heat_treat, y=spconc, color=severity))+geom_point()+facet_wrap('severity')

ggplot(data=data1, aes(x=heat_treat, y=spconc, color=Tub))+geom_point()
ggplot(data=data1, aes(x=heat_treat, y=spconc, color=severity))+geom_point()+facet_wrap('chip')

ggplot(data=data1, aes(x=heat_treat, y=p_mot, color=chytrid))+geom_point()
ggplot(data=data1, aes(x=heat_treat, y=p_fp, color=chytrid))+geom_point()
ggplot(data=data1, aes(x=heat_treat, y=p_vi, color=chytrid))+geom_point()

# Aim 1: Data Exploration
ggplot(data=data2, aes(x=chytrid, y=spconc, color=chip))+geom_point()
ggplot(data=data2, aes(x=chytrid, y=p_mot, color=chip))+geom_point()
ggplot(data=data2, aes(x=chytrid, y=p_fp, color=chip))+geom_point()
ggplot(data=data2, aes(x=chytrid, y=p_vi, color=chip))+geom_point()
 
# Aim 1: quantification results
ggplot(data=data2, aes(x=load, y=spconc, color=chip))+geom_point()+scale_x_log10()
ggplot(data=data2, aes(x=load, y=spconc, color=severity))+geom_point()+scale_x_log10()

# Aim 1: Categorised data exploration
ggplot(data=data2, aes(x=severity, y=spconc, color=chip))+geom_point()
ggplot(data=data2, aes(x=severity, y=p_mot, color=chip))+geom_point()
ggplot(data=data2, aes(x=severity, y=p_fp, color=chip))+geom_point()
ggplot(data=data2, aes(x=severity, y=p_vi, color=chip))+geom_point()

# Aim 2: Categorised data exploration
ggplot(data=data3, aes(x=chytrid, y=p_mot, color=chip))+geom_point()+facet_wrap('heat_treat')
ggplot(data=data3, aes(x=chytrid, y=p_fp, color=chip))+geom_point()+facet_wrap('heat_treat')
ggplot(data=data3, aes(x=chytrid, y=p_vi, color=chip))+geom_point()+facet_wrap('heat_treat')

# Aim 3 data exploration 
# All frogs
ggplot(data=data4, aes(x=heat_treat, y=spconc, color=chip))+geom_point()+facet_wrap('chip')
ggplot(data=data4, aes(x=heat_treat, y=spconc, color=chip))+geom_point()
# subset that has data in all timepoints
ggplot(data=data5, aes(x=heat_treat, y=spconc, color=chip))+geom_point()
ggplot(data=data5, aes(x=heat_treat, y=spconc, color=chip))+geom_point()+facet_wrap('chip')


#===========================================================================================================
# Model Development - sperm concentration
#===========================================================================================================

# Effect of Chytrid severity on sperm concentration - Aim 1
# pre heat-treat only (data2)

s1 <- glmmTMB(sperm ~ severity + offset(log_scale_ml_inv)  + (1|chip), data=data2, ziformula=~0, family=poisson)
summary(s1)
overdisp_fun(s1)
simulationOutput <- simulateResiduals(fittedModel = s1, plot = T ,n=1000)

s2 <- glmmTMB(sperm ~ severity + offset(log_scale_ml_inv)  + (1|bc), data=data2, ziformula=~0, family=poisson)
summary(s2)
overdisp_fun(s2)
simulationOutput <- simulateResiduals(fittedModel = s2, plot = T ,n=1000)

AIC(s1,s2)
drop1(s2,test="Chisq")

s2.g=ref_grid(s2)
emmeans(s2.g,"severity", type="response", offset=0)
emmip(s2, ~severity  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="Sperm Concentration (cells/mL)")

#Pairwise comparison of ratios
emmeans(s2.g,pairwise~severity,type="link")
pairs.sperm=emmeans(s2.g,pairwise~severity,type="response",adjust="none")
confint(pairs.sperm)

emmeans(s2.g,revpairwise~severity,type="link")
pairs.revsperm=emmeans(s2.g,revpairwise~severity,type="response",adjust="none")
confint(pairs.revsperm)

#==========================================================================================================================================
# Sperm concentration - heat treatment and chytrid status - Aim 2
# Without third time-point (data3)

P_ml <- glmmTMB(sperm~ chytrid + heat_treat + offset(log_scale_ml_inv)  + (1|chip), data=data3, ziformula=~0, family=poisson)
summary(P_ml)
overdisp_fun(P_ml)
simulationOutput <- simulateResiduals(fittedModel = P_ml, plot = T ,n=1000)

NB  <- glmmTMB(sperm~ chytrid + heat_treat+ offset(log_scale_ml_inv) + (1|chip), data=data3, ziformula=~0, family=nbinom2)
summary(NB)
overdisp_fun(NB)
simulationOutput <- simulateResiduals(fittedModel =NB , plot = T ,n=1000)

P_ml2 <- glmmTMB(sperm~ chytrid + heat_treat + offset(log_scale_ml_inv)  + (1|chip/rep), data=data3, ziformula=~0, family=poisson)
summary(P_ml2)
overdisp_fun(P_ml2)

P_ml3 <- glmmTMB(sperm~ chytrid + heat_treat + offset(log_scale_ml_inv)  + (1|bc), data=data3, ziformula=~0, family=poisson)
summary(P_ml3)
overdisp_fun(P_ml3)

AIC(P_ml, P_ml2, P_ml3)
drop1(P_ml3, test="Chisq")

# no interaction effect - final model = P_ml3
IP_ml3  <- glmmTMB(sperm~ chytrid + heat_treat + chytrid:heat_treat + offset(log_scale_ml_inv) + (1|bc), data=data3, ziformula=~0, family=poisson)

AIC(P_ml3, IP_ml3)
anova(P_ml3, IP_ml3)

emmip(P_ml3, chytrid~heat_treat  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="Sperm Concentration (cells/mL)")
simulationOutput <- simulateResiduals(fittedModel = P_ml3, plot = T ,n=1000)

# # Both look reasonable to me - Dharma has detected a problem, check traditional pearson residuals.
# # Extract Pearson residuals
pearson_resid <- residuals(P_ml3, type = "pearson")
#
# # Check the first few Pearson residuals
head(pearson_resid)

# # Plot Pearson residuals against fitted values
plot(fitted(P_ml3), pearson_resid,
     main = "Pearson Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Pearson Residuals")
abline(h = 0, col = "red", lty = 2)
#
# # Histogram of Pearson residuals
hist(pearson_resid, main = "Histogram of Pearson Residuals",
     xlab = "Pearson Residuals", col = "lightblue", border = "white")
# Extract Pearson residuals
pearson_resid <- residuals(NB2, type = "pearson")


ggplot(data=data3, aes(x=chytrid, y=spconc, color=chytrid))+geom_point()+facet_wrap('heat_treat')
emmip(P_ml3, chytrid~heat_treat  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="Sperm Concentration (cells/mL)")

# Odds Ratios and EMMs
P_ml3.g=ref_grid(P_ml3)
emmeans(P_ml3.g,"chytrid",by="heat_treat", type="response", offset=0)

confint(P_ml3.g,"chytrid",by="heat_treat",type="link")
confint(P_ml3.g,"chytrid",by="heat_treat",type="response")

#Pairwise comparison of ratios
emmeans(P_ml3.g,pairwise~chytrid,type="link")
pairs.sperm=emmeans(P_ml3.g,pairwise~chytrid,type="response",adjust="none")
confint(pairs.sperm)

emmeans(P_ml3.g,revpairwise~chytrid,type="link")
pairs.revsperm=emmeans(P_ml3.g,revpairwise~chytrid,type="response",adjust="none")
confint(pairs.revsperm)

emmeans(P_ml3.g,pairwise~heat_treat,type="link")
pairs.sperm=emmeans(P_ml3.g,pairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)

emmeans(P_ml3.g,revpairwise~heat_treat,type="link")
pairs.revsperm=emmeans(P_ml3.g,revpairwise~heat_treat,type="response",adjust="none")
confint(pairs.revsperm)

#=======================================================================================================================
# Effect of heat on chytrid postive frogs only - Aim 3
# (not enough replication to include chytrid negative frogs)
# positive frogs only (data 4)

t1 <- glmmTMB(sperm ~ heat_treat + offset(log_scale_ml_inv)  + (1|chip), data=data4, ziformula=~0, family=poisson)
summary(t1)
overdisp_fun(t1)
simulationOutput <- simulateResiduals(fittedModel = t1, plot = T ,n=1000)

t2 <- glmmTMB(sperm ~ heat_treat + offset(log_scale_ml_inv)  + (1|bc), data=data4, ziformula=~0, family=poisson)
summary(t2)
overdisp_fun(t2)
simulationOutput <- simulateResiduals(fittedModel = t2, plot = T ,n=1000)

AIC(t1,t2)
drop1(t2,test="Chisq")

t2.g=ref_grid(t2)
emmeans(t2.g,"heat_treat", type="response", offset=0)
emmip(t2, ~heat_treat  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="Sperm Concentration (cells/mL)")


#Pairwise comparison of ratios
emmeans(t2.g,pairwise~heat_treat,type="link")
pairs.sperm=emmeans(t2.g,pairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)

emmeans(t2.g,revpairwise~heat_treat,type="link")
pairs.revsperm=emmeans(t2.g,revpairwise~heat_treat,type="response",adjust="none")
confint(pairs.revsperm)

# Redo Analysis with data5 (only the four frogs in all three timepoints)
# no difference to overall trends and results

tt1 <- glmmTMB(sperm ~ heat_treat + offset(log_scale_ml_inv)  + (1|chip), data=data5, ziformula=~0, family=poisson)
summary(tt1)
overdisp_fun(tt1)
simulationOutput <- simulateResiduals(fittedModel = tt1, plot = T ,n=1000)

tt2 <- glmmTMB(sperm ~ heat_treat + offset(log_scale_ml_inv)  + (1|bc), data=data5, ziformula=~0, family=poisson)
summary(tt2)
overdisp_fun(tt2)
simulationOutput <- simulateResiduals(fittedModel = tt2, plot = T ,n=1000)

AIC(tt1,tt2)
drop1(tt2,test="Chisq")

tt2.g=ref_grid(tt2)
emmeans(tt2.g,"heat_treat", type="response", offset=0)
emmip(tt2, ~heat_treat  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="Sperm Concentration (cells/mL)")


#Pairwise comparison of ratios
emmeans(tt2.g,pairwise~heat_treat,type="link")
pairs.sperm=emmeans(tt2.g,pairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)

emmeans(tt2.g,revpairwise~heat_treat,type="link")
pairs.revsperm=emmeans(tt2.g,revpairwise~heat_treat,type="response",adjust="none")
confint(pairs.revsperm)



#====================================================================================================================
#Aim 1 - Effect of severity (data 2) forward progressive motility, total motility and vitality (membrane integrity)
#===============================

# forward progressive Motility
m1= glmer(p_fp ~ severity + (1|chip), weights = tot_mot, data=data2, family = binomial())
summary(m1)
plot(m1)
overdisp_fun((m1))
simulationOutput <- simulateResiduals(fittedModel = m1, plot = T ,n=1000)

m2= glmer(p_fp ~ severity + (1|chip/rep), weights = tot_mot, data=data2, family = binomial())
summary(m2)
plot(m2)
overdisp_fun((m2))

m3= glmmTMB(p_fp ~ severity + (1|bc), weights = tot_mot, data=data2, family = binomial())
summary(m3)
plot(m3)
overdisp_fun((m3))

AIC(m1,m3)
drop1(m3,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = m3, plot = T ,n=1000)

# Use m3
m3.g=ref_grid(m3)
emmeans(m3.g,"severity", type="response", offset=0)
emmip(m3, ~severity  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="Sperm Concentration (cells/mL)")
ggplot(data=data2, aes(x=severity, y=p_fp, color=chip))+geom_point()

#Pairwise comparison of ratios
emmeans(m3.g,revpairwise~severity,type="link")
pairs.revsperm=emmeans(m3.g,revpairwise~severity,type="response",adjust="none")
confint(pairs.revsperm)

# total Motility
m11= glmer(p_mot ~ severity + (1|chip), weights = tot_mot, data=data2, family = binomial())
summary(m11)
plot(m11)
overdisp_fun((m11))
simulationOutput <- simulateResiduals(fittedModel = m1, plot = T ,n=1000)

m12= glmer(p_mot ~ severity + (1|chip/rep), weights = tot_mot, data=data2, family = binomial())
summary(m12)
plot(m12)
overdisp_fun((m12))

m13= glmmTMB(p_mot ~ severity + (1|bc), weights = tot_mot, data=data2, family = binomial())
summary(m13)
plot(m13)
overdisp_fun((m13))

AIC(m11,m12,m13)
drop1(m13,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = m13, plot = T ,n=1000)

#Use m13
m13.g=ref_grid(m13)
emmeans(m13.g,"severity", type="response", offset=0)
emmip(m13, ~severity  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="Sperm Concentration (cells/mL)")
ggplot(data=data2, aes(x=severity, y=p_mot, color=chip))+geom_point()

#Pairwise comparison of ratios
emmeans(m13.g,revpairwise~severity,type="link")
pairs.revsperm=emmeans(m13.g,revpairwise~severity,type="response",adjust="none")
confint(pairs.revsperm)

# Vitality (MEMBRANE INTEGRITY)
m21= glmer(p_vi ~ severity + (1|chip), weights = vit_tot, data=data2, family = binomial())
summary(m21)
plot(m21)
overdisp_fun((m21))
simulationOutput <- simulateResiduals(fittedModel = m1, plot = T ,n=1000)

m22= glmer(p_vi ~ severity + (1|chip/rep), weights = vit_tot, data=data2, family = binomial())
summary(m22)
plot(m22)
overdisp_fun((m22))

m23= glmmTMB(p_vi ~ severity + (1|bc), weights = vit_tot, data=data2, family = binomial())
summary(m23)
plot(m23)
overdisp_fun((m23))

AIC(m21,m22,m23)
drop1(m23,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = m11, plot = T ,n=1000)

m23.g=ref_grid(m23)
emmeans(m23.g,"severity", type="response", offset=0)
emmip(m23, ~severity  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="Sperm Concentration (cells/mL)")
ggplot(data=data2, aes(x=severity, y=p_vi, color=chip))+geom_point()

#Pairwise comparison of ratios
emmeans(m23.g,revpairwise~severity,type="link")
pairs.revsperm=emmeans(m23.g,revpairwise~severity,type="response",adjust="none")
confint(pairs.revsperm)


#====================================================================================================================
#Aim 2 - Effect of pos/neg and heat_treat (data 3)
#===============================

# forward progressive Motility
mm1= glmer(p_fp ~ chytrid + heat_treat + (1|chip), weights = tot_mot, data=data3, family = binomial())
summary(mm1)
plot(m1)
overdisp_fun((mm1))

mm2= glmer(p_fp ~ chytrid + heat_treat 
           
           , weights = tot_mot, data=data3, family = binomial())
summary(mm2)
plot(mm2)
overdisp_fun((mm2))

mm3= glmmTMB(p_fp ~ chytrid + heat_treat + (1|bc), weights = tot_mot, data=data3, family = binomial())
summary(mm3)
plot(mm3)
overdisp_fun((mm3))

AIC(mm1,mm2,mm3)
drop1(mm3,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = mm3, plot = T ,n=1000)

#Check interaction - no interaction
mm4= glmmTMB(p_fp ~ chytrid + heat_treat+ chytrid:heat_treat + (1|bc), weights = tot_mot, data=data3, family = binomial())
drop1(mm4,test="Chisq")

# data check
ggplot(data=data3, aes(x=chytrid, y=p_fp, color=chytrid))+geom_point()+facet_grid(~heat_treat)
emmip(mm3, ~chytrid~heat_treat ,CIs=T,type="response",component="cond",offset=0) + labs(x="Chytrid",y="FP")

#EMMS
mm3.g=ref_grid(mm3)
emmeans(mm3.g,"heat_treat", by="chytrid", type="response", offset=0)
emmip(mm3, ~heat_treat~chytrid  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="FP")

#Pairwise comparison of ratios
emmeans(mm3.g,pairwise~heat_treat,type="link")
pairs.sperm=emmeans(mm3.g,pairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)

emmeans(mm3.g,pairwise~chytrid,type="link")
pairs.sperm=emmeans(mm3.g,pairwise~chytrid,type="response",adjust="none")
confint(pairs.sperm)


#MOTILITY
wm1= glmer(p_mot ~ chytrid + heat_treat + (1|chip), weights = tot_mot, data=data3, family = binomial())
summary(wm1)
plot(wm1)
overdisp_fun((wm1))

wm2= glmer(p_mot ~ chytrid + heat_treat + (1|chip/rep), weights = tot_mot, data=data3, family = binomial())
summary(wm2)
plot(wm2)
overdisp_fun((wm2))

wm3= glmmTMB(p_mot ~ chytrid + heat_treat + (1|bc), weights = tot_mot, data=data3, family = binomial())
summary(wm3)
plot(wm3)
overdisp_fun((wm3))

AIC(wm1,wm2,wm3)
drop1(wm3,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = wm3, plot = T ,n=1000)

#Check interaction - no interaction
wm4= glmmTMB(p_mot ~ chytrid + heat_treat+ chytrid:heat_treat + (1|bc), weights = tot_mot, data=data3, family = binomial())
drop1(wm4,test="Chisq")

# data check
ggplot(data=data3, aes(x=chytrid, y=p_mot, color=chytrid))+geom_point()+facet_grid(~heat_treat)
emmip(wm3, ~chytrid~heat_treat ,CIs=T,type="response",component="cond",offset=0) + labs(x="Chytrid",y="MOT")

#EMMS
wm3.g=ref_grid(wm3)
emmeans(wm3.g,"heat_treat", by="chytrid", type="response", offset=0)

#Pairwise comparison of ratios
emmeans(wm3.g,pairwise~heat_treat,type="link")
pairs.sperm=emmeans(wm3.g,pairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)

emmeans(wm3.g,pairwise~chytrid,type="link")
pairs.sperm=emmeans(wm3.g,pairwise~chytrid,type="response",adjust="none")
confint(pairs.sperm)


# MEMBRANE INTEGRITY
rm1= glmer(p_vi ~ chytrid + heat_treat + (1|chip), weights = vit_tot, data=data3, family = binomial())
summary(rm1)
plot(rm1)
overdisp_fun((rm1))

rm2= glmer(p_vi ~ chytrid + heat_treat + (1|chip/rep), weights = vit_tot, data=data3, family = binomial())
summary(rm2)
plot(rm2)
overdisp_fun((rm2))

rm3= glmmTMB(p_vi ~ chytrid + heat_treat + (1|bc), weights = vit_tot, data=data3, family = binomial())
summary(rm3)
plot(rm3)
overdisp_fun((wm3))

AIC(rm1,rm2,rm3)
drop1(rm3,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = rm3, plot = T ,n=1000)

#Check interaction - no interaction
rm4= glmmTMB(p_vi ~ chytrid + heat_treat+ chytrid:heat_treat + (1|bc), weights = vit_tot, data=data3, family = binomial())
drop1(rm4,test="Chisq")

# data check
ggplot(data=data3, aes(x=chytrid, y=p_vi, color=chytrid))+geom_point()+facet_grid(~heat_treat)
emmip(rm3, ~chytrid~heat_treat ,CIs=T,type="response",component="cond",offset=0) + labs(x="Chytrid",y="VIT")

#EMMS
rm3.g=ref_grid(rm3)
emmeans(rm3.g,"heat_treat", by="chytrid", type="response", offset=0)

#Pairwise comparison of ratios
emmeans(rm3.g,pairwise~heat_treat,type="link")
pairs.sperm=emmeans(rm3.g,pairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)

emmeans(rm3.g,pairwise~chytrid,type="link")
pairs.sperm=emmeans(rm3.g,pairwise~chytrid,type="response",adjust="none")
confint(pairs.sperm)

#====================================================================================================================
#Aim 3 - Effect of heat_treat (data 4)
#===============================

# forward progressive Motility
hh1= glmer(p_fp ~ heat_treat + (1|chip), weights = tot_mot, data=data4, family = binomial())
summary(hh1)
plot(hh1)
overdisp_fun((hh1))
simulationOutput <- simulateResiduals(fittedModel = m1, plot = T ,n=1000)

hh2= glmer(p_fp ~ heat_treat + (1|chip/rep), weights = tot_mot, data=data4, family = binomial())
summary(hh2)
plot(hh2)
overdisp_fun((hh2))

hh3= glmmTMB(p_fp ~ heat_treat + (1|bc), weights = tot_mot, data=data4, family = binomial())
summary(hh3)
plot(hh3)
overdisp_fun((hh3))

AIC(hh1,hh2,hh3)
drop1(hh3,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = hh3, plot = T ,n=1000)

#Use hh3
hh3.g=ref_grid(hh3)
emmeans(hh3.g,"heat_treat", type="response", offset=0)
emmip(hh3, ~heat_treat  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="FP")
ggplot(data=data4, aes(x=heat_treat, y=p_fp, color=chip))+geom_point()

#Pairwise comparison of ratios
emmeans(hh3.g,pairwise~heat_treat,type="link")
pairs.sperm=emmeans(hh3.g,pairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)

#Pairwise comparison of ratios
emmeans(hh3.g,revpairwise~heat_treat,type="link")
pairs.sperm=emmeans(hh3.g,revpairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)


# MOTILITY
th1= glmer(p_mot ~ heat_treat + (1|chip), weights = tot_mot, data=data4, family = binomial())
summary(th1)
plot(th1)
overdisp_fun((th1))
simulationOutput <- simulateResiduals(fittedModel = m1, plot = T ,n=1000)

th2= glmer(p_mot ~ heat_treat + (1|chip/rep), weights = tot_mot, data=data4, family = binomial())
summary(th2)
plot(th2)
overdisp_fun((th2))

th3= glmmTMB(p_mot ~ heat_treat + (1|bc), weights = tot_mot, data=data4, family = binomial())
summary(th3)
plot(th3)
overdisp_fun((th3))

AIC(th1,th2,th3)
drop1(th3,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = th3, plot = T ,n=1000)

th3.g=ref_grid(th3)
emmeans(th3.g,"heat_treat", type="response", offset=0)
emmip(th3, ~heat_treat  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="mot")
ggplot(data=data4, aes(x=heat_treat, y=p_mot, color=chip))+geom_point()

#Pairwise comparison of ratios
emmeans(th3.g,pairwise~heat_treat,type="link")
pairs.sperm=emmeans(th3.g,pairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)

# MEMBRANE INTEGRITY
vh1= glmer(p_vi ~ heat_treat + (1|chip), weights = vit_tot, data=data4, family = binomial())
summary(vh1)
plot(vh1)
overdisp_fun((vh1))
simulationOutput <- simulateResiduals(fittedModel = m1, plot = T ,n=1000)

vh2= glmer(p_vi ~ heat_treat + (1|chip/rep), weights = vit_tot, data=data4, family = binomial())
summary(vh2)
plot(vh2)
overdisp_fun((vh2))

vh3= glmmTMB(p_vi ~ heat_treat + (1|bc), weights = vit_tot, data=data4, family = binomial())
summary(vh3)
plot(vh3)
overdisp_fun((vh3))

AIC(vh1,vh2,vh3)
drop1(vh3,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = vh3, plot = T ,n=1000)

#Use vh3
vh3.g=ref_grid(vh3)
emmeans(vh3.g,"heat_treat", type="response", offset=0)
emmip(vh3, ~heat_treat  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="vi")
ggplot(data=data4, aes(x=heat_treat, y=p_vi, color=chip))+geom_point()

#Pairwise comparison of ratios
emmeans(vh3.g,revpairwise~heat_treat,type="link")
pairs.sperm=emmeans(vh3.g,revpairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)

#====================================================================================================================
#Aim 3 - Reanalyse Effect of heat_treat (data 5)
# No major change to trends and results - coauthor preferences are to use this analysis
#===============================

# forward progressive Motility
hh1= glmer(p_fp ~ heat_treat + (1|chip), weights = tot_mot, data=data5, family = binomial())
summary(hh1)
plot(hh1)
overdisp_fun((hh1))
#simulationOutput <- simulateResiduals(fittedModel = m1, plot = T ,n=1000)

hh2= glmer(p_fp ~ heat_treat + (1|chip/rep), weights = tot_mot, data=data5, family = binomial())
summary(hh2)
plot(hh2)
overdisp_fun((hh2))

hh3= glmmTMB(p_fp ~ heat_treat + (1|bc), weights = tot_mot, data=data5, family = binomial())
summary(hh3)
plot(hh3)
overdisp_fun((hh3))

AIC(hh1,hh2,hh3)
drop1(hh3,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = hh3, plot = T ,n=1000)

hh3.g=ref_grid(hh3)
emmeans(hh3.g,"heat_treat", type="response", offset=0)
emmip(hh3, ~heat_treat  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="FP")
ggplot(data=data4, aes(x=heat_treat, y=p_fp, color=chip))+geom_point()

#Pairwise comparison of ratios
emmeans(hh3.g,pairwise~heat_treat,type="link")
pairs.sperm=emmeans(hh3.g,pairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)


# MOTILITY
th1= glmer(p_mot ~ heat_treat + (1|chip), weights = tot_mot, data=data5, family = binomial())
summary(th1)
plot(th1)
overdisp_fun((th1))
#simulationOutput <- simulateResiduals(fittedModel = m1, plot = T ,n=1000)

th2= glmer(p_mot ~ heat_treat + (1|chip/rep), weights = tot_mot, data=data5, family = binomial())
summary(th2)
plot(th2)
overdisp_fun((th2))

th3= glmmTMB(p_mot ~ heat_treat + (1|bc), weights = tot_mot, data=data5, family = binomial())
summary(th3)
plot(th3)
overdisp_fun((th3))

AIC(th1,th2,th3)
drop1(th3,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = th3, plot = T ,n=1000)

th3.g=ref_grid(th3)
emmeans(th3.g,"heat_treat", type="response", offset=0)
emmip(th3, ~heat_treat  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="mot")
ggplot(data=data4, aes(x=heat_treat, y=p_mot, color=chip))+geom_point()

#Pairwise comparison of ratios
emmeans(th3.g,pairwise~heat_treat,type="link")
pairs.sperm=emmeans(th3.g,pairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)

# MEMBRANE INTEGRITY
vh1= glmer(p_vi ~ heat_treat + (1|chip), weights = vit_tot, data=data5, family = binomial())
summary(vh1)
plot(vh1)
overdisp_fun((vh1))
#simulationOutput <- simulateResiduals(fittedModel = m1, plot = T ,n=1000)

vh2= glmer(p_vi ~ heat_treat + (1|chip/rep), weights = vit_tot, data=data5, family = binomial())
summary(vh2)
plot(vh2)
overdisp_fun((vh2))

vh3= glmmTMB(p_vi ~ heat_treat + (1|bc), weights = vit_tot, data=data5, family = binomial())
summary(vh3)
plot(vh3)
overdisp_fun((vh3))

AIC(vh1,vh2,vh3)
drop1(vh3,test="Chisq")
simulationOutput <- simulateResiduals(fittedModel = vh3, plot = T ,n=1000)

vh3.g=ref_grid(vh3)
emmeans(vh3.g,"heat_treat", type="response", offset=0)
emmip(vh3, ~heat_treat  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="vi")
ggplot(data=data4, aes(x=heat_treat, y=p_vi, color=chip))+geom_point()

#Pairwise comparison of ratios
emmeans(vh3.g,revpairwise~heat_treat,type="link")
pairs.sperm=emmeans(vh3.g,revpairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)

emmeans(vh3.g,pairwise~heat_treat,type="link")
pairs.sperm=emmeans(vh3.g,pairwise~heat_treat,type="response",adjust="none")
confint(pairs.sperm)
#==========================================================================================================================
# Graphing for publication
#==========================================================================================================================
# Aim 1
#--------------------------------------------------------------------------------------------------------------------------
#spconc
s2 <- glmmTMB(sperm ~ severity + offset(log_scale_ml_inv)  + (1|bc), data=data2, ziformula=~0, family=poisson)

s2.g=ref_grid(s2)
emmeans(s2.g,"severity",type="response",offset=0)
spconc=emmeans(s2.g,"severity",type="response",offset=0)
df_spconc=as.data.frame(summary(spconc))[c('severity', 'rate', 'asymp.LCL', 'asymp.UCL')]

data2$severity <- factor(data2$severity, levels = c("Negative", "Moderate", "High"))
df_spconc$severity <- factor(df_spconc$severity, levels = c("Negative", "Moderate", "High"))

a=ggplot() +
  geom_col(data=df_spconc, aes(x=severity, y=rate), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=df_spconc, aes(x=severity, ymin=asymp.LCL, ymax=asymp.UCL), width=0.1,position= position_dodge(0.7))+
  geom_point(data=data3, aes(x=severity, y=spconc, fill=severity), position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Infection Severity")+ ylab("Sperm Concentration (x10⁸ Cells/mL)")+ 
  theme_classic()+
  scale_x_discrete(limits = c("Negative", "Moderate", "High"), labels = c("Negative", "Moderate", "High"))+
  scale_y_continuous(labels = scales::number_format(scale = 1e-8, accuracy = 0.1))+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank())+
  ggtitle("(a)")+ theme(plot.title = element_text(size = 11))

#position= position_dodge(0.7)

#fp
m3= glmmTMB(p_fp ~ severity + (1|bc), weights = tot_mot, data=data2, family = binomial())

m3.g=ref_grid(m3)
emmeans(m3.g,"severity",type="response",offset=0)
fp=emmeans(m3.g,"severity",type="response",offset=0)
df_fp=as.data.frame(summary(fp))[c('severity', 'prob', 'asymp.LCL', 'asymp.UCL')]

data2$severity <- factor(data2$severity, levels = c("Negative", "Moderate", "High"))
df_fp$severity <- factor(df_fp$severity, levels = c("Negative", "Moderate", "High"))

b=ggplot() +
  geom_col(data=df_fp, aes(x=severity, y=prob), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=df_fp, aes(x=severity, ymin=asymp.LCL, ymax=asymp.UCL), width=0.1,position=position_dodge(0.7))+
  geom_point(data=data2, aes(x=severity, y=p_fp, fill=severity), position= position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Infection Severity")+ ylab("Forward-progressive Sperm (%)")+ 
  theme_classic()+
  scale_fill_grey(start=0.4, end=0.8)+
  scale_x_discrete(limits = c("Negative", "Moderate", "High"), labels = c("Negative", "Moderate", "High"))+
  scale_y_continuous(labels = scales::number_format(scale = 100, accuracy = 1), limits = c(0, 1))+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank()) +
  ggtitle("(b)")+ theme(plot.title = element_text(size = 11))

#mot
m13= glmmTMB(p_mot ~ severity + (1|bc), weights = tot_mot, data=data2, family = binomial())

m13.g=ref_grid(m13)
emmeans(m13.g,"severity",type="response",offset=0)
mot=emmeans(m13.g,"severity",type="response",offset=0)
df_mot=as.data.frame(summary(mot))[c('severity', 'prob', 'asymp.LCL', 'asymp.UCL')]

data2$severity <- factor(data2$severity, levels = c("Negative", "Moderate", "High"))
df_mot$severity <- factor(df_mot$severity, levels = c("Negative", "Moderate", "High"))

c=ggplot() +
  geom_col(data=df_mot, aes(x=severity, y=prob), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=df_mot, aes(x=severity, ymin=asymp.LCL, ymax=asymp.UCL), width=0.1,position=position_dodge(0.7))+
  geom_point(data=data2, aes(x=severity, y=p_mot, fill=severity), position= position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Infection Severity")+ ylab("Total Sperm Motility (%)")+ 
  theme_classic()+
  scale_fill_grey(start=0.4, end=0.8)+
  scale_x_discrete(limits = c("Negative", "Moderate", "High"), labels = c("Negative", "Moderate", "High"))+
  scale_y_continuous(labels = scales::number_format(scale = 100, accuracy = 1), limits = c(0, 1)) +
  theme(legend.position="none")+
  ggtitle("(c)")+ theme(plot.title = element_text(size = 11))

#vit
m23= glmmTMB(p_vi ~ severity + (1|bc), weights = vit_tot, data=data2, family = binomial())

m23.g=ref_grid(m23)
emmeans(m23.g,"severity",type="response",offset=0)
vit=emmeans(m23.g,"severity",type="response",offset=0)
df_vit=as.data.frame(summary(vit))[c('severity', 'prob', 'asymp.LCL', 'asymp.UCL')]

data2$severity <- factor(data2$severity, levels = c("Negative", "Moderate", "High"))
df_vit$severity <- factor(df_vit$severity, levels = c("Negative", "Moderate", "High"))

d=ggplot() +
  geom_col(data=df_vit, aes(x=severity, y=prob), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=df_vit, aes(x=severity, ymin=asymp.LCL, ymax=asymp.UCL), width=0.1,position=position_dodge(0.7))+
  geom_point(data=data2, aes(x=severity, y=p_vi, fill=severity), position= position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Infection Severity")+ ylab("Intact Sperm Membranes (%)")+ 
  theme_classic()+
  scale_fill_grey(start=0.4, end=0.8)+
  scale_x_discrete(limits = c("Negative", "Moderate", "High"), labels = c("Negative", "Moderate", "High"))+
  scale_y_continuous(labels = scales::number_format(scale = 100, accuracy = 1), limits = c(0, 1))+
  theme(legend.position="none")+
  ggtitle("(d)")+ theme(plot.title = element_text(size = 11))

fig1=grid.arrange(a, b, c, d, ncol=2,nrow=2,widths = c(2,2))

ggsave("Fig1.tiff", fig1, units="cm", width=20, height=20, dpi=600, compression = 'lzw')

#--------------------------------------------------------------------------------------------------------------------------
# Aim 2
#--------------------------------------------------------------------------------------------------------------------------
#spconc
P_ml3 <- glmmTMB(sperm~ chytrid + heat_treat + offset(log_scale_ml_inv)  + (1|bc), data=data3, ziformula=~0, family=poisson)

P_ml3.g=ref_grid(P_ml3)
emmeans(P_ml3.g,"chytrid",by="heat_treat",type="response",offset=0)
spconc=emmeans(P_ml3.g,"chytrid",by="heat_treat",type="response",offset=0)
df_spconc=as.data.frame(summary(spconc))[c('heat_treat', 'chytrid', 'rate', 'asymp.LCL', 'asymp.UCL')]

data3$heat_treat <- factor(data3$heat_treat, levels = c("pre", "post"))
df_spconc$heat_treat <- factor(df_spconc$heat_treat, levels = c("pre", "post"))

a2=ggplot() +
  geom_col(data=df_spconc, aes(x=heat_treat, y=rate, fill=chytrid), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=df_spconc, aes(x=heat_treat, ymin=asymp.LCL, ymax=asymp.UCL,group=chytrid), width=0.1,position= position_dodge(0.7))+
  geom_point(data=data3, aes(x=heat_treat, y=spconc, fill=chytrid,group=chytrid), position= position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Time")+ ylab("Sperm Concentration (x10⁸ Cells/mL)")+ 
  theme_classic()+
  scale_fill_grey(start=0.4, end=0.8,name = "Initial Chytrid Status", labels = c("Negative", "Positive"))+
  scale_x_discrete(limits = c("pre", "post"), labels = c("Pre-treatment", "Post-treatment"))+
  scale_y_continuous(labels = scales::number_format(scale = 1e-8, accuracy = 0.1))+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank())+
  ggtitle("(a)")+ theme(plot.title = element_text(size = 11))


#fp
mm3= glmmTMB(p_fp ~ chytrid + heat_treat + (1|bc), weights = tot_mot, data=data3, family = binomial())

mm3.g=ref_grid(mm3)
emmeans(mm3.g,"chytrid",by="heat_treat",type="response")
fp=emmeans(mm3.g,"chytrid",by="heat_treat",type="response")
df_fp=as.data.frame(summary(fp))[c('heat_treat', 'chytrid', 'prob', 'asymp.LCL', 'asymp.UCL')]

data3$heat_treat <- factor(data3$heat_treat, levels = c("pre", "post"))
df_fp$heat_treat <- factor(df_fp$heat_treat, levels = c("pre", "post"))

b2=ggplot() +
  geom_col(data=df_fp, aes(x=heat_treat, y=prob, fill=chytrid), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=df_fp, aes(x=heat_treat, ymin=asymp.LCL, ymax=asymp.UCL,group=chytrid), width=0.1,position= position_dodge(0.7))+
  geom_point(data=data3, aes(x=heat_treat, y=p_fp, fill=chytrid,group=chytrid), position= position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Time")+ ylab("Forward-progressive Sperm (%)")+ 
  theme_classic()+
  scale_fill_grey(start=0.4, end=0.8,name = "Initial Chytrid Status", labels = c("Negative", "Positive"))+
  scale_x_discrete(limits = c("pre", "post"), labels = c("Pre-treatment", "Post-treatment"))+
  scale_y_continuous(labels = scales::number_format(scale = 100, accuracy = 1), limits = c(0, 1))+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank()) +
  ggtitle("(b)")+ theme(plot.title = element_text(size = 11))

#mot
wm3= glmmTMB(p_mot ~ chytrid + heat_treat + (1|bc), weights = tot_mot, data=data3, family = binomial())

wm3.g=ref_grid(wm3)
emmeans(wm3.g,"chytrid",by="heat_treat",type="response")
mot=emmeans(wm3.g,"chytrid",by="heat_treat",type="response")
df_mot=as.data.frame(summary(mot))[c('heat_treat', 'chytrid', 'prob', 'asymp.LCL', 'asymp.UCL')]

data3$heat_treat <- factor(data3$heat_treat, levels = c("pre", "post"))
df_mot$heat_treat <- factor(df_mot$heat_treat, levels = c("pre", "post"))

c2=ggplot() +
  geom_col(data=df_mot, aes(x=heat_treat, y=prob, fill=chytrid), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=df_mot, aes(x=heat_treat, ymin=asymp.LCL, ymax=asymp.UCL,group=chytrid), width=0.1,position= position_dodge(0.7))+
  geom_point(data=data3, aes(x=heat_treat, y=p_mot, fill=chytrid,group=chytrid), position= position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Time")+ ylab("Total Motility (%)")+ 
  theme_classic()+
  scale_fill_grey(start=0.4, end=0.8,name = "Initial Chytrid Status", labels = c("Negative", "Positive"))+
  scale_x_discrete(limits = c("pre", "post"), labels = c("Pre-treatment", "Post-treatment"))+
  scale_y_continuous(labels = scales::number_format(scale = 100, accuracy = 1), limits = c(0, 1))+
  theme(legend.position="none")+
  ggtitle("(c)")+ theme(plot.title = element_text(size = 11))

# membrane integrity
rm3= glmmTMB(p_vi ~ chytrid + heat_treat + (1|bc), weights = vit_tot, data=data3, family = binomial())

rm3.g=ref_grid(rm3)
emmeans(rm3.g,"chytrid",by="heat_treat",type="response")
vit=emmeans(rm3.g,"chytrid",by="heat_treat",type="response")
df_vit=as.data.frame(summary(vit))[c('heat_treat', 'chytrid', 'prob', 'asymp.LCL', 'asymp.UCL')]

data3$heat_treat <- factor(data3$heat_treat, levels = c("pre", "post"))
df_vit$heat_treat <- factor(df_vit$heat_treat, levels = c("pre", "post"))

d2_leg=ggplot() +
  geom_col(data=df_vit, aes(x=heat_treat, y=prob, fill=chytrid), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=df_vit, aes(x=heat_treat, ymin=asymp.LCL, ymax=asymp.UCL,group=chytrid), width=0.1,position= position_dodge(0.7))+
  geom_point(data=data3, aes(x=heat_treat, y=p_vi, fill=chytrid,group=chytrid), position= position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Time")+ ylab("Intact Sperm Membranes (%)")+ 
  theme_classic()+
  scale_fill_grey(start=0.4, end=0.8,name = "Initial Chytrid Status", labels = c("Negative", "Positive"))+
  scale_x_discrete(limits = c("pre", "post"), labels = c("Pre-treatment", "Post-treatment"))+
  scale_y_continuous(labels = scales::number_format(scale = 100, accuracy = 1), limits = c(0, 1))+
  ggtitle("(d)")+ theme(plot.title = element_text(size = 11))
  
get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
legend <- get_legend(d2_leg)
d2= d2_leg + theme(legend.position="none")

fig2=grid.arrange(a2, b2, legend, c2, d2, ncol=3,nrow=2,widths = c(2,2,1))

ggsave("Fig2.tiff", fig2, units="cm", width=25, height=20, dpi=600, compression = 'lzw')

#--------------------------------------------------------------------------------------------------------------------------
# Aim 3
#--------------------------------------------------------------------------------------------------------------------------
#spconc
tt2 <- glmmTMB(sperm ~ heat_treat + offset(log_scale_ml_inv)  + (1|bc), data=data5, ziformula=~0, family=poisson)

tt2.g=ref_grid(t2)
emmeans(tt2.g,"heat_treat",type="response",offset=0)
a3spconc=emmeans(tt2.g,"heat_treat",type="response",offset=0)
a3df_spconc=as.data.frame(summary(a3spconc))[c('heat_treat', 'rate', 'asymp.LCL', 'asymp.UCL')]

data5$heat_treat <- factor(data5$heat_treat, levels = c("pre", "post", "post2"))
a3df_spconc$heat_treat <- factor(a3df_spconc$heat_treat, levels = c("pre", "post", "post2"))

a3=ggplot() +
  geom_col(data=a3df_spconc, aes(x=heat_treat, y=rate), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=a3df_spconc, aes(x=heat_treat, ymin=asymp.LCL, ymax=asymp.UCL), width=0.1,position= position_dodge(0.7))+
  geom_point(data=data5, aes(x=heat_treat, y=spconc, fill=heat_treat), position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Time")+ ylab("Sperm Concentration (x10⁸ Cells/mL)")+ 
  theme_classic()+
  scale_x_discrete(limits = c("pre", "post", "post2"), labels = c("Pre-treatment", "Post-treatment", "Six months post"))+
  scale_y_continuous(labels = scales::number_format(scale = 1e-8, accuracy = 0.1))+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank())+
  ggtitle("(a)")+ theme(plot.title = element_text(size = 11))

#fp
hh3= glmmTMB(p_fp ~ heat_treat + (1|bc), weights = tot_mot, data=data5, family = binomial())

hh3.g=ref_grid(hh3)
emmeans(hh3.g,"heat_treat",type="response",offset=0)
a3fp=emmeans(hh3.g,"heat_treat",type="response",offset=0)
a3df_fp=as.data.frame(summary(a3fp))[c('heat_treat', 'prob', 'asymp.LCL', 'asymp.UCL')]

data5$heat_treat <- factor(data5$heat_treat, levels = c("pre", "post", "post2"))
a3df_fp$heat_treat <- factor(a3df_fp$heat_treat, levels = c("pre", "post", "post2"))

b3=ggplot() +
  geom_col(data=a3df_fp, aes(x=heat_treat, y=prob), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=a3df_fp, aes(x=heat_treat, ymin=asymp.LCL, ymax=asymp.UCL), width=0.1,position=position_dodge(0.7))+
  geom_point(data=data5, aes(x=heat_treat, y=p_fp, fill=heat_treat), position= position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Infection Severity")+ ylab("Forward-progressive Sperm (%)")+
  theme_classic()+
  scale_fill_grey(start=0.4, end=0.8)+
  scale_x_discrete(limits = c("pre", "post", "post2"), labels = c("Pre-treatment", "Post-treatment", "Six months post"))+
  scale_y_continuous(labels = scales::number_format(scale = 100, accuracy = 1), limits = c(0, 1))+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank()) +
  ggtitle("(b)")+ theme(plot.title = element_text(size = 11))

#mot
th3= glmmTMB(p_mot ~ heat_treat + (1|bc), weights = tot_mot, data=data5, family = binomial())

th3.g=ref_grid(th3)
emmeans(th3.g,"heat_treat",type="response",offset=0)
a3mot=emmeans(th3.g,"heat_treat",type="response",offset=0)
a3df_mot=as.data.frame(summary(a3mot))[c('heat_treat', 'prob', 'asymp.LCL', 'asymp.UCL')]

data5$heat_treat <- factor(data5$heat_treat, levels = c("pre", "post", "post2"))
a3df_mot$heat_treat <- factor(a3df_fp$heat_treat, levels = c("pre", "post", "post2"))

c3=ggplot() +
  geom_col(data=a3df_mot, aes(x=heat_treat, y=prob), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=a3df_mot, aes(x=heat_treat, ymin=asymp.LCL, ymax=asymp.UCL), width=0.1,position=position_dodge(0.7))+
  geom_point(data=data5, aes(x=heat_treat, y=p_mot, fill=heat_treat), position= position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Time")+ ylab("Total Sperm Motiliy (%)")+
  theme_classic()+
  scale_fill_grey(start=0.4, end=0.8)+
  scale_x_discrete(limits = c("pre", "post", "post2"), labels = c("Pre-treatment", "Post-treatment", "Six months post"))+
  scale_y_continuous(labels = scales::number_format(scale = 100, accuracy = 1), limits = c(0, 1))+
  theme(legend.position="none")+
  ggtitle("(c)")+ theme(plot.title = element_text(size = 11))

#vit
vh3= glmmTMB(p_vi ~ heat_treat + (1|bc), weights = vit_tot, data=data5, family = binomial())

vh3.g=ref_grid(vh3)
emmeans(vh3.g,"heat_treat",type="response",offset=0)
a3vit=emmeans(vh3.g,"heat_treat",type="response",offset=0)
a3df_vit=as.data.frame(summary(a3vit))[c('heat_treat', 'prob', 'asymp.LCL', 'asymp.UCL')]

data5$heat_treat <- factor(data5$heat_treat, levels = c("pre", "post", "post2"))
a3df_vit$heat_treat <- factor(a3df_vit$heat_treat, levels = c("pre", "post", "post2"))

d3_leg=ggplot() +
  geom_col(data=a3df_vit, aes(x=heat_treat, y=prob), position = position_dodge(0.7), width=0.5) +
  geom_errorbar(data=a3df_vit, aes(x=heat_treat, ymin=asymp.LCL, ymax=asymp.UCL), width=0.1,position=position_dodge(0.7))+
  geom_point(data=data5, aes(x=heat_treat, y=p_vi, fill=heat_treat), position= position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1)+
  xlab("Time")+ ylab("Intact Sperm Membranes (%)")+
  theme_classic()+
  scale_fill_grey(start=0.4, end=0.8)+
  theme(legend.position="none")+
  scale_x_discrete(limits = c("pre", "post", "post2"), labels = c("Pre-treatment", "Post-treatment", "Six months post"))+
  scale_y_continuous(labels = scales::number_format(scale = 100, accuracy = 1), limits = c(0, 1))+
  ggtitle("(d)")+ theme(plot.title = element_text(size = 11))



fig3_v2=grid.arrange(a3, b3, c3, d3_leg, ncol=2,nrow=2,widths = c(2,2))

ggsave("Fig3_v2.tiff", fig3_v2, units="cm", width=20, height=20, dpi=600, compression = 'lzw')