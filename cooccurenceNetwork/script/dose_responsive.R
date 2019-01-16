## This scripot is meant to draw the responsive curve with the facility of  Linear regression  Model(LM)


## Tank 

## Email:xnzhang@genetics.ac.cn

rm(list = ls())
setwd("/mnt/bai/xiaoning/xiaoxuan/180213/fun_New/")
print(paste("Your working directory is in",getwd()))

# install.packages("visreg")
# install.packages("sjPlot")
# devtools::install_github("strengejacke/strengejacke")
library(visreg)
library(MASS)
library(car)
library(sjPlot)
library(sjmisc)
########### to import the plotting theme() function ########################
source("script/plot_function.R")

### output directory assigned to include the pics & tables########################
# figures.dir <- paste("/mnt/bai/xiaoning/xiaoxuan/180213/bac_New/figure/",tem,'/',sep = '')
# table.dir <- paste("/mnt/bai/xiaoning/xiaoxuan/180213/bac_New/table/",tem,'/',sep = '')


design = read.table("doc/design.txt", header=T, row.names= 1, sep="\t") 

## the applicaiton of subset & !  for short without the $
sub_design <- subset(design,PlasmidID == "Scal-BI-12-4" & Other == "2:02:02"& !Description %in% "E00")



l1 = read.delim("data/usearch_map_L1/observation_table.txt", row.names= 1,  header=T, sep="\t")
l2 = read.delim("data/usearch_map_L2/observation_table.txt", row.names= 1,  header=T, sep="\t")
l3 = read.delim("data/usearch_map_L3/observation_table.txt", row.names= 1,  header=T, sep="\t")
l4 = read.delim("data/usearch_map_L4/observation_table.txt", row.names= 1,  header=T, sep="\t")

## merge and cbind difference 
## merge can only handle two data.frame 
l1$ID <- rownames(l1)
l2$ID <- rownames(l2)
l3$ID <- rownames(l3)
l4$ID <- rownames(l4)

merge_sp12 <- merge(l1,l2, by = "ID")
merge_sp34 <- merge(l3,l4, by = "ID")
merge_sp <- merge(merge_sp12,merge_sp34,by = "ID")
rownames(merge_sp) <- merge_sp$ID
merge_sp$ID <- NULL

## select usage in subset function  column select; row not to use select
otu_table <- subset(merge_sp, select = colnames(merge_sp)%in% rownames(sub_design))

sub_design <- subset(sub_design,rownames(sub_design) %in% colnames(otu_table) )
sub_design$ID <- row.names(sub_design)

otu_tr <- as.data.frame(t(otu_table))
otu_tr$ID <- row.names(otu_tr)
idx <- grep("BI-OS-12-4",colnames(otu_tr))
dat <- merge(sub_design,otu_tr,by = "ID")

dat$spike_read <- dat$`BI-OS-12-4`


## poission GLM

# dat=data.frame(
#   
#   cp_X=as.numeric(c(19.60204858,17.28012048,16.28012048,15.28012048)),
#   RA_Y=as.numeric(c(5.680165804,4.130875576,3.798621485,2.962170828))
#   
#   
# )

## dose responsive 

fit <- lm(spike_read ~  spike_concentration,data=dat)
visreg(fit = fit)

# summary(fit)





#offset 已知？
# fit<- dat.glm <- glm(`BI-OS-12-4` ~  log(spike_concentration),data=dat,family=poisson(link=log))
# fit<- dat.glm <- glm.nb(`BI-OS-12-4` ~  log(spike_concentration),data=dat,link=log)
#                      
#                      
# confint(fit, level = 0.95) # approximate 95% confidence intervals
# dose.p(fit, p = 0.50)      # LD50 for a dose-response curve
# 
# 
# 
# dat.pre=predict(dat.glm,type="response")
# 
# sjp.glm(fit = fit)
# sjp.frq(dat$`BI-OS-12-4`)
# 
# 
# sjp.glm(fit,type="pred")
# plot_model(fit)
# 
# ## I have used the predict() and lines() function to plot this model output:
# 
# plot(dat$spike_concentration,  dat$`BI-OS-12-4`,  xlab = "spike_concentration", ylab = "Spike_in_Reads", col=alpha("black",.35), font = 2, font.lab = 2)+abline(a=7.625171,b=1.986406e-06)
# # lines(dat$spike_concentration, dat.pre,lwd=4, col = "red")
# 
# layout(1)  #取消绘图区域分割
# plot(dat$spike_concentration,dat.pre,xlab='观测值',ylab='拟合值',main="poisson 拟合效果",pch="*")  #添加直线y=x，截距为0，斜率为1
# 
# plot(dat$`BI-OS-12-4`,exp(dat.pre),xlab='观测值',ylab='拟合值',main="poisson 拟合效果",pch="*")+abline(0,1)  #添加直线y=x，截距为0，斜率为1
# 
# 
# car::qqPlot(fit)
# 
# 
# 
# 
# 
# library(MASS)
# attach(dat)
# dat.glmnb <- glm.nb(RA_Y~log(cp_X),data=dat,link=log)  #负二项回归
# 
# 
# 
# 
# 
# 
# ####  Test
# 
# data(efc)
# # set basic theme options
# set_theme("forest",
#           axis.title.size = .85, 
#           axis.textsize = .85, 
#           legend.size = .8, 
#           geom.label.size = 3.5)
# 
# 
# fit2 <- glm(tot_sc_e ~ neg_c_7 + e42dep + c161sex,
#             data = efc, family = poisson(link = "log"))
# # fit incident rate ratios, we need three decimal points to see 
# # a difference to the negative binomial model...
# sjp.glm(fit2, digits = 3)
# 
# 
# fit3 <- glm.nb(tot_sc_e ~ neg_c_7 + e42dep + c161sex, data = efc)
# # fit incident rate ratios
# sjp.glm(fit3, digits = 3)
# 
# sjp.glm(fit, type = "eff")
# 

