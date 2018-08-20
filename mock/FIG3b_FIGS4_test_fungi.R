# originally by Tank

# Email:xnzhang@genetics.ac.cn

rm(list=ls())
# install.packages("DescTools")
library(DescTools)
library(dplyr)
library(vegan)
library(reshape2)
########### setting the working directory and print it ###################
tem = "final_test"
setwd("~/xiaoxuan/final/mock/")
print(paste("Your working directory is in",getwd()))


figures.dir <- paste("~/xiaoxuan/final/mock/fun/",tem,'/',sep = '')
table.dir <- paste("~/xiaoxuan/final/mock/table/",tem,'/',sep = '')

fig_flag <- dir.exists(figures.dir)
if( isTRUE(!fig_flag)){
  dir.create(figures.dir)
}

tab_flag <- dir.exists(table.dir)
if( isTRUE(!tab_flag)){
  dir.create(table.dir)
}
####################################


#####spike-in design in this batch #####################
spike <- c("BI-OS-11-3","BI-OS-12-4","BI-OS-10-2")


#1 design Mapping file and Sample id
design = read.table("fun/design.txt", header=T, row.names= 1, sep="\t") 
# design$SampleID <- row.names(design)
# sample_list <- as.matrix(row.names(design))

otu_table = read.table("fungi_RA_withoutSpike.txt", row.names= 1,  header=1, sep="\t")
# otu_table_num=as.data.frame(lapply(otu_table, function(x) as.numeric(gsub("\\%", "", x))/100))
# row.names(otu_table_num)=row.names(otu_table)
otu_table_t= as.data.frame(t(otu_table))
# otu_table_t<-data.frame(lapply(otu_table_t, function(x) as.numeric(gsub("\\%", "", x))/100))
otu_table_meta=merge(design,otu_table_t,by="row.names")
colnames(otu_table_meta)[1]= "SampleID"

dat=melt(otu_table_meta,id.vars = c("SampleID", "Other","Description"),measure.vars = row.names(otu_table))
# dat_s=dat[,c(3,5)]

# DunnettTest(value ~ Description, data = dat[,c(3:5)],control ="E00", conf.level=0.9)
# boxplot(Ozone ~ Month, data = airquality)

#################补充FIG 2b 和FIG S1中，3种比例条件下，去除spike-in的菌株RA 和E00的菌株RA 有无显著差别

vec1=row.names(otu_table)
vec2=unique(as.vector(dat$Other))
library("multcomp")

for (i in 1:length(vec1)){
  for (j in 1:length(vec2)){
    
    # data("recovery", package = "multcomp")
    RA.aov <- aov(value ~ Description, data = dat[dat$variable %in% row.names(otu_table)[i]& dat$Other %in% vec2[j],])
    RA.mc <- glht(RA.aov,
                  linfct = mcp(Description = "Dunnett"),
                  alternative = "less")
    
    print(paste0("dunnet-test for ",vec1[i]," ",vec2[j]," ratio  is constructing"))
    print(summary(RA.mc))
  }
}

#### to continue to rewrite with Rmarkdown 



#### 补充不同spike-in梯度下9株菌在1:1:1 VS 2:2:2时RA值有无统计学差异
dat_21=dat[dat$Other %in% c("1:1:1","2:2:2"), ]
vec3=row.names(otu_table)
vec4=unique(as.vector(dat_21$Description))

for (i in 1:length(vec3)){
  for (j in 1:length(vec4)){
    
    # data("recovery", package = "multcomp")
    RA.t <- t.test(value ~ Other, data = dat_21[dat_21$variable %in% vec3[i]& dat_21$Description%in% vec4[j],])
    # RA.mc <- glht(RA.t,
    #               linfct = mcp(Description = "t.test"),
    #               alternative = "less")
    
    
    print(paste0("t-test for ",vec3[i]," ",vec4[j]," ratio for 2a:2b:2c over 1a:1b:1c is constructing"))
    print(RA.t)
    
  }
}


######### 不同spike-in梯度下Asco的AA 在1:1:1 VS 2:2:1  时有无统计差异，理想情况是AA 无统计学差异

otu_table = read.table("fungi_AA_mock.txt", row.names= 1,  header=1, sep="\t")
# row.names(otu_table_num)=row.names(otu_table)
otu_table_t= as.data.frame(t(otu_table))
# otu_table_t<-data.frame(lapply(otu_table_t, function(x) as.numeric(gsub("\\%", "", x))/100))
otu_table_meta=merge(design,otu_table_t,by="row.names")
colnames(otu_table_meta)[1]= "SampleID"

dat=melt(otu_table_meta,id.vars = c("SampleID", "Other","Description"),measure.vars = row.names(otu_table))


dat_21=dat[dat$Other %in% c("1:1:1","2:2:1")&!dat$Description %in% "E00", ]
vec3=row.names(otu_table)
vec4=unique(as.vector(dat_21$Description))

for (i in 1:2){
  for (j in 1:length(vec4)){
    
    # data("recovery", package = "multcomp")
    AA.t <- t.test(value ~ Other, data = dat_21[dat_21$variable %in% vec3[i]& dat_21$Description%in% vec4[j],])
    # RA.mc <- glht(RA.t,
    #               linfct = mcp(Description = "t.test"),
    #               alternative = "less")
    
    
    print(paste0("t-test for ",vec3[i]," ",vec4[j]," AA ratio for 2a:2b:1c over 1a:1b:1c is constructing"))
    print(AA.t)
    
  }
}


