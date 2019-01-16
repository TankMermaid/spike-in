rm(list=ls())
options(warn=3)


#########The taxa are ordered by the significance of the correlation between their QMP abundance 

########### setting the working directory and print it ###################
tem <- "rootDisease"
setwd("~/xiaoxuan/180724/TARDF16/script/")
print(paste("Your working directory is in",getwd()))



########### to import the plotting theme() function ########################
source("plot_function.R")
# install.packages("ggpubr")
# install.packages("corrplot")
# install.packages("testthat")
library(testthat)
library(ggcorrplot)
library(corrplot)
library(pheatmap)
library(ggpubr)
library(vegan)
library(ggplot2)
library(reshape)
# library(multcomp)
library(ggsignif)
library("Biobase", quietly=T, warn.conflicts=F)
# library("edgeR", quietly=T, warn.conflicts=F)

library("ggplot2", quietly=T, warn.conflicts=F)
library("gplots", quietly=T, warn.conflicts=F)
library("grid", quietly=T, warn.conflicts=F)
library("RColorBrewer", quietly=T, warn.conflicts=F)
library("reshape2", quietly=T, warn.conflicts=F)
library(dplyr)
library(tidyr)
library(psych)
# library(tidyverse)
# library("VennDiagram", quietly=T, warn.conflicts=F)


######################## AA network construction


design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil"&design$Genotype %in% "Disease",]
id_rm= c("TAD1","TAD5","TAD8")
sub_design= sub_design[!rownames(sub_design) %in% id_rm,]

bac_aa <- read.delim("~/xiaoxuan/180724/TARDB28/spikein/otutab_norm_AA.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_aa <- read.delim("~/xiaoxuan/180724/TARDF16/spikein/otutab_norm_AA.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)

id <- colnames(fungi_aa)%in%colnames(bac_aa)
fungi_aa <- fungi_aa[,id]

id <- match(colnames(bac_aa),colnames(fungi_aa))
bac_aa <- bac_aa[,id]
bac_aa <- bac_aa[,colnames(bac_aa)%in% rownames(sub_design)]
fungi_aa <- fungi_aa[,colnames(fungi_aa)%in% rownames(sub_design)]


bac_aa$species <- rep("bac",nrow(bac_aa))
fungi_aa$species <- rep("fungi",nrow(fungi_aa))
merge_bac_fungi <-rbind(bac_aa,fungi_aa) 

###1. Filtering OTUs by occurrence frequency (i.e.,number of samples an OTU is Present 60% of the samples)
Abu <- merge_bac_fungi
Abu$species <- NULL
table <- Abu
table[table>0.0001] <- 1
table.generalist <- Abu[which(rowSums(table)>=5),]  ### 5 out of 7
Abu <- table.generalist

merge_bac_fungi <- merge_bac_fungi[rownames(merge_bac_fungi) %in% rownames(Abu),]


bac_tax <- read.delim("~/xiaoxuan/180724/TARDF16/result/rep_seq_tax_bac.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_tax <- read.delim("~/xiaoxuan/180724/TARDF16/result/rep_seqs_tax_filtered.txt", row.names= 1,header=F, sep="\t",stringsAsFactors = F)
fungi_tax$V10 <- NULL
fungi_tax$V9 <- NULL

colnames(fungi_tax) <- colnames(bac_tax)

tax_merge_bac_fungi <-rbind(fungi_tax,bac_tax) 

tax_merge_bac_fungi$ID <- rownames(tax_merge_bac_fungi)
merge_bac_fungi$ID <- rownames(merge_bac_fungi)
merge_abun_tax<- merge(merge_bac_fungi,tax_merge_bac_fungi,by = "ID")
write.table(merge_bac_fungi,file = "../table/aa_adjusted_bacAndFungiMerge_abundance_disease.txt",sep = '\t',row.names = T,quote = F)


######dplyr to data manipulation
# install.packages("psych")
library(psych)

aa_dis_a<- merge_abun_tax[,c(2:8,10,15)] %>%
  mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:8)])) %>%
  group_by(Kingdom,Genus) %>%
  summarize_all(sum) %>%
  arrange(desc(rowSumAbun))


aa_dis_a<- as.data.frame(aa_dis_a)
rownames(aa_dis_a)<- aa_dis_a$Genus
aa_dis_a$Kingdom <- NULL
aa_dis_a$Genus <- NULL
aa_dis_a <- as.data.frame(t(aa_dis_a))
aa_dis_a$g__unidentified <- NULL

#### note the most 210 threshold
aa_dis_occor_kendall = corr.test(aa_dis_a[1:222],use="pairwise",method="kendall",adjust="fdr",alpha=.05)
aa_dis_occor_kendall$r
# aa_occor_kendall$t
# aa_occor_kendall$n

aa_dis_occor_spearman <- corr.test(aa_dis_a[1:222],use="pairwise",method="spearman",adjust="fdr",alpha=.05)
aa_dis_occor_spearman$r

# ## to filter the top 30 genus according to AA profiling
# id <- rownames(aa_dis_occor_spearman$r)%in% rownames(ra_dis_sper_mat_top30)
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[id,id]
for (i in c(1:27,29:87,89:101,106:113,115:151,153:168,170:195,198:206,209:214,218:221) ){
  # idd <- grep(rownames(aa_dis_occor_spearman$r)[i],rownames(aa_dis_occor_spearman$r)) 
  # i=1
  print(i)
  # if (i==28) skip("bug")
  a_dis=aa_dis_occor_spearman$r
  a_dis_p=aa_dis_occor_spearman$p
  ids = a_dis_p[,colnames(a_dis_p) %in% colnames(a_dis_p)[i]] <0.01
  # ids
  

   a_dis[ids,ids]
  # a_dis_p[ids,ids]
  # ids = a_dis[,colnames(a_dis) %in% colnames(a_dis)[i]] <0
  # a_dis<- a_dis[ids,ids]
  # if (is.null(ids)) skip("no hit line ")
  write.table(a_dis[ids,ids],file = paste0("../table/global_association_genus/global_AAinteractionWIthBipolaris_dis",colnames(a_dis_p)[i],".txt"),sep = '\t',row.names = T,quote = F)
  # write.table(a_dis_p[ids,ids],file = "../table/global_AAinteractionWIthBipolaris_adjustedP_dis.txt",sep = '\t',row.names = T,quote = F)
  
}

# idd <- grep("g__Bipolaris",rownames(aa_dis_occor_spearman$r)) 
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:50,idd),c(1:50,idd)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:50,idd),c(1:50,idd)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[id,id]
# ids <- match(order(row.names(aa_dis_sper_mat_top30)),row.names(aa_dis_sper_mat_top30))
# 
# aa_dis_sper_mat_top30<- aa_dis_sper_mat_top30[ids,ids]
# aa_dis_sper_p_top30 <- aa_dis_sper_p_top30[ids,ids]


a_dis=aa_dis_occor_spearman$r
a_dis_p=aa_dis_occor_spearman$p
ids = a_dis_p[,colnames(a_dis_p) %in% "g__Bipolaris"] <0.01

a_dis[ids,ids]
a_dis_p[ids,ids]

write.table(a_dis[ids,ids],file = "../table/AAinteractionWIthBipolaris_dis.txt",sep = '\t',row.names = T,quote = F)
write.table(a_dis_p[ids,ids],file = "../table/AAinteractionWIthBipolaris_adjustedP_dis.txt",sep = '\t',row.names = T,quote = F)


#### global negative correlation
corrplot(aa_dis_occor_spearman$r,p.mat = aa_dis_occor_spearman$p, sig.level = .01,
         insig = "blank",  addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.65,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")






# aa_heal_occor_kendall$p[ids,ids]
# 
# pheatmap(aa_dis_sper_mat_top30)
# 
# name_vector <- c("Acholeplasma"   ,     "Promicromonospora" ,  "Phyllobacterium" ,   
#                  "Rhizobium"    ,       "Pseudomonas"           ,   
#                   "Sphingomonas"     ,   "Lentzea"       ,      "Flavobacterium" ,    
#                   "Variovorax"     ,     "Pseudoxanthomonas",     
#                  "Massilia"     ,       "Streptomyces"    ,    "Nocardioides"    ,   
#                   "Herbaspirillum"     , "Chryseobacterium"  ,  "Devosia"      ,       "Actinoplanes"     ,  
#                   "Kribbella"         , "Mesorhizobium"     ,  "Pantoea"           ,  "Mycobacterium" ,      "Microbacterium"    , 
#                   "Stenotrophomonas"  ,  "Steroidobacter"  ,    "Sphingopyxis"  ,
#                  "g__Lachnella" ,   "g__Scytalidium",   "g__Acremonium",
#                 "g__Ophiosphaerella",  "g__Magnaporthiopsis" ,  "g__Bipolaris"
# )
# ids <- match(name_vector,row.names(aa_dis_sper_mat_top30))
# 
# aa_dis_sper_mat_top30<- aa_dis_sper_mat_top30[ids,ids]
# aa_dis_sper_p_top30 <- aa_dis_sper_p_top30[ids,ids]


# corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .001
#          insig = "blank", addrect = 2,type = "upper",
#          tl.col = "black",tl.cex = 0.65,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")

pdf(file ="../plot/rootDisease/aa_dis_sper_mat_top50.pdf" )
corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .01,
         insig = "blank",  addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.65,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

dev.off()


corrplot(aa_dis_occor_spearman$r,p.mat = aa_dis_occor_spearman$p, sig.level = .001,
         insig = "blank",  addrect = 2,order="hclust",type = "upper",
         tl.col = "black",tl.cex = 0.15,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

# corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .001,
#          insig = "blank",  addrect = 2,type = "lower",
#          tl.col = "black",tl.cex = 0.65,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")


# bacVsFungiRatio <- sweep(tax_count_sum,,fungi,"/")


########## aa healthy
design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil"&design$Genotype %in% "Healthy",]
id_rm= c("TAH9")
sub_design= sub_design[!rownames(sub_design) %in% id_rm,]


bac_aa <- read.delim("~/xiaoxuan/180724/TARDB28/spikein/otutab_norm_AA.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_aa <- read.delim("~/xiaoxuan/180724/TARDF16/spikein/otutab_norm_AA.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)

id <- colnames(fungi_aa)%in%colnames(bac_aa)
fungi_aa <- fungi_aa[,id]
id <- match(colnames(bac_aa),colnames(fungi_aa))
bac_aa <- bac_aa[,id]
bac_aa <- bac_aa[,colnames(bac_aa)%in% rownames(sub_design)]
fungi_aa <- fungi_aa[,colnames(fungi_aa)%in% rownames(sub_design)]


bac_aa$species <- rep("bac",nrow(bac_aa))
fungi_aa$species <- rep("fungi",nrow(fungi_aa))
merge_bac_fungi <-rbind(bac_aa,fungi_aa) 


###1. Filtering OTUs by occurrence frequency (i.e.,number of samples an OTU is Present 60% of the samples)
Abu <- merge_bac_fungi
Abu$species <- NULL
table <- Abu
table[table>0.0001] <- 1
table.generalist <- Abu[which(rowSums(table)>=7),]  ### 7 out of 9
Abu <- table.generalist

merge_bac_fungi <- merge_bac_fungi[rownames(merge_bac_fungi) %in% rownames(Abu),]


bac_tax <- read.delim("~/xiaoxuan/180724/TARDF16/result/rep_seq_tax_bac.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_tax <- read.delim("~/xiaoxuan/180724/TARDF16/result/rep_seqs_tax_filtered.txt", row.names= 1,header=F, sep="\t",stringsAsFactors = F)
fungi_tax$V10 <- NULL
fungi_tax$V9 <- NULL

colnames(fungi_tax) <- colnames(bac_tax)

tax_merge_bac_fungi <-rbind(fungi_tax,bac_tax) 

tax_merge_bac_fungi$ID <- rownames(tax_merge_bac_fungi)
merge_bac_fungi$ID <- rownames(merge_bac_fungi)
merge_abun_tax<- merge(merge_bac_fungi,tax_merge_bac_fungi,by = "ID")
write.table(merge_bac_fungi,file = "../table/aa_adjusted_bacAndFungiMerge_abundance_healthy.txt",sep = '\t',row.names = T,quote = F)


######dplyr to data manipulation
# install.packages("psych")
library(psych)

aa_heal_a<- merge_abun_tax[,c(2:10,12,17)] %>%
  mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:10)])) %>%
  group_by(Kingdom,Genus) %>%
  summarize_all(sum) %>%
  arrange(desc(rowSumAbun))

aa_heal_a<- as.data.frame(aa_heal_a)
rownames(aa_heal_a)<- aa_heal_a$Genus
aa_heal_a$Kingdom <- NULL
aa_heal_a$Genus <- NULL
aa_heal_a <- as.data.frame(t(aa_heal_a))
aa_heal_a$g__unidentified <- NULL


#### note the most 300 threshold
aa_heal_occor_kendall = corr.test(aa_heal_a[1:220],use="pairwise",method="kendall",adjust="fdr",alpha=.05)
aa_heal_occor_kendall$r


aa_heal_occor_spearman <- corr.test(aa_heal_a[1:220],use="pairwise",method="spearman",adjust="fdr",alpha=.05)
aa_heal_occor_spearman$r

# id <- rownames(aa_heal_occor_spearman$r)%in% rownames(ra_heal_sper_mat_top30)
# aa_heal_sper_mat_top30 <- aa_heal_occor_spearman$r[id,id]
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:30),c(1:30)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:30),c(1:30)]
# aa_heal_sper_p_top30 <- aa_heal_occor_spearman$p[id,id]
idd <- grep("g__Bipolaris",rownames(aa_heal_occor_spearman$r)) 
aa_heal_sper_mat_top30 <- aa_heal_occor_spearman$r[c(1:50,idd),c(1:50,idd)]
aa_heal_sper_p_top30 <- aa_heal_occor_spearman$p[c(1:50,idd),c(1:50,idd)]

# 
# ids <- match(row.names(ra_heal_sper_mat_top30),row.names(aa_heal_sper_mat_top30))
# 
# aa_heal_sper_mat_top30<- aa_heal_sper_mat_top30[ids,ids]
# aa_heal_sper_p_top30 <- aa_heal_sper_p_top30[ids,ids]

# ggsave(pheatmap(aa_heal_sper_mat_top30),filename = "aa_heal_occor_spearman_top30.pdf")

a_heal=aa_heal_occor_spearman$r
a_heal_p=aa_heal_occor_spearman$p
ids = a_heal_p[,colnames(a_heal_p) %in% "g__Bipolaris"] <0.01

a_heal[ids,ids]
a_heal_p[ids,ids]

write.table(a_heal[ids,ids],file = "../table/AA_interactionWIthBipolaris_heal.txt",sep = '\t',row.names = T,quote = F)
write.table(a_heal_p[ids,ids],file = "../table/AA_interactionWIthBipolaris_adjustedP_heal.txt",sep = '\t',row.names = T,quote = F)


# name_vector <- c("Variovorax"   ,     "Flavobacterium" ,   "Actinoplanes"  ,    "Lentzea",          
#                   "Streptomyces"   ,   "Rhizobium"      ,   "Pseudoxanthomonas" ,"Promicromonospora",
#                   "Nocardioides"     ,   "Pseudomonas"  ,     "Massilia"       ,  
#                   "Sphingomonas"   ,   "Umezawaea"      ,   "Devosia"    ,       "Galbitalea" ,      
#                   "Mesorhizobium"    , "Herbaspirillum"  ,  "Caulobacter"  ,     "Phyllobacterium"  ,
#                   "Kribbella"      ,   "Arthrobacter",      "Shinella"  ,        "Aeromicrobium"  , "Kineosporia"     ,  "Chondromyces"   ,   "Rudaibacter"   ,   
#                   "Pantoea"          ,  "g__Podospora" ,"g__Cladosporium",   "g__Mycosphaerella" ,"g__Bipolaris"
# )
# ids <- match(name_vector,row.names(aa_heal_sper_mat_top30))
# 
# aa_heal_sper_mat_top30<- aa_heal_sper_mat_top30[ids,ids]
# aa_heal_sper_p_top30 <- aa_heal_sper_p_top30[ids,ids]

pdf(file ="../plot/rootDisease/aa_heal_sper_mat_top50.pdf" )

corrplot(aa_heal_sper_mat_top30,p.mat = aa_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")



dev.off()


# corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .01,
#          insig = "blank",  addrect = 2,type = "upper",
#          tl.col = "black",tl.cex = 0.65,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")


### RA network construction
#####disease 

design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 

sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil"& design$Genotype %in% "Disease",]
id_rm= c("TAD1","TAD5","TAD8")
sub_design= sub_design[!rownames(sub_design) %in% id_rm,]

bac_aa <- read.delim("~/xiaoxuan/180724/TARDB28/RA/otutab_norm_RA.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_aa <- read.delim("~/xiaoxuan/180724/TARDF16/RA/otutab_norm_RA.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)

id <- colnames(fungi_aa) %in% colnames(bac_aa)
fungi_aa <- fungi_aa[,id]

id <- match(colnames(bac_aa),colnames(fungi_aa))
bac_aa <- bac_aa[,id]
bac_aa <- bac_aa[,colnames(bac_aa)%in% rownames(sub_design)]
fungi_aa <- fungi_aa[,colnames(fungi_aa)%in% rownames(sub_design)]





bac_aa$species <- rep("bac",nrow(bac_aa))
fungi_aa$species <- rep("fungi",nrow(fungi_aa))
merge_bac_fungi <-rbind(bac_aa,fungi_aa) 


# 
# ###1. Filtering OTUs by occurrence frequency (i.e.,number of samples an OTU is Present 60% of the samples)
# Abu <- merge_bac_fungi
# Abu$species <- NULL
# table <- Abu
# table[table>0] <- 1
# table.generalist <- Abu[which(rowSums(table)>=6),]  ### 6 out of 10
# Abu <- table.generalist
# 
# merge_bac_fungi <- merge_bac_fungi[rownames(merge_bac_fungi) %in% rownames(Abu),]

bac_tax <- read.delim("~/xiaoxuan/180724/TARDF16/result/rep_seq_tax_bac.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_tax <- read.delim("~/xiaoxuan/180724/TARDF16/result/rep_seqs_tax_filtered.txt", row.names= 1,header=F, sep="\t",stringsAsFactors = F)
fungi_tax$V10 <- NULL
fungi_tax$V9 <- NULL



colnames(fungi_tax) <- colnames(bac_tax)

tax_merge_bac_fungi <-rbind(fungi_tax,bac_tax) 

tax_merge_bac_fungi$ID <- rownames(tax_merge_bac_fungi)
merge_bac_fungi$ID <- rownames(merge_bac_fungi)
merge_abun_tax<- merge(merge_bac_fungi,tax_merge_bac_fungi,by = "ID")
write.table(merge_bac_fungi,file = "../table/RA_adjusted_bacAndFungiMerge_abundance_disease.txt",sep = '\t',row.names = T,quote = F)


######dplyr to data manipulation
# install.packages("psych")


# ra_dis_a<- merge_abun_tax[,c(2:11,13,18)] %>%
#   mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:11)])) %>%
#   group_by(Kingdom,Genus) %>%
#   summarize_all(sum) %>%
#   arrange(desc(rowSumAbun))
ra_dis_a<- merge_abun_tax[,c(2:8,10,15)] %>%
  mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:8)])) %>%
  group_by(Kingdom,Genus) %>%
  summarize_all(sum) %>%
  arrange(desc(rowSumAbun))
# 
# ra_dis_a$g__unidentified <- NULL

ra_dis_a<- as.data.frame(ra_dis_a)
rownames(ra_dis_a)<- ra_dis_a$Genus
ra_dis_a$Kingdom <- NULL
ra_dis_a$Genus <- NULL
ra_dis_a <- as.data.frame(t(ra_dis_a))
ra_dis_a$g__unidentified <- NULL
ra_dis_a<- ra_dis_a[,colnames(ra_dis_a) %in% colnames(aa_dis_a)]

#### note the most 300 threshold using Kendall and spearman 
ra_dis_occor_kendall = corr.test(ra_dis_a[1:222],use="pairwise",method="kendall",adjust="fdr",alpha=.05)
ra_dis_occor_kendall$r

ra_dis_occor_spearman <- corr.test(ra_dis_a[1:222],use="pairwise",method="spearman",adjust="fdr",alpha=.05)
ra_dis_occor_spearman$r

# a_dis<- as.data.frame(ra_dis_occor_spearman$r)
# a_dis_p<- as.data.frame(ra_dis_occor_spearman$p)
## to filter the top 30 genus according to AA profiling
id <- rownames(ra_dis_occor_spearman$r)%in% rownames(aa_dis_sper_p_top30)
ra_dis_sper_mat_top30 <- ra_dis_occor_spearman$r[id,id]
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:30),c(1:30)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:30),c(1:30)]
ra_dis_sper_p_top30 <- ra_dis_occor_spearman$p[id,id]
ids <- match(row.names(aa_dis_sper_mat_top30),row.names(ra_dis_sper_mat_top30))

ra_dis_sper_mat_top30<- ra_dis_sper_mat_top30[ids,ids]
ra_dis_sper_p_top30 <- ra_dis_sper_p_top30[ids,ids]

a_dis=ra_dis_occor_spearman$r
a_dis_p=ra_dis_occor_spearman$p
ids = a_dis_p[,colnames(a_dis_p) %in% "g__Bipolaris"] <0.01

a_dis[ids,ids]
a_dis_p[ids,ids]

write.table(a_dis[ids,ids],file = "../table/RAinteractionWIthBipolaris_dis.txt",sep = '\t',row.names = T,quote = F)
write.table(a_dis_p[ids,ids],file = "../table/RAinteractionWIthBipolaris_adjustedP_dis.txt",sep = '\t',row.names = T,quote = F)


# ra_dis_sper_mat_top30 <- ra_dis_occor_spearman$r[c(1:40),c(1:40)]
# ra_dis_sper_p_top30 <- ra_dis_occor_spearman$p[c(1:40),c(1:40)]

pheatmap(ra_dis_sper_mat_top30)


pdf(file ="../plot/rootDisease/ra_dis_sper_mat_top50.pdf" )

corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "lower",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.length = 11,
         # cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation",
         tl.pos = "ld")

dev.off()


corrplot(ra_dis_occor_spearman$r,p.mat = ra_dis_occor_spearman$p, sig.level = .01,
         insig = "blank", addrect = 2,type = "lower",
         tl.col = "black",tl.cex = 0.15,
         cl.length = 11,
         # cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation",
         tl.pos = "ld")
# corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .01,
#          insig = "blank", addrect = 2,type = "lower",
#          tl.col = "black",tl.cex = 0.75,
#          cl.length = 11,
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation",
#          tl.pos = "ld")

ids = a_dis_p[,colnames(a_dis_p) %in% "g__Bipolaris"] <0.01

a_dis[ids,ids]
a_dis_p[ids,ids]

##### health
design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil"& design$Genotype %in% "Healthy",]
id_rm= c("TAH9")
sub_design= sub_design[!rownames(sub_design) %in% id_rm,]

bac_aa <- read.delim("~/xiaoxuan/180724/TARDB28/RA/otutab_norm_RA.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_aa <- read.delim("~/xiaoxuan/180724/TARDF16/RA/otutab_norm_RA.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)

id <- colnames(fungi_aa)%in%colnames(bac_aa)
fungi_aa <- fungi_aa[,id]

id <- match(colnames(bac_aa),colnames(fungi_aa))
bac_aa <- bac_aa[,id]
bac_aa <- bac_aa[,colnames(bac_aa)%in% rownames(sub_design)]
fungi_aa <- fungi_aa[,colnames(fungi_aa)%in% rownames(sub_design)]


bac_aa$species <- rep("bac",nrow(bac_aa))
fungi_aa$species <- rep("fungi",nrow(fungi_aa))
merge_bac_fungi <-rbind(bac_aa,fungi_aa) 
bac_tax <- read.delim("~/xiaoxuan/180724/TARDF16/result/rep_seq_tax_bac.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_tax <- read.delim("~/xiaoxuan/180724/TARDF16/result/rep_seqs_tax_filtered.txt", row.names= 1,header=F, sep="\t",stringsAsFactors = F)
fungi_tax$V10 <- NULL
fungi_tax$V9 <- NULL

colnames(fungi_tax) <- colnames(bac_tax)

tax_merge_bac_fungi <-rbind(fungi_tax,bac_tax) 

tax_merge_bac_fungi$ID <- rownames(tax_merge_bac_fungi)
merge_bac_fungi$ID <- rownames(merge_bac_fungi)
merge_abun_tax<- merge(merge_bac_fungi,tax_merge_bac_fungi,by = "ID")
write.table(merge_bac_fungi,file = "../table/RA_adjusted_bacAndFungiMerge_abundance_healthy.txt",sep = '\t',row.names = T,quote = F)


######dplyr to data manipulation
# install.packages("psych")
library(psych)

ra_heal_a<- merge_abun_tax[,c(2:10,12,17)] %>%
  mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:10)])) %>%
  group_by(Kingdom,Genus) %>%
  summarize_all(sum) %>%
  arrange(desc(rowSumAbun))

ra_heal_a<- as.data.frame(ra_heal_a)
rownames(ra_heal_a)<- ra_heal_a$Genus
ra_heal_a$Kingdom <- NULL
ra_heal_a$Genus <- NULL
ra_heal_a <- as.data.frame(t(ra_heal_a))
ra_heal_a$g__unidentified <- NULL
ra_heal_a<- ra_heal_a[,colnames(ra_heal_a) %in% colnames(aa_heal_a)]

#### note the most 300 threshold using Kendall and Spearman 
ra_heal_occor_kendall = corr.test(ra_heal_a[1:220],use="pairwise",method="kendall",adjust="fdr",alpha=.05)
ra_heal_occor_kendall$r

ra_heal_occor_spearman <- corr.test(ra_heal_a[1:220],use="pairwise",method="spearman",adjust="fdr",alpha=.05)
ra_heal_occor_spearman$r
# ra_heal_sper_mat_top30 <- ra_heal_occor_spearman$r[c(1:40),c(1:40)]
# ra_heal_sper_p_top30 <- ra_heal_occor_spearman$p[c(1:40),c(1:40)]

id <- rownames(ra_heal_occor_spearman$r)%in% rownames(aa_heal_occor_spearman$r)
ra_heal_sper_mat_top30 <- ra_heal_occor_spearman$r[id,id]
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:30),c(1:30)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:30),c(1:30)]
ra_heal_sper_p_top30 <- ra_heal_occor_spearman$p[id,id]
# aa_heal_sper_mat_top30 <- aa_heal_occor_spearman$r[c(1:30),c(1:30)]
# aa_heal_sper_p_top30 <- aa_heal_occor_spearman$p[c(1:30),c(1:30)]


ids <- match(row.names(aa_heal_sper_mat_top30),row.names(ra_heal_sper_mat_top30))

ra_heal_sper_mat_top30<- ra_heal_sper_mat_top30[ids,ids]
ra_heal_sper_p_top30 <- ra_heal_sper_p_top30[ids,ids]

ra_heal=ra_heal_occor_spearman$r
ra_heal_p=ra_heal_occor_spearman$p
ids = ra_heal_p[,colnames(ra_heal_p) %in% "g__Bipolaris"] <0.01
# ids = ra_heal_p[,colnames(ra_heal_p) %in% "g__Furasium"] <0.01

ra_heal[ids,ids]
ra_heal_p[ids,ids]


write.table(ra_heal[ids,ids],file = "../table/RA_interactionWIthBipolaris_heal.txt",sep = '\t',row.names = T,quote = F)
write.table(ra_heal_p[ids,ids],file = "../table/RA_interactionWIthBipolaris_adjustedP_heal.txt",sep = '\t',row.names = T,quote = F)



ggsave(pheatmap(ra_heal_sper_mat_top30),filename = "ra_heal_occor_spearman_top30.pdf")
a<- as.data.frame(ra_heal_occor_spearman$r)

pdf(file ="../plot/rootDisease/ra_heal_sper_mat_top50.pdf" )
corrplot(ra_heal_sper_mat_top30,p.mat = ra_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "lower",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")
# corrplot(ra_heal_sper_mat_top30,p.mat = ra_heal_sper_p_top30, sig.level = .01,
#          insig = "blank", addrect = 2,type = "lower",
#          tl.col = "black",tl.cex = 0.75,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")

dev.off()



