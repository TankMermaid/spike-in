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
# install.packages("ggcorrplot")
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

colnames(aa_dis_a)[1:100]

aa_dis_a<- as.data.frame(aa_dis_a)
rownames(aa_dis_a)<- aa_dis_a$Genus
aa_dis_a$Kingdom <- NULL
aa_dis_a$Genus <- NULL
aa_dis_a <- as.data.frame(t(aa_dis_a))
aa_dis_a$g__unidentified <- NULL


#### correlation of top 50 genus with microbial load
dis_micload<- as.data.frame(rowSums(aa_dis_a))
colnames(dis_micload) <- "microbial load"

top50_cor_micload_pearson<- corr.test(aa_dis_a[,c(1:50)],dis_micload,use = "complete",method="pearson",adjust = "fdr",alpha = .05)
top50_cor_micload_pearson$p
top50_cor_micload_pearson$r
write.table(top50_cor_micload_pearson$p,file = "../table/1028_AA_dis_top50_cor_micload_pearson_pvalue.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_pearson$r,file = "../table/1028_AA_dis_top50_cor_micload_pearson_coreffi.txt",sep = '\t',row.names = T,quote = F)

top50_cor_micload_spearman<- corr.test(aa_dis_a[,c(1:50)],dis_micload,use = "complete",method="spearman",adjust = "fdr",alpha = .05)
top50_cor_micload_spearman$p
top50_cor_micload_spearman$r
write.table(top50_cor_micload_spearman$p,file = "../table/1028_AA_dis_top50_cor_micload_spearman_pvalue.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$r,file = "../table/1028_AA_dis_top50_cor_micload_spearman_coreffi.txt",sep = '\t',row.names = T,quote = F)


idd <- grep("g__Bipolaris",colnames(aa_dis_a)) 
idd1 <- grep("Actinosynnema",colnames(aa_dis_a)) 

top50_cor_micload_pearson<- corr.test(aa_dis_a[,c(78)],dis_micload,use = "complete",method="pearson",adjust = "fdr",alpha = .05)
top50_cor_micload_pearson$p
top50_cor_micload_pearson$r
write.table(top50_cor_micload_pearson$p,file = "../table/1028_AA_dis_g__Bipolaris_cor_micload_pearson_pvalue.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_pearson$r,file = "../table/1028_AA_dis_g__Bipolaris_cor_micload_pearson_coreffi.txt",sep = '\t',row.names = T,quote = F)

top50_cor_micload_pearson<- corr.test(aa_dis_a[,c(78)],dis_micload,use = "complete",method="spearman",adjust = "fdr",alpha = .05)
top50_cor_micload_pearson$p
top50_cor_micload_pearson$r
write.table(top50_cor_micload_pearson$p,file = "../table/1028_AA_dis_g__Bipolaris_cor_micload_spearman_pvalue.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_pearson$r,file = "../table/1028_AA_dis_g__Bipolaris_cor_micload_spearman_coreffi.txt",sep = '\t',row.names = T,quote = F)




top50_cor_micload_spearman<- corr.test(aa_dis_a[,c(82)],dis_micload,use = "complete",method="pearson",adjust = "fdr",alpha = .05)
top50_cor_micload_spearman$p
top50_cor_micload_spearman$r
write.table(top50_cor_micload_spearman$p,file = "../table/1028_AA_dis_Actinosynnema_cor_micload_pearson_pvalue.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$r,file = "../table/1028_AA_dis_Actinosynnema_cor_micload_pearson_coreffi.txt",sep = '\t',row.names = T,quote = F)



top50_cor_micload_spearman<- corr.test(aa_dis_a[,c(82)],dis_micload,use = "complete",method="spearman",adjust = "fdr",alpha = .05)
top50_cor_micload_spearman$p
top50_cor_micload_spearman$r
write.table(top50_cor_micload_spearman$p,file = "../table/1028_AA_dis_Actinosynnema_cor_micload_spearman_pvalue.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$r,file = "../table/1028_AA_dis_Actinosynnema_cor_micload_spearman_coreffi.txt",sep = '\t',row.names = T,quote = F)


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
rownames(aa_dis_occor_spearman)
idd <- grep("g__Bipolaris",rownames(aa_dis_occor_spearman$r)) 
idd1 <- grep("g__Mycosphaerella",rownames(aa_dis_occor_spearman$r)) 
idd2 <- grep("g__Didymella",rownames(aa_dis_occor_spearman$r)) 
idd3 <- grep("Actinosynnema",rownames(aa_dis_occor_spearman$r)) 
idd4 <- grep("Reyranella",rownames(aa_dis_occor_spearman$r)) 
idd5 <- grep("Haliangium",rownames(aa_dis_occor_spearman$r)) 
# idd5 <- grep("Peredibacter",rownames(aa_dis_occor_spearman$r)) 
# idd6 <- grep("Polyangium",rownames(aa_dis_occor_spearman$r)) 

aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:50,idd,idd1,idd2,idd3,idd4,idd5),c(1:50,idd,idd1,idd2,idd3,idd4,idd5)]
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:50,idd,idd1,idd2,idd3,idd4,idd5,idd6),c(1:50,idd,idd1,idd2,idd3,idd4,idd5,idd6)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:50,idd,idd1,idd2,idd3,idd4,idd5,idd6),c(1:50,idd,idd1,idd2,idd3,idd4,idd5,idd6)]
aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:50,idd,idd1,idd2,idd3,idd4,idd5),c(1:50,idd,idd1,idd2,idd3,idd4,idd5)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[id,id]
# ids <- match(order(row.names(aa_dis_sper_mat_top30)),row.names(aa_dis_sper_mat_top30))

top50_plus_other <- aa_dis_a[,c(1:50,idd,idd1,idd2,idd3,idd4,idd5)]
# top50_plus_other <- aa_dis_a[,c(1:50,idd,idd1,idd2,idd3,idd4,idd5,idd6)]
# 
write.table(top50_plus_other,file = "../table/1030_AA_dis_top50_abundance.txt",sep = '\t',row.names = T,quote = F)
# write.table(aa_dis_sper_mat_top30,file = "../table/1024_AA_dis_top50_cor.txt",sep = '\t',row.names = T,quote = F)




write.table(aa_dis_sper_mat_top30,file = "../table/1024_AA_dis_top50_cor.txt",sep = '\t',row.names = T,quote = F)
write.table(aa_dis_sper_p_top30,file = "../table/1024_AA_adjustedP_dis_top50.txt",sep = '\t',row.names = T,quote = F)

# aa_dis_sper_mat_top30<- aa_dis_sper_mat_top30[ids,ids]
# aa_dis_sper_p_top30 <- aa_dis_sper_p_top30[ids,ids]


a_dis=aa_dis_occor_spearman$r
a_dis_p=aa_dis_occor_spearman$p
ids_p = as.data.frame(a_dis_p[a_dis_p[,colnames(a_dis_p) %in% "g__Bipolaris"] <0.01,colnames(a_dis_p) %in% "g__Bipolaris"])
colnames(ids_p)= "adjusted-P"


ids_m=a_dis[rownames(a_dis)%in% rownames(ids_p),colnames(a_dis)%in% rownames(ids_p)]



write.table(ids_m,file = "../table/1023_AAinteractionWIthBipolaris_dis.txt",sep = '\t',row.names = T,quote = F)
write.table(ids_p,file = "../table/1023_AAinteractionWIthBipolaris_adjustedP_dis.txt",sep = '\t',row.names = T,quote = F)


#### global negative correlation
# corrplot(aa_dis_occor_spearman$r,p.mat = aa_dis_occor_spearman$p, sig.level = .01,
#          insig = "blank",  addrect = 2,type = "upper",diag = F,
#          tl.col = "black",tl.cex = 0.65,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")






# aa_heal_occor_kendall$p[ids,ids]
# 
# pheatmap(aa_dis_sper_mat_top30)
# colname
name_vector <- c( "Promicromonospora" ,  "Phyllobacterium" ,    "Pseudomonas"  ,      
                  "Rhizobium"        ,   "Flavobacterium"    ,  "Sphingomonas"  ,     
                  "Variovorax"    ,      "Pseudoxanthomonas"  ,  "Streptomyces" ,   
                  "Chryseobacterium"  ,  "Lentzea"         ,    "Nocardioides"  ,     
                  "Actinoplanes"    ,    "Herbaspirillum"    ,  
                  "Mycobacterium"  ,     "Kribbella"        ,   "Massilia"   ,        
                    "Devosia"        ,     "Microbacterium"     ,
                  "Sphingopyxis"    ,    "Stenotrophomonas"  ,  "Mesorhizobium" ,     
                       "Caulobacter"   ,      "Aeromicrobium"   ,   
                  "Sorangium"     ,      "Agromyces"    ,       "Steroidobacter"   ,  
                  "Arthrobacter"  ,      "Chondromyces"   ,    
                  "Umezawaea"   ,        "Pantoea"    ,        "Rhodanobacter"    ,  
                     "Shinella"   ,         "Lysobacter"      ,   
                  "Solirubrobacter"    , "Bosea"          ,     "Glycomyces"    ,     
                  "Cellulomonas"     ,   "Bradyrhizobium"    ,  "Burkholderia"   ,    
                  "Rhodopseudomonas" ,   "Bacillus"    ,     "Actinosynnema",  "Reyranella","Haliangium",
                  "g__Acremonium" ,"g__Lachnella"      ,  "g__Scytalidium",
                  "g__Magnaporthiopsis" ,"g__Ophiosphaerella","g__Gibberella",   
                  "g__Pyrenochaetopsis", "g__Sarocladium", "g__Mycosphaerella","g__Didymella","g__Bipolaris" 
)
ids <- match(name_vector,row.names(aa_dis_sper_mat_top30))
# 
aa_dis_sper_mat_top30<- aa_dis_sper_mat_top30[ids,ids]
aa_dis_sper_p_top30 <- aa_dis_sper_p_top30[ids,ids]


# corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .001
#          insig = "blank", addrect = 2,type = "upper",
#          tl.col = "black",tl.cex = 0.65,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")

pdf(file ="../plot/rootDisease/aa_dis_sper_mat_top50_order_upper_1030.pdf" )
corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

dev.off()


### lower triangle not adjusted p-value
# pdf(file ="../plot/rootDisease/aa_dis_sper_mat_top50_order_lower.pdf" )
# corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .01,
#          insig = "blank",  addrect = 2, col = cm.colors(100),type = "lower",diag = F,
#          tl.col = "black",tl.cex = 0.65,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")
# 
# dev.off()

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
#          cl.lim = c(-1,1),plot.new()
#          title = "Genus-genus spearman correlation")

png(file ="../plot/rootDisease/aa_dis_sper_mat_top50_order_upper_1023.png" ,width = 11,height = 11)

corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .01,
         insig = "blank",  addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.65,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

dev.off()

png(file ="../plot/rootDisease/aa_dis_sper_mat_top50_order_lower.png" )
corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .01,
         insig = "blank",  addrect = 2,type = "lower",diag = F,
         tl.col = "black",tl.cex = 0.65,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

dev.off()
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

#### correlation of top 50 genus with microbial load
heal_micload<- as.data.frame(rowSums(aa_heal_a))
colnames(heal_micload) <- "microbial load"

top50_cor_micload_pearson<- corr.test(aa_heal_a[,c(1:50)],heal_micload,use = "complete",method="pearson",adjust = "fdr",alpha = .05)
top50_cor_micload_pearson$p
top50_cor_micload_pearson$r
write.table(top50_cor_micload_pearson$p,file = "../table/1028_AA_heal_top50_cor_micload_pearson_pvalue.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_pearson$r,file = "../table/1028_AA_heal_top50_cor_micload_pearson_coreffi.txt",sep = '\t',row.names = T,quote = F)

idd <- grep("g__Bipolaris",colnames(aa_heal_a)) 
idd1 <- grep("Actinosynnema",colnames(aa_heal_a)) 

top50_cor_micload_pearson<- corr.test(aa_heal_a[,c(133)],heal_micload,use = "complete",method="pearson",adjust = "fdr",alpha = .05)
top50_cor_micload_pearson$p
top50_cor_micload_pearson$r
write.table(top50_cor_micload_pearson$p,file = "../table/1028_AA_heal_g__Bipolaris_cor_micload_pearson_pvalue.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_pearson$r,file = "../table/1028_AA_heal_g__Bipolaris_cor_micload_pearson_coreffi.txt",sep = '\t',row.names = T,quote = F)

top50_cor_micload_pearson<- corr.test(aa_heal_a[,c(133)],heal_micload,use = "complete",method="spearman",adjust = "fdr",alpha = .05)
top50_cor_micload_pearson$p
top50_cor_micload_pearson$r
write.table(top50_cor_micload_pearson$p,file = "../table/1028_AA_heal_g__Bipolaris_cor_micload_spearman_pvalue.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_pearson$r,file = "../table/1028_AA_heal_g__Bipolaris_cor_micload_spearman_coreffi.txt",sep = '\t',row.names = T,quote = F)




top50_cor_micload_spearman<- corr.test(aa_heal_a[,c(88)],heal_micload,use = "complete",method="pearson",adjust = "fdr",alpha = .05)
top50_cor_micload_spearman$p
top50_cor_micload_spearman$r
write.table(top50_cor_micload_spearman$p,file = "../table/1028_AA_heal_Actinosynnema_cor_micload_pearson_pvalue.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$r,file = "../table/1028_AA_heal_Actinosynnema_cor_micload_pearson_coreffi.txt",sep = '\t',row.names = T,quote = F)


top50_cor_micload_spearman<- corr.test(aa_heal_a[,c(88)],heal_micload,use = "complete",method="spearman",adjust = "fdr",alpha = .05)
top50_cor_micload_spearman$p
top50_cor_micload_spearman$r
write.table(top50_cor_micload_spearman$p,file = "../table/1028_AA_heal_Actinosynnema_cor_micload_spearman_pvalue.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$r,file = "../table/1028_AA_heal_Actinosynnema_cor_micload_spearman_coreffi.txt",sep = '\t',row.names = T,quote = F)

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
idd1 <- grep("Haliangium",rownames(aa_heal_occor_spearman$r))
idd2 <- grep("g__Didymella",rownames(aa_heal_occor_spearman$r)) 
idd3 <- grep("Actinosynnema",rownames(aa_heal_occor_spearman$r)) 
idd4 <- grep("Reyranella",rownames(aa_heal_occor_spearman$r)) 

# idd5 <- grep("Peredibacter",rownames(aa_heal_occor_spearman$r))
# idd5 <- grep("Peredibacter",rownames(aa_heal_occor_spearman$r)) 
# idd6 <- grep("Polyangium",rownames(aa_heal_occor_spearman$r)) 

aa_heal_sper_mat_top30 <- aa_heal_occor_spearman$r[c(1:50,idd,idd1,idd2,idd3,idd4),c(1:50,idd,idd1,idd2,idd3,idd4)]
# aa_heal_sper_mat_top30 <- aa_heal_occor_spearman$r[c(1:50,idd,idd2,idd3,idd4,idd5,idd6),c(1:50,idd,idd2,idd3,idd4,idd5,idd6)]
aa_heal_sper_p_top30 <- aa_heal_occor_spearman$p[c(1:50,idd,idd1,idd2,idd3,idd4),c(1:50,idd,idd1,idd2,idd3,idd4)]
# aa_heal_sper_p_top30 <- aa_heal_occor_spearman$p[c(1:50,idd,idd1,idd2,idd3,idd4,idd5,idd6),c(1:50,idd,idd1,idd2,idd3,idd4,idd5,idd6)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[id,id]
# ids <- match(order(row.names(aa_dis_sper_mat_top30)),row.names(aa_dis_sper_mat_top30))

top50_plus_other <- aa_heal_a[,c(1:50,idd,idd1,idd2,idd3,idd4)]
# top50_plus_other <- aa_heal_a[,c(1:50,idd,idd2,idd3,idd4,idd5,idd6)]
# 
write.table(top50_plus_other,file = "../table/1030_AA_heal_top50_abundance.txt",sep = '\t',row.names = T,quote = F)

# aa_heal_sper_mat_top30 <- aa_heal_occor_spearman$r[c(1:50,idd),c(1:50,idd)]
# aa_heal_sper_p_top30 <- aa_heal_occor_spearman$p[c(1:50,idd),c(1:50,idd)]
write.table(aa_heal_sper_mat_top30,file = "../table/1105_AA_heal_top50_cor.txt",sep = '\t',row.names = T,quote = F)
write.table(aa_heal_sper_p_top30,file = "../table/1105_AA_adjustedP_heal_top50.txt",sep = '\t',row.names = T,quote = F)

# 
# ids <- match(row.names(ra_heal_sper_mat_top30),row.names(aa_heal_sper_mat_top30))
# 
# aa_heal_sper_mat_top30<- aa_heal_sper_mat_top30[ids,ids]
# aa_heal_sper_p_top30 <- aa_heal_sper_p_top30[ids,ids]

# ggsave(pheatmap(aa_heal_sper_mat_top30),filename = "aa_heal_occor_spearman_top30.pdf")

a_heal=aa_heal_occor_spearman$r
a_heal_p=aa_heal_occor_spearman$p
# ids = a_heal_p[,colnames(a_heal_p) %in% "g__Bipolaris"] <0.01

ids_p = as.data.frame(a_heal_p[a_heal_p[,colnames(a_heal_p) %in% "g__Bipolaris"] <0.01,colnames(a_heal_p) %in% "g__Bipolaris"])
colnames(ids_p)= "adjusted-P"


ids_m=a_heal[rownames(a_heal)%in% rownames(ids_p),colnames(a_heal)%in% rownames(ids_p)]



write.table(ids_m,file = "../table/1023_AAinteractionWIthBipolaris_heal.txt",sep = '\t',row.names = T,quote = F)
write.table(ids_p,file = "../table/1023_AAinteractionWIthBipolaris_adjustedP_heal.txt",sep = '\t',row.names = T,quote = F)





# name_vector <- c("Variovorax"    ,     "Flavobacterium"   ,  "Actinoplanes"    ,   "Lentzea"   ,         "Streptomyces"  ,    
#                  "Pseudoxanthomonas" , "Rhizobium"     ,     "Promicromonospora" , "Nocardioides"  ,         
#                   "Pseudomonas"    ,    "Massilia"      ,     "Umezawaea"    ,      "Devosia"   ,         "Galbitalea" ,       
#                  "Sphingomonas"    ,   "Mesorhizobium"   ,   "Herbaspirillum"   ,  "Caulobacter"  ,      "Phyllobacterium" ,  
#                  "Kribbella"   ,       "Arthrobacter"      , "Shinella"     ,      "Aeromicrobium"  ,    "Kineosporia"     ,  
#                     "Chondromyces" ,      "Rudaibacter"   ,       "Steroidobacter"  ,  
#                   "Stenotrophomonas"  , "Pantoea"    ,            "Ideonella"  ,       
#                   "Rhodanobacter"    ,  "Bradyrhizobium"    , "Phycicoccus"     ,   "Chryseobacterium"  , "Agromyces"  ,       
#                  "Microbacterium"   ,  "Solirubrobacter" ,   "Mycobacterium"    ,  "Bosea"       ,       "Altererythrobacter",
#                   "Aquabacterium"   ,     "Curtobacterium"   ,     "Dactylosporangium" ,"Actinosynnema",  "Reyranella","Peredibacter","Polyangium",
#                  "g__Podospora" , "g__Cladosporium" ,"g__Mycosphaerella", "g__Nigrospora"    , "g__Mortierella" ,"g__Edenia"  ,
#                  "g__Metacordyceps","g__Mycosphaerella","g__Didymella","g__Bipolaris" 
# )

name_vector <- c("Variovorax"    ,     "Flavobacterium"   ,  "Actinoplanes"    ,   "Lentzea"   ,         "Streptomyces"  ,    
                 "Pseudoxanthomonas" , "Rhizobium"     ,     "Promicromonospora" , "Nocardioides"  ,         
                  "Pseudomonas"    ,    "Massilia"      ,     "Umezawaea"    ,      "Devosia"   ,         "Galbitalea" ,       
                 "Sphingomonas"    ,   "Mesorhizobium"   ,   "Herbaspirillum"   ,  "Caulobacter"  ,      "Phyllobacterium" ,  
                 "Kribbella"   ,       "Arthrobacter"      , "Shinella"     ,      "Aeromicrobium"  ,    "Kineosporia"     ,  
                    "Chondromyces" ,      "Rudaibacter"   ,       "Steroidobacter"  ,  
                  "Stenotrophomonas"  , "Pantoea"    ,            "Ideonella"  ,       
                  "Rhodanobacter"    ,  "Bradyrhizobium"    , "Phycicoccus"     ,   "Chryseobacterium"  , "Agromyces"  ,       
                 "Microbacterium"   ,  "Solirubrobacter" ,   "Mycobacterium"    ,  "Bosea"       ,       "Altererythrobacter",
                  "Aquabacterium"   ,     "Curtobacterium"   ,     "Dactylosporangium" ,"Actinosynnema",  "Reyranella","Haliangium",
                 "g__Podospora" , "g__Cladosporium" , "g__Nigrospora"    , "g__Mortierella" ,"g__Edenia"  ,
                 "g__Metacordyceps","g__Mycosphaerella","g__Didymella","g__Bipolaris" 
)
ids <- match(name_vector,row.names(aa_heal_sper_mat_top30))

aa_heal_sper_mat_top30<- aa_heal_sper_mat_top30[ids,ids]
aa_heal_sper_p_top30 <- aa_heal_sper_p_top30[ids,ids]

pdf(file ="../plot/rootDisease/aa_heal_sper_mat_top50_sort_upper_1105.pdf" )

corrplot(aa_heal_sper_mat_top30,p.mat = aa_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")



dev.off()

pdf(file ="../plot/rootDisease/aa_heal_sper_mat_top50_sort_lower.pdf" )

corrplot(aa_heal_sper_mat_top30,p.mat = aa_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2, col = cm.colors(100),type = "lower",diag = F,
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

png(file ="../plot/rootDisease/aa_heal_sper_mat_top50_sort_upper_1023.png" )

corrplot(aa_heal_sper_mat_top30,p.mat = aa_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")



dev.off()

png(file ="../plot/rootDisease/aa_heal_sper_mat_top50_sort_lower.png" )

corrplot(aa_heal_sper_mat_top30,p.mat = aa_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "lower",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")



dev.off()


##### AA health without bipolaris

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
aa_heal_sper_mat_top30 <- aa_heal_occor_spearman$r[c(1:50),c(1:50)]
aa_heal_sper_p_top30 <- aa_heal_occor_spearman$p[c(1:50),c(1:50)]
write.table(aa_heal_sper_mat_top30,file = "../table/1025_AA_heal_top50_cor_without_Bipolaris.txt",sep = '\t',row.names = T,quote = F)
write.table(aa_heal_sper_p_top30,file = "../table/1025_AA_adjustedP_heal_top50_without_Bipolaris.txt",sep = '\t',row.names = T,quote = F)

# 
# ids <- match(row.names(ra_heal_sper_mat_top30),row.names(aa_heal_sper_mat_top30))
# 
# aa_heal_sper_mat_top30<- aa_heal_sper_mat_top30[ids,ids]
# aa_heal_sper_p_top30 <- aa_heal_sper_p_top30[ids,ids]

# ggsave(pheatmap(aa_heal_sper_mat_top30),filename = "aa_heal_occor_spearman_top30.pdf")

# a_heal=aa_heal_occor_spearman$r
# a_heal_p=aa_heal_occor_spearman$p
# # ids = a_heal_p[,colnames(a_heal_p) %in% "g__Bipolaris"] <0.01
# 
# ids_p = as.data.frame(a_heal_p[a_heal_p[,colnames(a_heal_p) %in% "g__Bipolaris"] <0.01,colnames(a_heal_p) %in% "g__Bipolaris"])
# colnames(ids_p)= "adjusted-P"
# 
# 
# ids_m=a_heal[rownames(a_heal)%in% rownames(ids_p),colnames(a_heal)%in% rownames(ids_p)]
# 
# 
# 
# write.table(ids_m,file = "../table/1023_AAinteractionWIthBipolaris_heal.txt",sep = '\t',row.names = T,quote = F)
# write.table(ids_p,file = "../table/1023_AAinteractionWIthBipolaris_adjustedP_heal.txt",sep = '\t',row.names = T,quote = F)





name_vector <- c("Variovorax"    ,     "Flavobacterium"   ,  "Actinoplanes"    ,   "Lentzea"   ,         "Streptomyces"  ,    
                 "Pseudoxanthomonas" , "Rhizobium"     ,     "Promicromonospora" , "Nocardioides"  ,         
                 "Pseudomonas"    ,    "Massilia"      ,     "Umezawaea"    ,      "Devosia"   ,         "Galbitalea" ,       
                 "Sphingomonas"    ,   "Mesorhizobium"   ,   "Herbaspirillum"   ,  "Caulobacter"  ,      "Phyllobacterium" ,  
                 "Kribbella"   ,       "Arthrobacter"      , "Shinella"     ,      "Aeromicrobium"  ,    "Kineosporia"     ,  
                 "Chondromyces" ,      "Rudaibacter"   ,       "Steroidobacter"  ,  
                 "Stenotrophomonas"  , "Pantoea"    ,            "Ideonella"  ,       
                 "Rhodanobacter"    ,  "Bradyrhizobium"    , "Phycicoccus"     ,   "Chryseobacterium"  , "Agromyces"  ,       
                 "Microbacterium"   ,  "Solirubrobacter" ,   "Mycobacterium"    ,  "Bosea"       ,       "Altererythrobacter",
                 "Aquabacterium"   ,     "Curtobacterium"   ,     "Dactylosporangium" ,
                 "g__Podospora" , "g__Cladosporium" ,"g__Mycosphaerella", "g__Nigrospora"    , "g__Mortierella" ,"g__Edenia"  ,
                 "g__Metacordyceps"
)
ids <- match(name_vector,row.names(aa_heal_sper_mat_top30))

aa_heal_sper_mat_top30<- aa_heal_sper_mat_top30[ids,ids]
aa_heal_sper_p_top30 <- aa_heal_sper_p_top30[ids,ids]

pdf(file ="../plot/rootDisease/aa_heal_sper_mat_top50_sort_upper_1025.pdf" )

corrplot(aa_heal_sper_mat_top30,p.mat = aa_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")



dev.off()

pdf(file ="../plot/rootDisease/aa_heal_sper_mat_top50_sort_lower.pdf" )

corrplot(aa_heal_sper_mat_top30,p.mat = aa_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2, col = cm.colors(100),type = "lower",diag = F,
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

png(file ="../plot/rootDisease/aa_heal_sper_mat_top50_sort_upper_1023.png" )

corrplot(aa_heal_sper_mat_top30,p.mat = aa_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")



dev.off()

png(file ="../plot/rootDisease/aa_heal_sper_mat_top50_sort_lower.png" )

corrplot(aa_heal_sper_mat_top30,p.mat = aa_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "lower",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")



dev.off()



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
top50_plus_other <- ra_dis_a[,colnames(ra_dis_a)%in% rownames(ra_dis_sper_mat_top30)]
# 
write.table(top50_plus_other,file = "../table/1030_RA_dis_top50_abundance.txt",sep = '\t',row.names = T,quote = F)

write.table(ra_dis_sper_mat_top30,file = "../table/1024_RA_dis_top50_cor.txt",sep = '\t',row.names = T,quote = F)
write.table(ra_dis_sper_p_top30,file = "../table/1024_RA_adjustedP_dis_top50.txt",sep = '\t',row.names = T,quote = F)


a_dis=ra_dis_occor_spearman$r
a_dis_p=ra_dis_occor_spearman$p

ids_p = as.data.frame(a_dis_p[a_dis_p[,colnames(a_dis_p) %in% "g__Bipolaris"] <0.01,colnames(a_dis_p) %in% "g__Bipolaris"])
colnames(ids_p)= "adjusted-P"


ids_m=a_dis[rownames(a_dis)%in% rownames(ids_p),colnames(a_dis)%in% rownames(ids_p)]



write.table(ids_m,file = "../table/1023_RAinteractionWIthBipolaris_dis.txt",sep = '\t',row.names = T,quote = F)
write.table(ids_p,file = "../table/1023_RAinteractionWIthBipolaris_adjustedP_dis.txt",sep = '\t',row.names = T,quote = F)



# write.table(a_dis[ids,ids],file = "../table/RAinteractionWIthBipolaris_dis.txt",sep = '\t',row.names = T,quote = F)
# write.table(a_dis_p[ids,ids],file = "../table/RAinteractionWIthBipolaris_adjustedP_dis.txt",sep = '\t',row.names = T,quote = F)






# ra_dis_sper_mat_top30 <- ra_dis_occor_spearman$r[c(1:40),c(1:40)]
# ra_dis_sper_p_top30 <- ra_dis_occor_spearman$p[c(1:40),c(1:40)]

pheatmap(ra_dis_sper_mat_top30)


pdf(file ="../plot/rootDisease/ra_dis_sper_mat_top50_upper_1030.pdf" )

corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

dev.off()

png(file ="../plot/rootDisease/ra_dis_sper_mat_top50_upper_1023.png" )

corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")


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

id <- rownames(ra_heal_occor_spearman$r)%in% rownames(aa_heal_sper_mat_top30)
ra_heal_sper_mat_top30 <- ra_heal_occor_spearman$r[id,id]
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:30),c(1:30)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:30),c(1:30)]
ra_heal_sper_p_top30 <- ra_heal_occor_spearman$p[id,id]
# aa_heal_sper_mat_top30 <- aa_heal_occor_spearman$r[c(1:30),c(1:30)] 
# aa_heal_sper_p_top30 <- aa_heal_occor_spearman$p[c(1:30),c(1:30)]


ids <- match(row.names(aa_heal_sper_mat_top30),row.names(ra_heal_sper_mat_top30))

ra_heal_sper_mat_top30<- ra_heal_sper_mat_top30[ids,ids]
ra_heal_sper_p_top30 <- ra_heal_sper_p_top30[ids,ids]

top50_plus_other <- ra_heal_a[,colnames(ra_heal_a)%in% rownames(ra_heal_sper_mat_top30)]
# 
write.table(top50_plus_other,file = "../table/1030_RA_heal_top50_abundance.txt",sep = '\t',row.names = T,quote = F)

write.table(ra_heal_sper_mat_top30,file = "../table/1024_RA_heal_top50_cor.txt",sep = '\t',row.names = T,quote = F)
write.table(ra_heal_sper_p_top30,file = "../table/1024_RA_adjustedP_heal_top50.txt",sep = '\t',row.names = T,quote = F)


ra_heal=ra_heal_occor_spearman$r
ra_heal_p=ra_heal_occor_spearman$p


ids_p = as.data.frame(ra_heal_p[ra_heal_p[,colnames(ra_heal_p) %in% "g__Bipolaris"] <0.01,colnames(ra_heal_p) %in% "g__Bipolaris"])
colnames(ids_p)= "adjusted-P"


ids_m=ra_heal[rownames(ra_heal)%in% rownames(ids_p),colnames(ra_heal)%in% rownames(ids_p)]



write.table(ids_m,file = "../table/1023_RAinteractionWIthBipolaris_heal.txt",sep = '\t',row.names = T,quote = F)
write.table(ids_p,file = "../table/1023_RAinteractionWIthBipolaris_adjustedP_heal.txt",sep = '\t',row.names = T,quote = F)

# ids = ra_heal_p[,colnames(ra_heal_p) %in% "g__Bipolaris"] <0.01
# # ids = ra_heal_p[,colnames(ra_heal_p) %in% "g__Furasium"] <0.01
# 
# ra_heal[ids,ids]
# ra_heal_p[ids,ids]
# 
# 
# write.table(ra_heal[ids,ids],file = "../table/RA_interactionWIthBipolaris_heal.txt",sep = '\t',row.names = T,quote = F)
# write.table(ra_heal_p[ids,ids],file = "../table/RA_interactionWIthBipolaris_adjustedP_heal.txt",sep = '\t',row.names = T,quote = F)



ggsave(pheatmap(ra_heal_sper_mat_top30),filename = "ra_heal_occor_spearman_top30.pdf")
a<- as.data.frame(ra_heal_occor_spearman$r)

pdf(file ="../plot/rootDisease/ra_heal_sper_mat_top50_upper_1105.pdf" )
corrplot(ra_heal_sper_mat_top30,p.mat = ra_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
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


png(file ="../plot/rootDisease/ra_heal_sper_mat_top50_upper_1023.png" )
corrplot(ra_heal_sper_mat_top30,p.mat = ra_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
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

##### RA  health without biolagris
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

write.table(ra_heal_sper_mat_top30,file = "../table/1025_RA_heal_top50_cor_without_Bipolaris.txt",sep = '\t',row.names = T,quote = F)
write.table(ra_heal_sper_p_top30,file = "../table/1025_RA_adjustedP_heal_top50_without_Bipolaris.txt",sep = '\t',row.names = T,quote = F)


ra_heal=ra_heal_occor_spearman$r
ra_heal_p=ra_heal_occor_spearman$p


# ids_p = as.data.frame(ra_heal_p[ra_heal_p[,colnames(ra_heal_p) %in% "g__Bipolaris"] <0.01,colnames(ra_heal_p) %in% "g__Bipolaris"])
# colnames(ids_p)= "adjusted-P"
# 
# 
# ids_m=ra_heal[rownames(ra_heal)%in% rownames(ids_p),colnames(ra_heal)%in% rownames(ids_p)]
# 
# 
# 
# write.table(ids_m,file = "../table/1023_RAinteractionWIthBipolaris_heal.txt",sep = '\t',row.names = T,quote = F)
# write.table(ids_p,file = "../table/1023_RAinteractionWIthBipolaris_adjustedP_heal.txt",sep = '\t',row.names = T,quote = F)

# ids = ra_heal_p[,colnames(ra_heal_p) %in% "g__Bipolaris"] <0.01
# # ids = ra_heal_p[,colnames(ra_heal_p) %in% "g__Furasium"] <0.01
# 
# ra_heal[ids,ids]
# ra_heal_p[ids,ids]
# 
# 
# write.table(ra_heal[ids,ids],file = "../table/RA_interactionWIthBipolaris_heal.txt",sep = '\t',row.names = T,quote = F)
# write.table(ra_heal_p[ids,ids],file = "../table/RA_interactionWIthBipolaris_adjustedP_heal.txt",sep = '\t',row.names = T,quote = F)



ggsave(pheatmap(ra_heal_sper_mat_top30),filename = "ra_heal_occor_spearman_top30.pdf")
a<- as.data.frame(ra_heal_occor_spearman$r)

pdf(file ="../plot/rootDisease/ra_heal_sper_mat_top50_upper_1025.pdf" )
corrplot(ra_heal_sper_mat_top30,p.mat = ra_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
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


png(file ="../plot/rootDisease/ra_heal_sper_mat_top50_upper_1023.png" )
corrplot(ra_heal_sper_mat_top30,p.mat = ra_heal_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
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

####Kendall test 

name_vector_dis <- c( "Promicromonospora" ,  "Phyllobacterium" ,    "Pseudomonas"  ,      
                  "Rhizobium"        ,   "Flavobacterium"    ,  "Sphingomonas"  ,     
                  "Variovorax"    ,      "Pseudoxanthomonas"  ,  "Streptomyces" ,   
                  "Chryseobacterium"  ,  "Lentzea"         ,    "Nocardioides"  ,     
                  "Actinoplanes"    ,    "Herbaspirillum"    ,  
                  "Mycobacterium"  ,     "Kribbella"        ,   "Massilia"   ,        
                  "Devosia"        ,     "Microbacterium"     ,
                  "Sphingopyxis"    ,    "Stenotrophomonas"  ,  "Mesorhizobium" ,     
                  "Caulobacter"   ,      "Aeromicrobium"   ,   
                  "Sorangium"     ,      "Agromyces"    ,       "Steroidobacter"   ,  
                  "Arthrobacter"  ,      "Chondromyces"   ,    
                  "Umezawaea"   ,        "Pantoea"    ,        "Rhodanobacter"    ,  
                  "Shinella"   ,         "Lysobacter"      ,   
                  "Solirubrobacter"    , "Bosea"          ,     "Glycomyces"    ,     
                  "Cellulomonas"     ,   "Bradyrhizobium"    ,  "Burkholderia"   ,    
                  "Rhodopseudomonas" ,   "Bacillus"    ,       
                  "g__Acremonium" ,"g__Lachnella"      ,  "g__Scytalidium",
                  "g__Magnaporthiopsis" ,"g__Ophiosphaerella","g__Gibberella",   
                  "g__Pyrenochaetopsis", "g__Sarocladium",  "g__Bipolaris" 
)

name_vector_heal <- c("Variovorax"    ,     "Flavobacterium"   ,  "Actinoplanes"    ,   "Lentzea"   ,         "Streptomyces"  ,    
                 "Pseudoxanthomonas" , "Rhizobium"     ,     "Promicromonospora" , "Nocardioides"  ,         
                 "Pseudomonas"    ,    "Massilia"      ,     "Umezawaea"    ,      "Devosia"   ,         "Galbitalea" ,       
                 "Sphingomonas"    ,   "Mesorhizobium"   ,   "Herbaspirillum"   ,  "Caulobacter"  ,      "Phyllobacterium" ,  
                 "Kribbella"   ,       "Arthrobacter"      , "Shinella"     ,      "Aeromicrobium"  ,    "Kineosporia"     ,  
                 "Chondromyces" ,      "Rudaibacter"   ,       "Steroidobacter"  ,  
                 "Stenotrophomonas"  , "Pantoea"    ,            "Ideonella"  ,       
                 "Rhodanobacter"    ,  "Bradyrhizobium"    , "Phycicoccus"     ,   "Chryseobacterium"  , "Agromyces"  ,       
                 "Microbacterium"   ,  "Solirubrobacter" ,   "Mycobacterium"    ,  "Bosea"       ,       "Altererythrobacter",
                 "Aquabacterium"   ,     "Curtobacterium"   ,     "Dactylosporangium" ,
                 "g__Podospora" , "g__Cladosporium" ,"g__Mycosphaerella", "g__Nigrospora"    , "g__Mortierella" ,"g__Edenia"  ,
                 "g__Metacordyceps"
)

dis_inter_heal<- intersect(name_vector_dis,name_vector_heal)
rank_dh <- seq(1,34,1)
heal_inter_dis<- intersect(name_vector_heal,name_vector_dis)
rank_hd <- seq(1,34,1)
dh <- data.frame(dis_inter_heal,rank_dh)
rownames(dh) <- dh$dis_inter_heal
hd <- data.frame(heal_inter_dis,rank_hd)
rownames(hd) <- hd$heal_inter_dis

dh
hd
id <- match(rownames(dh),rownames(hd))
hd<- hd[id,]
install.packages("Kendall")
library(Kendall)
kentall.test<- Kendall(hd$rank_hd,dh$rank_dh)

write.table(hd,file = "../table/dis_inter_heal.txt",sep = '\t',row.names = T,quote = F)
write.table(dh,file = "../table/heal_inter_dis.txt",sep = '\t',row.names = T,quote = F)
# write.csv2(kentall.test,file = "../table/kentall_test_34.txt",row.names = T,quote = F)
lapply(kentall.test, write, "../table/kendall_test_34.txt", append=TRUE)




