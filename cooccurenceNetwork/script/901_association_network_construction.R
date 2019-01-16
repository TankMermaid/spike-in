rm(list=ls())
options(warn=3)

########### setting the working directory and print it ###################
tem <- "rootDisease"
setwd("~/xiaoxuan/180724/TARDF16/script/")
print(paste("Your working directory is in",getwd()))


########### to import the plotting theme() function ########################
source("plot_function.R")
# install.packages("ggpubr")
install.packages("corrplot")
install.packages("ggcorrplot")
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


### RA network construction
#####disease 
design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 

sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil"& design$Genotype %in% "Disease",]

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
write.table(merge_bac_fungi,file = "../table/RA_adjusted_bacAndFungiMerge_abundance_disease.txt",sep = '\t',row.names = T,quote = F)


######dplyr to data manipulation
# install.packages("psych")


ra_dis_a<- merge_abun_tax[,c(2:11,13,18)] %>%
  mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:11)])) %>%
  group_by(Kingdom,Genus) %>%
  summarize_all(sum) %>%
  arrange(desc(rowSumAbun))

ra_dis_a<- as.data.frame(ra_dis_a)
rownames(ra_dis_a)<- ra_dis_a$Genus
ra_dis_a$Kingdom <- NULL
ra_dis_a$Genus <- NULL
ra_dis_a <- as.data.frame(t(ra_dis_a))


#### note the most 300 threshold using Kendall and spearman 
ra_dis_occor_kendall = corr.test(ra_dis_a[1:300],use="pairwise",method="kendall",adjust="fdr",alpha=.05)
ra_dis_occor_kendall$r

ra_dis_occor_spearman <- corr.test(ra_dis_a[1:300],use="pairwise",method="spearman",adjust="fdr",alpha=.05)
ra_dis_occor_spearman$r

ra_dis_sper_mat_top30 <- ra_dis_occor_spearman$r[c(1:30),c(1:30)]
ra_dis_sper_p_top30 <- ra_dis_occor_spearman$p[c(1:30),c(1:30)]

pheatmap(ra_dis_sper_mat_top30)

corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .05,
         insig = "blank",order = "hclust", addrect = 2,type = "upper",
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")




#####health
design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 

sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil"& design$Genotype %in% "Healthy",]

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

ra_heal_a<- merge_abun_tax[,c(2:11,13,18)] %>%
  mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:11)])) %>%
  group_by(Kingdom,Genus) %>%
  summarize_all(sum) %>%
  arrange(desc(rowSumAbun))

ra_heal_a<- as.data.frame(ra_heal_a)
rownames(ra_heal_a)<- ra_heal_a$Genus
ra_heal_a$Kingdom <- NULL
ra_heal_a$Genus <- NULL
ra_heal_a <- as.data.frame(t(ra_heal_a))


#### note the most 300 threshold using Kendall and Spearman 
ra_heal_occor_kendall = corr.test(ra_heal_a[1:300],use="pairwise",method="kendall",adjust="fdr",alpha=.05)
ra_heal_occor_kendall$r

ra_heal_occor_spearman <- corr.test(ra_heal_a[1:300],use="pairwise",method="spearman",adjust="fdr",alpha=.05)
ra_heal_occor_spearman$r
ra_heal_sper_mat_top30 <- ra_heal_occor_spearman$r[c(1:30),c(1:30)]
ra_heal_sper_p_top30 <- ra_heal_occor_spearman$p[c(1:30),c(1:30)]

ggsave(pheatmap(ra_heal_sper_mat_top30),filename = "ra_heal_occor_spearman_top30.pdf")


corrplot(ra_heal_sper_mat_top30,p.mat = ra_heal_sper_p_top30, sig.level = .01,
         insig = "blank",order = "hclust", addrect = 2,type = "upper",
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")





######################## AA network construction


design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil"&design$Genotype %in% "Disease",]


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

aa_dis_a<- merge_abun_tax[,c(2:11,13,18)] %>%
  mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:11)])) %>%
  group_by(Kingdom,Genus) %>%
  summarize_all(sum) %>%
  arrange(desc(rowSumAbun))

aa_dis_a<- as.data.frame(aa_dis_a)
rownames(aa_dis_a)<- aa_dis_a$Genus
aa_dis_a$Kingdom <- NULL
aa_dis_a$Genus <- NULL
aa_dis_a <- as.data.frame(t(aa_dis_a))


#### note the most 210 threshold
aa_dis_occor_kendall = corr.test(aa_dis_a[1:300],use="pairwise",method="kendall",adjust="fdr",alpha=.05)
aa_dis_occor_kendall$r
# aa_occor_kendall$t
# aa_occor_kendall$n

aa_dis_occor_spearman <- corr.test(aa_dis_a[1:300],use="pairwise",method="spearman",adjust="fdr",alpha=.05)
aa_dis_occor_spearman$r

## to filter the top 30 genus of 
aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:30),c(1:30)]
aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:30),c(1:30)]

pheatmap(aa_dis_sper_mat_top30)

corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .01,
         insig = "blank",order = "hclust", addrect = 2,type = "upper",
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")


# bacVsFungiRatio <- sweep(tax_count_sum,,fungi,"/")


########## aa healthy
design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil"&design$Genotype %in% "Healthy",]


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

aa_heal_a<- merge_abun_tax[,c(2:11,13,18)] %>%
  mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:11)])) %>%
  group_by(Kingdom,Genus) %>%
  summarize_all(sum) %>%
  arrange(desc(rowSumAbun))

aa_heal_a<- as.data.frame(aa_heal_a)
rownames(aa_heal_a)<- aa_heal_a$Genus
aa_heal_a$Kingdom <- NULL
aa_heal_a$Genus <- NULL
aa_heal_a <- as.data.frame(t(aa_heal_a))


#### note the most 300 threshold
aa_heal_occor_kendall = corr.test(aa_heal_a[1:300],use="pairwise",method="kendall",adjust="fdr",alpha=.05)
aa_heal_occor_kendall$r


aa_heal_occor_spearman <- corr.test(aa_heal_a[1:300],use="pairwise",method="spearman",adjust="fdr",alpha=.05)
aa_heal_occor_spearman$r
aa_heal_sper_mat_top30 <- aa_heal_occor_spearman$r[c(1:30),c(1:30)]
aa_heal_sper_p_top30 <- aa_heal_occor_spearman$p[c(1:30),c(1:30)]

ggsave(pheatmap(aa_heal_sper_mat_top30),filename = "aa_heal_occor_spearman_top30.pdf")


corrplot(aa_heal_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .001,
         insig = "blank",order = "hclust", addrect = 2,type = "upper",
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")



