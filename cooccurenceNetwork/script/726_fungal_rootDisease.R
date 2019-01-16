rm(list=ls())
options(warn=3)

########### setting the working directory and print it ###################
tem <- "rootDisease"
setwd("~/xiaoxuan/180724/TARDF16/script/")
print(paste("Your working directory is in",getwd()))


########### to import the plotting theme() function ########################
source("plot_function.R")
# install.packages("")
library(ggpubr)
library(vegan)
library(ggplot2)
library(reshape)
# library(multcomp)
library(ggsignif)
library("Biobase", quietly=T, warn.conflicts=F)
library("ggplot2", quietly=T, warn.conflicts=F)
library("gplots", quietly=T, warn.conflicts=F)
library("grid", quietly=T, warn.conflicts=F)
library("RColorBrewer", quietly=T, warn.conflicts=F)
library("reshape2", quietly=T, warn.conflicts=F)
library(dplyr)
library(tidyverse)

### output directory assigned to include the pics & tables########################
figures.dir <- paste("~/xiaoxuan/180724/TARDF16/plot/",tem,'/',sep = '')
table.dir <- paste("~/xiaoxuan/180724/TARDF16/table/",tem,'/',sep = '')


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
design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
sample_list <- as.matrix(row.names(design))

# ids <- design$Other %in% c("2:2:2","1:1:1","2:2:1")
# 
# design <- design[!ids,]
# design <- design[ids,]


#2 fungalteria_list to filter 
# fungl_li = read.delim("doc/Funterial_list.txt",row.names = 1, header=T, sep="") 
# rownames(fungl_li)

#3 OTUs file

# otu_table = read.table("result/usearch97/otu_silva_with_spike.txt", sep = '\t',row.names= 1,header = T)
otu_table = read.table("../result/otu_table_tax.txt", sep = '\t',row.names= 1,header = T)
# tax <- read.table("../result/rep_seqs_tax_filtered.txt", sep = '\t',row.names= 1,header = T)

tax_filtered = read.table("../result/rep_seqs_tax_filtered.txt", sep = '\t',row.names= 1,header = T)
# otu_table_filtered_with_spike <- rbind(otu_table_filtered,otu_table[26,])
# otu_table <- otu_table_filtered_with_spike 

# otu_table = read.delim("result/usearch97/otu_silva_without_spike.txt", row.names= 1, sep="\t")






############################################################subSet OTU for analysis nature sample feature otu distribution ########################
sub_table <- otu_table[,colnames(otu_table) %in% rownames(design)]
sub_table <- sub_table[rownames(sub_table)%in% rownames(tax_filtered),]
# sub_table <- otu_table[,colnames(otu_table) %in%  rownames(design)]




############rarefraction#####################
## a consistent random seed [set.seed(21336)] was used for reproducibility

set.seed(21336)
# rrarefy The sample can be a vector giving the sample sizes for each row.so you need the transpose
##rarefy(x, sample, se = FALSE, MARGIN = 1) rrarefy(x, sample)
sub_table<- as.data.frame(t(rrarefy(t(sub_table),sample = min(colSums(sub_table)))))
write.table(sub_table,file = "../table/filtered_rarefy.txt",sep = '\t',row.names = T,quote = F)

### to check the size whether to be same
colSums(sub_table)


#################reorder######################
ids <- match(rownames(design),colnames(sub_table))
sub_table_1 <- sub_table[,ids]
ids <- match(rownames(design),colnames(otu_table))
otu_table_1 <- otu_table[,ids]

######################################Order by ZH11.Fun.BI055.01 for  Scal-BI-12-4 E05/5
sub_OTU<- sub_table_1[order(sub_table_1[,1],decreasing = T),]
OTU<-otu_table_1[order(otu_table_1[,1],decreasing = T),]

otu_tax <- merge(sub_OTU,otu_table[,c(1,ncol(otu_table))], by="row.names")

write.table(otu_tax,file = "../table/filtered_OTU_rarefy_tax.txt",sep = '\t',row.names = F,quote = F)
write.table(sub_OTU,file = "../table/filtered_OTU_rarefy_ordered.txt",sep = '\t',row.names = T,quote = F)

# ids <- match(row.names(OTU),row.names(otu_table))
# otu_table<- otu_table[ids,]
# 
# pso <- OTU[(rownames(OTU)) %in% "OTU_9",]
# # 
# # 
# top_10_OTU_withSpike <- sub_OTU[1:11,]
# 
# OTU_top_11 <- OTU[1:11,]
# OTU_top_11_Tax <- cbind(otu_table[1:11,]$taxonomy,OTU_top_11)
# colnames(OTU_top_11_Tax)[1] <- "Taxa"
# 
# OTU_top_11 <- rbind(OTU_top_11,pso)
# row.names(top_10_OTU_withSpike) <- OTU_top_11


######E05######
# design_1 <- design[!design$PlasmidID %in% "Scal-BI-mixture",]


# pos <-  colnames(OTU_top_11)%in% rownames(design_1)
# sub_OTU <- OTU[,pos]
# # rownames(sub_OTU) <- OTU_top_11_Tax$Taxa

batch_group <- unique(design$Other)               # plasmid id 
# gra_group <- unique(design$Description)
condition_group <- unique(design$Spikein)
geno_group <- unique(design$Genotype)                 # Concentration


### OTU_1 spike and filter oTU2 and OTU4 
#-----------------------------------------------------------------------------------------------------------------------------------------------
# spike 
# id <- grep("OTU_2$",rownames(sub_OTU))
# sub_OTU <- sub_OTU[-id,]
# pos <- grep("OTU_4$",rownames(sub_OTU))
# sub_OTU <- sub_OTU[-pos,]
pos <- grep("OTU_1$",rownames(sub_OTU))
rownames(sub_OTU)[pos] <- "BI-OS-12-4"
idx <- match(spike[2],row.names(sub_OTU))
# sub_table_1 <- sub_table[,sub_table[idx,]>160]
sub_OTU$id <- rownames(sub_OTU)


# mapping spike 12-4
# idx1 <- match(spike[1],row.names(sub_OTU))
sub_table_2 <- sub_OTU
# sub_table_2 <- sub_OTU[-idx1,]

# idx1 <- match(spike[3],row.names(sub_table_2))
sub_table_2 <- sub_table_2



len <- length(sub_table_2[1,])
idx <- match(spike[2],row.names(sub_table_2))
col_sums <- colSums(sub_table_2[-idx,-len])


sub_table_nor_t <- t(sub_table_2[-idx,-len])/col_sums
# sub_table_nor <- t(sub_table_nor_t)

sub_table_nor <- as.data.frame(t(sub_table_nor_t)*10000)

# remove_sample_id <- c("Fun2DMH6309","Fun2DMH6310","Fun2DMH6311","Fun2DMH6303","Fun2WMH6304","Fun2WMH6308","Fun2WWYJ11","Fun2WWYJ14","Fun2WWYJ04","Fun2WWYJ10","FunWMH6301","FunWMH6307","FunDWYJ05","FunDWYJ07","FunDWYJ02","FunDWYJ10","FunWWYJ07","FunWWYJ14")
# rid <- match(remove_sample_id,colnames(sub_table_nor))
# sub_table_nor <- as.data.frame(sub_table_nor[,-rid])

write.table(sub_table_nor,file = "~/xiaoxuan/180724/TARDF16/RA/otutab_norm_RA.txt",sep = '\t',row.names = T,quote = F)


spikeAbun <- sub_table_2[idx,-len]/colSums(sub_table_2[,-len])
row.names(spikeAbun) <- "spike-in_Abundance"

delta <- col_sums/sub_table_2[idx,-len]
row.names(delta) <- "delta_value"

total_table<- t(rbind(sub_table_nor,spikeAbun,delta))
t_table <- as.data.frame(total_table)
t_table$id <- row.names(total_table)

idx2 <- match("spike-in_Abundance",colnames(t_table))
index <- cbind(total_table[, idx2], design[match(t_table[,length(t_table[1,])], design$SampleID), ])
colnames(index)[1] <- "value"




colors <- data.frame(group=geno_group,
                     color=c(c_red, c_dark_brown,c_black))

shapes <- data.frame(group=condition_group,
                     shape=c(19, 0))

# index$Condition <- factor(index$Condition, levels=shapes$group)
# index$Description <- factor(index$Description, levels=colors$group)
# index$Description <- NULL
# index$Mix_Ratio <- index$Other
# reorder boxplots

# l <- c("soil", "rhizosphere", "root", "pooled_nodules")
# index$compartment <- factor(index$compartment, levels=l)
# colors <- colors[match(l, colors$group), ]


detailName <- paste("rootDisease_fungal_Community_",spike[2],sep = '_')

p <- ggplot(index, aes(x=Genotype, y=value, color=Genotype)) +
  facet_wrap(~Spikein)+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=detailName, y="spike-in_Abundance") +
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme


p

# ggsave(paste(figures.dir, detailName,"spike_in_Abundance.pdf", sep="_"), p)
# ggsave("../Fun_New/figure/silva_spike_in_Abundance.pdf", p)

## to save in the new directory
ggsave("../plot/filtered_spike_in_Abundance.pdf",  p)

write.table(index,file = "../table/filtered_Facet_Natural_spikeAbun.txt",sep = '\t',row.names = T,quote = F)








############plot RA 

taxonomy = read.delim("../result/rep_seqs_tax_filtered.txt", row.names= 1,header=F, sep="\t",stringsAsFactors = F)
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","evalue","unknown")
taxonomy$full=taxonomy$kingdom


write.table(paste("../table/taxonomy\tSampAvsB\tPvalue",sep="\t"), file=paste("phylumpro.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)
# # select p__ProteoFunteria line
# idx=taxonomy$phylum=="p__ProteoFunteria"
# # 初始化full为门，并初化因子为字符方便修改
taxonomy$full=as.character(taxonomy$phylum)
# # 修改pro门为目
# taxonomy[idx,]$full=as.character(taxonomy[idx,]$class)


sub_OTU$id <- NULL
tax_count = merge(taxonomy, sub_OTU, by="row.names")


tax_count_sum = aggregate(tax_count[,-c(1:11)], by=tax_count[11],FUN=sum) # mean
# tax_count_sum$full[1] <- "other"
rownames(tax_count_sum) = tax_count_sum$full
tax_count_sum = tax_count_sum[,-1]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100

# 绘制样品组内各样品堆叠图
mean_sort = per[(order(-rowSums(per))), ] # decrease sort
colSums(mean_sort)

# Stackplot
mean_sort=as.data.frame(mean_sort)

other = colSums(mean_sort[15:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(15-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[15] = c("Low Abundance")
# ordered taxonomy
write.table(rownames(mean_sort), file="tax_phylumpro.topN", append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)

if (TRUE){
  tax_name = gsub("[\\w;_]+__","",rownames(mean_sort),perl=TRUE)
  j=1
  for (i in 1:length(tax_name)){
    if (tax_name[i]==""){
      tax_name[i]=paste("Noname",j,sep='')
      j=j+1
    }
  }	
  rownames(mean_sort) = tax_name # rowname unallowed same name
}

mean_sort$phylumpro = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("phylumpro")))
#data_all$phylumpro  = factor(data_all$phylumpro, levels=rownames(mean_sort))   # set taxonomy order by abundance, default alphabet
data_all = merge(data_all, design, by.x="variable", by.y = "row.names")
write.table(data_all, file="../table/phylum_RA.txt", append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)

if ("FALSE" == "TRUE") {
  data_all$group  = factor(data_all$group, levels=c("A50","A56","A58","D77","D86","I3703","IR24","soil","ZH11"))   # set group order
}

p = ggplot(data_all[data_all$Funtch %in% "NO.1",], aes(x=variable, y = value, fill = phylumpro,group=variable)) + 
  geom_bar(stat = "identity")+ 
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ Other, scales = "free_x", switch = "x") +  
  theme(strip.Funkground = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Conditions in Samples ")+ylab("Relative abundance (%)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme
p

ggsave("../plot/phylum_RA_sample_Hn.pdf", p)


# ggsave("tax_stack_phylumpro_sample.png", p, width = 5, height = 3)
p = ggplot(data_all[data_all$Funtch %in% "NO.2",], aes(x=variable, y = value, fill = phylumpro,group=variable)) + 
  # geom_bar(stat = "identity",position="fill", width=1)+ 
  geom_bar(stat = "identity")+ 
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ Other, scales = "free_x", switch = "x") +  
  theme(strip.Funkground = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Conditions in Samples ")+ylab("Relative abundance (%)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme
p

ggsave("../plot/phylum_RA_sample_Ah.pdf", p)




# p <- ggplot(data_all[data_all$Funtch %in% "NO.1"&data_all$phylumpro %in% c("Funteroidetes","ProteoFunteria"),], aes(x=Other, y=value,color=Other,full=Other)) +
#   geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
#   geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
#   labs(x="Growth Condition_batch_NO1", y="Relative phylum abundance") +
#   # geom_signif(comparisons =list(c("E05/5", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.8,tip_length = 0, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E04", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E03", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
# #             map_signif_level = T)+
#   facet_wrap(phylumpro~Genotype)+
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   stat_compare_means(method="wilcox.test", size=3)+
#   main_theme
# 
# 
# p
# 
# ggsave("../plot/FunteroidetesProteoFunteria_NO1_stat.pdf", p)
# 
# 
# p <- ggplot(data_all[data_all$Funtch %in% "NO.2"&data_all$phylumpro %in% c("Funteroidetes","ProteoFunteria"),], aes(x=Other, y=value,color=Other,full=Other)) +
#   geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
#   geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
#   labs(x="Growth Condition_batch_NO2", y="Relative phylum abundance") +
#   # geom_signif(comparisons =list(c("E05/5", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.8,tip_length = 0, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E04", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E03", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
#   #             map_signif_level = T)+
#   facet_wrap(phylumpro~Genotype)+
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   stat_compare_means(method="wilcox.test", size=3)+
#   main_theme
# 
# 
# p
# ggsave("../plot/FunteroidetesProteoFunteria_NO2_stat.pdf", p)
# 
# 
# p <- ggplot(data_all[data_all$Funtch %in% "NO.1"&data_all$phylumpro %in% c("ActinoFunteria","Firmicutes"),], aes(x=Other, y=value,color=Other,full=Other)) +
#   geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
#   geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
#   labs(x="Growth Condition_batch_NO1", y="Relative phylum abundance") +
#   # geom_signif(comparisons =list(c("E05/5", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.8,tip_length = 0, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E04", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E03", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
# #             map_signif_level = T)+
# facet_wrap(phylumpro~Genotype)+
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   stat_compare_means(method="wilcox.test", size=3)+
#   main_theme
# 
# 
# p
# ggsave("../plot/ActinoFunteriaFirmicutes_NO1_stat.pdf", p)
# 
# 
# p <- ggplot(data_all[data_all$Funtch %in% "NO.2"&data_all$phylumpro %in% c("ActinoFunteria","Firmicutes"),], aes(x=Other, y=value,color=Other,full=Other)) +
#   geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
#   geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
#   labs(x="Growth Condition_batch_NO2", y="Relative phylum abundance") +
#   # geom_signif(comparisons =list(c("E05/5", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.8,tip_length = 0, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E04", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E03", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
# #             map_signif_level = T)+
# facet_wrap(phylumpro~Genotype)+
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   stat_compare_means(method="wilcox.test", size=3)+
#   main_theme
# 
# 
# p
# ggsave("../plot/ActinoFunteriaFirmicutes_NO2_stat.pdf", p)
# 
# 
# p <- ggplot(data_all[data_all$Funtch %in% "NO.1"&data_all$phylumpro %in% c("AcidoFunteria","Low Abundance"),], aes(x=Other, y=value,color=Other,full=Other)) +
#   geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
#   geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
#   labs(x="Growth Condition_batch_NO1", y="Relative phylum abundance") +
#   # geom_signif(comparisons =list(c("E05/5", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.8,tip_length = 0, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E04", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E03", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
# #             map_signif_level = T)+
#   facet_wrap(phylumpro~Genotype)+
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   stat_compare_means(method="wilcox.test", size=3)+
#   main_theme
# 
# 
# p
# ggsave("../plot/AcidoFunteria_Abundance_NO1_stat.pdf", p)
# 
# p <- ggplot(data_all[data_all$Funtch %in% "NO.2"&data_all$phylumpro %in% c("AcidoFunteria","Low Abundance"),], aes(x=Other, y=value,color=Other,full=Other)) +
#   geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
#   geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
#   labs(x="Growth Condition_batch_NO2", y="Relative phylum abundance") +
#   # geom_signif(comparisons =list(c("E05/5", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.8,tip_length = 0, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E04", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E03", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
# #             map_signif_level = T)+
#   facet_wrap(phylumpro~Genotype)+
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   stat_compare_means(method="wilcox.test", size=3)+
#   main_theme
# 
# 
# p
# ggsave("../plot/AcidoFunteria_Abundance_NO2_stat.pdf", p)
# 
# p <- ggplot(data_all[data_all$Funtch %in% "NO.1"&data_all$phylumpro %in% c("Chloroflexi","Tenericutes"),], aes(x=Other, y=value,color=Other, fill=Other)) +
#   geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
#   geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
#   labs(x="Growth Condition_batch_NO1", y="Relative phylum abundance") +
#   # geom_signif(comparisons =list(c("E05/5", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.8,tip_length = 0, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E04", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E03", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
# #             map_signif_level = T)+
#   facet_wrap(phylumpro~Genotype)+
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   stat_compare_means(method="wilcox.test", size=3)+
#   main_theme
# 
# 
# p
# ggsave("../plot/ChloroflexiTenericutes_NO1_stat.pdf", p)
# 
# p <- ggplot(data_all[data_all$Funtch %in% "NO.2"&data_all$phylumpro %in% c("Chloroflexi","Tenericutes"),], aes(x=Other, y=value,color=Other, fill=Other)) +
#   geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
#   geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
#   labs(x="Growth Condition_batch_NO2", y="Relative phylum abundance") +
#   
#   # geom_signif(comparisons =list(c("E04", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
#   #             map_signif_level = T)+
#   # geom_signif(comparisons =list(c("E03", "E00")),
#   #             annotations = "No Sig.",
#   #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
# #             map_signif_level = T)+
#   facet_wrap(phylumpro~Genotype)+
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   # compare_means(value ~ Other)+
#   # stat_compare_means()+
#   stat_compare_means(method="wilcox.test", size=3)+
#   # geom_signif(comparisons =list(c("dry", "wet")),
#   #             # annotations = "No Sig.",
#   #             test = "wilcox.test",
#   #             y_position = 7.6,tip_length = 0, vjust=0.4,
#   #             map_signif_level = T)+
#   main_theme
# 
# 
# p
# ggsave("../plot/ChloroflexiTenericutes_NO2_stat.pdf", p)
# 


# plot absolute abundance 
qpcr <- read.delim("../doc/TArootdisease_qPCR.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)

# mapping spike 12-4
# idx1 <- match(spike[1],row.names(sub_OTU))
sub_table_2 <- sub_OTU

# idx1 <- match(spike[3],row.names(sub_table_2))
# sub_table_2 <- sub_table_2

len <- length(sub_table_2[1,])
idx <- match(spike[2],row.names(sub_table_2))
# col_sums <- colSums(sub_table_2[-idx,-len])

# idx <- match(spike[3],row.names(sub_table))
# len <- length(sub_table[1,])
internal_ref <- as.numeric(sub_table_2[idx,-len])
internal_ref=1+internal_ref
absoluteAbund <-as.data.frame(t(sub_table_2[-idx,-len])/internal_ref)

# 
# ####absAbundance relative to spike-in
# absAbundance <- sweep(sub_table_2,2,as.numeric(sub_table_f[idx,]),'/')
# 
# ####### reorder ##########
# 
# ord <- match(rownames(sub_design_filt),colnames(absAbundance))   ## column rearrange
# # ord1 <- match(rownames(Fun_li), rownames(sub_table))    ## row reorder
# absAbundance<- absAbundance[,ord]   ## parrellel
# write.table(absAbundance,file = "./table/absAbundance.xls",sep = '\t',row.names = T)
ord <- match(rownames(design),rownames(absoluteAbund))   ## column rearrange
# # ord1 <- match(rownames(Fun_li), rownames(sub_table))    ## row reorder
absoluteAbund<-absoluteAbund[ord,]   ## parrellel
Fun_reads_AA <- as.data.frame(t(absoluteAbund))
write.table(Fun_reads_AA,file = "../table/absAbundance.xls",sep = '\t',row.names = T)

micro_load <- as.data.frame(t(colSums(Fun_reads_AA)))
micro_load<- micro_load[,colnames(micro_load) %in% c("TAD1","TAD2","TAD3","TAD4","TAD5","TAD6","TAD7","TAD8","TAD9","TAD10","TAH1","TAH2","TAH3","TAH4","TAH5","TAH6","TAH7","TAH8","TAH9","TAH10")]


idp <- match(colnames(micro_load),row.names(qpcr))
qpcr_ratio<- qpcr[idp,]



micro_load <- as.data.frame(t(sweep(micro_load, 2, qpcr_ratio$Copynumber, "/")))
data_all = merge(micro_load, design, by = "row.names")
data_all$Genotype
colnames(data_all)[2] <- "Microbiome_Load"
write.table(data_all,file = "../table/modify_adjusted_fungi_microbiome_load.txt",sep = '\t',row.names = T,quote = F)



# micro_load <- as.data.frame(colSums(Fun_reads_AA))
# idp <- match(row.names(qpcr),rownames(micro_load))
# ratio<- as.data.frame(qpcr[idp,1])
# micro_load <- micro_load/ratio[1:197,1]
# data_all = merge(micro_load, design, by = "row.names")
# colnames(data_all)[2] <- "Microbiome_Load"
# write.table(data_all,file = "../table/adjusted_Fun_microbiome_load.txt",sep = '\t',row.names = T,quote = F)

# remove_sample_id <- c("Fun2WWYJ07","Fun2DMH6309","Fun2DMH6310","Fun2DMH6311","Fun2DMH6303","Fun2WMH6304","Fun2WMH6308","Fun2WWYJ11","Fun2WWYJ14","Fun2WWYJ04","Fun2WWYJ10","FunWMH6301","FunWMH6307","FunDWYJ05","FunDWYJ07","FunDWYJ02","FunDWYJ10","FunWWYJ07","FunWWYJ14")

# p = ggplot(data_all[data_all$site %in% "Hn"& data_all$Genotype %in% c("MH63","WYJ7DEP1")&!data_all$Row.names %in% remove_sample_id,], aes(x=Other, y = Microbiome_Load,color=Other, fill=Other)) + 
data_all<- data_all %>% drop_na()

remove_sample_id <- c("TAD1","TAD5","TAD8")

p = ggplot(data_all[!data_all$Row.names %in% remove_sample_id,], aes(x=Genotype, y = Microbiome_Load,color=Genotype, fill=Genotype)) + 
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Genotype), position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~Other , scales = "free_x", switch = "x") +  
  # theme(strip.Funkground = element_blank())+
  # 关闭x轴刻度和标签
  # theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups_Batch_Hainan")+ylab("Microbiome_Fungi Load (Total Reads / Spike-in)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  stat_compare_means(method="wilcox.test", size=3)+
  # geom_signif(comparisons =list(c("dry", "wet")),
  #             # annotations = "No Sig.",
  #             test = "wilcox.test",
  #             y_position = 7.6,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  main_theme
p

ggsave("../plot/microbiome_load_DiseaseVsHealthy_fungi_removedTAD158.pdf", p)


# p = ggplot(data_all[data_all$Site %in% "AnHui"& data_all$Genotype %in% c("MH63","WYJ7DEP.1")&!data_all$Row.names %in% remove_sample_id,], aes(x=Condition, y = Microbiome_Load,color=Condition, fill=Condition)) + 
#   geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
#   geom_jitter(aes(shape=Condition), position=position_jitter(0.17), size=1, alpha=0.7) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   # scale_y_continuous(labels = scales::percent) + 
#   # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
#   facet_grid( ~Genotype , scales = "free_x", switch = "x") +  
#   # theme(strip.Funkground = element_blank())+
#   # 关闭x轴刻度和标签
#   # theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
#   xlab("Groups_Batch_Ah")+ylab("Microbiome_Fungi Load (Total Reads / Spike-in)")+
#   theme(axis.text.x = element_text(size=7.5,angle = 90))+
#   stat_compare_means(method="wilcox.test", size=3)+
#   # geom_signif(comparisons =list(c("dry", "wet")),
#   #             # annotations = "No Sig.",
#   #             test = "wilcox.test",
#   #             y_position = 7.6,tip_length = 0, vjust=0.4,
#   #             map_signif_level = T)+
#   main_theme
# p
# 
# 
# ggsave("../plot_719_remove_sample/microbiome_load_Ah.pdf", p)


###AA taxa plot
qpcr <- read.delim("../doc/TArootdisease_qPCR.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
# idp <- match(colnames(Fun_reads_AA),row.names(qpcr))
# ratio<- as.data.frame(qpcr)

Fun_reads_AA<- Fun_reads_AA[,colnames(Fun_reads_AA) %in% c("TAD1","TAD2","TAD3","TAD4","TAD5","TAD6","TAD7","TAD8","TAD9","TAD10","TAH1","TAH2","TAH3","TAH4","TAH5","TAH6","TAH7","TAH8","TAH9","TAH10")]

idp <- match(colnames(Fun_reads_AA),row.names(qpcr))
qpcr_ratio<- qpcr[idp,]
Fun_reads_AA_cor <- sweep(Fun_reads_AA, 2, qpcr_ratio$Copynumber, "/")

Fun_reads_AA_cor <- Fun_reads_AA_cor * 1000000+0.0001

colSums(Fun_reads_AA_cor)

write.table(Fun_reads_AA_cor,file = "~/xiaoxuan/180724/TARDF16/table/adjusted_absAbundance.txt",sep = '\t',row.names = T,quote = F)
write.table(Fun_reads_AA_cor,file = "~/xiaoxuan/180724/TARDF16/spikein/otutab_norm_AA.txt",sep = '\t',row.names = T,quote = F)




# write.table(Fun_reads_AA_cor,file = "../table/adjusted_Fun_reads_AA_needToChangeLabelOfMh63andZH.txt",sep = '\t',row.names = T,quote = F)
tax_count = merge(taxonomy, Fun_reads_AA_cor, by="row.names")


tax_count_sum = aggregate(tax_count[,-c(1:11)], by=tax_count[11],FUN=sum) # mean
# tax_count_sum$full[1] <- "other"
rownames(tax_count_sum) = tax_count_sum$full
tax_count_sum = tax_count_sum[,-1]

# per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100

# 绘制样品组内各样品堆叠图
mean_sort = tax_count_sum[(order(-rowSums(tax_count_sum))), ] # decrease sort
colSums(mean_sort)

# Stackplot
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[15:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(15-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[15] = c("Low Abundance")
# ordered taxonomy
write.table(rownames(mean_sort), file="AA_phylumpro.txt", append = FALSE, sep="\t", quote=F, row.names=T, col.names=T)

if (TRUE){
  tax_name = gsub("[\\w;_]+__","",rownames(mean_sort),perl=TRUE)
  j=1
  for (i in 1:length(tax_name)){
    if (tax_name[i]==""){
      tax_name[i]=paste("Noname",j,sep='')
      j=j+1
    }
  }	
  rownames(mean_sort) = tax_name # rowname unallowed same name
}

mean_sort$phylumpro = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("phylumpro")))
#data_all$phylumpro  = factor(data_all$phylumpro, levels=rownames(mean_sort))   # set taxonomy order by abundance, default alphabet
data_all = merge(data_all, design, by.x="variable", by.y = "row.names")


p = ggplot(data_all[data_all$Funtch %in% "NO.1",], aes(x=variable, y = value, fill = phylumpro,group=variable)) + 
  geom_bar(stat = "identity")+ 
  # geom_bar(stat = "identity",position="fill")+ 
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ Other, scales = "free_x", switch = "x") +  
  theme(strip.Funkground = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups")+ylab("Absolute phylum abundance (reads per unit host DNA)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme
p

ggsave("../plot/AA_phylum_stack_phylumpro_sample_Hn.pdf", p)

p = ggplot(data_all[data_all$Funtch %in% "NO.2",], aes(x=variable, y = value, fill = phylumpro,group=variable)) + 
  geom_bar(stat = "identity")+ 
  # geom_bar(stat = "identity",position="fill")+ 
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ Other, scales = "free_x", switch = "x") +  
  theme(strip.Funkground = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups")+ylab("Absolute phylum abundance (reads per unit host DNA)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme
p
ggsave("../plot/AA_phylum_stack_phylumpro_sample_Ah.pdf", p)


######AA taxa statistics

p <- ggplot(data_all[data_all$Funtch %in% "NO.1"&data_all$phylumpro %in% c("Funteroidetes","ProteoFunteria"),], aes(x=Other, y=value,color=Other,full=Other)) +
  geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  labs(x="Growth Condition_batch_NO1", y="Absolute phylum abundance") +
  # geom_signif(comparisons =list(c("E05/5", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.8,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E04", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E03", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
#             map_signif_level = T)+
facet_wrap(phylumpro~Genotype)+
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  stat_compare_means(method="wilcox.test", size=3)+
  main_theme


p

ggsave("../plot/AA_FunteroidetesProteoFunteria_NO1_stat.pdf", p)


p <- ggplot(data_all[data_all$Funtch %in% "NO.2"&data_all$phylumpro %in% c("Funteroidetes","ProteoFunteria"),], aes(x=Other, y=value,color=Other,full=Other)) +
  geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  labs(x="Growth Condition_batch_NO2", y="Relative phylum abundance") +
  # geom_signif(comparisons =list(c("E05/5", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.8,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E04", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E03", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
#             map_signif_level = T)+
facet_wrap(phylumpro~Genotype)+
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  stat_compare_means(method="wilcox.test", size=3)+
  main_theme


p
ggsave("../plot/AA_FunteroidetesProteoFunteria_NO2_stat.pdf", p)


p <- ggplot(data_all[data_all$Funtch %in% "NO.1"&data_all$phylumpro %in% c("ActinoFunteria","Firmicutes"),], aes(x=Other, y=value,color=Other,full=Other)) +
  geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  labs(x="Growth Condition_batch_NO1", y="Relative phylum abundance") +
  # geom_signif(comparisons =list(c("E05/5", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.8,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E04", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E03", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
#             map_signif_level = T)+
facet_wrap(phylumpro~Genotype)+
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  stat_compare_means(method="wilcox.test", size=3)+
  main_theme


p
ggsave("../plot/AA_ActinoFunteriaFirmicutes_NO1_stat.pdf", p)


p <- ggplot(data_all[data_all$Funtch %in% "NO.2"&data_all$phylumpro %in% c("ActinoFunteria","Firmicutes"),], aes(x=Other, y=value,color=Other,full=Other)) +
  geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  labs(x="Growth Condition_batch_NO2", y="Relative phylum abundance") +
  # geom_signif(comparisons =list(c("E05/5", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.8,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E04", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E03", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
#             map_signif_level = T)+
facet_wrap(phylumpro~Genotype)+
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  stat_compare_means(method="wilcox.test", size=3)+
  main_theme


p
ggsave("../plot/AA_ActinoFunteriaFirmicutes_NO2_stat.pdf", p)


p <- ggplot(data_all[data_all$Funtch %in% "NO.1"&data_all$phylumpro %in% c("AcidoFunteria","Low Abundance"),], aes(x=Other, y=value,color=Other,full=Other)) +
  geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  labs(x="Growth Condition_batch_NO1", y="Relative phylum abundance") +
  # geom_signif(comparisons =list(c("E05/5", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.8,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E04", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E03", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
#             map_signif_level = T)+
facet_wrap(phylumpro~Genotype)+
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  stat_compare_means(method="wilcox.test", size=3)+
  main_theme


p
ggsave("../plot/AA_AcidoFunteria_Abundance_NO1_stat.pdf", p)

p <- ggplot(data_all[data_all$Funtch %in% "NO.2"&data_all$phylumpro %in% c("AcidoFunteria","Low Abundance"),], aes(x=Other, y=value,color=Other,full=Other)) +
  geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  labs(x="Growth Condition_batch_NO2", y="Relative phylum abundance") +
  # geom_signif(comparisons =list(c("E05/5", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.8,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E04", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E03", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
#             map_signif_level = T)+
facet_wrap(phylumpro~Genotype)+
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  stat_compare_means(method="wilcox.test", size=3)+
  main_theme


p
ggsave("../plot/AA_AcidoFunteriaLow_Abundance_NO2_stat.pdf", p)

p <- ggplot(data_all[data_all$Funtch %in% "NO.1"&data_all$phylumpro %in% c("Chloroflexi","Tenericutes"),], aes(x=Other, y=value,color=Other, fill=Other)) +
  geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  labs(x="Growth Condition_batch_NO1", y="Relative phylum abundance") +
  # geom_signif(comparisons =list(c("E05/5", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.8,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E04", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E03", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
#             map_signif_level = T)+
facet_wrap(phylumpro~Genotype)+
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  stat_compare_means(method="wilcox.test", size=3)+
  main_theme


p
ggsave("../plot/AA_ChloroflexiTenericutes_NO1_stat.pdf", p)

p <- ggplot(data_all[data_all$Funtch %in% "NO.2"&data_all$phylumpro %in% c("Chloroflexi","Tenericutes"),], aes(x=Other, y=value,color=Other, fill=Other)) +
  geom_boxplot(alpha=0.51, outlier.size=0, size=0.7, width=0.1) +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  labs(x="Growth Condition_batch_NO2", y="Relative phylum abundance") +
  
  # geom_signif(comparisons =list(c("E04", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.73,tip_length = 0.01, vjust=0.4,
  #             map_signif_level = T)+
  # geom_signif(comparisons =list(c("E03", "E00")),
  #             annotations = "No Sig.",
  #             y_position = 4.65,tip_length = 0.01, vjust=0.4,
  #             map_signif_level = T)+
  facet_wrap(phylumpro~Genotype)+
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  # compare_means(value ~ Other)+
  # stat_compare_means()+
  stat_compare_means(method="wilcox.test", size=3)+
  # geom_signif(comparisons =list(c("dry", "wet")),
  #             # annotations = "No Sig.",
  #             test = "wilcox.test",
  #             y_position = 7.6,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  main_theme


p
ggsave("../plot/AA_ChloroflexiTenericutes_NO2_stat.pdf", p)


# #### beta diversity #######
# # source("https://bioconductor.org/biocLite.R")
# # biocLite("metagenomeSeq")
# # 
# # library(metagenomeSeq)
# # library(biomformat)
# 
# 
# ## Distances between samples was calculated as Bray-Curtis dissimilarities of subsampled and square-root transformed OTU counts

##filter spike 
otu_table_for_beta <- sub_table_2[-1,-ncol(sub_table_2)]
otu_table_for_beta <- sweep(otu_table_for_beta,2,colSums(otu_table_for_beta),FUN = '/')

otu_table_for_beta<- cbind(sub_table_2[-1,]$id,otu_table_for_beta)
colnames(otu_table_for_beta)[1] <- "OTU_id"
# otu_choose <- otu_table_for_beta[1:20,]


write.table(otu_table_for_beta,file = "../table/otu_table_for_beta.txt",sep = '\t',quote = F,row.names = F,col.names = T)

## 
# biom convert -i otu_table_for_beta.txt -o otu_table_for_beta.biom --table-type="OTU table" --to-json
# 
# beta_diversity.py  -i otu_table_for_beta.biom -m bray_curtis -o beta


#### CSS normalization
# 
# loadMeta("../Fun_New/table/otu_table.biom")
# 
# "CSS" <-function(input_path, out_path, output_CSS_statistics=NULL) {
#   obj = load_biom(input_path)
#   p = cumNormStatFast(obj)
#   obj = cumNorm(obj, p = p)
#   if (!is.null(output_CSS_statistics)) {
#     exportStats(obj, p=p, file = file.path(output_CSS_statistics))
#   }
#   write_biom(MRexperiment2biom(obj, norm=TRUE, log=TRUE), out_path)
# }
# 
# CSS("../Fun_New/table/otu_table.biom","../Fun_New/table/otu_table_css.biom" , 1)
# 
# "CSS" <-function(input_path, out_path, output_CSS_statistics=NULL) {
#   obj = load_biom("../Fun_New/table/otu_table.biom")
#   p = cumNormStatFast(obj)
#   obj = cumNorm(obj, p = p)
#   if (!is.null(output_CSS_statistics)) {
#     exportStats(obj, p=p, file = file.path(output_CSS_statistics))
#   }
#   write_biom(MRexperiment2biom(obj, norm=TRUE, log=TRUE), out_path)
# }


# bray_curtis<-as.data.frame(vegdist(t(sub_table_2[,-ncol(sub_table_2)]+1),method="bray"))
# pcoa<-capscale(dis~1)

###### after beta_div calculation by qiime 
###### RA  but not the CSS or DESEQ2 method better 

bray_curtis <- read.table("../spikein/beta_716/bray_curtis_otutab_1.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../Fun_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../Fun_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 

design = design[rownames(design) %in% rownames(bray_curtis),]
id <- match(row.names(design[design$Site %in% "AnHui",]),rownames(bray_curtis))



ah_bray_curtis <- bray_curtis[id,id]
design <- design[rownames(design) %in% rownames(ah_bray_curtis),]
##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~1
k <- 3
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y","z")
# points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Fun_spike_perturbation"),each=3,4)

points <- cbind(points, design[match(rownames(points), design$SampleID), ])



colors <- data.frame(group=geno_group,
                     color=c(c_red, c_dark_brown,c_black,c_orange))

shapes <- data.frame(group=condition_group,
                     shape=c(19, 0))





# points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # points$genotype <- factor(points$genotype, levels=shapes$group)
# points$GroupID <- factor(points$GroupID, levels=colors$group)

# plot PCo 1 and 2

# p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Site)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   main_theme +
#   theme(legend.position="top")
# p
# 
# ggsave("../spikein/PCoA12_BC_Site.pdf", p)

p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Condition)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  # facet_wrap(~ Funtch)+
  main_theme +
  theme(legend.position="top")
p


ggsave("../spikein/PCoA12_BC_Condition_Ah.pdf", p)



# p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Site)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
#   main_theme +
#   theme(legend.position="top")
# p
# 
# ggsave("../spikein/PCoA13_BC_Site.pdf", p)

p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Condition)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p


ggsave("../spikein/PCoA13_BC_Condition_Ah.pdf", p)
ggsave("../Fun_New/figure/eulcian_PCoA_BC.pdf", p)
write.table(points,file = "../table/WD_Natural_fungi_betaDiv_points.txt",sep = '\t',row.names = F,quote = F)


###AA HaiNan
bray_curtis <- read.table("../spikein/beta_716/bray_curtis_otutab_1.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../Fun_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../Fun_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 
design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
design = design[rownames(design) %in% rownames(bray_curtis),]
id <- match(row.names(design[design$Site %in% "HaiNan",]),rownames(bray_curtis))



ah_bray_curtis <- bray_curtis[id,id]
##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~1
k <- 3
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y","z")
# points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Fun_spike_perturbation"),each=3,4)

points <- cbind(points, design[match(rownames(points), design$SampleID), ])



colors <- data.frame(group=geno_group,
                     color=c(c_red, c_dark_brown,c_black,c_orange))

shapes <- data.frame(group=condition_group,
                     shape=c(19, 0))





# points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # points$genotype <- factor(points$genotype, levels=shapes$group)
# points$GroupID <- factor(points$GroupID, levels=colors$group)

# plot PCo 1 and 2

# p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Site)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   main_theme +
#   theme(legend.position="top")
# p
# 
# ggsave("../spikein/PCoA12_BC_Site.pdf", p)

p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Condition)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  # facet_wrap(~ Funtch)+
  main_theme +
  theme(legend.position="top")
p


ggsave("../spikein/PCoA12_BC_Condition_Hainan.pdf", p)



# p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Site)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   main_theme +
#   theme(legend.position="top")
# p
# 
# ggsave("../spikein/PCoA13_BC_Site.pdf", p)

p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Condition)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p


ggsave("../spikein/PCoA13_BC_Condition_HaiNan.pdf", p)
ggsave("../Fun_New/figure/eulcian_PCoA_BC.pdf", p)
write.table(points,file = "../table/WD_Natural_fungi_betaDiv_points.txt",sep = '\t',row.names = F,quote = F)


######RA method
design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
bray_curtis <- read.table("../RA/beta_div/bray_curtis_otutab_norm.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../Fun_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../Fun_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 
design = design[rownames(design) %in% rownames(bray_curtis),]
# id <- match(row.names(design[design$Spikein %in% "spike",]),rownames(bray_curtis))
id <- match(row.names(design[design$Spikein %in% "spike"&!design$SampleID %in%c("TAD1","TAD5","TAD8"),]),rownames(bray_curtis))



ah_bray_curtis <- bray_curtis[id,id]
design <- design[rownames(design) %in% rownames(ah_bray_curtis),]
##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~10
k <- 3
pcoa <- cmdscale(ah_bray_curtis, k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y","z")
# points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Fun_spike_perturbation"),each=3,4)

points <- cbind(points, design[match(rownames(points), design$SampleID), ])
points <- points[-nrow(points),]



colors <- data.frame(group=geno_group,
                     color=c(c_red, c_dark_brown,c_black))

shapes <- data.frame(group=condition_group,
                     shape=c(19,0))





# points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # points$genotype <- factor(points$genotype, levels=shapes$group)
# points$GroupID <- factor(points$GroupID, levels=colors$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Genotype)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("../RA/730_fungi_PCoA12_BC_removedTAD158.pdf", p)

# p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Other)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   facet_wrap(~ Funtch)+
#   main_theme +
#   theme(legend.position="top")
# p
# 
# 
# ggsave("../plot/PCoA12_BC_WD.pdf", p)



p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Genotype)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("../RA/730_fungi_PCoA13_BCremovedTAD158.pdf", p)

p <- ggplot(points, aes(x=y, y=z, color=Genotype,shape=Genotype)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("../RA/730_fungi_PCoA23_BCremovedTAD158.pdf", p)




design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
bray_curtis <- read.table("../RA_rem_sample_720/beta/bray_curtis_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
#
# bray_curtis <- read.table("../RA/beta/bray_curtis_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis <- read.table("../RA/beta/bray_curtis_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../Fun_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../Fun_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 
design = design[rownames(design) %in% rownames(bray_curtis),]
id <- match(row.names(design[design$Site %in% "HaiNan"&!design$Genotype %in% "MH63ZH",]),rownames(bray_curtis))



ah_bray_curtis <- bray_curtis[id,id]
##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~10
k <- 3
pcoa <- cmdscale(ah_bray_curtis, k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y","z")
# points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Fun_spike_perturbation"),each=3,4)

points <- cbind(points, design[match(rownames(points), design$SampleID), ])
points <- points[-nrow(points),]



colors <- data.frame(group=geno_group,
                     color=c(c_red, c_dark_brown,c_black))

shapes <- data.frame(group=condition_group,
                     shape=c(19, 0))





# points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # points$genotype <- factor(points$genotype, levels=shapes$group)
# points$GroupID <- factor(points$GroupID, levels=colors$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Condition)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("../RA_rem_sample_720/718_fungi_PCoA12_BC_Conditon_HaiNan.pdf", p)

# p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Other)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   facet_wrap(~ Funtch)+
#   main_theme +
#   theme(legend.position="top")
# p
# 
# 
# ggsave("../plot/PCoA12_BC_WD.pdf", p)



p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Condition)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("../RA_rem_sample_720/718_fungi_PCoA13_BC_Condition_Hainan.pdf", p)


p <- ggplot(points, aes(x=y, y=z, color=Genotype,shape=Condition)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("../RA_rem_sample_720/718_fungi_PCoA23_BC_Condition_Hainan.pdf", p)
# p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Other)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   main_theme +
#   theme(legend.position="top")
# p
# 
# 
# ggsave("../plot/PCoA13_BC_WD.pdf", p)

####### Cluster dendrogram ##############
bray_curtis <- read.table("../table/beta/bray_curtis_otu_table_for_beta.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../Fun_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eu


# dat_rare <- read.delim("PyroTagger_counts_rarefied_A.txt", header=T, row.names=1, sep="\t")
# # application of threshold for lower presence limit of relative abundance
# threshold <- 5    # 0.5% presence: use integer 5 considering x1000 normalization
# dat_rare_thres <- dat_rare[ apply(dat_rare,1,max) > threshold, ]
# # log2-transformation
# dat_rare_thres_log <- log2(dat_rare_thres + 1)
# dim(dat_rare_thres_log)
# calculate Bray-Curtis dissimilarity 
otu_table_for_beta$OTU_id <- NULL

# log2-transformation
# dat_rare_thres_log <- log2(dat_rare_thres + 1)
dat_dist <- vegdist(t(otu_table_for_beta), method="bray")
# clustering
dat_dist_clu <- hclust(dat_dist , "average")
# pdf("Figure_4a.pdf", width=15, height=5)
plot(dat_dist_clu, main = "Clutering Dentrogram of  Beta-diversity / Bray Curtis", lwd=1.5, cex=.5)





weighed_unifrac <- read.table("../Fun_New/table/beta_div_unifrac/weighted_unifrac_table.from_txt_json.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../Fun_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../Fun_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 

##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~1
k <- 2
pcoa <- cmdscale(weighed_unifrac , k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")
points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Fun_spike_perturbation"),each=3,4)

points <- cbind(points, design_1[match(rownames(points), design_1$SampleID), ])



colors <- data.frame(group=unique(points$Description),
                     color=c(c_red, c_very_dark_green, c_black ,c_dark_brown))


shapes <- data.frame(group=unique(points$Concentration),
                     shape=c(19, 10, 5))





# points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # points$genotype <- factor(points$genotype, levels=shapes$group)
# points$GroupID <- factor(points$GroupID, levels=colors$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=Description,shape=Concentration)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("Weighted_Unifrac_PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("Weighted_Unifrac_PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("../Fun_New/figure/PCoA_W_unifrac.pdf", p)

#####unweighted Unifrac

un_weighed_unifrac <- read.table("../Fun_New/table/beta_div_unifrac/unweighted_unifrac_table.from_txt_json.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../Fun_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../Fun_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 

##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~1
k <- 2
pcoa <- cmdscale(un_weighed_unifrac, k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")
points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Fun_spike_perturbation"),each=3,4)

points <- cbind(points, design_1[match(rownames(points), design_1$SampleID), ])



colors <- data.frame(group=unique(points$Description),
                     color=c(c_red, c_very_dark_green, c_black ,c_dark_brown))


shapes <- data.frame(group=unique(points$Concentration),
                     shape=c(19, 10, 5))





# points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # points$genotype <- factor(points$genotype, levels=shapes$group)
# points$GroupID <- factor(points$GroupID, levels=colors$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=Description,shape=Concentration)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("Weighted_Unifrac_PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("Weighted_Unifrac_PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("../Fun_New/figure/PCoA_W_unifrac.pdf", p
