rm(list=ls())
options(warn=3)

########### setting the working directory and print it ###################
tem <- "nature_perturbation"
setwd("~/xiaoxuan/180213/180731_bac_pert/script")
print(paste("Your working directory is in",getwd()))


########### to import the plotting theme() function ########################
source("plot_function.R")
# install.packages("ggpubr")
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
# library("VennDiagram", quietly=T, warn.conflicts=F)

### output directory assigned to include the pics & tables########################
figures.dir <- paste("~/xiaoxuan/180213/180731_bac_pert/plot/",tem,'/',sep = '')
table.dir <- paste("~/xiaoxuan/180213/180731_bac_pert/table/",tem,'/',sep = '')


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
design = read.table("~/xiaoxuan/180213/180731_bac_pert/doc/zh11_design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
sample_list <- as.matrix(row.names(design))

#3 OTUs file

otu_table = read.table("~/xiaoxuan/180213/180731_bac_pert/result/otutab_header_not_comment.txt", sep = '\t',row.names= 1,header = T)
tax = read.table("~/xiaoxuan/180213/180731_bac_pert/result/taxonomy_8.txt", sep = '\t',row.names= 1,header = T)
# otu_table_filtered = read.table("../result/filtered_sorted_otutab_tax_siftOTU1.txt", sep = '\t',row.names= 1,header = T)
# otu_table_filtered_with_spike <- rbind(otu_table_filtered,otu_table[26,])
# otu_table <- otu_table_filtered_with_spike 

# otu_table = read.delim("result/usearch97/otu_silva_without_spike.txt", row.names= 1, sep="\t")






############################################################subSet OTU for analysis nature sample feature otu distribution ########################
sub_table <- otu_table[,colnames(otu_table) %in% rownames(design)]
# sub_table <- otu_table[,colnames(otu_table) %in%  rownames(design)]




############rarefraction#####################
## a consistent random seed [set.seed(21336)] was used for reproducibility

set.seed(21336)
# rrarefy The sample can be a vector giving the sample sizes for each row.so you need the transpose
##rarefy(x, sample, se = FALSE, MARGIN = 1) rrarefy(x, sample)
sub_table<- as.data.frame(t(rrarefy(t(sub_table),sample = min(colSums(sub_table)))))
write.table(sub_table,file = "~/xiaoxuan/180213/180731_bac_pert/table/filtered_rarefy.txt",sep = '\t',row.names = T,quote = F)

### to check the size whether to be same
colSums(sub_table)


#################reorder######################
ids <- match(rownames(design),colnames(sub_table))
sub_table_1 <- sub_table[,ids]
ids <- match(rownames(design),colnames(otu_table))
otu_table_1 <- otu_table[,ids]

######################################Order by ZH11.Bac.BI055.01 for  Scal-BI-12-4 E05/5
sub_OTU<- sub_table_1[order(sub_table_1[,1],decreasing = T),]
OTU<-otu_table_1[order(otu_table_1[,1],decreasing = T),]


otu_tax <- merge(sub_OTU,tax, by="row.names")

write.table(otu_tax,file = "~/xiaoxuan/180213/180731_bac_pert/table/filtered_OTU_rarefy_tax.txt",sep = '\t',row.names = F,quote = F)     # 无意修改了一下 无视就可以了
write.table(sub_OTU,file = "~/xiaoxuan/180213/180731_bac_pert/table/filtered_OTU_rarefy_ordered.txt",sep = '\t',row.names = T,quote = F)



batch_group <- unique(design$spike_concentration)               # plasmid id 
# gra_group <- unique(design$Description)
condition_group <- unique(design$Description)
geno_group <- unique(design$Genotype)                 # Concentration


### OTU_2 spike 
#-----------------------------------------------------------------------------------------------------------------------------------------------
# spike 
# # id <- grep("OTU_1$",rownames(sub_OTU))
# sub_OTU <- sub_OTU[-id,]
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
sub_table_nor <- t(sub_table_nor_t)


# remove_sample_id <- c("Bac2DMH6309","Bac2DMH6310","Bac2DMH6311","Bac2DMH6303","Bac2WMH6304","Bac2WMH6308","Bac2WWYJ11","Bac2WWYJ14","Bac2WWYJ04","Bac2WWYJ10","BacWMH6301","BacWMH6307","BacDWYJ05","BacDWYJ07","BacDWYJ02","BacDWYJ10","BacWWYJ07","BacWWYJ14")
# rid <- match(remove_sample_id,colnames(sub_table_nor))
# sub_table_nor <- as.data.frame(sub_table_nor[,-rid])

write.table(sub_table_nor,file = "~/xiaoxuan/180213/180731_bac_pert/RA/otutab_norm_RA.txt",sep = '\t',row.names = T,quote = F)


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
                     shape=c(19, 0,12, 7))

# index$Other <- factor(index$Other, levels=shapes$group)
# index$Description <- factor(index$Description, levels=colors$group)
# index$Description <- NULL
# index$Mix_Ratio <- index$Other
# reorder boxplots

# l <- c("soil", "rhizosphere", "root", "pooled_nodules")
# index$compartment <- factor(index$compartment, levels=l)
# colors <- colors[match(l, colors$group), ]


detailName <- paste("naturalPerturbation_bacteria_Community_",spike[2],sep = '_')

p <- ggplot(index, aes(x=Genotype, y=value, color=Genotype)) +
  facet_wrap(~Description)+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=detailName, y="spike-in_Abundance") +
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme


p

# ggsave(paste(figures.dir, detailName,"spike_in_Abundance.pdf", sep="_"), p)
# ggsave("../bac_New/figure/silva_spike_in_Abundance.pdf", p)

## to save in the new directory
ggsave("~/xiaoxuan/180213/180731_bac_pert/plot/filtered_spike_in_Abundance.pdf",  p)

write.table(index,file = "~/xiaoxuan/180213/180731_bac_pert/table/filtered_Facet_Natural_spikeAbun.txt",sep = '\t',row.names = T,quote = F)


#### alpha diversity comparison 

shannon <- as.data.frame(vegan::diversity(sub_table_2[-1,-ncol(sub_table_2)],index = "shannon",MARGIN = 2))
colnames(shannon)[1] <- "shannon_value"

index <- cbind(design,shannon)


##rep usage each and times
# index$Other <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)
df <- melt(index,id.vars = colnames(index)[c(4:6,8:10)],measure.vars = "shannon_value")


  
  
p <- ggplot(df, aes(x=Description, y=value)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5) +
  geom_jitter(aes(shape=Other,color=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  labs(x="Perturbation Concentration", y="shannon index") +
  geom_signif(comparisons =list(c("E05/5", "E00")),
              annotations = "No Sig.",
              y_position = 4.8,tip_length = 0, vjust=0.4,
              map_signif_level = T)+
  geom_signif(comparisons =list(c("E04", "E00")),
              annotations = "No Sig.",
              y_position = 4.73,tip_length = 0.01, vjust=0.4,
              map_signif_level = T)+
  geom_signif(comparisons =list(c("E03", "E00")),
              annotations = "No Sig.",
              y_position = 4.65,tip_length = 0.01, vjust=0.4,
              map_signif_level = T)+
  facet_wrap(~Genotype)+
  # scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  main_theme


p

# ggsave(paste(figures.dir, "alpha_shannon_boxplot.pdf", sep=""), p)
ggsave("~/xiaoxuan/180213/180731_bac_pert/plot/alpha_shannon_boxplot.pdf", p)

write.table(df,file = "~/xiaoxuan/180213/180731_bac_pert/table/alpha_shannon.txt",sep = '\t',row.names = F,quote = F)
# # p <- ggviolin(index,x="Description", y="shannon_value",fill = "Description",
# #               palette = c(c_red, c_dark_brown,c_black,c_orange),
# #               add = "boxplot", add.params = list(fill = "white"))+
# #               stat_compare_means(comparisons = list(c("E05/5", "E00"), c("E4", "E00"),c("E03", "E00")), label = "p.signif")
# #               # stat_compare_means(comparisons = my_comparisons, label = "p.signif")+  # Add significance levels
# #               # stat_compare_means(label.y = 50)                                       # Add global the p-value
# 
# 
# cl_group <- unique(sub_design_1$Description)
# mix_ratio_group <- unique(index$Other)
# colors <- data.frame(group=cl_group,
#                      color=c(c_red, c_dark_brown,c_black,c_orange))
# 
# 
# 
# shapes <- data.frame(group=mix_ratio_group,
#                      shape=c(19, 0, 24))
# index$Other <- factor(index$Other, levels=shapes$group)
# index$Description <- factor(index$Description, levels=colors$group)
# index$Mix_Ratio <- index$Other
# # reorder boxplots
# 
# 
# 
# p <- ggplot(index, aes(x=Description, y=shannon_value)) +
#   geom_violin(trim=FALSE, fill="red")+
#   geom_boxplot(width=0.1, fill="white")+
#   # geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill=rep(c(c_red, c_dark_brown,c_black,c_orange),3)) +
#   geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
#   labs(x="Perturbation Concentration", y="shannon index") +
#   geom_signif(comparisons =list(c("E05/5", "E00")),
#               annotations = "No Sig.",
#               y_position = 4.8,tip_length = 0, vjust=0.4,
#               map_signif_level = T)+
#   geom_signif(comparisons =list(c("E04", "E00")),
#               annotations = "No Sig.",
#               y_position = 4.73,tip_length = 0.01, vjust=0.4,
#               map_signif_level = T)+
#   geom_signif(comparisons =list(c("E03", "E00")),
#               annotations = "No Sig.",
#               y_position = 4.65,tip_length = 0.01, vjust=0.4,
#               map_signif_level = T)+
#   facet_wrap(~Other)+
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   main_theme
# p
# 
# # ggsave(paste(figures.dir, "alpha_shannon_violin.pdf", sep=""), p)
# ggsave("../bac_New/figure/perturbation_unoise_bac/silva_alpha_shannon_violin.pdf", p)
# write.table(index,file = "../bac_New/table/silva_Facet_Natural_bacteria_AA_alpha.txt",sep = '\t',row.names = F,quote = F)








############plot RA 

taxonomy = read.delim("../result/rep_seqs_tax.txt", row.names= 1,header=F, sep="\t",stringsAsFactors = F)
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","evalue","unknown")
taxonomy$full=taxonomy$kingdom


write.table(paste("../table/taxonomy\tSampAvsB\tPvalue",sep="\t"), file=paste("phylumpro.txt", sep=""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=F)
# # select p__Proteobacteria line
# idx=taxonomy$phylum=="p__Proteobacteria"
# # 初始化full为门，并初化因子为字符方便修改
taxonomy$full=as.character(taxonomy$phylum)
# # 修改pro门为目
# taxonomy[idx,]$full=as.character(taxonomy[idx,]$class)


sub_OTU$id <- NULL
tax_count = merge(taxonomy, sub_OTU, by="row.names")


tax_count_sum = aggregate(tax_count[,-c(1:11)], by=tax_count[11],FUN=sum) # mean
tax_count_sum$full[1] <- "other"
rownames(tax_count_sum) = tax_count_sum$full
tax_count_sum = tax_count_sum[,-1]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100

# 绘制样品组内各样品堆叠图
mean_sort = per[(order(-rowSums(per))), ] # decrease sort
colSums(mean_sort)

# Stackplot
mean_sort=as.data.frame(mean_sort)

other = colSums(mean_sort[21:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(21-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[21] = c("Low Abundance")
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

p = ggplot(data_all[data_all$site %in% "NO.1",], aes(x=variable, y = value, fill = phylumpro,group=variable)) + 
  geom_bar(stat = "identity")+ 
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ Other, scales = "free_x", switch = "x") +  
  theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Conditions in Samples ")+ylab("Relative abundance (%)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme
p

ggsave("../plot/phylum_RA_sample_NO1.pdf", p)


# ggsave("tax_stack_phylumpro_sample.png", p, width = 5, height = 3)
p = ggplot(data_all[data_all$site %in% "NO.2",], aes(x=variable, y = value, fill = phylumpro,group=variable)) + 
  # geom_bar(stat = "identity",position="fill", width=1)+ 
  geom_bar(stat = "identity")+ 
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ Other, scales = "free_x", switch = "x") +  
  theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Conditions in Samples ")+ylab("Relative abundance (%)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme
p

ggsave("../plot/phylum_RA_sample_NO2.pdf", p)




p <- ggplot(data_all[data_all$site %in% "NO.1"&data_all$phylumpro %in% c("Bacteroidetes","Proteobacteria"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p

ggsave("../plot/BacteroidetesProteobacteria_NO1_stat.pdf", p)


p <- ggplot(data_all[data_all$site %in% "NO.2"&data_all$phylumpro %in% c("Bacteroidetes","Proteobacteria"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/BacteroidetesProteobacteria_NO2_stat.pdf", p)


p <- ggplot(data_all[data_all$site %in% "NO.1"&data_all$phylumpro %in% c("Actinobacteria","Firmicutes"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/ActinobacteriaFirmicutes_NO1_stat.pdf", p)


p <- ggplot(data_all[data_all$site %in% "NO.2"&data_all$phylumpro %in% c("Actinobacteria","Firmicutes"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/ActinobacteriaFirmicutes_NO2_stat.pdf", p)


p <- ggplot(data_all[data_all$site %in% "NO.1"&data_all$phylumpro %in% c("Acidobacteria","Low Abundance"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/Acidobacteria_Abundance_NO1_stat.pdf", p)

p <- ggplot(data_all[data_all$site %in% "NO.2"&data_all$phylumpro %in% c("Acidobacteria","Low Abundance"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/Acidobacteria_Abundance_NO2_stat.pdf", p)

p <- ggplot(data_all[data_all$site %in% "NO.1"&data_all$phylumpro %in% c("Chloroflexi","Tenericutes"),], aes(x=Other, y=value,color=Other, fill=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/ChloroflexiTenericutes_NO1_stat.pdf", p)

p <- ggplot(data_all[data_all$site %in% "NO.2"&data_all$phylumpro %in% c("Chloroflexi","Tenericutes"),], aes(x=Other, y=value,color=Other, fill=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  # geom_signif(comparisons =list(c("dry", "wet")),
  #             # annotations = "No Sig.",
  #             test = "kruskal.test",
  #             y_position = 7.6,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  main_theme


p
ggsave("../plot/ChloroflexiTenericutes_NO2_stat.pdf", p)



# plot absolute abundance 
qpcr <- read.delim("~/xiaoxuan/180213/180731_bac_pert/doc/qPCR.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)

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


####absAbundance relative to spike-in
# absAbundance <- sweep(sub_table_2,2,as.numeric(sub_table_f[idx,]),'/')

####### reorder ##########

ord <- match(rownames(design),rownames(absoluteAbund))   ## column rearrange
# # ord1 <- match(rownames(bac_li), rownames(sub_table))    ## row reorder
absoluteAbund<-absoluteAbund[ord,]   ## parrellel
bac_reads_AA <- as.data.frame(t(absoluteAbund))

# remove_sample_id <- c("Bac2DMH6309","Bac2DMH6310","Bac2DMH6311","Bac2DMH6303","Bac2WMH6304","Bac2WMH6308","Bac2WWYJ11","Bac2WWYJ14","Bac2WWYJ04","Bac2WWYJ10","BacWMH6301","BacWMH6307","BacDWYJ05","BacDWYJ07","BacDWYJ02","BacDWYJ10","BacWWYJ07","BacWWYJ14")
# rid <- match(remove_sample_id,colnames(bac_reads_AA))
# bac_reads_AA <- as.data.frame(bac_reads_AA[,-rid])

write.table(bac_reads_AA,file = "~/xiaoxuan/180213/180731_bac_pert/spikein/unadjusted_absAbundance.xls",sep = '\t',row.names = T)


# bac_reads_AA <- t(absoluteAbund)

micro_load <- as.data.frame(t(colSums(bac_reads_AA)))
# merge(micro_load, design, by = "row.names")

# micro_load <- micro_load[,colnames(micro_load) %in% row.names(qpcr)]
idp <- match(colnames(micro_load),row.names(qpcr))
qpcr_ratio<- qpcr[idp,]
micro_load <- as.data.frame(t(sweep(micro_load, 2, qpcr_ratio$Copynumber*3/qpcr_ratio$spike_host_ratio, "/")))
data_all = merge(micro_load, design, by = "row.names")
colnames(data_all)[2] <- "Microbiome_Load"
data_all <- data_all %>% drop_na()
write.table(data_all,file = "~/xiaoxuan/180213/180731_bac_pert/table/modify_adjusted_bac_microbiome_load.txt",sep = '\t',row.names = T,quote = F)


# remove_sample_id <- c("Bac2DMH6309","Bac2DMH6310","Bac2DMH6311","Bac2DMH6303","Bac2WMH6304","Bac2WMH6308","Bac2WWYJ11","Bac2WWYJ14","Bac2WWYJ04","Bac2WWYJ10","BacWMH6301","BacWMH6307","BacDWYJ05","BacDWYJ07","BacDWYJ02","BacDWYJ10","BacWWYJ07","BacWWYJ14")

p = ggplot(data_all[!data_all$Description %in% c("E00")&!data_all$Genotype %in% "ZH11_Plas_Bac",], aes(x=Genotype, y = Microbiome_Load,color=Genotype, fill=Genotype)) + 
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Genotype), position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~Description , scales = "free_x", switch = "x") +  
  # theme(strip.Funkground = element_blank())+
  # 关闭x轴刻度和标签
  # theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups_Batch_Hainan")+ylab("Microbiome_Bacteria Load (Total Reads / Spike-in)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  stat_compare_means(method="wilcox.test", size=3)+
  # geom_signif(comparisons =list(c("dry", "wet")),
  #             # annotations = "No Sig.",
  #             test = "wilcox.test",
  #             y_position = 7.6,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  main_theme
p

ggsave("../plot/microbiome_load_OneVstwofold.pdf", p)






# p = ggplot(data_all[data_all$site %in% "Ah"& data_all$Genotype %in% c("MH63","WYJ7DEP1")&!data_all$Row.names %in% remove_sample_id,], aes(x=Other, y = Microbiome_Load,color=Other, fill=Other)) + 
#   geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
#   geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   # scale_y_continuous(labels = scales::percent) + 
#   # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
#   facet_grid( ~Genotype , scales = "free_x", switch = "x") +  
#   # theme(strip.background = element_blank())+
#   # 关闭x轴刻度和标签
#   # theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
#   xlab("Groups_Batch_AnHui")+ylab("Microbiome Load (Total Reads / Spike-in)")+
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
# ggsave("~/xiaoxuan/180528/180627_AQ/bac/plot_719_remove_sample/microbiome_load_Anhui.pdf", p)


# ######## soil Type comparison
# 
# micro_load <- as.data.frame(t(colSums(bac_reads_AA)))
# merge(micro_load, design, by = "row.names")
# colnames(data_all)[2] <- "Microbiome_Load" 
# 
# ###need correction




###AA taxa plot
qpcr <- read.delim("~/xiaoxuan/180213/180731_bac_pert/doc/qPCR.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
# idp <- match(colnames(bac_reads_AA),row.names(qpcr))
# ratio<- as.data.frame(qpcr)

idp <- match(colnames(bac_reads_AA),row.names(qpcr))

qpcr_ratio<- qpcr[idp,]

idp1 <- match(colnames(bac_reads_AA),row.names(qpcr_ratio))
idp1

colnames(bac_reads_AA)
rownames(qpcr_ratio)
bac_reads_AA_cor <- sweep(bac_reads_AA, 2, as.vector(qpcr_ratio$Copynumber), "/")
bac_reads_AA_cor <- sweep(bac_reads_AA_cor, 2, as.vector(qpcr_ratio$spike_host_ratio/3), "*")


# remove_sample_id <- c("Bac2DMH6309","Bac2DMH6310","Bac2DMH6311","Bac2DMH6303","Bac2WMH6304","Bac2WMH6308","Bac2WWYJ11","Bac2WWYJ14","Bac2WWYJ04","Bac2WWYJ10","BacWMH6301","BacWMH6307","BacDWYJ05","BacDWYJ07","BacDWYJ02","BacDWYJ10","BacWWYJ07","BacWWYJ14")
# rid <- match(remove_sample_id,colnames(bac_reads_AA_cor))
# bac_reads_AA_cor <- as.data.frame(bac_reads_AA_cor[,-rid])

# head(bac_reads_AA)
# head(bac_reads_AA_cor)


colSums(bac_reads_AA_cor)
write.table(bac_reads_AA_cor,file = "~/xiaoxuan/180213/180731_bac_pert/table/adjusted_absAbundance.txt",sep = '\t',row.names = T,quote = F)
# write.table(bac_reads_AA_cor,file = "~/xiaoxuan/180528/180627_AQ/bac/table/adjusted_absAbundance_1.txt",sep = '\t',row.names = T,quote = F)
write.table(bac_reads_AA_cor,file = "~/xiaoxuan/180213/180731_bac_pert/spikein/pert_otutab_norm_AA.txt",sep = '\t',row.names = T,quote = F)
# write.table(bac_reads_AA_cor,file = "~/xiaoxuan/180528/180627_AQ/bac/spike_rem_sample_720/otutab.txt",sep = '\t',row.names = T,quote = F)
# write.table(bac_reads_AA_cor,file = "~/xiaoxuan/180528/180627_AQ/bac/spikein/otutab_1.txt",sep = '\t',row.names = T,quote = F)
write.table(bac_reads_AA,file = "~/xiaoxuan/180528/180627_AQ/bac/spikein/AA.txt",sep = '\t',row.names = T,quote = F)


# bac_reads_AA_cor <- bac_reads_AA/qpcr_ratio$CopyNumber[1:197]
# ratio<- as.data.frame(qpcr[idp,1])
# bac_reads_AA_cor <- bac_reads_AA/ratio[1:197,1]

# bac_reads_AA[,1]/qpcr_ratio$CopyNumber[1]
# 
# bac_reads_AA_cor[,1]


write.table(bac_reads_AA_cor,file = "../table/adjusted_bac_reads_AA.txt",sep = '\t',row.names = T,quote = F)
tax_count = merge(taxonomy, bac_reads_AA_cor, by="row.names")



tax_count_sum = aggregate(tax_count[,-c(1:11)], by=tax_count[11],FUN=sum) # mean
tax_count_sum$full[1] <- "other"
rownames(tax_count_sum) = tax_count_sum$full
tax_count_sum = tax_count_sum[,-1]

# per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100

# 绘制样品组内各样品堆叠图
mean_sort = tax_count_sum[(order(-rowSums(tax_count_sum))), ] # decrease sort
colSums(mean_sort)

# Stackplot
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[41:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(41-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[41] = c("Low Abundance")
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


p = ggplot(data_all[data_all$site %in% "NO.1",], aes(x=variable, y = value, fill = phylumpro,group=variable)) + 
  geom_bar(stat = "identity")+ 
  # geom_bar(stat = "identity",position="fill")+ 
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ Other, scales = "free_x", switch = "x") +  
  theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups")+ylab("Absolute phylum abundance (reads per unit host DNA)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme
p

ggsave("../plot/AA_phylum_stack_phylumpro_sample_NO1.pdf", p)

p = ggplot(data_all[data_all$site %in% "NO.2",], aes(x=variable, y = value, fill = phylumpro,group=variable)) + 
  geom_bar(stat = "identity")+ 
  # geom_bar(stat = "identity",position="fill")+ 
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ Other, scales = "free_x", switch = "x") +  
  theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups")+ylab("Absolute phylum abundance (reads per unit host DNA)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme
p
ggsave("../plot/AA_phylum_stack_phylumpro_sample_NO2.pdf", p)


######AA taxa statistics

p <- ggplot(data_all[data_all$site %in% "NO.1"&data_all$phylumpro %in% c("Bacteroidetes","Proteobacteria"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p

ggsave("../plot/AA_BacteroidetesProteobacteria_NO1_stat.pdf", p)


p <- ggplot(data_all[data_all$site %in% "NO.2"&data_all$phylumpro %in% c("Bacteroidetes","Proteobacteria"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/AA_BacteroidetesProteobacteria_NO2_stat.pdf", p)


p <- ggplot(data_all[data_all$site %in% "NO.1"&data_all$phylumpro %in% c("Actinobacteria","Firmicutes"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/AA_ActinobacteriaFirmicutes_NO1_stat.pdf", p)


p <- ggplot(data_all[data_all$site %in% "NO.2"&data_all$phylumpro %in% c("Actinobacteria","Firmicutes"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/AA_ActinobacteriaFirmicutes_NO2_stat.pdf", p)


p <- ggplot(data_all[data_all$site %in% "NO.1"&data_all$phylumpro %in% c("Acidobacteria","Low Abundance"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/AA_Acidobacteria_Abundance_NO1_stat.pdf", p)

p <- ggplot(data_all[data_all$site %in% "NO.2"&data_all$phylumpro %in% c("Acidobacteria","Low Abundance"),], aes(x=Other, y=value,color=Other,full=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/AA_AcidobacteriaLow_Abundance_NO2_stat.pdf", p)

p <- ggplot(data_all[data_all$site %in% "NO.1"&data_all$phylumpro %in% c("Chloroflexi","Tenericutes"),], aes(x=Other, y=value,color=Other, fill=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  main_theme


p
ggsave("../plot/AA_ChloroflexiTenericutes_NO1_stat.pdf", p)

p <- ggplot(data_all[data_all$site %in% "NO.2"&data_all$phylumpro %in% c("Chloroflexi","Tenericutes"),], aes(x=Other, y=value,color=Other, fill=Other)) +
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
  stat_compare_means(method="kruskal.test", size=3)+
  # geom_signif(comparisons =list(c("dry", "wet")),
  #             # annotations = "No Sig.",
  #             test = "kruskal.test",
  #             y_position = 7.6,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  main_theme


p
ggsave("../plot/AA_ChloroflexiTenericutes_NO2_stat.pdf", p)






# p
# 
# 
# # ggsave(paste(figures.dir,"Facet_Natural_bacteria_RA_mean_box_line_chart.pdf", sep=""), p)
# ggsave("../bac_New/figure/perturbation_unoise_bac/silva_Facet_Natural_bacteria_RA_mean_box_line_chart.pdf", p)
# write.table(mean_table,file = "../bac_New/table/silva_Facet_Natural_bacteria_RA_mean_box_line.txt",sep = '\t',row.names = F,quote = F)
# 
# #### alpha diversity comparison 
# 
# shannon <- as.data.frame(diversity(sub_table_2[-1,-ncol(sub_table_2)],index = "shannon",MARGIN = 2))
# colnames(shannon)[1] <- "shannon_value"
# 
# index <- cbind(design,shannon)
# 


# #               y_position = 4.8,tip_length = 0, vjust=0.4,
# #               map_signif_level = T)+
# #   geom_signif(comparisons =list(c("E04", "E00")),
# #               annotations = "No Sig.",
# #               y_position = 4.73,tip_length = 0.01, vjust=0.4,
# #               map_signif_level = T)+
# #   geom_signif(comparisons =list(c("E03", "E00")),
# #               annotations = "No Sig.",
# #               y_position = 4.65,tip_length = 0.01, vjust=0.4,
# #               map_signif_level = T)+
# #   facet_wrap(~Other)+
# #   scale_colour_manual(values=as.character(colors$color)) +
# #   scale_shape_manual(values=shapes$shape) +
# #   main_theme
# # p
# # 
# # # ggsave(paste(figures.dir, "alpha_shannon_violin.pdf", sep=""), p)
# # ggsave("../bac_New/figure/perturbation_unoise_bac/silva_alpha_shannon_violin.pdf", p)
# # write.table(index,file = "../bac_New/table/silva_Facet_Natural_bacteria_AA_alpha.txt",sep = '\t',row.names = F,quote = F)
# 
# 
# 
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

sub_table_2$id <- row.names(sub_table_2) 
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
# loadMeta("../bac_New/table/otu_table.biom")
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
# CSS("../bac_New/table/otu_table.biom","../bac_New/table/otu_table_css.biom" , 1)
# 
# "CSS" <-function(input_path, out_path, output_CSS_statistics=NULL) {
#   obj = load_biom("../bac_New/table/otu_table.biom")
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

bray_curtis <- read.table("~/xiaoxuan/180528/180627_AQ/bac/spikein/beta/bray_curtis_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../bac_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../bac_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 
design = design[rownames(design) %in% rownames(bray_curtis),]
id <- match(row.names(design[design$site %in% "Ah",]),rownames(bray_curtis))



ah_bray_curtis <- bray_curtis[id,id]

##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~1
k <- 3
pcoa <- cmdscale(ah_bray_curtis, k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y","z")
# points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)

points <- cbind(points, design[match(rownames(points), design$SampleID), ])



colors <- data.frame(group=geno_group,
                     color=c(c_red, c_dark_brown,c_black,c_orange))

shapes <- data.frame(group=condition_group,
                     shape=c(19, 0))





# points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # points$genotype <- factor(points$genotype, levels=shapes$group)
# points$GroupID <- factor(points$GroupID, levels=colors$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Other)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("~/xiaoxuan/180528/180627_AQ/bac/spikein/beta/AA_bray_curtis_pcoa12_Ah_wetDry.pdf", p)


p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Other)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("~/xiaoxuan/180528/180627_AQ/bac/spikein/beta/AA_bray_curtis_pcoa13_Ah_wetDry.pdf", p)




#####Hainan pcoa 
bray_curtis <- read.table("~/xiaoxuan/180528/180627_AQ/bac/spikein/beta/bray_curtis_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../bac_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../bac_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 
design = design[rownames(design) %in% rownames(bray_curtis),]
id <- match(row.names(design[design$site %in% "Hn",]),rownames(bray_curtis))



hn_bray_curtis <- bray_curtis[id,id]

##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~1
k <- 3
pcoa <- cmdscale(hn_bray_curtis, k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y","z")
# points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)

points <- cbind(points, design[match(rownames(points), design$SampleID), ])



colors <- data.frame(group=geno_group,
                     color=c(c_red, c_dark_brown,c_black,c_orange))

shapes <- data.frame(group=condition_group,
                     shape=c(19, 0))





# points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # points$genotype <- factor(points$genotype, levels=shapes$group)
# points$GroupID <- factor(points$GroupID, levels=colors$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Other)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("~/xiaoxuan/180528/180627_AQ/bac/spikein/beta/AA_bray_curtis_pcoa12_Hn_wetDry.pdf", p)


p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Other)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("~/xiaoxuan/180528/180627_AQ/bac/spikein/beta/AA_bray_curtis_pcoa13_Hn_wetDry.pdf", p)

# p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Other)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   facet_wrap(~ site)+
#   main_theme +
#   theme(legend.position="top")
# p
# 
# 
# ggsave("../plot/PCoA12_BC_WD.pdf", p)



# p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=site)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   main_theme +
#   theme(legend.position="top")
# p
# 
# ggsave("~/xiaoxuan/180528/180627_AQ/bac/spikein/beta/AA_bray_curtis_pcoa13.pdf", p)


#####RA AH
bray_curtis <- read.table("~/xiaoxuan/180213/180731_bac_pert/RA/beta_div/bray_curtis_otutab_norm.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../bac_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../bac_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 
design = read.table("~/xiaoxuan/180213/180731_bac_pert/doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
sample_list <- as.matrix(row.names(design))

design = design[rownames(design) %in% rownames(bray_curtis),]
# id <- match(row.names(design[design$site %in% "Ah"&!design$Genotype %in% "MH63ZH",]),rownames(bray_curtis))
# design = design[rownames(design) %in% rownames(bray_curtis),]
id <- match(row.names(design),rownames(bray_curtis))


# ah_bray_curtis <- bray_curtis
ah_bray_curtis <- bray_curtis[id,id]


##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~1
k <- 3
pcoa <- cmdscale(ah_bray_curtis, k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y","z")
# points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)

points <- cbind(points, design[match(rownames(points), design$SampleID), ])



colors <- data.frame(group=geno_group,
                     color=c(c_red, c_dark_brown,c_black))

shapes <- data.frame(group=condition_group,
                     shape=c(19, 6))





# points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # points$genotype <- factor(points$genotype, levels=shapes$group)
# points$GroupID <- factor(points$GroupID, levels=colors$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Genotype)) +
  geom_point(alpha=.7, size=5) +
  facet_grid(~Description)+
  scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("~/xiaoxuan/180213/180731_bac_pert/RA/beta_div/802_RA_bray_curtis_pcoa12.pdf", p)

# p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Other)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   facet_wrap(~ site)+
#   main_theme +
#   theme(legend.position="top")
# p
# 
# 
# ggsave("../plot/PCoA12_BC_WD.pdf", p)



p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Genotype)) +
  geom_point(alpha=.7, size=5) +
  facet_grid(~Description)+
  scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("~/xiaoxuan/180213/180731_bac_pert/RA/beta_div/802_RA_bray_curtis_pcoa13.pdf", p)

p <- ggplot(points, aes(x=y, y=z, color=Genotype,shape=Genotype)) +
  geom_point(alpha=.7, size=5) +
  facet_grid(~Description)+
  scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(Genotypevalues=shapes$shape) +
  labs(x=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("~/xiaoxuan/180213/180731_bac_pert/RA/beta_div/802_RA_bray_curtis_pcoa23.pdf", p)
#### Hainan

bray_curtis <- read.table("~/xiaoxuan/180528/180627_AQ/bac/RA_rem_sample_720/beta/bray_curtis_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
#
# bray_curtis <- read.table("~/xiaoxuan/180528/180627_AQ/bac/RA/beta/bray_curtis_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis <- read.table("~/xiaoxuan/180528/180627_AQ/bac/RA/beta/bray_curtis_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../bac_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../bac_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 
design = design[rownames(design) %in% rownames(bray_curtis),]
id <- match(row.names(design[design$site %in% "Hn"&!design$Genotype %in% "MH63ZH",]),rownames(bray_curtis))



hn_bray_curtis <- bray_curtis[id,id]


##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~1
k <- 3
pcoa <- cmdscale(hn_bray_curtis, k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y","z")
# points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)

points <- cbind(points, design[match(rownames(points), design$SampleID), ])



colors <- data.frame(group=geno_group,
                     color=c(c_red, c_dark_brown,c_black,c_orange))

shapes <- data.frame(group=condition_group,
                     shape=c(19, 0))





# points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # points$genotype <- factor(points$genotype, levels=shapes$group)
# points$GroupID <- factor(points$GroupID, levels=colors$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Other)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("~/xiaoxuan/180528/180627_AQ/bac/RA_rem_sample_720/beta/721_RA_bray_curtis_pcoa12_Hn_Wetdry.pdf", p)


p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Other)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("~/xiaoxuan/180528/180627_AQ/bac/RA_rem_sample_720/beta/721_RA_bray_curtis_pcoa13_Hn_Wetdry.pdf", p)

p <- ggplot(points, aes(x=y, y=z, color=Genotype,shape=Other)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("~/xiaoxuan/180528/180627_AQ/bac/RA_rem_sample_720/beta/721_RA_bray_curtis_pcoa23_Hn_Wetdry.pdf", p)

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
# ggsave("../bac_New/figure/eulcian_PCoA_BC.pdf", p)
# write.table(points,file = "../table/WD_Natural_bacteria_betaDiv_points.txt",sep = '\t',row.names = F,quote = F)

# bray_curtis <- read.table("../table/beta/bray_curtis_otu_table_for_beta.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# # bray_curtis_1 <- read.table("../bac_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# # eulian <- read.table("../bac_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# # 
# 
# ##  blank \t  need sed to filter  
# 
# #### Unconstrained ordination
# # pcoa<-capscale(dis~1
# k <- 3
# pcoa <- cmdscale(bray_curtis, k=k, eig=T)
# # pcoa <- cmdscale(eulian, k=k, eig=T)
# points <- pcoa$points
# eig <- pcoa$eig
# points <- as.data.frame(points)
# colnames(points) <- c("x", "y","z")
# # points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)
# 
# points <- cbind(points, design[match(rownames(points), design$SampleID), ])
# 
# 
# 
# colors <- data.frame(group=geno_group,
#                      color=c(c_red, c_dark_brown,c_black,c_orange))
# 
# shapes <- data.frame(group=condition_group,
#                      shape=c(19, 0))
# 
# 
# 
# 
# 
# # points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # # points$genotype <- factor(points$genotype, levels=shapes$group)
# # points$GroupID <- factor(points$GroupID, levels=colors$group)
# 
# # plot PCo 1 and 2
# 
# p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=site)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   main_theme +
#   theme(legend.position="top")
# p
# 
# ggsave("../plot/PCoA12_BC_Batch.pdf", p)
# 
# p <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Other)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   facet_wrap(~ site)+
#   main_theme +
#   theme(legend.position="top")
# p
# 
# 
# ggsave("../plot/PCoA12_BC_WD.pdf", p)
# 
# 
# 
# p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=site)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   main_theme +
#   theme(legend.position="top")
# p
# 
# ggsave("../plot/PCoA13_BC_Batch.pdf", p)
# 
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
# 
# ####### Cluster dendrogram ##############
# bray_curtis <- read.table("../table/beta/bray_curtis_otu_table_for_beta.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# # bray_curtis_1 <- read.table("../bac_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# # eu
# 
# 
# # dat_rare <- read.delim("PyroTagger_counts_rarefied_A.txt", header=T, row.names=1, sep="\t")
# # # application of threshold for lower presence limit of relative abundance
# # threshold <- 5    # 0.5% presence: use integer 5 considering x1000 normalization
# # dat_rare_thres <- dat_rare[ apply(dat_rare,1,max) > threshold, ]
# # # log2-transformation
# # dat_rare_thres_log <- log2(dat_rare_thres + 1)
# # dim(dat_rare_thres_log)
# # calculate Bray-Curtis dissimilarity 
# otu_table_for_beta$OTU_id <- NULL
# 
# # log2-transformation
# # dat_rare_thres_log <- log2(dat_rare_thres + 1)
# dat_dist <- vegdist(t(otu_table_for_beta), method="bray")
# # clustering
# dat_dist_clu <- hclust(dat_dist , "average")
# # pdf("Figure_4a.pdf", width=15, height=5)
# plot(dat_dist_clu, main = "Clutering Dentrogram of  Beta-diversity / Bray Curtis", lwd=1.5, cex=.5)
# 
# 
# 
# 
# 
# weighed_unifrac <- read.table("../bac_New/table/beta_div_unifrac/weighted_unifrac_table.from_txt_json.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# # bray_curtis_1 <- read.table("../bac_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# # eulian <- read.table("../bac_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# # 
# 
# ##  blank \t  need sed to filter  
# 
# #### Unconstrained ordination
# # pcoa<-capscale(dis~1
# k <- 2
# pcoa <- cmdscale(weighed_unifrac , k=k, eig=T)
# # pcoa <- cmdscale(eulian, k=k, eig=T)
# points <- pcoa$points
# eig <- pcoa$eig
# points <- as.data.frame(points)
# colnames(points) <- c("x", "y")
# points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)
# 
# points <- cbind(points, design_1[match(rownames(points), design_1$SampleID), ])
# 
# 
# 
# colors <- data.frame(group=unique(points$Description),
#                      color=c(c_red, c_very_dark_green, c_black ,c_dark_brown))
# 
# 
# shapes <- data.frame(group=unique(points$Concentration),
#                      shape=c(19, 10, 5))
# 
# 
# 
# 
# 
# # points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # # points$genotype <- factor(points$genotype, levels=shapes$group)
# # points$GroupID <- factor(points$GroupID, levels=colors$group)
# 
# # plot PCo 1 and 2
# 
# p <- ggplot(points, aes(x=x, y=y, color=Description,shape=Concentration)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("Weighted_Unifrac_PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("Weighted_Unifrac_PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   main_theme +
#   theme(legend.position="top")
# p
# 
# ggsave("../bac_New/figure/PCoA_W_unifrac.pdf", p)
# 
# #####unweighted Unifrac
# 
# un_weighed_unifrac <- read.table("../bac_New/table/beta_div_unifrac/unweighted_unifrac_table.from_txt_json.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# # bray_curtis_1 <- read.table("../bac_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# # eulian <- read.table("../bac_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# # 
# 
# ##  blank \t  need sed to filter  
# 
# #### Unconstrained ordination
# # pcoa<-capscale(dis~1
# k <- 2
# pcoa <- cmdscale(un_weighed_unifrac, k=k, eig=T)
# # pcoa <- cmdscale(eulian, k=k, eig=T)
# points <- pcoa$points
# eig <- pcoa$eig
# points <- as.data.frame(points)
# colnames(points) <- c("x", "y")
# points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)
# 
# points <- cbind(points, design_1[match(rownames(points), design_1$SampleID), ])
# 
# 
# 
# colors <- data.frame(group=unique(points$Description),
#                      color=c(c_red, c_very_dark_green, c_black ,c_dark_brown))
# 
# 
# shapes <- data.frame(group=unique(points$Concentration),
#                      shape=c(19, 10, 5))
# 
# 
# 
# 
# 
# # points$soiltype <- factor(points$soiltype, levels=shapes$group)
# # # points$genotype <- factor(points$genotype, levels=shapes$group)
# # points$GroupID <- factor(points$GroupID, levels=colors$group)
# 
# # plot PCo 1 and 2
# 
# p <- ggplot(points, aes(x=x, y=y, color=Description,shape=Concentration)) +
#   geom_point(alpha=.7, size=5) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("Weighted_Unifrac_PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("Weighted_Unifrac_PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
#   main_theme +
#   theme(legend.position="top")
# p
# 
# ggsave("../bac_New/figure/PCoA_W_unifrac.pdf", p)
# 
# 
# 
# 
# 
# 
# 
# 
# ############################### RA accuracy ####################
# # ####A bac_perturbation B onefold  C twofold
# # 
# # 
# # mean_BA <- mean_table[!c(1:nrow(mean_table))%%3==0,]
# # mean_C <- mean_table[c(1:nrow(mean_table))%%3==0,]
# # mean_A <- mean_BA[!c(1:nrow(mean_BA))%%2==0,]
# # mean_B <- mean_BA[c(1:nrow(mean_BA))%%2==0,]
# # 
# # ratio_BA <- mean_B$value/mean_A$value
# # ratio_BC <- mean_B$value/mean_C$value
# # 
# # mean_A$ration <- ratio_BA
# # mean_C$ration <- ratio_BC
# # 
# # 
# # detailName <- "Feature OTUs in Natural Community"
# # 
# # ######BA########
# # p <- ggplot(mean_A[!mean_A$Group.3 %in% "OTU_9",], aes(x=Group.3, y=ration,shape=Group.2,color=Group.2,group=Group.2)) +
# #   facet_wrap(~Group.2)+
# #   geom_line() +
# #   geom_point(size=3, fill="white") +
# #   scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
# #   scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
# #   labs(x=detailName, y="Relative Abundance Ration GroupB(one fold host) vs GroupA(bacteria spike perturbation)") +
# #   geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
# #   theme(axis.text.x = element_text(size=10,angle = 90))+
# #   # geom_hline(mean(mean_A$ration),colour='black',lwd=0.36,linetype="dashed")+
# #   main_theme
# # 
# # p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # # p <- p+geom_hline(yintercept = mean(mean_A[!mean_A$Group.3 %in% "OTU_9",]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # 
# # 
# # p
# # 
# # # ggsave(paste(figures.dir,"RA_ratio_B_oneFold_vs_A_BacSpikePertur_line_chart.pdf", sep=""), p)
# # ggsave("../bac_New/figure/perturbation_unoise_bac/silva_RA_ratio_B_oneFold_vs_A_BacSpikePertur_line_chart.pdf", p)
# # write.table(mean_A,file = "../bac_New/table/silva_RA_ratio_B_oneFold_vs_A_BacSpikePertur_line_chart_accuracy.txt",sep = '\t',row.names = F,quote = F)
# # 
# # ####BC###
# # p <- ggplot(mean_C[!mean_C$Group.3 %in% "OTU_9",], aes(x=Group.3, y=ration,shape=Group.2,color=Group.2,group=Group.2)) +
# #   facet_wrap(~Group.2)+
# #   geom_line() +
# #   geom_point(size=3, fill="white") +
# #   scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
# #   scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
# #   labs(x=detailName, y="Relative Abundance Ratio GroupB(one fold host) vs GroupC (two fold)") +
# #   geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
# #   theme(axis.text.x = element_text(size=10,angle = 90))+
# #   main_theme
# # 
# # p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # # p <- p+geom_hline(yintercept = mean(ratio_AC[-c(16)]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # 
# # p
# # 
# # # ggsave(paste(figures.dir,"RA_ratio_B_oneFold_vsC_twoFold_line_chart.pdf", sep=""), p)
# # ggsave("../bac_New/figure/perturbation_unoise_bac/silva_RA_ratio_B_oneFold_vsC_twoFold_line_chart.pdf", p)
# # write.table(mean_C,file = "../bac_New/table/silva_RA_ratio_B_oneFold_vsC_twoFold_line_chart_accuracy.txt",sep = '\t',row.names = F,quote = F)
# # 
# # 
# # # plot absolute abundance 
# # 
# # 
# # # mapping spike 12-4
# # # idx1 <- match(spike[1],row.names(sub_OTU))
# # sub_table_2 <- sub_OTU
# # 
# # # idx1 <- match(spike[3],row.names(sub_table_2))
# # # sub_table_2 <- sub_table_2
# # 
# # len <- length(sub_table_2[1,])
# # idx <- match(spike[2],row.names(sub_table_2))
# # # col_sums <- colSums(sub_table_2[-idx,-len])
# # 
# # # idx <- match(spike[3],row.names(sub_table))
# # # len <- length(sub_table[1,])
# # internal_ref <- as.numeric(sub_table_2[idx,-len])
# # internal_ref=1+internal_ref
# # absoluteAbund <-as.data.frame(t(sub_table_2[-idx,-len])/internal_ref)
# # 
# # bac_reads_AA <- t(absoluteAbund)
# # 
# # micro_load <- as.data.frame(colSums(bac_reads_AA))
# # write.table(bac_reads_AA,file = "../table/bac_reads_AA.txt",sep = '\t',row.names = F,quote = F)
# # 
# # 
# # 
# # 
# # # index <- cbind(design[match(row.names(absoluteAbund), design$SampleID), ],absoluteAbund )
# # # 
# # # index<- index[,c(1:25)]
# # # index$Other <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)
# # # 
# # # # design <- subset(design,design$Description %in% "E06")
# # # # idx5 <- match(rownames(design),colnames(absoluteAbund)) 
# # # # absoluteAbund <- absoluteAbund[,idx5]
# # # colors <- data.frame(group=cl_group,
# # #                      color=c(c_red, c_dark_brown,c_blue,c_green))
# # # 
# # # 
# # # 
# # # shapes <- data.frame(group=mix_ratio_group,
# # #                      shape=c(19, 0, 24))
# # # 
# # # index$Other <- factor(index$Other, levels=shapes$group)
# # # index$Mix_Ratio <- index$Other
# # # index$Description <- factor(index$Description,levels =colors$group )
# # # 
# # # index
# # # 
# # # df <- melt(index)
# # # 
# # # mean_table <- aggregate(df$value, by = list(df$Mix_Ratio,df$variable,df$Description),FUN = mean)
# # # mean_table$variable <- rownames(mean_table)
# # # colnames(mean_table)[4] <- "value"
# # # detailName <- "Facet Natural Feature bacteria "
# # # p <- ggplot(mean_table, aes(x=Group.2,y=value,color=Group.2)) +
# # #   facet_wrap(Group.1~Group.3,scales = "free")+
# # #   geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
# # #   geom_jitter(aes(shape=Group.3), position=position_jitter(0.17), size=1, alpha=0.7) +
# # #   # scale_colour_manual(values=as.character(colors$color)) +
# # #   # scale_shape_manual(values=shapes$shape) +
# # #   labs(x=detailName, y="AA_mean") +
# # #   theme(axis.text.x = element_text(size=10,angle = 90))+
# # #   main_theme
# # # 
# # # p
# # # 
# # # 
# # # # ggsave(paste(figures.dir,"Facet_Natural_bacteria_AA_mean_box_line_chart.pdf", sep=""), p)
# # # 
# # # ggsave("../bac_New/figure/perturbation_unoise_bac/silva_Facet_Natural_bacteria_AA_mean_box_line_chart.pdf", p)
# # # write.table(mean_table,file = "../bac_New/table/silva_Facet_Natural_bacteria_AA_mean_box_line.txt",sep = '\t',row.names = F,quote = F)
# # # 
# # # ####AA Ratio
# # # 
# # # 
# # # mean_BA <- mean_table[!c(1:nrow(mean_table))%%3==0,]
# # # mean_C <- mean_table[c(1:nrow(mean_table))%%3==0,] # disturbation
# # # mean_A <- mean_BA[!c(1:nrow(mean_BA))%%2==0,]  # one fold
# # # mean_B <- mean_BA[c(1:nrow(mean_BA))%%2==0,]  # two fold
# # # 
# # # ratio_AB <- mean_A$value/mean_B$value
# # # ratio_AC <- mean_A$value/mean_C$value
# # # 
# # # mean_A$ration <- ratio_AB
# # # mean_C$ration <- ratio_AC
# # # 
# # # 
# # # detailName <- "Feature OTUs in Natural Community"
# # # 
# # # ######AB########
# # # p <- ggplot(mean_A, aes(x=Group.2, y=ration,shape=Group.2,color=Group.3,group=Group.3)) +
# # #   facet_wrap(~Group.3)+
# # #   geom_line() +
# # #   geom_point(size=3, fill="white") +
# # #   scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
# # #   scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
# # #   labs(x=detailName, y="Absolute Abundance Ration GroupA(one fold host) vs GroupB(two fold host)") +
# # #   geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
# # #   theme(axis.text.x = element_text(size=10,angle = 90))+
# # #   main_theme
# # # 
# # # p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # # # p <- p+geom_hline(yintercept = mean(mean_A$ration[-c(16)]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # # 
# # # 
# # # p
# # # 
# # # # ggsave(paste(figures.dir,"modified_AA_ratio_B_twoFold_vs_A_onefold_line_chart.pdf", sep=""), p)
# # # ggsave("../bac_New/figure/perturbation_unoise_bac/silva_modified_AA_ratio_A_onefold_vsB_twoFold__line_chart.pdf", p)
# # # 
# # # write.table(mean_A,file = "../bac_New/table/silva_modified_Facet_Natural_bacteria_AA_chart_line_AB.txt",sep = '\t',row.names = F,quote = F)
# # # 
# # # 
# # # 
# # # ####AC###
# # # p <- ggplot(mean_C, aes(x=Group.2, y=ration,shape=Group.2,color=Group.3,group=Group.3)) +
# # #   facet_wrap(~Group.3)+
# # #   geom_line() +
# # #   geom_point(size=3, fill="white") +
# # #   scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
# # #   scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
# # #   labs(x=detailName, y="Absolute Abundance Ratio GroupA(one fold host) vs GroupC (perturbation)") +
# # #   geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
# # #   theme(axis.text.x = element_text(size=10,angle = 90))+
# # #   main_theme
# # # 
# # # p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # # # p <- p+geom_hline(yintercept = mean(ratio_AC[-c(16)]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # # 
# # # p
# # # 
# # # # ggsave(paste(figures.dir,"AA_ratio_A_oneFold_vsC_perturbation_line_chart.pdf", sep=""), p)
# # # ggsave("../bac_New/figure/perturbation_unoise_bac/silva_AA_ratio_A_oneFold_vsC_perturbation_line_chart.pdf", p)
# # # write.table(mean_C,file = "../bac_New/table/silva_modified_Facet_Natural_bacteria_AA_mean_chart_line_AC.txt",sep = '\t',row.names = F,quote = F)
# # # 
