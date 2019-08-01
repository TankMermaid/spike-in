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


detailName <- paste("Perturbation_bacteria_Community_",spike[2],sep = '_')


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
# ggsave("~/xiaoxuan/180213/180731_bac_pert/plot/alpha_shannon_boxplot.pdf", p)
ggsave("~/xiaoxuan/180213/180731_bac_pert/plot/perturbationBacteriaShannonBoxplot.pdf", p)
