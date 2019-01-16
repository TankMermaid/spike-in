rm(list=ls())
options(warn=3)

########### setting the working directory and print it ###################
tem <- "wetDryRice"
setwd("~/xiaoxuan/180528/180706_fungal/script/")
print(paste("Your working directory is in",getwd()))


########### to import the plotting theme() function ########################
source("plot_function.R")
# install.packages("ggpubr")
library(ggpubr)
library(vegan)
library(ggplot2)
library(reshape)
library(multcomp)
library(ggsignif)
library("Biobase", quietly=T, warn.conflicts=F)
library("edgeR", quietly=T, warn.conflicts=F)

library("ggplot2", quietly=T, warn.conflicts=F)
library("gplots", quietly=T, warn.conflicts=F)
library("grid", quietly=T, warn.conflicts=F)
library("RColorBrewer", quietly=T, warn.conflicts=F)
library("reshape2", quietly=T, warn.conflicts=F)
library("VennDiagram", quietly=T, warn.conflicts=F)

design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 

bac_aa <- read.delim("~/xiaoxuan/180528/180627_AQ/bac/table/adjusted_bac_reads_AA_for_merge.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_aa <- read.delim("../table/adjusted_fungi_reads_AA_for_merge.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)

id <- match(colnames(bac_aa),colnames(fungi_aa))
bac_aa <- bac_aa[,id]

bac_aa$species <- rep("bac",nrow(bac_aa))
fungi_aa$species <- rep("fungi",nrow(fungi_aa))
merge_bac_fungi <-rbind(bac_aa,fungi_aa) 

write.table(merge_bac_fungi,file = "../table/adjusted_bacAndFungiMerge_abundance.txt",sep = '\t',row.names = T,quote = F)

# bacVsFungiRatio <- sweep(tax_count_sum,,fungi,"/")



# bac_sum<- as.data.frame(colSums(bac_aa))
# fun_sum <- as.data.frame(colSums(fungi_aa))
# merge_sum <- colSums(merge_bac_fungi)


micro_load <- as.data.frame(colSums(merge_bac_fungi[,-ncol(merge_bac_fungi)]))
# idp <- match(row.names(qpcr),rownames(micro_load))
# ratio<- as.data.frame(qpcr[idp,1])
# micro_load <- micro_load/ratio[1:197,1]
data_all = merge(micro_load, design, by = "row.names")
colnames(data_all)[2] <- "Microbiome_Load"
write.table(data_all,file = "../table/adjusted_bacAndFungiMerge_microbiome_load.txt",sep = '\t',row.names = T,quote = F)

p = ggplot(data_all[data_all$Site %in% "HaiNan"& data_all$Genotype %in% c("MH63","WYJ7DEP.1"),], aes(x=Condition, y = Microbiome_Load,color=Condition, fill=Condition)) + 
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Condition), position=position_jitter(0.17), size=1, alpha=0.7) +
  # scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~Genotype , scales = "free_x", switch = "x") +  
  # theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  # theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups_Batch_HaiNan")+ylab("Microbiome(Bac&Fungi) Load (Total Reads / Spike-in)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  stat_compare_means(method="wilcox.test", size=3)+
  # geom_signif(comparisons =list(c("dry", "wet")),
  #             # annotations = "No Sig.",
  #             test = "wilcox.test",
  #             y_position = 7.6,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  main_theme
p

ggsave("../plot/mergeBacFungi_microbiome_load_Hn.pdf", p)

p = ggplot(data_all[data_all$Site %in% "AnHui"& data_all$Genotype %in% c("MH63","WYJ7DEP.1") &!data_all$Row.names %in% c("Fun2WWYJ04","Fun2WWYJ07","Fun2WWYJ14"),], aes(x=Condition, y = Microbiome_Load,color=Condition, fill=Condition)) + 
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Condition), position=position_jitter(0.17), size=1, alpha=0.7) +
  # scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~Genotype , scales = "free_x", switch = "x") +  
  # theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  # theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups_Batch_AnHui")+ylab("Microbiome(Bac&Fungi) Load (Total Reads / Spike-in)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  stat_compare_means(method="wilcox.test", size=3)+
  # geom_signif(comparisons =list(c("dry", "wet")),
  #             # annotations = "No Sig.",
  #             test = "wilcox.test",
  #             y_position = 7.6,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  main_theme
p

ggsave("../plot/mergeBacFungi_microbiome_load_AnHui.pdf", p)

##### Bac vs Fungi Ratio 
tax_count_sum = aggregate(merge_bac_fungi[,-ncol(merge_bac_fungi)], by=merge_bac_fungi[ncol(merge_bac_fungi)],FUN=sum)
rownames(tax_count_sum) <- tax_count_sum$species

df <- melt(tax_count_sum)

data_all = merge(df,design,by.x="variable", by.y="row.names")
data_all$species <- factor(data_all$species,levels = c("bac","fungi"))


dat <- data_all[data_all$Site %in% "HaiNan" & data_all$Genotype %in% c("MH63","WYJ7DEP.1"),]
dat$GroupID <- NULL
mean_table <- aggregate(dat$value,by = list(dat$species,dat$Condition,dat$Genotype),FUN = mean)
colnames(mean_table) <- c("species","Condition","Genotype","value")
st.err <- function(x) {
  sd(x)/sqrt(length(x))
}

std_table <- aggregate(dat$value,by = list(dat$species,dat$Condition,dat$Genotype),FUN =st.err)

colnames(std_table) <- c("species","Condition","Genotype","stderr")
mean_table <- table_merge <- cbind(std_table,mean_table)
# mean_table$variable <- rownames(mean_table)


p = ggplot(mean_table, aes(x=Condition, y = value, fill = species,group=Condition)) + 
  geom_bar(stat = "identity",width = 0.2,position = position_dodge(0.4))+
  geom_errorbar(aes(ymin=value-stderr, ymax=value+stderr),
                width=.1,                    # Width of the error bars
                position=position_dodge(.1))+
  # geom_bar(stat = "identity",position="fill")+
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ Genotype, scales = "free_x", switch = "x") +
  theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups")+ylab("Absolute kingdom abundance (reads per unit host DNA)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  theme_bw(base_size = 12) +
  # scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values = c("darkred", "purple"))+
  main_theme
p
ggsave("../plot/bacVsFunIndividual_HaiNan.pdf", p)


p = ggplot(data_all[data_all$Site %in% "HaiNan"& data_all$Genotype %in% c("MH63","WYJ7DEP.1"),], aes(x=variable, y = value, fill = species,group=variable)) + 
  geom_bar(stat = "identity")+ 
  # geom_bar(stat = "identity",position="fill")+ 
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ Condition, scales = "free_x", switch = "x") +
  theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups")+ylab("Absolute kingdom abundance (reads per unit host DNA)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme
p
ggsave("../plot/mergeBacFungi_microbiome_load_individual_HaiNan.pdf", p)

p = ggplot(data_all[data_all$Site %in% "AnHui"& data_all$Genotype %in% c("MH63","WYJ7DEP.1"),], aes(x=variable, y = value, fill = species,group=variable)) + 
  geom_bar(stat = "identity")+ 
  # geom_bar(stat = "identity",position="fill")+ 
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ Condition, scales = "free_x", switch = "x") +
  theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups")+ylab("Absolute kingdom abundance (reads per unit host DNA)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme
p
ggsave("../plot/mergeBacFungi_microbiome_load_individual_AnHui.pdf", p)

write.table(data_all,file = "../table/bacFunInd.txt",sep = '\t',row.names = T,quote = F)

tax_count_sum$species <- NULL
bacVsFungiRatio <- tax_count_sum[1,]/tax_count_sum[2,]
write.table(bacVsFungiRatio,file = "../table/bacVsFungiRatio.txt",sep = '\t',row.names = T,quote = F)


t.tax_count <- t(tax_count_sum)

data_all = merge( t.tax_count,design, by="row.names")
colnames(data_all)[1] <- "SampleID"
write.table(data_all,file = "../table/bacAndFungiInd_microbiome_load.txt",sep = '\t',row.names = T,quote = F)


p = ggplot(data_all[data_all$Bactch %in% "NO1",], aes(x=Other, y = Microbiome_Load,color=Other, fill=Other)) + 
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
  # scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  # scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~Genotype , scales = "free_x", switch = "x") +  
  # theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  # theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups_Batch_NO.1")+ylab("Microbiome(Bac&Fungi) Load (Total Reads / Spike-in)")+
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  stat_compare_means(method="wilcox.test", size=3)+
  # geom_signif(comparisons =list(c("dry", "wet")),
  #             # annotations = "No Sig.",
  #             test = "wilcox.test",
  #             y_position = 7.6,tip_length = 0, vjust=0.4,
  #             map_signif_level = T)+
  main_theme
p

ggsave("../plot/mergeBacFungi_microbiome_load_NO1.pdf", p)




library(BiocInstaller)
source("http://www.bioconductor.org/biocLite.R")
useDevel()
biocLite("microbiome")

library(microbiome)  
