rm(list=ls())
options(warn=3)


#########The taxa are ordered by the significance of the correlation between their QMP abundance 

########### setting the working directory and print it ###################
tem <- "wetDryRice"
setwd("~/xiaoxuan/180528/180627_AQ/bac/script/")
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


######################## AA phylum Plotting construction


map= design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
map$SampleID <- rownames(map)
# sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil"&design$Genotype %in% "Disease",]
# id_rm= c("TAD1","TAD5","TAD8")
# sub_design= sub_design[!rownames(sub_design) %in% id_rm,]

bac_aa <- read.delim("../spike_rem_sample_720/otutab.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
# fungi_aa <- read.delim("~/xiaoxuan/180724/TARDF16/spikein/otutab_norm_AA.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)

# 
# id <- colnames(fungi_aa)%in%colnames(bac_aa)
# fungi_aa <- fungi_aa[,id]

# id <- match(colnames(bac_aa),colnames(fungi_aa))
# bac_aa <- bac_aa[,id]
bac_aa <- bac_aa[,colnames(bac_aa)%in% rownames(design)]
# bac_aa <- bac_aa[,colnames(bac_aa)%in% rownames(sub_design)]
# fungi_aa <- fungi_aa[,colnames(fungi_aa)%in% rownames(sub_design)]


# bac_aa$species <- rep("bac",nrow(bac_aa))
# fungi_aa$species <- rep("fungi",nrow(fungi_aa))
merge_bac_fungi <-bac_aa

# ###1. Filtering OTUs by occurrence frequency (i.e.,number of samples an OTU is Present 60% of the samples)
# Abu <- merge_bac_fungi
# Abu$species <- NULL
# table <- Abu
# table[table>0.0001] <- 1
# table.generalist <- Abu[which(rowSums(table)>=5),]  ### 5 out of 7
# Abu <- table.generalist
# 
# merge_bac_fungi <- merge_bac_fungi[rownames(merge_bac_fungi) %in% rownames(Abu),]


tax_merge_bac_fungi <- bac_tax <- read.delim("../result/rep_seq_tax_bac.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
# fungi_tax <- read.delim("~/xiaoxuan/180724/TARDF16/result/rep_seqs_tax_filtered.txt", row.names= 1,header=F, sep="\t",stringsAsFactors = F)
# fungi_tax$V10 <- NULL
# fungi_tax$V9 <- NULL

# colnames(fungi_tax) <- colnames(bac_tax)

# tax_merge_bac_fungi <-rbind(fungi_tax,bac_tax) 

tax_merge_bac_fungi$ID <- rownames(tax_merge_bac_fungi)
merge_bac_fungi$ID <- rownames(merge_bac_fungi)
merge_abun_tax<- merge(merge_bac_fungi,tax_merge_bac_fungi,by = "ID")
write.table(merge_bac_fungi,file = "../table/aa_abundance_dryWet.txt",sep = '\t',row.names = T,quote = F)


######dplyr to data manipulation
# install.packages("psych")
library(psych)

phy.aa<- aa_dis_a<- merge_abun_tax[,c(2:181,183)] %>%
  mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:181)])) %>%
  group_by(Phylum) %>%
  summarize_all(sum) %>%
  arrange(desc(rowSumAbun))%>%
  gather(key = "SampleID", value = "AbsAb",-Phylum,-rowSumAbun)


# Get the top taxa in the whole dataset
phy.ord <- phy.aa

phy.top <- c(as.character(head(phy.ord, n = 25)$Phylum))

# Generate the df that will be used for plotting
prof.plot <- phy.aa %>%
  filter(Phylum %in% phy.ord$Phylum[1:12]) %>%
  inner_join(map, by = "SampleID")%>%
  arrange(desc(AbsAb))


order_de<- prof.plot[,c(3:4)]%>%
  group_by(SampleID) %>%
  summarize_all(sum) %>%
  arrange(desc(AbsAb))


prof.plot$SampleID <- factor(prof.plot$SampleID,levels =unique(order_de$SampleID) )

# ordr_de <- as.data.frame(order_de)

# ids <- match(rownames(order_de),)
# prof.plot <- prof.plot[ids,]
# my.labs <- list("Acidobacteria",
#                 "Actinobacteria",
#                 "Bacteroidetes",
#                 "Chloroflexi",
#                 "Firmicutes",
#                 "Gemmatimonadetes",
#                 "Planctomycetes",
#                 "Verrucomicrobia",
#                 expression(paste(alpha, " - Proteobacteria")),
#                 expression(paste(beta, " - Proteobacteria")),
#                 expression(paste(delta, " - Proteobacteria")),
#                 expression(paste(gamma, " - Proteobacteria")))


write.table(prof.plot,file = "../table/aa_abundance_dryWet.txt",sep = '\t',row.names = T,quote = F)

# Plot
phyloplot <-ggplot() +
  # ggplot() +
  geom_bar(data = prof.plot, aes(x=SampleID, y=AbsAb, fill = Phylum), stat = "identity", position = "stack") +
  facet_grid( Other~ site, scales = "free", space = "free") +
  # scale_fill_manual(name = "Phylum",
  #                   values = c(brewer.pal(8, "Set2")[1:8], brewer.pal(11, "RdYlBu")[7:10])) +
  # geom_point(data = filter(prof.plot, Genotype!= "Bulksoil"),
  #            aes(paste(Genotype, Other),-0.08, color = Genotype),
  #            size = 5, shape = 15) +
  # scale_color_manual(values = c(brewer.pal(12, "Paired")[7:10]),
  #                    guide = FALSE) +
  labs(y = "Absolute Abundance") +
  theme_light()+
  theme(text = element_text(size = 17),
        legend.position = "right",
        legend.text = element_text(size =8),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_legend(ncol = 1))
  
  ggsave(phyloplot,filename = "../plot_719_remove_sample/aa_abundance.pdf")
  
######################## RA phylum Plotting construction
  
  
  rm(list=ls())
  options(warn=3)
  
  
  #########The taxa are ordered by the significance of the correlation between their QMP abundance 
  
  ########### setting the working directory and print it ###################
  tem <- "wetDryRice"
  setwd("~/xiaoxuan/180528/180627_AQ/bac/script/")
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
  
  
  map= design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
  map$SampleID <- rownames(map)
  # sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil"&design$Genotype %in% "Disease",]
  # id_rm= c("TAD1","TAD5","TAD8")
  # sub_design= sub_design[!rownames(sub_design) %in% id_rm,]
  
  bac_aa <- read.delim("../RA_rem_sample_720/otutab_norm_RA.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
  # fungi_aa <- read.delim("~/xiaoxuan/180724/TARDF16/spikein/otutab_norm_AA.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
  
  # 
  # id <- colnames(fungi_aa)%in%colnames(bac_aa)
  # fungi_aa <- fungi_aa[,id]
  
  # id <- match(colnames(bac_aa),colnames(fungi_aa))
  # bac_aa <- bac_aa[,id]
  bac_aa <- bac_aa[,colnames(bac_aa)%in% rownames(design)]
  # bac_aa <- bac_aa[,colnames(bac_aa)%in% rownames(sub_design)]
  # fungi_aa <- fungi_aa[,colnames(fungi_aa)%in% rownames(sub_design)]
  
  
  # bac_aa$species <- rep("bac",nrow(bac_aa))
  # fungi_aa$species <- rep("fungi",nrow(fungi_aa))
  merge_bac_fungi <-bac_aa
  
  # ###1. Filtering OTUs by occurrence frequency (i.e.,number of samples an OTU is Present 60% of the samples)
  # Abu <- merge_bac_fungi
  # Abu$species <- NULL
  # table <- Abu
  # table[table>0.0001] <- 1
  # table.generalist <- Abu[which(rowSums(table)>=5),]  ### 5 out of 7
  # Abu <- table.generalist
  # 
  # merge_bac_fungi <- merge_bac_fungi[rownames(merge_bac_fungi) %in% rownames(Abu),]
  
  
  tax_merge_bac_fungi <- bac_tax <- read.delim("../result/rep_seq_tax_bac.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
  # fungi_tax <- read.delim("~/xiaoxuan/180724/TARDF16/result/rep_seqs_tax_filtered.txt", row.names= 1,header=F, sep="\t",stringsAsFactors = F)
  # fungi_tax$V10 <- NULL
  # fungi_tax$V9 <- NULL
  
  # colnames(fungi_tax) <- colnames(bac_tax)
  
  # tax_merge_bac_fungi <-rbind(fungi_tax,bac_tax) 
  
  tax_merge_bac_fungi$ID <- rownames(tax_merge_bac_fungi)
  merge_bac_fungi$ID <- rownames(merge_bac_fungi)
  merge_abun_tax<- merge(merge_bac_fungi,tax_merge_bac_fungi,by = "ID")
  write.table(merge_bac_fungi,file = "../table/ra_abundance_dryWet.txt",sep = '\t',row.names = T,quote = F)
  
  
  ######dplyr to data manipulation
  # install.packages("psych")
  library(psych)
  
  phy.aa<- aa_dis_a<- merge_abun_tax[,c(2:181,183)] %>%
    mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:181)])) %>%
    group_by(Phylum) %>%
    summarize_all(sum) %>%
    arrange(desc(rowSumAbun))%>%
    gather(key = "SampleID", value = "RelAb",-Phylum,-rowSumAbun)
  
  
  # Get the top taxa in the whole dataset
  phy.ord <- phy.aa
  
  phy.top <- c(as.character(head(phy.ord, n = 25)$Phylum))
  
  # Generate the df that will be used for plotting
  prof.plot <- phy.aa %>%
    filter(Phylum %in% phy.ord$Phylum[1:12]) %>%
    inner_join(map, by = "SampleID")%>%
    arrange(desc(RelAb))
  
  
  order_de<- prof.plot[,c(3:4)]%>%
    group_by(SampleID) %>%
    summarize_all(sum) %>%
    arrange(desc(RelAb))
  
  
  prof.plot$SampleID <- factor(prof.plot$SampleID,levels =unique(order_de$SampleID) )
  
  # ordr_de <- as.data.frame(order_de)
  
  # ids <- match(rownames(order_de),)
  # prof.plot <- prof.plot[ids,]
  # my.labs <- list("Acidobacteria",
  #                 "Actinobacteria",
  #                 "Bacteroidetes",
  #                 "Chloroflexi",
  #                 "Firmicutes",
  #                 "Gemmatimonadetes",
  #                 "Planctomycetes",
  #                 "Verrucomicrobia",
  #                 expression(paste(alpha, " - Proteobacteria")),
  #                 expression(paste(beta, " - Proteobacteria")),
  #                 expression(paste(delta, " - Proteobacteria")),
  #                 expression(paste(gamma, " - Proteobacteria")))
  
  # Plot
  # phyloplot <-ggplot() +
  phyloplot <-ggplot() +
    geom_bar(data = prof.plot, aes(x=SampleID, y=RelAb, fill = Phylum), stat = "identity", position = "stack") +
    facet_grid( Other~ site, scales = "free", space = "free") +
    # scale_fill_manual(name = "Phylum",
    #                   values = c(brewer.pal(8, "Set2")[1:8], brewer.pal(11, "RdYlBu")[7:10])) +
    # geom_point(data = filter(prof.plot, Genotype!= "Bulksoil"),
    #            aes(paste(Genotype, Other),-0.08, color = Genotype),
    #            size = 5, shape = 15) +
    # scale_color_manual(values = c(brewer.pal(12, "Paired")[7:10]),
    #                    guide = FALSE) +
    labs(y = "Relative Abundance") +
    theme_light()+
    theme(text = element_text(size = 17),
          legend.position = "right",
          legend.text = element_text(size =8),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    guides(fill = guide_legend(ncol = 1))
  
  ggsave(phyloplot,filename = "../plot_719_remove_sample/ra_abundance.pdf")
  write.table(prof.plot,file = "../table/ra_abundance_dryWet.txt",sep = '\t',row.names = T,quote = F)