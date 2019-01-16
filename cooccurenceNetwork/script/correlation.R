rm(list = ls())
setwd("~/xiaoxuan/180213/bac_New/")


design <- read.table("doc/design.txt",header=T, row.names= 1, sep="\t")
sub_design <- subset(design, design$Other %in% "2:2:2")
sub_design <- subset(sub_design, sub_design$Description  %in% c("E05","E05/5","E05/10","E05/20") & sub_design$PlasmidID %in% "Scal-BI-12-4")

concent <- read.table("doc/concentr.txt",header=T, row.names= 1, sep="\t")

l1 <- read.table("data/usearch_map_L1/observation_table.txt",header=T, row.names= 1,sep='\t')
l2 <- read.table("data/usearch_map_L2/observation_table.txt",header=T, row.names= 1,sep='\t')
l3 <- read.table("data/usearch_map_L3/observation_table.txt",header=T, row.names= 1,sep='\t')
l5 <- read.table("data/usearch_map_L5/observation_table.txt",header=T, row.names= 1,sep='\t')



