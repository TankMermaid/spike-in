# scripts for 16S spike-in analysis for mock-Bac-E05-10 

rm(list=ls())
tem <- "E05-10_mock_bac"
setwd("/mnt/bai/xiaoxuan/20180402spike-in/bac_New")
print(paste("Your working directory is in",getwd()))


########### to import the plotting theme() function 
source("script/plot_function.R")

### output directory assigned to include the pics & tables########################
figures.dir <- paste("/mnt/bai/xiaoxuan/20180402spike-in/bac_New/figure/",tem,'/',sep = '')
table.dir <- paste("/mnt/bai/xiaoxuan/20180402spike-in/bac_New/table/",tem,'/',sep = '')


fig_flag <- dir.exists(figures.dir)
if( isTRUE(!fig_flag)){
  dir.create(figures.dir)
}

tab_flag <- dir.exists(table.dir)
if( isTRUE(!tab_flag)){
  dir.create(table.dir)
}


##### spike-in design in this batch #####################
spike <- c("BI-OS-11-3","BI-OS-12-4","BI-OS-10-2")

#1 design Mapping file 

design = read.table("doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)  #在design 表最后一列加了一列SampleID的信息，该信息就是row names

#调取mock-EO5/10的design 信息， 筛选条件为2个，分别是Plasmid ID中的Scal-BI-12-4和Description中的E05/10
sub_design <- design[design$PlasmidID %in% "Scal-BI-12-4" & design$Description %in% "E05/10",]

###2 bacteria_list to filter 
bac_li = read.delim("doc/bacterial_list.txt",row.names = 1, header=T, sep="") 

rownames(bac_li)

###3 OTUs file 
otu_table = read.delim("data/usearch_map_L3/observation_table.txt", row.names= 1,  header=T, sep="\t")

otu_table$id <- row.names(otu_table) #在otu_table 表最后一列加了一列id的信息，该信息就是row names

sub_table <- subset(otu_table,  otu_table$id %in%  rownames(bac_li))#调取otu_table中的otu，筛选条件为otu_table 中id和bac_li的rownames相同的otu

####### reorder ##########

ord <- match(rownames(sub_design),colnames(sub_table))   ## column rearrange
ord1 <- match(rownames(bac_li), rownames(sub_table))    ## row reorder
sub_table<- sub_table[ord1,ord]   ## parrellel

## ----Pseudocount to  escape the error of division of 0---------------------------------------------------------
sub_table <- sub_table+1

# spike BI-12-4
idx <- match(spike[2],row.names(sub_table))
sub_table_2 <- sub_table[-idx,]  #从sub_table里过滤掉BI-12-4的信息

####relAbundance without spike-in
relAbundance <- sweep(sub_table_2,2,colSums(sub_table_2),'/')

####relAbundance with spike-in
relAbun_withSpike <- sweep(sub_table,2,colSums(sub_table),'/')

write.table(relAbundance,"/mnt/bai/xiaoxuan/20180402spike-in/bac_New/table/E05-10_mock_bac/relAbundance.txt",sep='\t')
write.table(relAbun_withSpike,"/mnt/bai/xiaoxuan/20180402spike-in/bac_New/table/E05-10_mock_bac/relAbun_withSpike.txt",sep='\t')
write.table(sub_table,"/mnt/bai/xiaoxuan/20180402spike-in/bac_New/table/E05-10_mock_bac/sub_table.txt",sep='\t')
write.table(sub_table_2,"/mnt/bai/xiaoxuan/20180402spike-in/bac_New/table/E05-10_mock_bac/sub_table_2.txt",sep='\t')
