# scripts for 16S spike-in analysis for E05-10 
#
# originally by Tank

#xnzhang@genetics.ac.cn

rm(list=ls())
########### setting the working directory and print it ###################

tem <- "E05-10_mock_bac"
setwd("/mnt/bai/xiaoning/xiaoxuan/180213/bac_New/")
print(paste("Your working directory is in",getwd()))


########### to import the plotting theme() function ########################
source("script/plot_function.R")

### output directory assigned to include the pics & tables########################
figures.dir <- paste("/mnt/bai/xiaoning/xiaoxuan/180213/bac_New/figure/",tem,'/',sep = '')
table.dir <- paste("/mnt/bai/xiaoning/xiaoxuan/180213/bac_New/table/",tem,'/',sep = '')

fig_flag <- dir.exists(figures.dir)
if( isTRUE(!fig_flag)){
  dir.create(figures.dir)
}

tab_flag <- dir.exists(table.dir)
if( isTRUE(!tab_flag)){
  dir.create(table.dir)
}
####################################


##### spike-in design in this batch #####################
 spike <- c("BI-OS-11-3","BI-OS-12-4","BI-OS-10-2")



#1 design Mapping file 
design = read.table("doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)

sub_design <- design[design$PlasmidID %in% "Scal-BI-12-4" & design$Description %in% "E05/10",]
# sample_list <- as.matrix(row.names(design))


#2 bacteria_list to filter 
bac_li = read.delim("doc/bacterial_list.txt",row.names = 1, header=T, sep="") 
rownames(bac_li)

#3 OTUs file 
otu_table = read.delim("data/usearch_map_L3/observation_table.txt", row.names= 1,  header=T, sep="\t")
# otu_table = read.delim("usearch/observation_table.txt", row.names= 1,  header=T, sep="\t")
otu_table$id <- row.names(otu_table)


sub_table <- subset(otu_table,  otu_table$id %in%  rownames(bac_li))

####### reorder ##########

ord <- match(rownames(sub_design),colnames(sub_table))   ## column rearrange
ord1 <- match(rownames(bac_li), rownames(sub_table))    ## row reorder
sub_table<- sub_table[ord1,ord]   ## parrellel

## ----Pseudocount to  escape the error of division of 0---------------------------------------------------------
sub_table <- sub_table+1

#-----------------------------------------------------------------------------------------------------------------------------------------------
# spike BI-12-4
idx <- match(spike[2],row.names(sub_table))
sub_table_2 <- sub_table[-idx,]
####relAbundance without spike-in
relAbundance <- sweep(sub_table_2,2,colSums(sub_table_2),'/')

####relAbundance with spike-in
relAbun_withSpike <- sweep(sub_table,2,colSums(sub_table),'/')

###### scaling factor ###########
## extract the spike-in 
Pos_12_4 <- grep("BI-OS-12-4",rownames(sub_table))

Spike_12_4 <- as.matrix(sub_table[Pos_12_4,])

#Calculate size factor:
sizeFactor_Spike<-Spike_12_4/mean(Spike_12_4)

#Apply size factor to each sample
sub.counts<- sub_table

for(i in seq(ncol(sub.counts))){
  sub.counts[,i]<-sub.counts[,i]/sizeFactor_Spike[i]
}




# idx <- match(spike[2],row.names(sub_table_2))
col_sums <- colSums(sub_table_2[-idx,-len])


sub_table_nor_t <- t(sub_table_2[-idx,-len])/col_sums
sub_table_nor <- t(sub_table_nor_t)





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




gra_group <- unique(design$PlasmidID)               # plasmid id 
# gra_group <- unique(design$Description)
mix_ratio_group <- unique(design$Other)
cl_group <- unique(design$Plasmid)

colors <- data.frame(group=cl_group,
                     color=c(c_red, c_dark_brown))

shapes <- data.frame(group=mix_ratio_group,
                     shape=c(19, 0, 24))

index$Other <- factor(index$Other, levels=shapes$group)
index$Plasmid <- factor(index$Plasmid, levels=colors$group)

index$PlasType <- index$Plasmid
index$Mix_Ratio <- index$Other
# reorder boxplots

# l <- c("soil", "rhizosphere", "root", "pooled_nodules")
# index$compartment <- factor(index$compartment, levels=l)
# colors <- colors[match(l, colors$group), ]


detailName <- paste("Synthetic_Bacteria_Community_ratio",spike[2],sep = '_')

p <- ggplot(index, aes(x=Mix_Ratio, y=value, color=PlasType)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Mix_Ratio), position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=detailName, y="E05/10-spike-in_Abundance") +
  main_theme

p

ggsave(paste(figures.dir, detailName,"spike_in_Abundance.pdf", sep="_"), p)


#plot delta



idx2 <- match("delta_value",colnames(t_table))
index <- cbind(total_table[, idx2], design[match(t_table[,length(t_table[1,])], design$SampleID), ])
colnames(index)[1] <- "value"


# E05 E04 bad result
# index <- index[index$Description %in% "E06",]

# 

colors <- data.frame(group=cl_group,
                     color=c(c_red, c_dark_brown))

shapes <- data.frame(group=mix_ratio_group,
                     shape=c(19, 0, 24))

index$Other <- factor(index$Other, levels=shapes$group)
index$Plasmid <- factor(index$Plasmid, levels=colors$group)

index$PlasType <- index$Plasmid
index$Mix_Ratio <- index$Other
# index <- subset(index, index$Plasmid %in% c("Circular"))


p <- ggplot(index, aes(x=Mix_Ratio, y=value, color=PlasType)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Mix_Ratio), position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=detailName, y="delta_value") +
  main_theme

p


ggsave(paste(figures.dir,detailName, "delta_value.pdf", sep="_"), p)









#plot RA 


idx1 <- match("spike-in_Abundance",colnames(t_table))
idx2 <- match("delta_value",colnames(t_table))
index <- cbind(total_table[, -c(idx2,idx1)], design[match(t_table[,length(t_table[1,])], design$SampleID), ])

# index <- index[index$Description %in% "E06",]



colors <- data.frame(group=cl_group,
                     color=c(c_red, c_dark_brown))



shapes <- data.frame(group=mix_ratio_group,
                     shape=c(19, 0, 24))

index$Other <- factor(index$Other, levels=shapes$group)
index$Plasmid <- factor(index$Plasmid, levels=colors$group)

index$PlasType <- index$Plasmid
index$Mix_Ratio <- index$Other

index
df <- melt(index)
# reorder boxplots

# l <- c("soil", "rhizosphere", "root", "pooled_nodules")
# index$compartment <- factor(index$compartment, levels=l)
# colors <- colors[match(l, colors$group), ]


# index_line <- as.data.frame(t(index[c(1:10)]))
# sub_design <- subset(design[c(55:57),],design[c(55:57),]$Description %in% "E05")
# # sub_design <- subset(design[c(55:57,82:84),],design[c(55:57,82:84),]$Description %in% "E05")
# idx4 <- match(rownames(sub_design),colnames(index_line))
# index_line <- index_line[,idx4]
# 
# index_line$name <-factor(row.names(index_line))
# str(index_line)

# df <- melt(index_line,id.vars = "name")

# for (i in 1:length(spike)){
#   for (j in 1:length(cl_group)){
#     for (m in 1:length(mix_ratio_group)){
#       df1<-df[df$PlasmidID %in% spike[i] & df$Plasmid %in% cl_group[j] & df$Other %in% mix_ratio_group[m],]
#       detailName <- paste("Synthetic_Bacteria_RA",spike[i],cl_group[j], mix_ratio_group[m],sep ='_')
#       
#       p <- ggplot(df1, aes(x=variable, y=value,color=SampleID,group=SampleID,shape=variable)) +
#         geom_line() +
#         geom_point(size=3, fill="white") +
#         scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black))) +
#         scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11)) +
#         labs(x=detailName, y="Relative Abundance within Synthetic Bacteria") +
#         main_theme
#       
#       # p <- ggplot(index, aes(x=, y=index[,1], group=Mix_Ratio,shape=Mix_Ratio)) +
#       #   geom_line( size=0.7) +
#       #   geom_point(size=3, fill="white") +
#       #   scale_colour_manual(values=as.character(colors$color)) +
#       #   scale_shape_manual(values=shapes$shape) +
#       #   labs(x="gradient_of_Bacteria_mix_ratio", y="Relative Abundance within Synthetic Bacteria") +
#       #   main_theme
#       
#       p
#       
#       ggsave(paste(figures.dir,detailName,"_line_chart.pdf", sep=""), p)
#     }}}


#facet 
######### Facet according MIx ratio and PlasmidType
groupID <- unique(df$PlasmidID)

df1<-df[df$PlasmidID %in% groupID[1] , ]
# df1<-df[df$PlasmidID %in% groupID[1] & df$Plasmid %in% cl_group[2], ]
detailName <- paste("Facet_Synthetic_Bacteria_RA",groupID[1],sep ='_')
write.table(df[,c(10:14)],file = paste(figures.dir,detailName,"_aaValue.xls"),sep = '\t',row.names = F)
# p <- ggplot(df1, aes(x=variable, y=value,color=SampleID,shape=variable)) +
err.mean <- tapply(df1$value,df1[,12:13],mean)
err.me<- err.mean[1:3,] 
data.frame(err.me)

p <- ggplot(df1, aes(x=variable, y=value,color=SampleID,group=SampleID,shape=variable)) +
  facet_wrap(PlasType~Mix_Ratio,scales = "free")+
  # geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1) +
  geom_line() +
  geom_point(size=3, fill="white") +
  # scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black))) +
  scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3)) +
  labs(x=detailName, y="Relative Abundance within Synthetic Bacteria") +
  main_theme +
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14))+
  theme(axis.text.x = element_text(size=10,angle = 90))

# p <- ggplot(index, aes(x=, y=index[,1], group=Mix_Ratio,shape=Mix_Ratio)) +
#   geom_line( size=0.7) +
#   geom_point(size=3, fill="white") +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x="gradient_of_Bacteria_mix_ratio", y="Relative Abundance within Synthetic Bacteria") +
#   main_theme

# p <- p+geom_line(aes(x=,y=allerr$err.mean),colour='black')
p
ggsave(paste(figures.dir,detailName,"_line_chart.pdf", sep=""), p)










# plot absolute abundance 
idx <- match(spike[2],row.names(sub_table))
sub_table_1 <- sub_table[,sub_table[idx,]>160]

# mapping spike 12-4
idx1 <- match(spike[1],row.names(sub_table_1))
sub_table_2 <- sub_table_1[-idx1,]


# idx1 <- match(spike[3],row.names(sub_table_2))
# sub_table_2 <- sub_table_2[-idx1,]

len <- length(sub_table_2[1,])
idx <- match(spike[2],row.names(sub_table_2))
# col_sums <- colSums(sub_table_2[-idx,-len])

# idx <- match(spike[3],row.names(sub_table))
# len <- length(sub_table[1,])
internal_ref <- as.numeric(sub_table_2[idx,-len])

absoluteAbund <-as.data.frame(t(sub_table_2[-idx,-len])/internal_ref)
index <- cbind(absoluteAbund, design[match(row.names(absoluteAbund), design$SampleID), ])

# sub_design <- subset(design,design$Description %in% "E06")
# idx5 <- match(rownames(sub_design),colnames(absoluteAbund)) 
# absoluteAbund <- absoluteAbund[,idx5]
colors <- data.frame(group=cl_group,
                     color=c(c_red, c_dark_brown))



shapes <- data.frame(group=mix_ratio_group,
                     shape=c(19, 0, 24))

index$Other <- factor(index$Other, levels=shapes$group)
index$Plasmid <- factor(index$Plasmid, levels=colors$group)

index$PlasType <- index$Plasmid
index$Mix_Ratio <- index$Other

index
df <- melt(index)

# for (i in 1:length(spike)){
#   for (j in 1:length(cl_group)){
#     for (m in 1:length(mix_ratio_group)){
#       df1<-df[df$PlasmidID %in% spike[i] & df$Plasmid %in% cl_group[j] & df$Other %in% mix_ratio_group[m],]
#       detailName <- paste("Synthetic_Bacteria_AA",spike[i],cl_group[j], mix_ratio_group[m],sep ='_')
#       
#       p <- ggplot(df1, aes(x=variable, y=value,color=SampleID,group=SampleID,shape=variable)) +
#         geom_line() +
#         geom_point(size=3, fill="white") +
#         scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black))) +
#         scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11)) +
#         labs(x=detailName, y="Absolute Abundance within Synthetic Bacteria") +
#         main_theme
#       
#       # p <- ggplot(index, aes(x=, y=index[,1], group=Mix_Ratio,shape=Mix_Ratio)) +
#       #   geom_line( size=0.7) +
#       #   geom_point(size=3, fill="white") +
#       #   scale_colour_manual(values=as.character(colors$color)) +
#       #   scale_shape_manual(values=shapes$shape) +
#       #   labs(x="gradient_of_Bacteria_mix_ratio", y="Relative Abundance within Synthetic Bacteria") +
#       #   main_theme
#       
#       p
#       
#       ggsave(paste(figures.dir,detailName,"_line_chart.pdf", sep=""), p)
#     }}}
#**************************************************************************************************************************************
#facet 
######### Facet according MIx ratio and PlasmidType
groupID <- unique(df$PlasmidID)

df1<-df[df$PlasmidID %in% groupID[1] , ]
# df1<-df[df$PlasmidID %in% groupID[1] & df$Plasmid %in% cl_group[2], ]
detailName <- paste("Facet_Synthetic_Bacteria_AA",groupID[1],sep ='_')
write.table(df[,c(10:14)],file = paste(figures.dir,detailName,"_aaValue.xls"),sep = '\t',row.names = F)
# p <- ggplot(df1, aes(x=variable, y=value,color=SampleID,shape=variable)) +
p <- ggplot(df1, aes(x=variable, y=value,color=SampleID,group=SampleID,shape=variable)) +
  facet_wrap(PlasType~Mix_Ratio,scales = "free")+
  geom_line() +
  geom_point(size=3, fill="white") +
  # scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black))) +
  scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3)) +
  labs(x=detailName, y="Absolute Abundance within Synthetic Bacteria") +
  main_theme +
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14))+
  theme(axis.text.x = element_text(size=10,angle = 90))

# p <- ggplot(index, aes(x=, y=index[,1], group=Mix_Ratio,shape=Mix_Ratio)) +
#   geom_line( size=0.7) +
#   geom_point(size=3, fill="white") +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x="gradient_of_Bacteria_mix_ratio", y="Relative Abundance within Synthetic Bacteria") +
#   main_theme

p

ggsave(paste(figures.dir,detailName,"_line_chart.pdf", sep=""), p)


#***********************************************************************************************************************************************

# absoluteAbund$name <-factor(row.names(absoluteAbund))
# str(absoluteAbund)

# df <- melt(absoluteAbund,id.vars = "name")

# p <- ggplot(df, aes(x=name, y=value,color=variable,group=variable)) +
#   geom_line() +
#   geom_point(size=3, fill="white") +
#   # scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   scale_y_continuous(breaks=c(0.5,0.6,0.7,0.8,0.9,1,2,3,4))+
#   theme_bw()+
#   labs(x="gradient_of_Fungal_mix_ratio", y="Absolute Abundance within Synthetic Bacteria") +
#   main_theme
# 
# p
# 
# ggsave(paste(figures.dir, "AA_line_chart.pdf", sep=""), p)


#plot RA with internal ref 

# spike 12-4
idx <- match(spike[2],row.names(sub_table))
sub_table_1 <- sub_table[,sub_table[idx,]>160]

# mapping spike 12-4
idx1 <- match(spike[1],row.names(sub_table_1))
sub_table_2 <- sub_table_1[-idx1,]

# idx1 <- match(spike[2],row.names(sub_table_2))
# sub_table_2 <- sub_table_2[-idx1,]

len <- length(sub_table_2[1,])
col_sums <- colSums(sub_table_2[,-len])


sub_table_nor <- t(sub_table_2[,-len])/col_sums
# sub_table_nor_t <- t(sub_table_2[,-len])/col_sums
# sub_table_nor <- t(sub_table_nor_t)

index <- cbind(sub_table_nor, design[match(row.names(sub_table_nor), design$SampleID), ])

# sub_design <- subset(design,design$Description %in% "E06")
# idx5 <- match(rownames(sub_design),colnames(absoluteAbund)) 
# absoluteAbund <- absoluteAbund[,idx5]
colors <- data.frame(group=cl_group,
                     color=c(c_red, c_dark_brown))



shapes <- data.frame(group=mix_ratio_group,
                     shape=c(19, 0, 24))

index$Other <- factor(index$Other, levels=shapes$group)
index$Plasmid <- factor(index$Plasmid, levels=colors$group)

index$PlasType <- index$Plasmid
index$Mix_Ratio <- index$Other

index
df <- melt(index)
df

# for (i in 1:length(spike)){
#   for (j in 1:length(cl_group)){
#     for (m in 1:length(mix_ratio_group)){
#       df1<-df[df$PlasmidID %in% spike[i] & df$Plasmid %in% cl_group[j] & df$Other %in% mix_ratio_group[m],]
#       detailName <- paste("Synthetic_Bacteria_RA_withSpike",spike[i],cl_group[j], mix_ratio_group[m],sep ='_')
#       
#       p <- ggplot(df1, aes(x=variable, y=value,color=SampleID,group=SampleID,shape=variable)) +
#         geom_line() +
#         geom_point(size=3, fill="white") +
#         scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black))) +
#         scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3)) +
#         labs(x=detailName, y="Absolute Abundance within Synthetic Bacteria") +
#         main_theme
#       
#       # p <- ggplot(index, aes(x=, y=index[,1], group=Mix_Ratio,shape=Mix_Ratio)) +
#       #   geom_line( size=0.7) +
#       #   geom_point(size=3, fill="white") +
#       #   scale_colour_manual(values=as.character(colors$color)) +
#       #   scale_shape_manual(values=shapes$shape) +
#       #   labs(x="gradient_of_Bacteria_mix_ratio", y="Relative Abundance within Synthetic Bacteria") +
#       #   main_theme
#       
#       p
#       
#       ggsave(paste(figures.dir,detailName,"_line_chart.pdf", sep=""), p)
#     }}}

#**************************************************************************************************************************************
#facet 
######### Facet according MIx ratio and PlasmidType
groupID <- unique(df$PlasmidID)

df1<-df[df$PlasmidID %in% groupID[1] , ]
# df1<-df[df$PlasmidID %in% groupID[1] & df$Plasmid %in% cl_group[2], ]
detailName <- paste("Facet_Synthetic_Bacteria_RA_withSpike",groupID[1],sep ='_')
write.table(df[,c(10:14)],file = paste(figures.dir,detailName,"_aaValue.xls"),sep = '\t',row.names = F)
# p <- ggplot(df1, aes(x=variable, y=value,color=SampleID,shape=variable)) +
p <- ggplot(df1, aes(x=variable, y=value,color=SampleID,group=SampleID,shape=variable)) +
  facet_wrap(PlasType~Mix_Ratio,scales = "free")+
  geom_line() +
  geom_point(size=3, fill="white") +
  # scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black))) +
  scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3)) +
  labs(x=detailName, y="Relative Abundance within Synthetic Bacteria + Spike-in") +
  main_theme +
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14))+
  theme(axis.text.x = element_text(size=10,angle = 90))

# p <- ggplot(index, aes(x=, y=index[,1], group=Mix_Ratio,shape=Mix_Ratio)) +
#   geom_line( size=0.7) +
#   geom_point(size=3, fill="white") +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x="gradient_of_Bacteria_mix_ratio", y="Relative Abundance within Synthetic Bacteria") +
#   main_theme

p

ggsave(paste(figures.dir,detailName,"_line_chart.pdf", sep=""), p)

