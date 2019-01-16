rm(list=ls())

########### setting the working directory and print it ###################
tem <- "Nature_disturbation_OTU_E05-5"
setwd("~/xiaoxuan/180213/fung/")
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

### output directory assigned to include the pics & tables########################
figures.dir <- paste("~/xiaoxuan/180213/fung//plot/Unoise/",tem,'/',sep = '')
table.dir <- paste("~/xiaoxuan/180213/fung/table/",tem,'/',sep = '')


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
design = read.table("doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
sample_list <- as.matrix(row.names(design))

ids <- design$Other %in% c("2:2:2","1:1:1","2:2:1")

sub_design <- design[!ids,]
# sub_design <- design[ids,]


#2 bacteria_list to filter 
# fungl_li = read.delim("doc/bacterial_list.txt",row.names = 1, header=T, sep="") 
# rownames(fungl_li)

#3 OTUs file

otu_table = read.table("result/unoise/otu_table_tax.txt", row.names= 1,  header=1, sep="\t")


id <- grep("OTU_4$",row.names(otu_table))
otu_table <- otu_table[-id,]
id <- grep("OTU_39$",row.names(otu_table))

otu_table[id,]
# mock_otu = read.table("result/unoise/candidate_otu_12.txt",  header=F, sep="\t")
# otu_table = read.delim("usearch/observation_table.txt", row.names= 1,  header=T, sep="\t")
# otu_table$id <- row.names(otu_table)
# otu_table = read.delim("usearch_map_L4/observation_table.txt", row.names= 1,  header=T, sep="\t")
# otu_table = read.delim("usearch_map_L4/observation_table.txt", row.names= 1,  header=T, sep="\t")
# otu_table = read.delim("usearch_map_L4/observation_table.txt", row.names= 1,  header=T, sep="\t")
# otu_table1 = read.delim("bwa-short_map_L4_1.0/observation_table.txt",  header=T, sep="\t")


# sub_table <- otu_table[,c(ncol(otu_table),colnames(otu_table) %in%  rownames(sub_design))]
# otu_table[otu_table$taxonomy %in% "Unclassified",]




############################################################subSet OTU for analysis nature sample feature otu distribution ########################
sub_table <- otu_table[,colnames(otu_table) %in% rownames(sub_design)]
# sub_table <- otu_table[,colnames(otu_table) %in%  rownames(sub_design)]




############rarefraction#####################
## a consistent random seed [set.seed(21336)] was used for reproducibility

set.seed(21336)
# rrarefy The sample can be a vector giving the sample sizes for each row.so you need the transpose
##rarefy(x, sample, se = FALSE, MARGIN = 1) rrarefy(x, sample)
sub_table<- as.data.frame(t(rrarefy(t(sub_table),sample = min(colSums(sub_table)))))
### to check the size whether to be same
colSums(sub_table)


#################reorder######################
ids <- match(rownames(sub_design),colnames(sub_table))
sub_table_1 <- sub_table[,ids]
ids <- match(rownames(sub_design),colnames(otu_table))
otu_table_1 <- otu_table[,ids]

######################################Order by ZH11.Bac.BI055.01 for  Scal-BI-12-4 E05/5
sub_OTU<- sub_table_1[order(sub_table_1[,1],decreasing = T),]
sub_OTU$id <- rownames(sub_OTU)
otu_table_ca <- as.data.frame(otu_table[rownames(otu_table) %in% rownames(sub_OTU),c(1,ncol(otu_table))])
otu_table_ca$GOS.Fun.BI05.08 <- NULL
otu_table_ca$id <- rownames(otu_table_ca)
colnames(otu_table_ca)[1] <- "Taxa"

merge_OTU_taxa <- merge(sub_OTU,otu_table_ca,by ="id" )
write.table(merge_OTU_taxa,file = "../fun_New/table/rarefied_fungal_OTU_table_taxa.txt",sep = '\t',row.names = F)

OTU<-otu_table_1[order(otu_table_1[,1],decreasing = T),]

ids <- match(row.names(OTU),row.names(otu_table))
otu_table<- otu_table[ids,]

pso <- OTU[(rownames(OTU)) %in% "OTU_6",]
# 
# 
top_10_OTU_withSpike <- sub_OTU[1:11,]

OTU_top_11 <- OTU[1:11,]
OTU_top_11_Tax <- cbind(otu_table[1:11,]$taxonomy,OTU_top_11)
colnames(OTU_top_11_Tax)[1] <- "Taxa"

OTU_top_11 <- rbind(OTU_top_11,pso)
# row.names(top_10_OTU_withSpike) <- OTU_top_11


######E05######
sub_design_1 <- sub_design[!sub_design$PlasmidID %in% "Scal-BI-mixture",]


pos <-  colnames(OTU_top_11)%in% rownames(sub_design_1)
E5_mo <- OTU[,pos]
# rownames(E5_mo) <- OTU_top_11_Tax$Taxa

gra_group <- unique(design$PlasmidID)               # plasmid id 
# gra_group <- unique(design$Description)
mix_ratio_group <- unique(design$Other)
cl_group <- unique(design$Description)                 # Concentration



#-----------------------------------------------------------------------------------------------------------------------------------------------
# spike 
rownames(E5_mo)[1] <- "BI-OS-12-4"
idx <- match(spike[2],row.names(E5_mo))
# sub_table_1 <- sub_table[,sub_table[idx,]>160]
E5_mo$id <- rownames(E5_mo)


# mapping spike 12-4
# idx1 <- match(spike[1],row.names(E5_mo))
sub_table_2 <- E5_mo
# sub_table_2 <- E5_mo[-idx1,]

# idx1 <- match(spike[3],row.names(sub_table_2))
sub_table_2 <- sub_table_2



len <- length(sub_table_2[1,])
idx <- match(spike[2],row.names(sub_table_2))
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




colors <- data.frame(group=cl_group,
                     color=c(c_red, c_dark_brown,c_black,c_orange,c_blue,c_sea_green,c_yellow))

shapes <- data.frame(group=mix_ratio_group,
                     shape=c(19, 0, 24))

index$Other <- factor(index$Other, levels=shapes$group)
index$Description <- factor(index$Description, levels=colors$group)
index$Mix_Ratio <- index$Other
# reorder boxplots

# l <- c("soil", "rhizosphere", "root", "pooled_nodules")
# index$compartment <- factor(index$compartment, levels=l)
# colors <- colors[match(l, colors$group), ]


detailName <- paste("Natural_fungal_Community_Concentration",spike[2],sep = '_')

p <- ggplot(index, aes(x=Mix_Ratio, y=value, color=Description)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Mix_Ratio), position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=detailName, y="spike-in_Abundance") +
  theme(axis.text.x = element_text(size=7.5,angle = 90))+
  main_theme


p

ggsave(paste(figures.dir, detailName,"spike_in_Abundance.pdf", sep="_"), p)

## to save in the new directory
ggsave("../fun_New/figure/perturbation_unoise_fung/spike_in_Abundance.pdf",  p)
write.table(index,file = "./table/spike_in_Abundance_fun.xls",sep = '\t',row.names = T)


#plot RA 

idx1 <- match("spike-in_Abundance",colnames(t_table))
idx2 <- match("delta_value",colnames(t_table))
index <- cbind(design[match(t_table[,length(t_table[1,])], design$SampleID), ], total_table[, -c(idx2,idx1)])
index<- index[,c(1:25,grep("OTU_2$",colnames(index)))]

colors <- data.frame(group=cl_group,
                     color=c(c_red, c_dark_brown,c_black,c_orange,c_sea_green,c_yellow,c_green))
  

shapes <- data.frame(group=mix_ratio_group,
                     shape=c(19, 0, 24))

#,c_sea_green,c_yellow,c_green,c_very_dark_green,"purple","pink",c_dark_red

index$Other <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)
# index$Other <- factor(index$Other, levels=shapes$group)
index$Description <- factor(index$Description, levels=colors$group)
index$Mix_Ratio <- index$Other

df <- melt(index)

mean_table <- aggregate(df$value,by = list(df$Mix_Ratio,df$Description,df$variable),FUN = mean)
mean_table$variable <- rownames(mean_table)
colnames(mean_table)[4] <- 'value'
p <- ggplot(mean_table, aes(x=Group.3,y=value,color=Group.2)) +
  facet_wrap(Group.1~Group.2,scales = "free")+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Group.3), position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=detailName, y="RA_mean") +
  theme(axis.text.x = element_text(size=10,angle = 90))+
  main_theme

p


ggsave(paste(figures.dir,"Facet_Natural_fung_RA_mean_box_line_chart.pdf", sep=""), p)
ggsave("../fun_New/figure/perturbation_unoise_fung/Facet_Natural_fung_RA_mean_box_line_chart.pdf", p)
write.table(mean_table,file = "./table/Facet_Natural_fung_RA_mean_box_line_chart.xls",sep = '\t',row.names = T)

#### alpha diversity comparison 

shannon <- as.data.frame(diversity(sub_table_2[-1,-ncol(sub_table_2)],index = "shannon",MARGIN = 2))
colnames(shannon)[1] <- "shannon_value"

index <- cbind(sub_design_1,shannon)


##rep usage each and times
index$Other <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)
df <- melt(index)

p <- ggplot(df, aes(x=Description, y=value)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill=rep(c(c_red, c_dark_brown,c_black,c_orange),3)) +
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
  facet_wrap(~Other)+
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  main_theme


p

ggsave(paste(figures.dir, "alpha_shannon_boxplot.pdf", sep=""), p)
ggsave("../fun_New/figure/perturbation_unoise_fung/alpha_shannon_boxplot.pdf", p)
write.table(df,file = "./table/alpha_fun.xls",sep = '\t',row.names = T)

# p <- ggviolin(index,x="Description", y="shannon_value",fill = "Description",
#               palette = c(c_red, c_dark_brown,c_black,c_orange),
#               add = "boxplot", add.params = list(fill = "white"))+
#               stat_compare_means(comparisons = list(c("E05/5", "E00"), c("E4", "E00"),c("E03", "E00")), label = "p.signif")
#               # stat_compare_means(comparisons = my_comparisons, label = "p.signif")+  # Add significance levels
#               # stat_compare_means(label.y = 50)                                       # Add global the p-value


cl_group <- unique(sub_design_1$Description)
mix_ratio_group <- unique(index$Other)
colors <- data.frame(group=cl_group,
                     color=c(c_red, c_dark_brown,c_black,c_orange))



shapes <- data.frame(group=mix_ratio_group,
                     shape=c(19, 0, 24))
index$Other <- factor(index$Other, levels=shapes$group)
index$Description <- factor(index$Description, levels=colors$group)
index$Mix_Ratio <- index$Other
# reorder boxplots



p <- ggplot(index, aes(x=Description, y=shannon_value)) +
  geom_violin(trim=FALSE, fill="red")+
  geom_boxplot(width=0.1, fill="white")+
  # geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill=rep(c(c_red, c_dark_brown,c_black,c_orange),3)) +
  geom_jitter(aes(shape=Other), position=position_jitter(0.17), size=1, alpha=0.7) +
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
  facet_wrap(~Other)+
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  main_theme
p

ggsave(paste(figures.dir, "alpha_shannon_violin.pdf", sep=""), p)
ggsave("../fun_New/figure/perturbation_unoise_fung/alpha_shannon_violin.pdf", p)
write.table(index,file = "./table/fun_shannon_violin.xls",sep = '\t',row.names = T)



#### beta diversity #######
# source("https://bioconductor.org/biocLite.R")
# biocLite("metagenomeSeq")
# 
# library(metagenomeSeq)
# library(biomformat)
# 
#   
# n## Distances between samples was calculated as Bray-Curtis dissimilarities of subsampled and square-root transformed OTU counts
# 
# ##filter spike 
# otu_table_for_beta <- sub_table_2[-1,-ncol(sub_table_2)]
# # otu_table_for_beta <- sweep(otu_table_for_beta,2,colSums(otu_table_for_beta),FUN = '/')
# 
# otu_table_for_beta<- cbind(sub_table_2[-1,]$id,otu_table_for_beta)
# colnames(otu_table_for_beta)[1] <- "OTU_id"
# otu_choose <- otu_table_for_beta[1:20,]
# 
# 
# write.table(otu_choose,file = "../bac_New/table/otu_table_for_beta.txt",sep = '\t',quote = F,row.names = F,col.names = T)


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

bray_curtis <- read.table("../bac_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, check.names=F,stringsAsFactors = F)
# bray_curtis <- read.table("../bac_New/table/beta_div/euclidean_otu_table.txt", sep="\t", header=T, check.names=F,stringsAsFactors = F)
# 

##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~1
k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")
points$Concentration <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)

points <- cbind(points, sub_design_1[match(rownames(points), sub_design_1$SampleID), ])



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
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("../bac_New/figure/PCoA_BC.pdf", p)
ggsave("../bac_New/figure/eulcian_PCoA_BC.pdf", p)



############################### RA accuracy ####################
####A bac_perturbation B onefold  C twofold


mean_BA <- mean_table[!c(1:nrow(mean_table))%%3==0,]
mean_C <- mean_table[c(1:nrow(mean_table))%%3==0,]
mean_A <- mean_BA[!c(1:nrow(mean_BA))%%2==0,]
mean_B <- mean_BA[c(1:nrow(mean_BA))%%2==0,]

ratio_BA <- mean_B$value/mean_A$value
ratio_BC <- mean_B$value/mean_C$value

mean_A$ration <- ratio_BA
mean_C$ration <- ratio_BC


detailName <- "Feature OTUs in Natural Community"

######BA########
p <- ggplot(mean_A[!mean_A$Group.3 %in% "OTU_6",], aes(x=Group.3, y=ration,shape=Group.2,color=Group.2,group=Group.2)) +
  facet_wrap(~Group.2)+
  geom_line() +
  geom_point(size=3, fill="white") +
  scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
  scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
  labs(x=detailName, y="Relative Abundance Ration GroupB(one fold host) vs GroupA(bacteria spike perturbation)") +
  geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
  theme(axis.text.x = element_text(size=10,angle = 90))+
  # geom_hline(mean(mean_A$ration),colour='black',lwd=0.36,linetype="dashed")+
  main_theme

p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# p <- p+geom_hline(yintercept = mean(mean_A[!mean_A$Group.3 %in% "OTU_9",]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))


p

ggsave(paste(figures.dir,"RA_ratio_B_oneFold_vs_A_funSpikePertur_line_chart.pdf", sep=""), p)
ggsave("../fun_New/figure/perturbation_unoise_fung/RA_ratio_B_oneFold_vs_A_FunSpikePertur_line_chart.pdf", p)
write.table(mean_A,file = "./table/RA_ratio_B_oneFold_vs_A_FunSpikePertur_line_chart.xls",sep = '\t',row.names = T)


####BC###
p <- ggplot(mean_C[!mean_C$Group.3 %in% "OTU_6",], aes(x=Group.3, y=ration,shape=Group.2,color=Group.2,group=Group.2)) +
  facet_wrap(~Group.2)+
  geom_line() +
  geom_point(size=3, fill="white") +
  scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
  scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
  labs(x=detailName, y="Relative Abundance Ratio GroupB(one fold host) vs GroupC (two fold)") +
  geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
  theme(axis.text.x = element_text(size=10,angle = 90))+
  main_theme

p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# p <- p+geom_hline(yintercept = mean(ratio_AC[-c(16)]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))

p

ggsave(paste(figures.dir,"RA_ratio_B_oneFold_vsC_twoFold_line_chart.pdf", sep=""), p)
ggsave("../fun_New/figure/perturbation_unoise_fung/RA_ratio_B_oneFold_vsC_twoFold_line_chart.pdf", p)
write.table(mean_C,file = "./table/RA_ratio_B_oneFold_vsC_twoFold_line_chart.xls",sep = '\t',row.names = T)



# plot absolute abundance 


# mapping spike 12-4
# idx1 <- match(spike[1],row.names(E5_mo))
sub_table_2 <- E5_mo

# idx1 <- match(spike[3],row.names(sub_table_2))
sub_table_2 <- sub_table_2

len <- length(sub_table_2[1,])
idx <- match(spike[2],row.names(sub_table_2))
# col_sums <- colSums(sub_table_2[-idx,-len])

# idx <- match(spike[3],row.names(sub_table))
# len <- length(sub_table[1,])
internal_ref <- as.numeric(sub_table_2[idx,-len])
internal_ref=1+internal_ref
absoluteAbund <-as.data.frame(t(sub_table_2[-idx,-len])/internal_ref)


fun_reads_AA <- t(absoluteAbund)
write.table(fun_reads_AA,file = "../fun_New/table/fun_reads_AA.txt",sep = '\t',row.names = F,quote = F)




index <- cbind(design[match(row.names(absoluteAbund), design$SampleID), ],absoluteAbund )

index<- index[,c(1:25)]
index$Other <- rep(c("One_fold_plant","Two_fold_plant","Bac_spike_perturbation"),each=3,4)

# sub_design <- subset(design,design$Description %in% "E06")
# idx5 <- match(rownames(sub_design),colnames(absoluteAbund)) 
# absoluteAbund <- absoluteAbund[,idx5]
colors <- data.frame(group=cl_group,
                     color=c(c_red, c_dark_brown,c_blue,c_green))



shapes <- data.frame(group=mix_ratio_group,
                     shape=c(19, 0, 24))

index$Other <- factor(index$Other, levels=shapes$group)
index$Mix_Ratio <- index$Other
index$Description <- factor(index$Description,levels =colors$group )

index

df <- melt(index)

mean_table <- aggregate(df,by = list(df$Mix_Ratio,df$variable,df$Description),FUN = mean)
mean_table$variable <- rownames(mean_table)

detailName <- "Facet Natural Feature fungal"
p <- ggplot(mean_table, aes(x=Group.2,y=value,color=Group.2)) +
  facet_wrap(Group.1~Group.3,scales = "free")+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  geom_jitter(aes(shape=Group.3), position=position_jitter(0.17), size=1, alpha=0.7) +
  # scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x=detailName, y="AA_mean") +
  theme(axis.text.x = element_text(size=10,angle = 90))+
  main_theme

p


ggsave(paste(figures.dir,"Facet_Natural_fungal_AA_mean_box_line_chart.pdf", sep=""), p)
ggsave("../fun_New/figure/perturbation_unoise_fung/Facet_Natural_fungal_AA_mean_box_line_chart.pdf", p)
write.table(mean_table,file = "./table/Facet_Natural_fungal_AA_mean_box_line.xls",sep = '\t',row.names = T)


####AA Ratio


mean_BA <- mean_table[!c(1:nrow(mean_table))%%3==0,]
mean_C <- mean_table[c(1:nrow(mean_table))%%3==0,]
mean_A <- mean_BA[!c(1:nrow(mean_BA))%%2==0,]
mean_B <- mean_BA[c(1:nrow(mean_BA))%%2==0,]


ratio_AB <- mean_A$value/mean_B$value
ratio_AC <- mean_A$value/mean_C$value

mean_A$ration <- ratio_AB
mean_C$ration <- ratio_AC




detailName <- "Feature OTUs in Natural Community"

######AB########
p <- ggplot(mean_A, aes(x=Group.2, y=ration,shape=Group.2,color=Group.3,group=Group.3)) +
  facet_wrap(~Group.3)+
  geom_line() +
  geom_point(size=3, fill="white") +
  scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
  scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
  labs(x=detailName, y="Absolute Abundance Ration GroupA(one fold host) vs GroupB(two fold host)") +
  geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
  theme(axis.text.x = element_text(size=10,angle = 90))+
  main_theme

p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# p <- p+geom_hline(yintercept = mean(mean_A$ration[-c(16)]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))


p

ggsave(paste(figures.dir,"modi_AA_ratio_A_oneFold_vs_B_twofold_line_chart.pdf", sep=""), p)
ggsave("../fun_New/figure/perturbation_unoise_fung/modi_AA_ratio_A_oneFold_vs_B_twofold_line_chart.pdf", p)
write.table(mean_A,file = "./table/modi_AA_ratio_A_oneFold_vs_B_twofold_line_chart.xls",sep = '\t',row.names = T)






####AC###
p <- ggplot(mean_C, aes(x=Group.2, y=ration,shape=Group.2,color=Group.3,group=Group.3)) +
  facet_wrap(~Group.3)+
  geom_line() +
  geom_point(size=3, fill="white") +
  scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
  scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
  labs(x=detailName, y="Absolute Abundance Ratio GroupA(one fold host) vs GroupC (perturbation)") +
  geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
  theme(axis.text.x = element_text(size=10,angle = 90))+
  main_theme

p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# p <- p+geom_hline(yintercept = mean(ratio_AC[-c(16)]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))

p

ggsave(paste(figures.dir,"modiAA_ratio_A_oneFold_vsC_perturbation_line_chart.pdf", sep=""), p)
ggsave("../fun_New/figure/perturbation_unoise_fung/modi—AA_ratio_A_oneFold_vsC_perturbation_line_chart.pdf", p)
write.table(mean_C,file = "./table/modi_AA_ratio_A_oneFold_vsC_perturbation_line_chart.pdf.xls",sep = '\t',row.names = T)

