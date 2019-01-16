#### this script mainly to perform the PCA based on the euclian matrix calcucated from AA table of fungi


rm(list=ls())
options(warn=3)

########### setting the working directory and print it ###################
tem <- "wetDryRice"
setwd("~/xiaoxuan/180528/180706_fungal/script/")
print(paste("Your working directory is in",getwd()))


########### to import the plotting theme() function ########################
source("plot_function.R")



design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
sample_list <- as.matrix(row.names(design))
batch_group <- unique(design$Site)               # plasmid id 
# gra_group <- unique(design$Description)
condition_group <- unique(design$Condition)
geno_group <- unique(design$Genotype)                 # Concentration



bray_curtis <- read.table("../spike_rem_sample_720/beta_720/euclidean_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../Fun_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../Fun_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 

design = design[rownames(design) %in% rownames(bray_curtis),]
id <- match(row.names(design[design$Site %in% "AnHui"&!design$Genotype %in% "MH63ZH",]),rownames(bray_curtis))



ah_bray_curtis <- bray_curtis[id,id]
design <- design[rownames(design) %in% rownames(ah_bray_curtis),]
##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~1
pcoa <-prcomp(ah_bray_curtis)



# pcoa <- cmdscale(ah_bray_curtis, k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$rotation
eig <- pcoa$sdev
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


ggsave("../spike_rem_sample_720/PCA12_BC_Condition_Ah.pdf", p)



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
# ggsave("../spikein/PCA13_BC_Site.pdf", p)

p <- ggplot(points, aes(x=x, y=z, color=Genotype,shape=Condition)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p


ggsave("../spike_rem_sample_720//PCA13_BC_Condition_Ah.pdf", p)



ggsave("../Fun_New/figure/eulcian_PCoA_BC.pdf", p)
write.table(points,file = "../table/WD_Natural_fungi_betaDiv_points.txt",sep = '\t',row.names = F,quote = F)


###AA HaiNan
bray_curtis <- read.table("../spike_rem_sample_720/beta_720/euclidean_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../Fun_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../Fun_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 
design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
design = design[rownames(design) %in% rownames(bray_curtis),]
id <- match(row.names(design[design$Site %in% "HaiNan"&!design$Genotype %in% "MH63ZH",]),rownames(bray_curtis))



ah_bray_curtis <- bray_curtis[id,id]
##  blank \t  need sed to filter  

#### Unconstrained ordination
# pcoa<-capscale(dis~1
pcoa <-prcomp(ah_bray_curtis)



# pcoa <- cmdscale(ah_bray_curtis, k=k, eig=T)
# pcoa <- cmdscale(eulian, k=k, eig=T)
points <- pcoa$rotation
eig <- pcoa$sdev
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
  labs(x=paste("PCA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  # facet_wrap(~ Funtch)+
  main_theme +
  theme(legend.position="top")
p


ggsave("../spike_rem_sample_720/PCA12_BC_Condition_Hainan.pdf", p)



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
  labs(x=paste("PCA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p


ggsave("../spike_rem_sample_720/PCA13_BC_Condition_HaiNan.pdf", p)


p <- ggplot(points, aes(x=y, y=z, color=Genotype,shape=Condition)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p


ggsave("../spike_rem_sample_720/PCA23_BC_Condition_HaiNan.pdf", p)
ggsave("../Fun_New/figure/eulcian_PCoA_BC.pdf", p)
write.table(points,file = "../table/WD_Natural_fungi_betaDiv_points.txt",sep = '\t',row.names = F,quote = F)


######RA method
design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
bray_curtis <- read.table("../RA/beta/bray_curtis_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../Fun_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../Fun_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 
design = design[rownames(design) %in% rownames(bray_curtis),]
id <- match(row.names(design[design$Site %in% "AnHui"&!design$Genotype %in% "MH63ZH",]),rownames(bray_curtis))



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
                     color=c(c_red, c_dark_brown,c_black,c_orange))

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

ggsave("../RA/718_fungi_PCoA12_BC_Conditon_ah.pdf", p)

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

ggsave("../RA/718_fungi_PCoA13_BC_Condition_ah.pdf", p)

p <- ggplot(points, aes(x=y, y=z, color=Genotype,shape=Condition)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("../RA/718_fungi_PCoA23_BC_Condition_ah.pdf", p)




design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
bray_curtis <- read.table("../RA/beta/bray_curtis_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
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
                     color=c(c_red, c_dark_brown,c_black,c_orange))

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

ggsave("../RA/718_fungi_PCoA12_BC_Conditon_HaiNan.pdf", p)

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

ggsave("../RA/718_fungi_PCoA13_BC_Condition_Hainan.pdf", p)


p <- ggplot(points, aes(x=y, y=z, color=Genotype,shape=Condition)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top")
p

ggsave("../RA/718_fungi_PCoA23_BC_Condition_Hainan.pdf", p)
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

ggsave("../Fun_New/figure/PCoA_W_unifrac.pdf", p)








############################### RA accuracy ####################
# ####A Fun_perturbation B onefold  C twofold
# 
# 
# mean_BA <- mean_table[!c(1:nrow(mean_table))%%3==0,]
# mean_C <- mean_table[c(1:nrow(mean_table))%%3==0,]
# mean_A <- mean_BA[!c(1:nrow(mean_BA))%%2==0,]
# mean_B <- mean_BA[c(1:nrow(mean_BA))%%2==0,]
# 
# ratio_BA <- mean_B$value/mean_A$value
# ratio_BC <- mean_B$value/mean_C$value
# 
# mean_A$ration <- ratio_BA
# mean_C$ration <- ratio_BC
# 
# 
# detailName <- "Feature OTUs in Natural Community"
# 
# ######BA########
# p <- ggplot(mean_A[!mean_A$Group.3 %in% "OTU_9",], aes(x=Group.3, y=ration,shape=Group.2,color=Group.2,group=Group.2)) +
#   facet_wrap(~Group.2)+
#   geom_line() +
#   geom_point(size=3, fill="white") +
#   scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
#   scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
#   labs(x=detailName, y="Relative Abundance Ration GroupB(one fold host) vs GroupA(Funteria spike perturbation)") +
#   geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
#   theme(axis.text.x = element_text(size=10,angle = 90))+
#   # geom_hline(mean(mean_A$ration),colour='black',lwd=0.36,linetype="dashed")+
#   main_theme
# 
# p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # p <- p+geom_hline(yintercept = mean(mean_A[!mean_A$Group.3 %in% "OTU_9",]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# 
# 
# p
# 
# # ggsave(paste(figures.dir,"RA_ratio_B_oneFold_vs_A_FunSpikePertur_line_chart.pdf", sep=""), p)
# ggsave("../Fun_New/figure/perturbation_unoise_Fun/silva_RA_ratio_B_oneFold_vs_A_FunSpikePertur_line_chart.pdf", p)
# write.table(mean_A,file = "../Fun_New/table/silva_RA_ratio_B_oneFold_vs_A_FunSpikePertur_line_chart_accuracy.txt",sep = '\t',row.names = F,quote = F)
# 
# ####BC###
# p <- ggplot(mean_C[!mean_C$Group.3 %in% "OTU_9",], aes(x=Group.3, y=ration,shape=Group.2,color=Group.2,group=Group.2)) +
#   facet_wrap(~Group.2)+
#   geom_line() +
#   geom_point(size=3, fill="white") +
#   scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
#   scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
#   labs(x=detailName, y="Relative Abundance Ratio GroupB(one fold host) vs GroupC (two fold)") +
#   geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
#   theme(axis.text.x = element_text(size=10,angle = 90))+
#   main_theme
# 
# p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # p <- p+geom_hline(yintercept = mean(ratio_AC[-c(16)]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# 
# p
# 
# # ggsave(paste(figures.dir,"RA_ratio_B_oneFold_vsC_twoFold_line_chart.pdf", sep=""), p)
# ggsave("../Fun_New/figure/perturbation_unoise_Fun/silva_RA_ratio_B_oneFold_vsC_twoFold_line_chart.pdf", p)
# write.table(mean_C,file = "../Fun_New/table/silva_RA_ratio_B_oneFold_vsC_twoFold_line_chart_accuracy.txt",sep = '\t',row.names = F,quote = F)
# 
# 
# # plot absolute abundance 
# 
# 
# # mapping spike 12-4
# # idx1 <- match(spike[1],row.names(sub_OTU))
# sub_table_2 <- sub_OTU
# 
# # idx1 <- match(spike[3],row.names(sub_table_2))
# # sub_table_2 <- sub_table_2
# 
# len <- length(sub_table_2[1,])
# idx <- match(spike[2],row.names(sub_table_2))
# # col_sums <- colSums(sub_table_2[-idx,-len])
# 
# # idx <- match(spike[3],row.names(sub_table))
# # len <- length(sub_table[1,])
# internal_ref <- as.numeric(sub_table_2[idx,-len])
# internal_ref=1+internal_ref
# absoluteAbund <-as.data.frame(t(sub_table_2[-idx,-len])/internal_ref)
# 
# Fun_reads_AA <- t(absoluteAbund)
# 
# micro_load <- as.data.frame(colSums(Fun_reads_AA))
# write.table(Fun_reads_AA,file = "../table/Fun_reads_AA.txt",sep = '\t',row.names = F,quote = F)
# 
# 
# 
# 
# # index <- cbind(design[match(row.names(absoluteAbund), design$SampleID), ],absoluteAbund )
# # 
# # index<- index[,c(1:25)]
# # index$Other <- rep(c("One_fold_plant","Two_fold_plant","Fun_spike_perturbation"),each=3,4)
# # 
# # # design <- subset(design,design$Description %in% "E06")
# # # idx5 <- match(rownames(design),colnames(absoluteAbund)) 
# # # absoluteAbund <- absoluteAbund[,idx5]
# # colors <- data.frame(group=cl_group,
# #                      color=c(c_red, c_dark_brown,c_blue,c_green))
# # 
# # 
# # 
# # shapes <- data.frame(group=mix_ratio_group,
# #                      shape=c(19, 0, 24))
# # 
# # index$Other <- factor(index$Other, levels=shapes$group)
# # index$Mix_Ratio <- index$Other
# # index$Description <- factor(index$Description,levels =colors$group )
# # 
# # index
# # 
# # df <- melt(index)
# # 
# # mean_table <- aggregate(df$value, by = list(df$Mix_Ratio,df$variable,df$Description),FUN = mean)
# # mean_table$variable <- rownames(mean_table)
# # colnames(mean_table)[4] <- "value"
# # detailName <- "Facet Natural Feature Funteria "
# # p <- ggplot(mean_table, aes(x=Group.2,y=value,color=Group.2)) +
# #   facet_wrap(Group.1~Group.3,scales = "free")+
# #   geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
# #   geom_jitter(aes(shape=Group.3), position=position_jitter(0.17), size=1, alpha=0.7) +
# #   # scale_colour_manual(values=as.character(colors$color)) +
# #   # scale_shape_manual(values=shapes$shape) +
# #   labs(x=detailName, y="AA_mean") +
# #   theme(axis.text.x = element_text(size=10,angle = 90))+
# #   main_theme
# # 
# # p
# # 
# # 
# # # ggsave(paste(figures.dir,"Facet_Natural_Funteria_AA_mean_box_line_chart.pdf", sep=""), p)
# # 
# # ggsave("../Fun_New/figure/perturbation_unoise_Fun/silva_Facet_Natural_Funteria_AA_mean_box_line_chart.pdf", p)
# # write.table(mean_table,file = "../Fun_New/table/silva_Facet_Natural_Funteria_AA_mean_box_line.txt",sep = '\t',row.names = F,quote = F)
# # 
# # ####AA Ratio
# # 
# # 
# # mean_BA <- mean_table[!c(1:nrow(mean_table))%%3==0,]
# # mean_C <- mean_table[c(1:nrow(mean_table))%%3==0,] # disturbation
# # mean_A <- mean_BA[!c(1:nrow(mean_BA))%%2==0,]  # one fold
# # mean_B <- mean_BA[c(1:nrow(mean_BA))%%2==0,]  # two fold
# # 
# # ratio_AB <- mean_A$value/mean_B$value
# # ratio_AC <- mean_A$value/mean_C$value
# # 
# # mean_A$ration <- ratio_AB
# # mean_C$ration <- ratio_AC
# # 
# # 
# # detailName <- "Feature OTUs in Natural Community"
# # 
# # ######AB########
# # p <- ggplot(mean_A, aes(x=Group.2, y=ration,shape=Group.2,color=Group.3,group=Group.3)) +
# #   facet_wrap(~Group.3)+
# #   geom_line() +
# #   geom_point(size=3, fill="white") +
# #   scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
# #   scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
# #   labs(x=detailName, y="Absolute Abundance Ration GroupA(one fold host) vs GroupB(two fold host)") +
# #   geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
# #   theme(axis.text.x = element_text(size=10,angle = 90))+
# #   main_theme
# # 
# # p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # # p <- p+geom_hline(yintercept = mean(mean_A$ration[-c(16)]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # 
# # 
# # p
# # 
# # # ggsave(paste(figures.dir,"modified_AA_ratio_B_twoFold_vs_A_onefold_line_chart.pdf", sep=""), p)
# # ggsave("../Fun_New/figure/perturbation_unoise_Fun/silva_modified_AA_ratio_A_onefold_vsB_twoFold__line_chart.pdf", p)
# # 
# # write.table(mean_A,file = "../Fun_New/table/silva_modified_Facet_Natural_Funteria_AA_chart_line_AB.txt",sep = '\t',row.names = F,quote = F)
# # 
# # 
# # 
# # ####AC###
# # p <- ggplot(mean_C, aes(x=Group.2, y=ration,shape=Group.2,color=Group.3,group=Group.3)) +
# #   facet_wrap(~Group.3)+
# #   geom_line() +
# #   geom_point(size=3, fill="white") +
# #   scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_green))) +
# #   scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
# #   labs(x=detailName, y="Absolute Abundance Ratio GroupA(one fold host) vs GroupC (perturbation)") +
# #   geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
# #   theme(axis.text.x = element_text(size=10,angle = 90))+
# #   main_theme
# # 
# # p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # # p <- p+geom_hline(yintercept = mean(ratio_AC[-c(16)]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# # 
# # p
# # 
# # # ggsave(paste(figures.dir,"AA_ratio_A_oneFold_vsC_perturbation_line_chart.pdf", sep=""), p)
# # ggsave("../Fun_New/figure/perturbation_unoise_Fun/silva_AA_ratio_A_oneFold_vsC_perturbation_line_chart.pdf", p)
# # write.table(mean_C,file = "../Fun_New/table/silva_modified_Facet_Natural_Funteria_AA_mean_chart_line_AC.txt",sep = '\t',row.names = F,quote = F)
# # 
