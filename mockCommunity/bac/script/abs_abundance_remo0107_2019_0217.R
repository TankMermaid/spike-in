## This scripot is meant to draw the relative abundance with SpikeIn and without SpikeIn 
## these scripts are based on the single_rarefraction result
## Tank 

## Email:xnzhang@genetics.ac.cn

rm(list = ls())
setwd("/mnt/bai/xiaoning/xiaoxuan/180213/bac_New/")
print(paste("Your working directory is in",getwd()))

library(ggplot2)
library(reshape)
library(multcomp)
library(ggsignif)
library(vegan)
source("script/plot_function.R")

### output directory assigned to include the pics & tables########################
# figures.dir <- paste("/mnt/bai/xiaoning/xiaoxuan/180213/bac_New/figure/",tem,'/',sep = '')
# table.dir <- paste("/mnt/bai/xiaoning/xiaoxuan/180213/bac_New/table/",tem,'/',sep = '')


design = read.table("doc/design.txt", header=T, row.names= 1, sep="\t") 

## the applicaiton of subset & !  for short without the $
## if row not to use select in subset function
sub_design <- subset(design,PlasmidID == "Scal-BI-12-4"& Other %in% c("2:02:02","2:02:01","1:01:01" ))
# sub_design <- subset(design,PlasmidID == "Scal-BI-12-4" & Other == "2:02:02")


## remove outlier sample but 2019 0724 modified back 
# pos1 <- grep("*01$|*07|*08$",rownames(sub_design))
# sub_design_filt <- sub_design[-pos1,]
sub_design_filt <- sub_design
# sub_design_filt <- sub_design

l1 = read.delim("data/usearch_map_L1/observation_table.txt", row.names= 1,  header=T, sep="\t")
l2 = read.delim("data/usearch_map_L2/observation_table.txt", row.names= 1,  header=T, sep="\t")
l3 = read.delim("data/usearch_map_L3/observation_table.txt", row.names= 1,  header=T, sep="\t")
l4 = read.delim("data/usearch_map_L4/observation_table.txt", row.names= 1,  header=T, sep="\t")
l5 = read.delim("data/usearch_map_L5/observation_table.txt", row.names= 1,  header=T, sep="\t")

## modification of control with spike-in with dummy variable to facilitate the merge operation
dumm <- replicate(15,0)
l5 <- rbind(l5,dumm)
rownames(l5)[12] <- "BI-OS-12-4"

l1$ID <- rownames(l1)
l2$ID <- rownames(l2)
l3$ID <- rownames(l3)
l4$ID <- rownames(l4)
l5$ID <- rownames(l5)


merge_sp12 <- merge(l1,l2, by = "ID")
merge_sp34 <- merge(l3,l4, by = "ID")
merge_sp <- merge(merge_sp12,merge_sp34,by = "ID")

## add l5 as control E00
merge_sp <- merge(merge_sp,l5)
rownames(merge_sp) <- merge_sp$ID
merge_sp$ID <- NULL


## ----Pseudocount to  escape the error of division of 0---------------------------------------------------------
sub_table <- merge_sp

##------------------filter samples not mapping to sub_design-----------------------
sub_table_f <- sub_table[,colnames(sub_table) %in% rownames(sub_design_filt)]
sub_design_filt <- sub_design_filt[rownames(sub_design_filt) %in% colnames(sub_table_f),]

## remove "BI-OS-10-2"
sub_table_f <- sub_table_f[-c(4,6),]
set.seed(21336)
# rrarefy The sample can be a vector giving the sample sizes for each row.so you need the transpose
##rarefy(x, sample, se = FALSE, MARGIN = 1) rrarefy(x, sample)
sub_table_f<- as.data.frame(t(rrarefy(t(sub_table_f),sample = min(colSums(sub_table_f)))))
colSums(sub_table_f)
write.table(sub_table,file = "../bac_New/table/bac_OTU_rarefy.txt",sep = '\t',row.names = T,quote = F)

# sub_table_f<- read.table( "../bac_New/table/20190218_bac_OTU_rarefy.txt",sep = '\t',row.names = 1,header = T)
# sub_table_f <- sub_table_f[,match(rownames(sub_design_filt),colnames(sub_table_f))]

#-----------------------------------------------------------------------------------------------------------------------------------------------
# spike BI-12-4
##### spike-in design in this batch #####################
spike <- c("BI-OS-11-3","BI-OS-12-4","BI-OS-10-2")
idx <- match(spike[2],row.names(sub_table_f))
sub_table_2 <- sub_table_f[-idx,]


####absAbundance relative to spike-in
absAbundance <- sweep(sub_table_2,2,as.numeric(sub_table_f[idx,]),'/')

####### reorder ##########

spike_host_ratio = read.table("doc/spike_host_ratio.txt", header=T, row.names= 1, sep="\t") 

spike_host_ratio <- spike_host_ratio[rownames(spike_host_ratio) %in% rownames(sub_design_filt),]

ord <- match(rownames(sub_design_filt),colnames(absAbundance))   ## column rearrange
# ord1 <- match(rownames(bac_li), rownames(sub_table))    ## row reorder
absAbundance<- absAbundance[,ord]   ## parrellel
spike_host_ratio<- spike_host_ratio[ord,]

absAbundance <- sweep(absAbundance,2,as.numeric(spike_host_ratio$spike_host_ratio),'*')



write.table(absAbundance,file = "./table/20190725_mockBacAAremo010708.xls",sep = '\t',row.names = T)



###### join design and melt for mean and ggplot 
melt_design <- sub_design_filt[rownames(sub_design_filt) %in% colnames(absAbundance),]
melt_design$spike_concentration <-NULL 
idx <- match(rownames(melt_design),rownames(t(absAbundance)))
t_ra <- as.data.frame(t(absAbundance)[idx,])

df <- cbind(melt_design,t_ra)


##### Mean and variance calculation 
index <- melt(df)
mean_table <- aggregate(index$value,by = list(index$Other,index$Description,index$variable),FUN = mean)
mean_table$variable <- rownames(mean_table)
mean_table$abs_value <- mean_table$x


mean_AB <- mean_table[!c(1:nrow(mean_table))%%3==0,]
mean_C <- mean_table[c(1:nrow(mean_table))%%3==0,]  ##### C 2:2:2
mean_A <- mean_AB[!c(1:nrow(mean_AB))%%2==0,]       ##### A 1:1:1
mean_B <- mean_AB[c(1:nrow(mean_AB))%%2==0,]        ##### B 2:2:1

ratio_CA <- mean_C$abs_value/mean_A$abs_value
ratio_CB <- mean_C$abs_value/mean_B$abs_value
ratio_BA <- mean_B$abs_value/mean_A$abs_value

mean_A$ration_CA <- ratio_CA
mean_C$ration_CB <- ratio_CB
mean_B$ration_CA <- ratio_BA

mean_A$x <- NULL


pre_BA <- mean_B[c(1:30),]
pre_BA$x <- NULL

# colnames(pre_BA)[7] <- "ration_CA"

aft_BA <- mean_B[c(31:50),]
aft_BA$x <- NULL
# colnames(aft_BA)[6] <- "ration_CA"


mean_A <- rbind(mean_A,pre_BA)


# write.table(mean_A,file = "./table/mean_A_111_and_CA_ratio.xls",sep = '\t',row.names = T)
write.table(mean_A,file = "./table/mean_A_111_and_CA_ratio_merge.xls",sep = '\t',row.names = T)
write.table(mean_C,file = "./table/mean_C_222_and_CB_ratio.xls",sep = '\t',row.names = T)
write.table(mean_B,file = "./table/mean_B_221.xls",sep = '\t',row.names = T)
write.table(mean_B,file = "./table/mean_B_221_and_BA_ratio.xls",sep = '\t',row.names = T)



# p <- ggplot(mean_table, aes(x=Group.2,y=value,color=Group.1)) +
#   facet_wrap(~Group.1,scales = "free")+
#   geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
#   geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x="Combined spike Feature OTUs in Natural Community", y="AA_mean") +
#   theme(axis.text.x = element_text(size=10,angle = 90))+
#   main_theme
# 
# p
# 
# 
# ggsave(paste(figures.dir,"Facet_Natural_bacteria_Combinedspike_AA_mean_box_line_chart.pdf", sep=""), p)

######CA########
p <- ggplot(mean_A[!mean_A$Group.3 %in% "Act-322" & !mean_A$Group.2 %in% "E00",], aes(x=Group.3, y=ration_CA,shape=Group.3,group=Group.1)) +
  geom_line() +
  facet_wrap(Group.1~Group.2)+
  geom_point(size=3, fill="white") +
  scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black))) +
  scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
  labs(x="Mock Taxa", y="Absolute Abundance Ratio GroupC(2:2:2) vs GroupA(1:1:1)") +
  geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
  theme(axis.text.x = element_text(size=10,angle = 90))+
  main_theme

p <- p+geom_hline(yintercept = 2,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
p <- p+geom_hline(yintercept = 1.5,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# p <- p+geom_hline(yintercept = mean(ratio_AB[-c(2,16)]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))

p

ggsave("./figure/AA_ratio_CvsA_line_chart.pdf", p)


# #### error
# expectedRatio <- 2
# AA_ratio <- mean_A[!mean_A$Group.3 %in% "Act-322" & mean_A$Group.2 %in% "E05/20",6]
# 
# error_AA <- log2(expectedRatio)-log2(AA_ratio)
# df <- data.frame(x=rep("bacMock_E5/20",9),
#                  y=error_AA)
# 
# p <- ggplot(df,aes(x=x,y=y))+
#   geom_boxplot()+xlab("Absolute abundance approach")+
#   theme(axis.text.x = element_text(size=10))+
#   geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) +
#   ylab(expression(paste(atop(paste("Error",sep=""), "Expected vs Observed"))))+
#   main_theme
#  
# 
# p




#### error
expectedRatio <- 2
AA_ratio <- mean_A[!mean_A$Group.3 %in% "Act-322" ,]

AA_ratio$error_AA <- log2(expectedRatio)-log2(AA_ratio$ration_CA)

AA_ratio<- AA_ratio[!AA_ratio$Group.2 %in% "E00",]
p <- ggplot(AA_ratio,aes(x=AA_ratio$Group.2,y=error_AA,color=Group.2))+
  geom_boxplot()+xlab("Absolute abundance approach")+
  theme(axis.text.x = element_text(size=10))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) +
  ylab(expression(paste(atop(paste("Error",sep=""), "Expected vs Observed"))))+
  
  main_theme


p

ggsave("./figure/AA_ratio_error_boxplot.pdf", p)

### errror aft_BA
expectedRatio <- 1
aft_BA <- aft_BA[!aft_BA$Group.3 %in% "Act-322" ,]

aft_BA$error_AA <- log2(expectedRatio)-log2(aft_BA$ration_CA)

aft_BA<- aft_BA[!aft_BA$Group.2 %in% "E00",]
p <- ggplot(aft_BA,aes(x=aft_BA$Group.2,y=error_AA,color=Group.2))+
  geom_boxplot()+xlab("Absolute abundance approach")+
  theme(axis.text.x = element_text(size=10))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) +
  ylab(expression(paste(atop(paste("Error",sep=""), "Expected vs Observed"))))+
  
  main_theme


p

ggsave("./figure/aftAB_AA_ratio_error_boxplot.pdf", p)


## AA merge
AA_ar<- AA_ratio <- rbind(AA_ratio,aft_BA)


##### error AA and RA together
RA_ar <- ra_error<- read.table("table/RA_ratio_error_merge.txt", header=T,sep="\t") 


AA_ratio$approach <- rep("AA",length(AA_ratio$Group.1))
ra_error$approach <- rep("RA",length(ra_error$Group.1))

colnames(ra_error)[7] <- "error"
colnames(AA_ratio)[7] <- "error"

df_error <- rbind(AA_ratio,ra_error)

df_error<- df_error[!is.na(df_error$ration_CA),]
p <- ggplot(df_error,aes(x=df_error$Group.2,y=error,color=Group.2))+
  facet_wrap(~approach)+
  geom_boxplot()+xlab("Absolute abundance approach               vs               Relative abundance approach")+
  theme(axis.text.x = element_text(size=10))+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7) +
  ylab(expression(paste(atop(paste("Error",sep=""), "Expected vs Observed"))))+
  
  main_theme


p



ggsave("./figure/AA_ratio_vsRA_ratio_error_boxplot_merge.pdf", p)

write.table(df_error,file = "./table/AAvsRA_ratio_error_merge.txt",sep = '\t',row.names = F,quote = F)


## AA and RA column merge

# r2<-sprintf("italic(R^2)==%.3f",1)
ar_data <- cbind(AA_ar,RA_ar)
p <- ggplot(ar_data, aes(x=abs(error_RA), y=abs(error_AA),shape=Group.2,color=Group.2) )+
  # facet_wrap(Group.1~Group.2)+
  # geom_point(size=3, fill="white") +
  scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black,c_orange))) +
  scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
  labs(x="RA_error", y="AA_error") +
  geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
  theme(axis.text.x = element_text(size=10,angle = 90))+
  geom_abline(intercept = 0,slope=1,)+
  
  # p<-p+geom_text(x=1,aes(label=r2,y=y),data=labels,parse=TRUE,hjust=0,size=8)  
  main_theme

p 

ggsave("./figure/AA_ratio_vsRA_ratio_error_line—plot_merge.pdf", p)

write.table(ar_data,file = "./table/AAvsRA_ratio_error_col——merge.txt",sep = '\t',row.names = F,quote = F)




####CB###
p <- ggplot(mean_C[!mean_C$Group.3 %in% "Act-322" & !mean_C$Group.2 %in% "E00",], aes(x=Group.3, y=ration_CB,shape=Group.3,group=Group.1)) +
  geom_line() +
  facet_wrap(Group.1~Group.2)+
  geom_point(size=3, fill="white") +
  scale_colour_manual(values=as.character(c(c_red, c_dark_brown,c_black))) +
  scale_shape_manual(values=c(19,0,24,5,10,16,13,2,20,11,3,4,7,8,9,1)) +
  labs(x="Mock Taxa", y="Absolute Abundance Ratio GroupC(2:2:2) vs GroupB(2:2:1)") +
  geom_jitter(aes(shape=Group.2), position=position_jitter(0.17), size=1, alpha=0.7) +
  theme(axis.text.x = element_text(size=10,angle = 90))+
  main_theme

p <- p+geom_hline(yintercept = 2,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
p <- p+geom_hline(yintercept = 1.5,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
p <- p+geom_hline(yintercept = 1,colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))
# p <- p+geom_hline(yintercept = mean(ratio_AB[-c(2,16)]),colour='black',lwd=0.36,linetype="dashed",ylab("Expected Value"))

p

ggsave("./figure/AA_ratio_CvsB_line_chart.pdf", p)

#### setting parametres of plotting
Mock_ratio <- unique(melt_design$Other)
Concentration <- unique(melt_design$Description)

index <- melt(df)
Species <- unique(index$variable)

colors <- data.frame(group=Mock_ratio,
                     color=c(c_red, c_dark_brown,c_black))

colors <- data.frame(group=Species,
                     color=c(c_red, c_dark_brown,c_black,c_dark_brown,c_dark_red,c_green,c_orange,c_sea_green,c_very_dark_green,c_yellow))



shapes <- data.frame(group=Species,
                     shape=c(19, 0, 24,11,6,18,7,16,2,15))

color_con <- data.frame(group=Concentration,
                        color=c(c_red, c_dark_brown,c_black,c_green,c_yellow))

index$Species <- factor(index$variable, levels=shapes$group)
index$Mock_ratio <- factor(index$Other, levels=colors$group)
index$Concentration <- factor(index$Description,levels =color_con$group )


p <- ggplot(index[Mock_ratio %in% "2:02:02",], aes(x=Concentration, y=value, shape=Species,color=Species)) +
  # facet_wrap(~,scales = "free")+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
  # geom_point()+
  geom_jitter(aes(shape=Species), position=position_jitter(0.07), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x="Mock Concentration Ratio", y="spike-in_Abundance") +
  theme(axis.text.x = element_text(size=10,angle = 90))+
  main_theme

p

ggsave("./figure/mock_concentrationVSspike_in_Abundance.pdf", p)

# p <- ggplot(index, aes(x=Species, y=value, color=Mock_ratio)) +
#   facet_wrap(~Concentration,scales = "free")+
#   geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
#   geom_jitter(aes(shape=Concentration), position=position_jitter(0.07), size=1, alpha=0.7) +
#   scale_colour_manual(values=as.character(colors$color)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x="Mock ratio", y="spike-in_Abundance") +
#   theme(axis.text.x = element_text(size=10,angle = 90))+
#   main_theme
# 
# p



#### Dunnett's Test with Control

for( i in c(10:ncol(df))){
  
  df_1 <- df[df$Other %in% "2:02:02",c(9,i)]
  df_1$Concentration <- factor(df_1$Description)
  
  aov <- aov( df_1[,2]~ Concentration, df_1)
  dunnett_test <- summary(glht(aov, linfct=mcp(Concentration="Dunnett")))
  
  # shapes <- data.frame(group=Species,
  #                      shape=c(19, 0, 24,11,6,18,7,16,2,15))
  
  color_con <- data.frame(group=Concentration,
                          color=c(c_red, c_dark_brown,c_black,c_green,c_yellow))
  
  
  p<- ggplot(df_1, aes(x=Description, y=df_1[,2])) +
    geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.4, fill="transparent") +
    geom_jitter(aes(color=Concentration), position=position_jitter(0.17), size=1, alpha=0.7) +
    scale_colour_manual(values=as.character(color_con$color)) +
    # scale_shape_manual(values=shapes$shape) +
    labs(x="Mock ratio", y="spike-in_Abundance") +
    geom_signif(comparisons = list(c("E05", "E00"), c("E05/5", "E00"),c("E05/10", "E00"),c("E05/20", "E00")),
                map_signif_level = TRUE, textsize=6)
  
  p
  ggsave(paste("./figure/dunnett_test",i,".pdf",sep = ''), p)
  
  # write.csv(colnames(df_1)[2],file = "table/dunnett_test.xls",append = TRUE)
  # write.csv(dunnett_test,file = "table/dunnett_test.xls",append = TRUE)
  
  
}



## ----ErrorBackground-----------------------------------------------------
#Calculate error between expected and observed log2 ratios:














