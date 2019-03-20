rm(list=ls())
options(warn=3)


#### in case of sink() function

#########The taxa are ordered by the significance of the correlation between their QMP abundance 

########### setting the working directory and print it ###################
setwd("~/xiaoxuan/180724/CoNet/script/")
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

# library("ggplot2", quietly=T, warn.conflicts=F)
# library("gplots", quietly=T, warn.conflicts=F)
# library("grid", quietly=T, warn.conflicts=F)
library("RColorBrewer", quietly=T, warn.conflicts=F)
library("reshape2", quietly=T, warn.conflicts=F)
library(dplyr)
library(tidyr)
library(psych)
# library(tidyverse)
# library("VennDiagram", quietly=T, warn.conflicts=F)


## to set the date as the version
########flag==TRUE,usearch ; FALSE, unoise; project_name:.project; .version: as output version

.flag <- TRUE
# .flag <- FALSE
## to set the date as the version
.date <- Sys.Date()

### sub version swich
.subflag <- TRUE
# .subv <- FALSE

if(.subflag == TRUE){
.sub_ver <- "Top300"
threashold=300
}else
{
  .sub_ver <- ''
}


if(.flag==TRUE)
{.meth <- "Usearch" 
}else
{.meth <- "Unoise"}


.project <- "ConetRootDis"
.version <- paste0(.date,.project,.sub_ver,.meth)



### output directory assigned to include the pics & tables########################
figures.dir <- paste0("~/xiaoxuan/180724/CoNet/",.meth,'/',"plot/",sep = '')
table.dir   <- paste0("~/xiaoxuan/180724/CoNet/",.meth,'/',"table/",sep = '')


if(.flag==TRUE){
  
### OTU AA analysis
.botu <- paste0("~/xiaoxuan/180724/CoNet/",.meth,"/AA/otutab_norm_AA.txt")
.fotu <- paste0("~/xiaoxuan/180724/CoNet/",.meth,"/fAA/otutab_norm_AA.txt")

#### Tax
.btax <- paste0("~/xiaoxuan/180724/CoNet/",.meth,"/Tax/Tax.txt")
.ftax <- paste0("~/xiaoxuan/180724/CoNet/",.meth,"/fTax/fTax.txt")

 
}else
{
  ### OTU AA analysis
  .botu <- paste0("~/xiaoxuan/180724/CoNet/",.meth,"/AA/otutab_norm_AA.txt")
  .fotu <- paste0("~/xiaoxuan/180724/CoNet/",.meth,"/fAA/otutab_norm_AA.txt")
  
  #### Tax
  .btax <- paste0("~/xiaoxuan/180724/CoNet/",.meth,"/Tax/Tax.txt")
  .ftax <- paste0("~/xiaoxuan/180724/CoNet/",.meth,"/fTax/fTax.txt")
  
  }
######################## AA network construction removal of sample TAD1,TAD5 TAD8 TAH2 TAH7 TAH9

design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil",]
id_rm= c("TAD1","TAD5","TAD8","TAH2","TAH7","TAH9")
sub_design= sub_design[!rownames(sub_design) %in% id_rm,]


bac_aa <- read.delim(.botu, row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_aa <- read.delim(.fotu, row.names= 1,header=T, sep="\t",stringsAsFactors = F)

id <- colnames(fungi_aa)%in%colnames(bac_aa)

fungi_aa <- fungi_aa[,id]

id <- match(colnames(bac_aa),colnames(fungi_aa))
bac_aa <- bac_aa[,id]
bac_aa <- bac_aa[,colnames(bac_aa)%in% rownames(sub_design)]
fungi_aa <- fungi_aa[,colnames(fungi_aa)%in% rownames(sub_design)]


bac_aa$species <- rep("bac",nrow(bac_aa))
fungi_aa$species <- rep("fungi",nrow(fungi_aa))
merge_bac_fungi <-rbind(bac_aa,fungi_aa) 

###1. Filtering OTUs by occurrence frequency (i.e.,number of samples an OTU is Present 60% of the samples)
Abu <- merge_bac_fungi
Abu$species <- NULL
table <- Abu
table[table>0.0001] <- 1
table.generalist <- Abu[which(rowSums(table)>=5),]  ### 5 out of 14
Abu <- table.generalist

merge_bac_fungi <- merge_bac_fungi[rownames(merge_bac_fungi) %in% rownames(Abu),]



######### Tax #############
bac_tax <- read.delim(.btax, row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_tax <- read.delim(.ftax, row.names= 1,header=F, sep="\t",stringsAsFactors = F)

colnames(fungi_tax)<- colnames(bac_tax)
# fungi_tax$V10 <- NULL
# fungi_tax$V9 <- NULL

colnames(fungi_tax) <- colnames(bac_tax)

tax_merge_bac_fungi <-rbind(fungi_tax,bac_tax) 
tax_merge_bac_fungi <- tax_merge_bac_fungi[-1,]

tax_merge_bac_fungi$ID <- rownames(tax_merge_bac_fungi)
merge_bac_fungi$ID <- rownames(merge_bac_fungi)
merge_abun_tax<- merge(merge_bac_fungi,tax_merge_bac_fungi,by = "ID")
# merge_abun_tax<- merge_abun_tax[-1,]

write.table(merge_bac_fungi,file = paste0(table.dir,.version,"BacAndFunMerDisHealShf.txt"),sep = '\t',row.names = T,quote = F,append = F)



  
######dplyr to data manipulation
# install.packages("psych")
library(psych)

aa_dis_a<- merge_abun_tax[,c(2:15,17,22)] %>%
  mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:15)])) %>%
  group_by(Kingdom,Genus) %>%
  summarize_all(sum) %>%
  arrange(desc(rowSumAbun))

colnames(aa_dis_a)

aa_dis_a<- as.data.frame(aa_dis_a)
rownames(aa_dis_a)<- paste(aa_dis_a$Kingdom,aa_dis_a$Genus,sep = '_')
aa_dis_a$Kingdom <- NULL
aa_dis_a$Genus <- NULL
aa_dis_a <- as.data.frame(t(aa_dis_a))
aa_dis_a$Fungi_Unassigned <- NULL


#### correlation of top 50 genus with microbial load
dis_micload<- as.data.frame(rowSums(aa_dis_a))
colnames(dis_micload) <- "microbial load"

top50_cor_micload_pearson<- corr.test(aa_dis_a[,c(1:50)],dis_micload,use = "complete",method="pearson",adjust = "fdr",alpha = .05)
top50_cor_micload_pearson$p
top50_cor_micload_pearson$r
write.table(top50_cor_micload_pearson$p,file = paste(table.dir,.version,"AAShfTop50CorMicloadPearsonPvalue.txt",sep = ''),sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_pearson$r,file = paste(table.dir,.version,"AAShfTop50CorMicloadPearsonCoeffi.txt",sep = ''),sep = '\t',row.names = T,quote = F)

top50_cor_micload_spearman<- corr.test(aa_dis_a[,c(1:ncol(aa_dis_a))],dis_micload,use = "complete",method="spearman",adjust = "fdr",alpha = .05)
top50_cor_micload_spearman$p
top50_cor_micload_spearman$r
# write.table(top50_cor_micload_spearman$p,file = "../table/190318_AA_shuffle_allgenus_cor_micload_spearman_pvalue.txt",sep = '\t',row.names = T,quote = F)
# write.table(top50_cor_micload_spearman$r,file = "../table/190318_AA_shuffle_allgenus_cor_micload_spearman_coreffi.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$p,file = paste(table.dir,.version,"AAShfAllGenusCorMicloadPearsonPvalue.txt",sep = ''),sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$r,file = paste(table.dir,.version,"AAShfAllGenusCorMicloadPearsonCoeffi.txt",sep = ''),sep = '\t',row.names = T,quote = F)


idd <- grep("Bipolaris",colnames(aa_dis_a)) 
idd1 <- grep("Fusarium",colnames(aa_dis_a)) 

top50_cor_micload_pearson<- corr.test(aa_dis_a[,c(idd)],dis_micload,use = "complete",method="pearson",adjust = "fdr",alpha = .05)
top50_cor_micload_pearson$p
top50_cor_micload_pearson$r
# write.table(top50_cor_micload_pearson$p,file = "../table/190318_AA_shuffle_Bipolaris_cor_micload_pearson_pvalue.txt",sep = '\t',row.names = T,quote = F)
# write.table(top50_cor_micload_pearson$r,file = "../table/190318_AA_shuffle_Bipolaris_cor_micload_pearson_coreffi.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_pearson$p,file = paste(table.dir,.version,"AAShfBipolarisCorMicloadPearsonPvalue.txt",sep = ''),sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_pearson$r,file = paste(table.dir,.version,"AAShfBipolarisCorMicloadPearsonCoeffi.txt",sep = ''),sep = '\t',row.names = T,quote = F)




top50_cor_micload_spearman<- corr.test(aa_dis_a[,c(idd)],dis_micload,use = "complete",method="spearman",adjust = "fdr",alpha = .05)
top50_cor_micload_spearman$p
top50_cor_micload_spearman$r
# write.table(top50_cor_micload_pearson$p,file = "../table/190318_AA_shuffle_Bipolaris_cor_micload_spearman_pvalue.txt",sep = '\t',row.names = T,quote = F)
# write.table(top50_cor_micload_pearson$r,file = "../table/190318_AA_shuffle_Bipolaris_cor_micload_spearman_coreffi.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$p,file = paste(table.dir,.version,"AAShfBipolarisCorMicloadSpearmanPvalue.txt",sep = ''),sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$r,file = paste(table.dir,.version,"AAShfBipolarisCorMicloadSpearmanCoeffi.txt",sep = ''),sep = '\t',row.names = T,quote = F)




top50_cor_micload_pearson<- corr.test(aa_dis_a[,c(idd1)],dis_micload,use = "complete",method="pearson",adjust = "fdr",alpha = .05)
top50_cor_micload_pearson$p
top50_cor_micload_pearson$r
# write.table(top50_cor_micload_spearman$p,file = "../table/190318_AA_shuffle_Fusarium_cor_micload_pearson_pvalue.txt",sep = '\t',row.names = T,quote = F)
# write.table(top50_cor_micload_spearman$r,file = "../table/190318_AA_shuffle_Fusarium_cor_micload_pearson_coreffi.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$p,file = paste(table.dir,.version,"AAShfFusariumCorMicloadPearsonPvalue.txt",sep = ''),sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$r,file = paste(table.dir,.version,"AAShfFusariumsCorMicloadPearsonCoeffi.txt",sep = ''),sep = '\t',row.names = T,quote = F)



top50_cor_micload_spearman<- corr.test(aa_dis_a[,c(idd1)],dis_micload,use = "complete",method="spearman",adjust = "fdr",alpha = .05)
top50_cor_micload_spearman$p
top50_cor_micload_spearman$r
# write.table(top50_cor_micload_spearman$p,file = "../table/190318_AA_shuffle_Fusarium_cor_micload_spearman_pvalue.txt",sep = '\t',row.names = T,quote = F)
# write.table(top50_cor_micload_spearman$r,file = "../table/190318_AA_shuffle_Fusarium_cor_micload_spearman_coreffi.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$p,file = paste(table.dir,.version,"AAShfFusariumCorMicloadSpearmanPvalue.txt",sep = ''),sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$r,file = paste(table.dir,.version,"AAShfFusariumCorMicloadSpearmanCoeffi.txt",sep = ''),sep = '\t',row.names = T,quote = F)



#### note the most 210 threshold
# aa_dis_occor_kendall = corr.test(aa_dis_a[1:222],use="pairwise",method="kendall",adjust="fdr",alpha=.05)
# aa_dis_occor_kendall$r
# aa_occor_kendall$t
# aa_occor_kendall$n
if (.subflag==TRUE){
  col_length=threashold
}else{
  col_length <- ncol(aa_dis_a)
}

aa_dis_occor_spearman <- corr.test(aa_dis_a[1:col_length],use="pairwise",method="spearman",adjust="fdr",alpha=.05)
aa_dis_occor_spearman$r

# ## to filter the top 30 genus according to AA profiling
# id <- rownames(aa_dis_occor_spearman$r)%in% rownames(ra_dis_sper_mat_top30)
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[id,id]
rownames(aa_dis_occor_spearman$r)

idd <- grep("Fungi_Bipolaris",rownames(aa_dis_occor_spearman$r))

idd1 <- grep("Fungi_Fusarium",rownames(aa_dis_occor_spearman$r))
# idd2 <- grep("g__Mycosphaerella",rownames(aa_dis_occor_spearman$r))
# idd3 <- grep("Kaistia",rownames(aa_dis_occor_spearman$r))
# idd4 <- grep("g__Didymella",rownames(aa_dis_occor_spearman$r))
# idd5 <- grep("Kinneretia",rownames(aa_dis_occor_spearman$r))
# idd6 <- grep("Methylotenera",rownames(aa_dis_occor_spearman$r))
# idd7 <- grep("Aquimonas",rownames(aa_dis_occor_spearman$r))
# idd8 <- grep("Reyranella",rownames(aa_dis_occor_spearman$r))
# idd9 <- grep("Cryobacterium",rownames(aa_dis_occor_spearman$r))
# idd10 <- grep("Ferruginibacter",rownames(aa_dis_occor_spearman$r))
# idd11 <- grep("Desulfuromonas",rownames(aa_dis_occor_spearman$r))
# idd12<- grep("Mucilaginibacter",rownames(aa_dis_occor_spearman$r))

# ind<- c(idd2)
# # idd3 <- grep("g__Mycosphaerella",rownames(aa_dis_occor_spearman$r))
# 
# 


# ids <- match(row.names(ids_m),rownames(aa_dis_occor_spearman$r))
# # aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:50,idd),c(1:50,idd)]
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(ids,idd2),c(ids,idd2)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(ids,idd2),c(ids,idd2)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:50,idd),c(1:50,idd)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[id,id]
# ids <- match(order(row.names(aa_dis_sper_mat_top30)),row.names(aa_dis_sper_mat_top30))

# top50_plus_other <- aa_dis_a[,c(1:50,idd)]
# top50_plus_other <- aa_dis_a[,c(1:50,idd,idd1,idd2)]
# 


# aa_dis_sper_mat_top30<- aa_dis_sper_mat_top30[ids,ids]
# aa_dis_sper_p_top30 <- aa_dis_sper_p_top30[ids,ids]


a_dis=aa_dis_occor_spearman$r
a_dis_p=aa_dis_occor_spearman$p
rr1 <- a_dis_p
rr1[lower.tri(rr1)] <- t(rr1)[lower.tri(rr1)]
a_dis_p<- rr1

#### p-value 0.05
# ids_p = as.data.frame(a_dis_p[a_dis_p[,colnames(a_dis_p) %in% "Fungi_Bipolaris"] <0.05,colnames(a_dis_p) %in% "Fungi_Bipolaris"])
# colnames(ids_p)= "adjusted-P"
# 
# 
# ids_m=a_dis[rownames(a_dis)%in% rownames(ids_p),colnames(a_dis)%in% rownames(ids_p)]
# ids_p <- a_dis_p[rownames(a_dis_p) %in% rownames(ids_m),colnames(a_dis_p)%in% rownames(ids_m)]

# idd2 <- grep("Mucilaginibacter",rownames(aa_dis_occor_spearman$r))
# idd3 <- grep("g__Mycosphaerella",rownames(aa_dis_occor_spearman$r))


# ids <- match(row.names(ids_m),rownames(aa_dis_occor_spearman$r))[c(1:55)]
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:55,idd),c(1:55,idd)]
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:55,idd1,idd),c(1:55,idd1,idd)]
# # aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:55,idd),c(1:55,idd)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:55,idd1,idd),c(1:55,idd1,idd)]
# # aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(ids,idd),c(ids,idd)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(ids,idd),c(ids,idd)]




# write.table(ids_m,file = "../table/190318_AA_shuffle_AAinteractionWIthBipolarisFusarium_top55_p0.05.txt",sep = '\t',row.names = T,quote = F)
# write.table(ids_p,file = "../table/190318_AA_shuffle_AAinteractionWIthBipolarisFusarium_top55adjustedP_0.05.txt",sep = '\t',row.names = T,quote = F)
# write.table(top50_plus_other,file = "../table/190318_AA_shuffle_bipolaris-associated_abundance.txt",sep = '\t',row.names = T,quote = F)
# write.table(aa_dis_sper_mat_top30,file = "../table/190318_AA_shuffle_bipolaris-associated_cor.txt",sep = '\t',row.names = T,quote = F)




# write.table(aa_dis_sper_mat_top30,file = "../table/190318_AA_shuffle_Bipolaris-associated_corp0.05.txt",sep = '\t',row.names = T,quote = F)
# write.table(aa_dis_sper_p_top30,file = "../table/190318_AA_shuffle_Bipolaris-associated_adjustedPp0.05.txt",sep = '\t',row.names = T,quote = F)



# gsub(pattern =' |\t',',',colnames(aa_dis_sper_p_top30),perl = T)
# name_vector <- c( "Promicromonospora" , "Acholeplasma"   ,    "Phyllobacterium" ,   "Rhizobium"  ,        "Variovorax"   ,      "Lentzea"     ,       "Flavobacterium" ,    "Sphingomonas",      
#                    "Pseudoxanthomonas" , "Nocardioides"  ,     "Streptomyces"    ,   "Actinoplanes" ,      "Massilia"     ,      "Herbaspirillum"  ,   "Devosia"   ,         "Scytalidium",       
#                    "Chryseobacterium" ,  "Acremonium" ,        "Mesorhizobium"  ,    "Kribbella"    ,      "Mycobacterium"   ,   "Stenotrophomonas" ,  "Microbacterium" ,    "Steroidobacter"  ,  
#                    "Aeromicrobium"    ,  "Sorangium"     ,     "Caulobacter"    ,    "Sphingopyxis"    ,   "Arthrobacter" ,      "Agromyces",       "Umezawaea"  ,        "Shinella" ,         
#                    "Rhodanobacter"  ,    "Solirubrobacter" ,   "Bradyrhizobium"  ,   "Burkholderia"   ,    "Lysobacter"  ,       "Bosea"       ,       "Cellulomonas"   ,  "Bacillus"  ,     
#                    "Glycomyces"    ,     "Podospora"      ,    "Altererythrobacter", "Amycolatopsis"  ,    "Fusarium"      ,"Phycicoccus"  ,  "Sphingobium"   ,     "Gaiella" ,          
#                    "Blastococcus"     ,  "Paenibacillus"   ,   "Rhodopseudomonas"  , "Conexibacter"   ,    "Pseudonocardia"   ,  "Nitrobacter" , "Piscinibacter", "Bipolaris"  )


if(.flag == TRUE){
  
  aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:55,idd),c(1:55,idd)]
  # aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:55,idd),c(1:55,idd)]
  aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:55,idd),c(1:55,idd)]
  colnames(aa_dis_sper_mat_top30)

name_vector <- c("Bacteria_Promicromonospora","Bacteria_Pseudomonas","Bacteria_Phyllobacterium","Bacteria_Variovorax",
                 "Bacteria_Flavobacterium","Bacteria_Rhizobium","Bacteria_Lentzea","Bacteria_Pseudoxanthomonas",
                 "Bacteria_Sphingomonas","Bacteria_Nocardioides","Bacteria_Actinoplanes","Bacteria_Streptomyces",
                 "Bacteria_Chryseobacterium","Bacteria_Herbaspirillum","Bacteria_Massilia",
                 "Bacteria_Devosia","Bacteria_Kribbella","Bacteria_Mycobacterium",
                 "Bacteria_Microbacterium","Bacteria_Mesorhizobium","Bacteria_Stenotrophomonas","Bacteria_Sphingopyxis",
                 "Bacteria_Aeromicrobium","Bacteria_Caulobacter","Bacteria_Steroidobacter","Bacteria_Arthrobacter",
                 "Bacteria_Umezawaea","Bacteria_Chondromyces","Bacteria_Sorangium","Bacteria_Agromyces",
                 "Bacteria_Rhodanobacter","Bacteria_Galbitalea","Bacteria_Pantoea","Bacteria_Shinella",
                 "Bacteria_Bosea","Bacteria_Solirubrobacter","Bacteria_Cellulomonas",
                 "Bacteria_Lysobacter","Bacteria_Bradyrhizobium","Bacteria_Rudaibacter","Bacteria_Glycomyces",
                 "Bacteria_Burkholderia","Bacteria_Altererythrobacter","Bacteria_Bacillus","Bacteria_Amycolatopsis",
                 "Bacteria_Rhodopseudomonas","Bacteria_Phycicoccus","Bacteria_Sphingobium","Bacteria_Gaiella",
                 "Bacteria_Blastococcus","Bacteria_Ideonella",
                 "Fungi_Acremonium","Fungi_Scytalidium","Fungi_Gibberella","Fungi_Holtermanniella" ,
                 
                 "Fungi_Bipolaris" )

ids <- match(name_vector,row.names(aa_dis_sper_mat_top30))
# 
aa_dis_sper_mat_top30<- aa_dis_sper_mat_top30[ids,ids]
aa_dis_sper_p_top30 <- aa_dis_sper_p_top30[ids,ids]

write.table(aa_dis_sper_p_top30,file =paste0(table.dir,.version,"OrdShfAASperMatTop55FuBiUppTriQv0.05.txt"),sep = '\t',row.names = T,quote = F)
write.table(aa_dis_sper_mat_top30,file =paste0(table.dir,.version,"OrdShfAASperMatTop55FuBiUppTriCoeffi.txt"),sep = '\t',row.names = T,quote = F)

# pdf(file ="../plot/rootDisease/190318_namenotorder_shuffle_AA_sper_mat_top50_order_upper_p_value_0.05.pdf" )
# corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .05,
#          insig = "blank", addrect = 2,type = "upper",diag = F,
#          tl.col = "black",tl.cex = 0.75,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")
# 
# dev.off()
# 
# ordered_name_vector <- c("Bacteria_Promicromonospora"              ,        "Bacteria_Acholeplasma"                 ,          "Bacteria_Phyllobacterium"            ,           
#                          "Bacteria_Rhizobium"                      ,        "Bacteria_Variovorax"                 ,            "Bacteria_Lentzea"                ,               
#                          "Bacteria_Flavobacterium"              ,           "Bacteria_Sphingomonas"                  ,         "Bacteria_Nocardioides"        ,                  
#                          "Bacteria_Streptomyces"                 ,          "Bacteria_Pseudoxanthomonas"      ,                "Bacteria_Devosia"           ,                    
#                                         "Bacteria_Kribbella"         ,                    
#                          "Bacteria_Massilia"                     ,          "Bacteria_Mesorhizobium"         ,                 "Bacteria_Chryseobacterium"   ,                   
#                          "Bacteria_Herbaspirillum"                 ,        "Bacteria_Lysobacter"                        ,     "Bacteria_Mycobacterium"       ,                  
#                          "Bacteria_Microbacterium"                      ,   "Bacteria_Stenotrophomonas"      ,                 "Bacteria_Steroidobacter"     ,                   
#                          "Bacteria_Aeromicrobium"                  ,        "Bacteria_Caulobacter"                        ,    "Bacteria_Sphingopyxis"        ,                  
#                          "Bacteria_Shinella"                        ,       "Bacteria_Agromyces"                    ,          "Bacteria_Umezawaea"         ,                    
#                          "Bacteria_Bradyrhizobium"                 ,        "Bacteria_Solirubrobacter"             ,           "Bacteria_Bosea"            ,                     
#                          "Bacteria_Burkholderia"                    ,       "Bacteria_Cellulomonas"              ,             "Bacteria_Polaromonas"         ,                  
#                          "Bacteria_Bacillus"                           ,    "Bacteria_Glycomyces"              ,                               
#                          "Bacteria_Rhodanobacter"                   ,                  "Bacteria_Xanthomonas"         ,                  
#                          "Bacteria_Altererythrobacter"                 ,    "Bacteria_Amycolatopsis"        ,                  "Bacteria_Phycicoccus"          ,                 
#                          "Bacteria_Rhodopseudomonas"                 ,      "Bacteria_Phenylobacterium"           ,            "Bacteria_Saccharibacteria_genera_incertae_sedis",
#                          "Bacteria_Paenibacillus"                    ,      "Bacteria_Blastococcus"               ,            "Bacteria_Pseudonocardia"                ,        
#                          "Bacteria_Sphingobium"                        ,    "Bacteria_Microvirga"                       ,      "Bacteria_Dyella"                         ,      
#                          "Fungi_Scytalidium"                      ,         "Fungi_Acremonium"               ,   "Fungi_Podospora"               ,  
#                          "Fungi_Fusarium"                     , "Fungi_Holtermanniella" , "Fungi_Bipolaris" )
# 
# ids <- match(ordered_name_vector,row.names(aa_dis_sper_mat_top30))
# # 
# aa_dis_sper_mat_top30<- aa_dis_sper_mat_top30[ids,ids]
# aa_dis_sper_p_top30 <- aa_dis_sper_p_top30[ids,ids]


pdf(file = paste0(figures.dir,.version,"OrdShfAASperMatTop55FuBiUppTriP0.05.pdf") )
corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .05,
         insig = "blank", addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

dev.off()

# pdf(file ="../plot/rootDisease/190318_usearch_ordername_shuffle_AA_sper_mat_top50_order_upper_p_value_0.01.pdf" )
pdf(file = paste0(figures.dir,.version,"OrdShfAASperMatTop55FuBiUppTriP0.01.pdf") )
corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

dev.off()


# pdf(file ="../plot/rootDisease/190318_usearch_ordername_shuffle_AA_sper_mat_top50_order_upper_p_value_0.0001.pdf" )
pdf(file = paste0(figures.dir,.version,"OrdShfAASperMatTop55FuBiUppTriP0.001.pdf") )
corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .001,
         insig = "blank", addrect = 2,type = "upper",diag = F,
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

dev.off()

}else
{
  aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:55,idd),c(1:55,idd)]
  # aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:55,idd),c(1:55,idd)]
  aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:55,idd),c(1:55,idd)]
  
  fh <- paste0(table.dir,.version,"AssociationGenusNameSort.txt")
  sink(fh)
  colnames(aa_dis_sper_mat_top30)
  # sort(colnames(aa_dis_sper_mat_top30))
  sink()
  
  system(paste0("sed -i 's/\\s\\+/,/g' ",fh))
  system(paste0("sed -i 's/^,//g' ",fh))
  system(paste0("sed -i 's/^\\[[1-9]\\+\\],//g' ",fh))
  
  
  name_vector <- c("Bacteria_Promicromonospora","Bacteria_Pseudomonas","Bacteria_Phyllobacterium","Bacteria_Rhizobium","Bacteria_Flavobacterium",
                   "Bacteria_Variovorax","Bacteria_Lentzea","Bacteria_Sphingomonas","Bacteria_Nocardioides","Bacteria_Streptomyces",
                   "Bacteria_Chryseobacterium","Bacteria_Actinoplanes","Bacteria_Pseudoxanthomonas",
                   "Bacteria_Devosia","Bacteria_Kribbella","Bacteria_Mycobacterium","Bacteria_Lysobacter","Bacteria_Microbacterium",
                   "Bacteria_Massilia","Bacteria_Herbaspirillum","Bacteria_Mesorhizobium","Bacteria_Stenotrophomonas","Bacteria_Aeromicrobium",
                   "Bacteria_Sphingopyxis","Bacteria_Steroidobacter","Bacteria_Caulobacter","Bacteria_Arthrobacter","Bacteria_Umezawaea",
                   "Bacteria_Agromyces","Bacteria_Shinella","Bacteria_Chondromyces","Bacteria_Sorangium","Bacteria_Pantoea",
                   "Bacteria_Solirubrobacter","Bacteria_Bradyrhizobium","Bacteria_Bosea","Bacteria_Couchioplanes","Bacteria_Cellulomonas",
                   "Bacteria_Polaromonas","Bacteria_Galbitalea","Bacteria_Burkholderia","Bacteria_Rhodanobacter",
                   "Bacteria_Glycomyces","Bacteria_Xanthomonas","Bacteria_Bacillus",
                   "Bacteria_Rudaibacter","Bacteria_Rhodopseudomonas","Bacteria_Altererythrobacter","Bacteria_Pseudonocardia","Bacteria_Phycicoccus",
                   "Fungi_Acremonium","Fungi_Scytalidium","Fungi_Podospora","Fungi_Myrmecridium","Fungi_Fusarium",
                   "Fungi_Bipolaris" )

  ids <- match(name_vector,row.names(aa_dis_sper_mat_top30))
  # 
  aa_dis_sper_mat_top30<- aa_dis_sper_mat_top30[ids,ids]
  aa_dis_sper_p_top30 <- aa_dis_sper_p_top30[ids,ids]
  
  write.table(aa_dis_sper_p_top30,file =paste0(table.dir,.version,"OrdShfAASperMatTop55FuBiUppTriQv0.05.txt"),sep = '\t',row.names = T,quote = F)
  write.table(aa_dis_sper_mat_top30,file =paste0(table.dir,.version,"OrdShfAASperMatTop55FuBiUppTriCoeffi.txt"),sep = '\t',row.names = T,quote = F)
  
  
  pdf(file = paste0(figures.dir,.version,"OrdShfAASperMatTop55FuBiUppTriP0.05.pdf") )
  corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .05,
           insig = "blank", addrect = 2,type = "upper",diag = F,
           tl.col = "black",tl.cex = 0.75,
           cl.pos = "b",cl.length = 11,cl.align.text = "l",
           cl.lim = c(-1,1),
           title = "Genus-genus spearman correlation")
  
  dev.off()
  
  # pdf(file ="../plot/rootDisease/190318_usearch_ordername_shuffle_AA_sper_mat_top50_order_upper_p_value_0.01.pdf" )
  pdf(file = paste0(figures.dir,.version,"OrdShfAASperMatTop55FuBiUppTriP0.01.pdf") )
  corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .01,
           insig = "blank", addrect = 2,type = "upper",diag = F,
           tl.col = "black",tl.cex = 0.75,
           cl.pos = "b",cl.length = 11,cl.align.text = "l",
           cl.lim = c(-1,1),
           title = "Genus-genus spearman correlation")
  
  dev.off()
  
  
  # pdf(file ="../plot/rootDisease/190318_usearch_ordername_shuffle_AA_sper_mat_top50_order_upper_p_value_0.0001.pdf" )
  pdf(file = paste0(figures.dir,.version,"OrdShfAASperMatTop55FuBiUppTriP0.001.pdf") )
  corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .001,
           insig = "blank", addrect = 2,type = "upper",diag = F,
           tl.col = "black",tl.cex = 0.75,
           cl.pos = "b",cl.length = 11,cl.align.text = "l",
           cl.lim = c(-1,1),
           title = "Genus-genus spearman correlation")
  
  dev.off()
  
}

# pdf(file ="../plot/rootDisease/190318_shuffle_AA_sper_mat_bipolarisAssociated_order_upper.pdf" )
# corrplot(ids_m,p.mat = ids_p, sig.level = .01,
#          insig = "blank", addrect = 2,type = "upper",diag = F,
#          tl.col = "black",tl.cex = 0.75,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")
# 
# dev.off()


#### global negative correlation
# corrplot(aa_dis_occor_spearman$r,p.mat = aa_dis_occor_spearman$p, sig.level = .01,
#          insig = "blank",  addrect = 2,type = "upper",diag = F,
#          tl.col = "black",tl.cex = 0.65,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")




# #### p-value 0.001
# ids_p = as.data.frame(a_dis_p[a_dis_p[,colnames(a_dis_p) %in% "Fungi_Bipolaris"] <0.001,colnames(a_dis_p) %in% "Fungi_Bipolaris"])
# colnames(ids_p)= "adjusted-P"
# 
# 
# ids_m=a_dis[rownames(a_dis)%in% rownames(ids_p),colnames(a_dis)%in% rownames(ids_p)]
# ids_p <- a_dis_p[rownames(a_dis_p) %in% rownames(ids_m),colnames(a_dis_p)%in% rownames(ids_m)]
# 
# # idd2 <- grep("Mucilaginibacter",rownames(aa_dis_occor_spearman$r))
# # idd3 <- grep("g__Mycosphaerella",rownames(aa_dis_occor_spearman$r))
# 
# 
# ids <- match(row.names(ids_m),rownames(aa_dis_occor_spearman$r))[c(1:26)]
# # aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:50,idd),c(1:50,idd)]
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(ids,idd),c(ids,idd)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(ids,idd),c(ids,idd)]
# 
# 
# write.table(ids_m,file = "../table/190318_AA_shuffle_AAinteractionWIthBipolaris_dis_p0.001.txt",sep = '\t',row.names = T,quote = F)
# write.table(ids_p,file = "../table/190318_AA_shuffle_AAinteractionWIthBipolaris_adjustedP_0.001.txt",sep = '\t',row.names = T,quote = F)
# # write.table(top50_plus_other,file = "../table/190318_AA_shuffle_bipolaris-associated_abundance.txt",sep = '\t',row.names = T,quote = F)
# # write.table(aa_dis_sper_mat_top30,file = "../table/190318_AA_shuffle_bipolaris-associated_cor.txt",sep = '\t',row.names = T,quote = F)
# 
# 
# 
# 
# write.table(aa_dis_sper_mat_top30,file = "../table/190318_AA_shuffle_bipolaris-associated_corp0.001.txt",sep = '\t',row.names = T,quote = F)
# write.table(aa_dis_sper_p_top30,file = "../table/190318_AA_shuffle_bipolaris-associated_adjustedPp0.001.txt",sep = '\t',row.names = T,quote = F)
# 
# # aa_heal_occor_kendall$p[ids,ids]
# # 
# # pheatmap(aa_dis_sper_mat_top30)
# colnames(ids_m)
# name_vector <- c( "Promicromonospora" ,  "Flavobacterium"   ,   "Rhizobium"      ,     "Sphingomonas"   ,     "Chryseobacterium" ,      
#                   "Mycobacterium"    ,   "Microbacterium"   ,   "Lysobacter"     ,     
#                   "Burkholderia"     ,   "Bacillus"       ,     "Amycolatopsis" ,      "Sphingobium"   ,      "Luteimonas"    ,      "Achromobacter"   ,   
#                   "Microvirga"       ,   "Novosphingobium"  ,   "Kaistia"     ,            "Chitinophaga"   ,     
#                   "Hyphomicrobium"   ,   "Iamia"     ,          "Pseudorhodoferax"  ,  "Dissulfuribacter" ,   "Nocardia"   ,         "Methylophilus"    ,  
#                   "Reyranella"      ,    "Afipia"       ,       "Rhizocola"     ,     "Nonomuraea"      ,         "Haloactinopolyspora",
#                   "Panacagrimonas"   ,   "Povalibacter"   ,     "Pedomicrobium"  ,     "Variibacter"   ,      "Pedobacter"    ,      "Chelativorans"  ,    
#                   "Chelatococcus"     ,"Mucilaginibacter", "g__Acremonium"  ,   "g__Ophiosphaerella" , "g__Gibberella"    , "g__Pyrenochaetopsis"  ,"g__Pseudogymnoascus", 
#                   "g__Didymella"   ,"g__Bipolaris" 
# )
# 
# 
# # name_vector <- c("Promicromonospora",   "Flavobacterium" ,     "Variovorax"    ,      "Rhizobium"     ,      "Pseudomonas"      ,   "Phyllobacterium" ,   
# #                  "Lentzea"       ,      "Pseudoxanthomonas" ,  "Actinoplanes"  ,      "Sphingomonas",        "Nocardioides" ,       "Streptomyces",       
# #                  "Chryseobacterium",    "Massilia"      ,      "Herbaspirillum"  ,    "Devosia"      ,     "Kribbella"   ,       
# #                  "Mesorhizobium"    ,   "Mycobacterium"    , "Caulobacter"  ,      
# #                  "Umezawaea"         ,       "Microbacterium"    , "Aeromicrobium"   ,    "Stenotrophomonas"  , 
# #                  "Galbitalea"    ,      "Arthrobacter"   ,     "Steroidobacter" ,     "Chondromyces",        "Sphingopyxis" ,       "Shinella"        ,   
# #                  "Pantoea"        ,     "Agromyces"   ,        "Sorangium"     ,      "Rhodanobacter"  ,      "Bradyrhizobium"  ,   
# #                  "Solirubrobacter"   ,  "Rudaibacter"     ,    "Bosea"          ,  "Lysobacter"  ,
# #                  "Kineosporia"     ,    "Cellulomonas"      ,"Reyranella" ,"Mucilaginibacter", "g__Acremonium"  ,"g__Lachnella"   ,     "g__Scytalidium"   ,   "g__Magnaporthiopsis", 
# #                  "g__Podospora" ,  "g__Ophiosphaerella",      "g__Gibberella",        "g__Cladosporium",        "g__Pyrenochaetopsis",  "g__Didymella", "g__Bipolaris" 
# #                  
# # )
# 
# 
# ids <- match(name_vector,row.names(aa_dis_sper_mat_top30))
# # 
# aa_dis_sper_mat_top30<- aa_dis_sper_mat_top30[ids,ids]
# aa_dis_sper_p_top30 <- aa_dis_sper_p_top30[ids,ids]
# 
# ids1 <- match(name_vector,colnames(aa_dis_a))
# aa_44_genus <- aa_dis_a[,ids1]
# write.table(aa_44_genus,file = "../table/190318_AA_shuffle_disease-associated_44Genus.txt",sep = '\t',row.names = T,quote = F)
# 
# # corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .001
# #          insig = "blank", addrect = 2,type = "upper",
# #          tl.col = "black",tl.cex = 0.65,
# #          cl.pos = "b",cl.length = 11,cl.align.text = "l",
# #          cl.lim = c(-1,1),
# #          title = "Genus-genus spearman correlation")
# 
# pdf(file ="../plot/rootDisease/190318_shuffle_AA_sper_mat_top50_order_upper_1030_p0.001.pdf" )
# corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .001,
#          insig = "blank", addrect = 2,type = "upper",diag = F,
#          tl.col = "black",tl.cex = 0.75,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")
# 
# dev.off()
# 
# #          title = "Genus-genus spearman correlation")
# 
# pdf(file ="../plot/rootDisease/190318_shuffle_AA_sper_mat_top50_order_upper_0.01.pdf" )
# corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .01,
#          insig = "blank", addrect = 2,type = "upper",diag = F,
#          tl.col = "black",tl.cex = 0.75,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# pdf(file ="../plot/rootDisease/190318_shuffle_AA_sper_mat_top50_order_upper_p_value_0.001.pdf" )
# corrplot(aa_dis_sper_mat_top30,p.mat = aa_dis_sper_p_top30, sig.level = .001,
#          insig = "blank", addrect = 2,type = "upper",diag = F,
#          tl.col = "black",tl.cex = 0.75,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")
# 
# dev.off()



#### correlation of disease-associated genus with microbial load
dis_micload<- as.data.frame(rowSums(aa_dis_a))
colnames(dis_micload) <- "microbial load"


top50_cor_micload_pearson<- corr.test(aa_dis_a[,ids],dis_micload,use = "complete",method="pearson",adjust = "fdr",alpha = .05)
top50_cor_micload_pearson$p
top50_cor_micload_pearson$r
write.table(top50_cor_micload_pearson$p,file = paste0(table.dir,.version,"AATop55DisAssoCorMicloadPearsonPv.txt"),sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_pearson$r,file = paste0(table.dir,.version,"AATop55DisAssoCorMicloadPearsonCoeffi.txt"),sep = '\t',row.names = T,quote = F)
# write.table(top50_cor_micload_pearson$r,file = "../table/190318_AA_shuffle_disease-associated_cor_micload_pearson_coreffi.txt",sep = '\t',row.names = T,quote = F)

top50_cor_micload_spearman<- corr.test(aa_dis_a[,ids],dis_micload,use = "complete",method="spearman",adjust = "fdr",alpha = .05)
top50_cor_micload_spearman$p
top50_cor_micload_spearman$r
# write.table(top50_cor_micload_spearman$p,file = "../table/190318_AA_shuffle_disease-associated_cor_micload_spearman_pvalue.txt",sep = '\t',row.names = T,quote = F)
# write.table(top50_cor_micload_spearman$r,file = "../table/190318_AA_shuffle_disease-associated_cor_micload_spearman_coreffi.txt",sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$p,file = paste0(table.dir,.version,"AATop55DisAssoCorMicloadSpearmanPv.txt"),sep = '\t',row.names = T,quote = F)
write.table(top50_cor_micload_spearman$r,file = paste0(table.dir,.version,"AATop55DisAssoCorMicloadSpearmanCoeffi.txt"),sep = '\t',row.names = T,quote = F)






### RA network construction
#####disease 

### OTU RA analysis
.botu <- paste0("~/xiaoxuan/180724/CoNet/",.meth,"/RA/otutab_norm_RA.txt")
.fotu <- paste0("~/xiaoxuan/180724/CoNet/",.meth,"/fRA/otutab_norm_RA.txt")


design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 

sub_design <- design[design$Spikein %in% "spike"& !design$Other %in% "Bulksoil",]
id_rm= c("TAD1","TAD5","TAD8","TAH2","TAH7","TAH9")
sub_design= sub_design[!rownames(sub_design) %in% id_rm,]

bac_aa <- read.delim(.botu, row.names= 1,header=T, sep="\t",stringsAsFactors = F)
fungi_aa <- read.delim(.fotu, row.names= 1,header=T, sep="\t",stringsAsFactors = F)

id <- colnames(fungi_aa) %in% colnames(bac_aa)
fungi_aa <- fungi_aa[,id]

id <- match(colnames(bac_aa),colnames(fungi_aa))
bac_aa <- bac_aa[,id]
bac_aa <- bac_aa[,colnames(bac_aa)%in% rownames(sub_design)]
fungi_aa <- fungi_aa[,colnames(fungi_aa)%in% rownames(sub_design)]





bac_aa$species <- rep("bac",nrow(bac_aa))
fungi_aa$species <- rep("fungi",nrow(fungi_aa))
merge_bac_fungi <-rbind(bac_aa,fungi_aa) 


# 
# ###1. Filtering OTUs by occurrence frequency (i.e.,number of samples an OTU is Present 60% of the samples)
# Abu <- merge_bac_fungi
# Abu$species <- NULL
# table <- Abu
# table[table>0] <- 1
# table.generalist <- Abu[which(rowSums(table)>=6),]  ### 6 out of 10
# Abu <- table.generalist
# 
# merge_bac_fungi <- merge_bac_fungi[rownames(merge_bac_fungi) %in% rownames(Abu),]

# bac_tax <- read.delim("/mnt/bai/xiaoning/xiaoxuan/180724/TARDB28/unoise/result/taxonomy_8_match.txt", row.names= 1,header=T, sep="\t",stringsAsFactors = F)
# fungi_tax <- read.delim("/mnt/bai/qinyuan/xiaoxuan/TAD/190313/AA/result/taxonomy_8.txt", row.names= 1,header=F, sep="\t",stringsAsFactors = F)
# fungi_tax$V10 <- NULL
# fungi_tax$V9 <- NULL



colnames(fungi_tax) <- colnames(bac_tax)

tax_merge_bac_fungi <-rbind(fungi_tax,bac_tax) 
tax_merge_bac_fungi <- tax_merge_bac_fungi[-1,]
tax_merge_bac_fungi$ID <- rownames(tax_merge_bac_fungi)
merge_bac_fungi$ID <- rownames(merge_bac_fungi)
merge_abun_tax<- merge(merge_bac_fungi,tax_merge_bac_fungi,by = "ID")
# write.table(merge_bac_fungi,file = "../table/190318_usearch_RA_adjusted_bacAndFungiMerge_abundance_disease.txt",sep = '\t',row.names = T,quote = F)
write.table(merge_bac_fungi,file = paste0(table.dir,.version,"RABacAndFunMerDisHealShf.txt"),sep = '\t',row.names = T,quote = F,append = F)

######dplyr to data manipulation
# install.packages("psych")


# ra_dis_a<- merge_abun_tax[,c(2:11,13,18)] %>%
#   mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:11)])) %>%
#   group_by(Kingdom,Genus) %>%
#   summarize_all(sum) %>%
#   arrange(desc(rowSumAbun))
ra_dis_a<- merge_abun_tax[,c(2:15,17,22)] %>%
  mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:15)])) %>%
  group_by(Kingdom,Genus) %>%
  summarize_all(sum) %>%
  arrange(desc(rowSumAbun))
# # ra_dis_a<- merge_abun_tax[,c(2:17,19,24)] %>%
# mutate(rowSumAbun=rowSums(merge_abun_tax[,c(2:17)])) %>%
#   group_by(Kingdom,Genus) %>%
#   summarize_all(sum) %>%
#   arrange(desc(rowSumAbun))
# ra_dis_a$g__unidentified <- NULL

ra_dis_a<- as.data.frame(ra_dis_a)
rownames(ra_dis_a)<- paste(ra_dis_a$Kingdom,ra_dis_a$Genus,sep = '_')
ra_dis_a$Kingdom <- NULL
ra_dis_a$Genus <- NULL
ra_dis_a <- as.data.frame(t(ra_dis_a))
ra_dis_a$Fungi_Unassigned <- NULL
ra_dis_a<- ra_dis_a[,colnames(ra_dis_a) %in% colnames(aa_dis_a)]


# aa_dis_a<- as.data.frame(aa_dis_a)
# rownames(aa_dis_a)<- paste(aa_dis_a$Kingdom,aa_dis_a$Genus,sep = '_')
# aa_dis_a$Kingdom <- NULL
# aa_dis_a$Genus <- NULL
# aa_dis_a <- as.data.frame(t(aa_dis_a))
# aa_dis_a$Fungi_Unassigned <- NULL











#### note the most 300 threshold using Kendall and spearman 
# ra_dis_occor_kendall = corr.test(ra_dis_a[1:300],use="pairwise",method="kendall",adjust="fdr",alpha=.05)
# ra_dis_occor_kendall$r

if (.subflag==TRUE){
  col_length <- threashold
  
}else
{
  col_length <- 350
}

ra_dis_occor_spearman <- corr.test(ra_dis_a[1:col_length],use="pairwise",method="spearman",adjust="fdr",alpha=.05)
ra_dis_occor_spearman$r

# a_dis<- as.data.frame(ra_dis_occor_spearman$r)
# a_dis_p<- as.data.frame(ra_dis_occor_spearman$p)
## to filter the top 30 genus according to AA profiling
id <- rownames(ra_dis_occor_spearman$r)%in% rownames(aa_dis_sper_p_top30)
ra_dis_sper_mat_top30 <- ra_dis_occor_spearman$r[id,id]
# aa_dis_sper_mat_top30 <- aa_dis_occor_spearman$r[c(1:30),c(1:30)]
# aa_dis_sper_p_top30 <- aa_dis_occor_spearman$p[c(1:30),c(1:30)]
ra_dis_sper_p_top30 <- ra_dis_occor_spearman$p[id,id]
ids <- match(row.names(aa_dis_sper_mat_top30),row.names(ra_dis_sper_mat_top30))

ra_dis_sper_mat_top30<- ra_dis_sper_mat_top30[ids,ids]
ra_dis_sper_p_top30 <- ra_dis_sper_p_top30[ids,ids]
top50_plus_other <- ra_dis_a[,colnames(ra_dis_a)%in% rownames(ra_dis_sper_mat_top30)]
# 
write.table(top50_plus_other,file = paste0(table.dir,.version,"RAShfTop55Abun.txt"),sep = '\t',row.names = T,quote = F)

write.table(ra_dis_sper_mat_top30,file = paste0(table.dir,.version,"RAShfTop55AbunCoeffi.txt"),sep = '\t',row.names = T,quote = F)
write.table(ra_dis_sper_p_top30,file =paste0(table.dir,.version,"RAShfTop55AbunPv.txt"),sep = '\t',row.names = T,quote = F)

# 
# a_dis=ra_dis_occor_spearman$r
# a_dis_p=ra_dis_occor_spearman$p
# 
# rr1 <- a_dis_p
# rr1[lower.tri(rr1)] <- t(rr1)[lower.tri(rr1)]
# a_dis_p<- rr1


# 
# ###### p=0.05
# 
# ids_p_RA = as.data.frame(a_dis_p[a_dis_p[,colnames(a_dis_p) %in% "Fungi_Bipolaris"] <0.05,colnames(a_dis_p) %in% "Fungi_Bipolaris"])
# colnames(ids_p_RA)= "adjusted-P"
# 
# 
# ids_m_RA=a_dis[rownames(a_dis)%in% rownames(ids_p_RA),colnames(a_dis)%in% rownames(ids_p_RA)]
# colnames(ids_m_RA)
# 
# 
# write.table(ids_m_RA,file = "../table/190318_RA_shuffleinteractionWIthBipolaris_dis_p0.05.txt",sep = '\t',row.names = T,quote = F)
# write.table(ids_p_RA,file = "../table/190318_RA_shuffleinteractionWIthBipolaris_adjustedP_p0.05.txt",sep = '\t',row.names = T,quote = F)



# pdf(file ="../plot/rootDisease/190318_namenoorderRA_shuffle_sper_mat_top50_upper_pvalue0.05.pdf" )
# 
# corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .05,
#          insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
#          tl.col = "black",tl.cex = 0.75,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")
# 
# dev.off()

# pdf(file ="../plot/rootDisease/190318_usearch_nameorderRA_shuffle_sper_mat_top50_upper_pvalue0.05.pdf" )
pdf(file = paste0(figures.dir,.version,"OrdShfRASperMatTop55FuBiLowerTriP0.05.pdf") )
corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .05,
         insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

dev.off()

# pdf(file ="../plot/rootDisease/190318_usearch_nameorderRA_shuffle_sper_mat_top50_upper_pvalue0.01.pdf" )
pdf(file = paste0(figures.dir,.version,"OrdShfRASperMatTop55FuBiLowerTriP0.01.pdf") )
corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .01,
         insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

dev.off()

# pdf(file ="../plot/rootDisease/190318_usearch_nameorderRA_shuffle_sper_mat_top50_upper_pvalue0.001.pdf" )
pdf(file = paste0(figures.dir,.version,"OrdShfRASperMatTop55FuBiLowerTriP0.001.pdf") )
corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .001,
         insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
         tl.col = "black",tl.cex = 0.75,
         cl.pos = "b",cl.length = 11,cl.align.text = "l",
         cl.lim = c(-1,1),
         title = "Genus-genus spearman correlation")

dev.off()

# 
# unique(colnames(ra_dis_sper_mat_top30))
# unique(colnames(aa_dis_sper_mat_top30))
# 
# 
# ###### p=0.001
# 
# ids_p_RA = as.data.frame(a_dis_p[a_dis_p[,colnames(a_dis_p) %in% "g__Bipolaris"] <0.001,colnames(a_dis_p) %in% "g__Bipolaris"])
# colnames(ids_p_RA)= "adjusted-P"
# 
# 
# ids_m_RA=a_dis[rownames(a_dis)%in% rownames(ids_p_RA),colnames(a_dis)%in% rownames(ids_p_RA)]
# colnames(ids_m_RA)
# 
# 
# write.table(ids_m_RA,file = "../table/190318_RA_shuffleinteractionWIthBipolaris_dis_p0.001.txt",sep = '\t',row.names = T,quote = F)
# write.table(ids_p_RA,file = "../table/190318_RA_shuffleinteractionWIthBipolaris_adjustedP_p0.001.txt",sep = '\t',row.names = T,quote = F)
# 
# 
# 
# # write.table(a_dis[ids,ids],file = "../table/RAinteractionWIthBipolaris_dis.txt",sep = '\t',row.names = T,quote = F)
# # write.table(a_dis_p[ids,ids],file = "../table/RAinteractionWIthBipolaris_adjustedP_dis.txt",sep = '\t',row.names = T,quote = F)
# 
# 
# 
# 
# 
# 
# # ra_dis_sper_mat_top30 <- ra_dis_occor_spearman$r[c(1:40),c(1:40)]
# # ra_dis_sper_p_top30 <- ra_dis_occor_spearman$p[c(1:40),c(1:40)]
# 
# # pheatmap(ra_dis_sper_mat_top30)
# 
# 
# pdf(file ="../plot/rootDisease/190318_RA_shuffle_sper_mat_top50_upper_1030p0.001.pdf" )
# 
# corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .001,
#          insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
#          tl.col = "black",tl.cex = 0.75,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")
# 
# dev.off()
# 
# 
# pdf(file ="../plot/rootDisease/190318_RA_shuffle_sper_mat_top50_upper_pvalue0.01.pdf" )
# 
# corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .01,
#          insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
#          tl.col = "black",tl.cex = 0.75,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")
# 
# dev.off()
# 
# 
# 
# pdf(file ="../plot/rootDisease/190318_RA_shuffle_sper_mat_top50_upper_pvalue0.05.pdf" )
# 
# corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .05,
#          insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
#          tl.col = "black",tl.cex = 0.75,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")
# 
# dev.off()
# 
# 
# 
# 
# pdf(file ="../plot/rootDisease/190318_RA_shuffle_sper_mat_top50_upper_pvalue0.001.pdf" )
# 
# corrplot(ra_dis_sper_mat_top30,p.mat = ra_dis_sper_p_top30, sig.level = .001,
#          insig = "blank", addrect = 2,type = "upper",diag = F,col = cm.colors(100),
#          tl.col = "black",tl.cex = 0.75,
#          cl.pos = "b",cl.length = 11,cl.align.text = "l",
#          cl.lim = c(-1,1),
#          title = "Genus-genus spearman correlation")
# 
# dev.off()
# 
# ids2 <- match(name_vector,colnames(ra_dis_a))
# ra_44_genus <- ra_dis_a[,ids2]
# write.table(ra_44_genus,file = "../table/190318_RA_shuffle_disease-associated_44Genus.txt",sep = '\t',row.names = T,quote = F)
# 


