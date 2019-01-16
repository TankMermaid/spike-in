setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(ggplot2)
library(dplyr)

count <- read.csv("table/20180816-FIGS12.csv", header = T)
MH <- filter(count, genotype == "MH63")
WYJ <- filter(count, genotype == "WYJ")

## MH63 - RA

MHRA <- filter(MH, method == "RA")

MR1 <- ggplot(MHRA , aes(x= phylum, y = count, fill = site )) + 
  geom_bar(stat = "identity", width=0.7) + 
  labs(x="Phylum", y="Relative abundance (RA)") +
  facet_grid(trent ~ ., scales = "free", space = "free")  + scale_fill_brewer(palette = "Set1")+
  scale_x_discrete(limits=c("Chloroflexi","Firmicutes","Alphaproteobacteria",
                            "Bacteroidetes","Acidobacteria","Actinobacteria",
                            "Gammaproteobacteria","Deltaproteobacteria","Betaproteobacteria"))+
  coord_flip()+
  theme_bw()
MR1
ggsave(paste("figure/20180816-FIGS12-MH-RA-1.pdf", sep=""), MR1, width = 9.5, height = 5.5)


MR2 <- ggplot(MHRA , aes(x= phylum, y = count, fill = site )) + 
  geom_bar(stat = "identity", width=0.7) + 
  labs(x="Phylum", y="Relative abundance (RA)") +
  facet_grid(trent ~ ., scales = "free", space = "free")  + scale_fill_brewer(palette = "Set1")+
  scale_x_discrete(limits=c("Chloroflexi","Firmicutes","Bacteroidetes","Acidobacteria",
                            "Deltaproteobacteria","Gammaproteobacteria","Betaproteobacteria",
                            "Actinobacteria","Alphaproteobacteria"))+
  coord_flip()+
  theme_bw()
MR2
ggsave(paste("figure/20180816-FIGS12-MH-RA-2.pdf", sep=""), MR2, width = 9.5, height = 5.5)


## MH63 - QA

MHQA <- filter(MH, method == "QA")

MQ1 <- ggplot(MHQA , aes(x= phylum, y = count, fill = site )) + 
  geom_bar(stat = "identity", width=0.7) + 
  labs(x="Phylum", y="Quantitative abundance (QA)") +
  facet_grid(trent ~ ., scales = "free", space = "free")  + scale_fill_brewer(palette = "Set1")+
  scale_x_discrete(limits=c("Chloroflexi","Firmicutes","Alphaproteobacteria","Acidobacteria",
                            "Gammaproteobacteria","Actinobacteria","Bacteroidetes","Deltaproteobacteria",
                            "Betaproteobacteria"))+
  coord_flip()+
  theme_bw()
MQ1
ggsave(paste("figure/20180816-FIGS12-MH-QA-1.pdf", sep=""), MQ1, width = 9.5, height = 5.5)


MQ2 <- ggplot(MHQA , aes(x= phylum, y = count, fill = site )) + 
  geom_bar(stat = "identity", width=0.7) + 
  labs(x="Phylum", y="Quantitative abundance (QA)") +
  facet_grid(trent ~ ., scales = "free", space = "free")  + scale_fill_brewer(palette = "Set1")+
  scale_x_discrete(limits=c("Firmicutes","Bacteroidetes","Gammaproteobacteria","Chloroflexi",
                            "Acidobacteria","Deltaproteobacteria","Betaproteobacteria",
                            "Alphaproteobacteria","Actinobacteria"))+
  coord_flip()+
  theme_bw()
MQ2
ggsave(paste("figure/20180816-FIGS12-MH-QA-2.pdf", sep=""), MQ2, width = 9.5, height = 5.5)

## WYJ - RA

WRA <- filter(WYJ, method == "RA")

WR1 <- ggplot(WRA , aes(x= phylum, y = count, fill = site )) + 
  geom_bar(stat = "identity", width=0.7) + 
  labs(x="Phylum", y="Relative abundance (RA)") +
  facet_grid(trent ~ ., scales = "free", space = "free")  + scale_fill_brewer(palette = "Set1")+
  scale_x_discrete(limits=c("Chloroflexi","Bacteroidetes","Gammaproteobacteria","Acidobacteria",
                            "Firmicutes","Alphaproteobacteria","Betaproteobacteria","Actinobacteria",
                            "Deltaproteobacteria"))+
  coord_flip()+
  theme_bw()
WR1
ggsave(paste("figure/20180816-FIGS12-WYJ-RA-1.pdf", sep=""), WR1, width = 9.5, height = 5.5)


WR2 <- ggplot(WRA , aes(x= phylum, y = count, fill = site )) + 
  geom_bar(stat = "identity", width=0.7) + 
  labs(x="Phylum", y="Relative abundance (RA)") +
  facet_grid(trent ~ ., scales = "free", space = "free")  + scale_fill_brewer(palette = "Set1")+
  scale_x_discrete(limits=c("Chloroflexi","Gammaproteobacteria","Firmicutes",
                            "Acidobacteria","Bacteroidetes","Actinobacteria","Alphaproteobacteria",
                            "Betaproteobacteria","Deltaproteobacteria"))+
  coord_flip()+
  theme_bw()
WR2
ggsave(paste("figure/20180816-FIGS12-WYJ-RA-2.pdf", sep=""), WR2, width = 9.5, height = 5.5)


## WYJ - QA

WQA <- filter(WYJ, method == "QA")

WQ1 <- ggplot(WQA , aes(x= phylum, y = count, fill = site )) + 
  geom_bar(stat = "identity", width=0.7) + 
  labs(x="Phylum", y="Quantitative abundance (QA)") +
  facet_grid(trent ~ ., scales = "free", space = "free")  + scale_fill_brewer(palette = "Set1")+
  scale_x_discrete(limits=c("Chloroflexi","Firmicutes","Acidobacteria","Bacteroidetes",
                            "Gammaproteobacteria","Alphaproteobacteria","Actinobacteria",
                            "Betaproteobacteria","Deltaproteobacteria"))+
  coord_flip()+
  theme_bw()
WQ1
ggsave(paste("figure/20180816-FIGS12-WYJ-QA-1.pdf", sep=""), WQ1, width = 9.5, height = 5.5)


WQ2 <- ggplot(WQA , aes(x= phylum, y = count, fill = site )) + 
  geom_bar(stat = "identity", width=0.7) + 
  labs(x="Phylum", y="Quantitative abundance (QA)") +
  facet_grid(trent ~ ., scales = "free", space = "free")  + scale_fill_brewer(palette = "Set1")+
  scale_x_discrete(limits=c("Chloroflexi","Bacteroidetes","Firmicutes","Acidobacteria",
                            "Gammaproteobacteria","Deltaproteobacteria","Betaproteobacteria",
                            "Alphaproteobacteria","Actinobacteria"))+
  coord_flip()+
  theme_bw()
WQ2
ggsave(paste("figure/20180816-FIGS12-WYJ-QA-2.pdf", sep=""), WQ2, width = 9.5, height = 5.5)













