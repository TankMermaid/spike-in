setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))

library(ggplot2)
library(plyr)
library(scales)
library(reshape2)

## All

data <- read.csv("data/Fig5c-heatmap-p.csv", header = T)

## HN 12个门

data_hn <- dplyr::select(data, Bacteria, starts_with("HN"))

all <- gather(data_hn, key = group , value = logFC, `HN_MH_RA`:`HN_WYJ_AA`)

p_hn <- ggplot(all, aes(group, Bacteria)) + 
  geom_tile(aes(fill = logFC),colour = "grey50") + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_y_discrete(limits=c("Gammaproteobacteria","Deltaproteobacteria","Betaproteobacteria",
                            "Alphaproteobacteria","Verrucomicrobia","Spirochaetes","Nitrospirae",
                            "Ignavibacteriae","Firmicutes","Chloroflexi","Bacteroidetes","Actinobacteria",
                            "Acidobacteria"))+
  scale_x_discrete(limits=c("HN_MH_RA","HN_MH_AA","HN_WYJ_RA","HN_WYJ_AA"))
p_hn

ggsave(paste("result/Fig5c-Heatmap-HN.pdf", sep=""), p_hn, width = 3.5, height = 4)


## AH 12个门

data_ah <- select(data, Bacteria, starts_with("AH"))

all2 <- gather(data_ah, key = group , value = logFC, `AH_MH_RA`:`AH_WYJ_AA`)

p_ah <- ggplot(all2, aes(group, Bacteria)) + 
  geom_tile(aes(fill = logFC),colour = "grey50") + 
  scale_fill_gradient2(low = "#6D9EC1",high = "sandybrown", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_y_discrete(limits=c("Gammaproteobacteria","Deltaproteobacteria","Betaproteobacteria",
                            "Alphaproteobacteria","Verrucomicrobia","Spirochaetes","Nitrospirae",
                            "Firmicutes","Chloroflexi","Bacteroidetes","Actinobacteria",
                            "Acidobacteria"))+
  scale_x_discrete(limits=c("AH_MH_RA","AH_MH_AA","AH_WYJ_RA","AH_WYJ_AA"))
p_ah
ggsave(paste("result/Fig5c-Heatmap-AH.pdf", sep=""), p_ah , width = 3.5, height = 4)


