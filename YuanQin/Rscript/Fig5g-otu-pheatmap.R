setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(pheatmap)


data_AHWYJ <- read.csv("data/Fig5g-otu-heatmap-3-AHWYJ.csv", header = T,row.names = 1)
pheatmap(data_AHWYJ, cluster_cols = F, cluster_rows = T, border_color = NA)
pheatmap(data_AHWYJ, cluster_cols = F, cluster_rows = T, border_color = NA,  filename = "result/Fig5g-HM-AHWYJ.pdf", width = 7, height = 9)

data_AHMH63 <- read.csv("data/Fig5g-otu-heatmap-3-AHMH63.csv", header = T,row.names = 1)
pheatmap(data_AHMH63, cluster_cols = F, cluster_rows = T, border_color = NA)
pheatmap(data_AHMH63, cluster_cols = F, cluster_rows = T, border_color = NA, filename = "result/Fig5g-HM-AHMH63.pdf", width = 7, height = 9)



data_HNWYJ <- read.csv("data/Fig5g-otu-heatmap-3-HNWYJ.csv", header = T,row.names = 1)
pheatmap(data_HNWYJ, cluster_cols = F, cluster_rows = T, border_color = NA)
pheatmap(data_HNWYJ, cluster_cols = F, cluster_rows = T, border_color = NA,  filename = "result/Fig5g-HM-HNWYJ.pdf", width = 7, height = 9)



data_HNMH63 <- read.csv("data/Fig5g-otu-heatmap-3-HNMH63.csv", header = T,row.names = 1)
pheatmap(data_HNMH63, cluster_cols = F, cluster_rows = T, border_color = NA)
pheatmap(data_HNMH63, cluster_cols = F, cluster_rows = T, border_color = NA, filename = "result/Fig5g-HM-HNMH63.pdf", width = 7, height = 9)
