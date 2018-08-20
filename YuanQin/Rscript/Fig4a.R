setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))
library(dplyr)
library(ggplot2)
library(broom)

bac <- read.csv("data/Fig4a-lm-bac.csv", header = T)
fun <- read.csv("data/Fig4a-lm-fun.csv", header = T)

## Bacteria E05-E00

data_e05_bac <- dplyr::select(bac, OTUID, E05_bac, E00_bac)

x1 <- data_e05_bac$E05_bac
y1 <- data_e05_bac$E00_bac

cor.test(x1,y1,method="pearson")
lm.fit1 <- lm(y1~0+x1)
coef(lm.fit1)

plot(data_e05_bac$E05_bac,data_e05_bac$E00_bac,
     xlab ="OTU reads (%) Unamended root DNA", 
     ylab ="OTU reads (%) spiked root DNA")

abline(lm.fit1)

## Bacteria E04-E00

data_e04_bac <- dplyr::select(bac, OTUID, E04_bac, E00_bac)

x2 <- data_e04_bac$E04_bac
y2 <- data_e04_bac$E00_bac

cor.test(x2,y2,method="pearson")
lm.fit2 <-lm(y2~0+x2)
coef(lm.fit2)

plot(data_e04_bac$E04_bac,data_e04_bac$E00_bac,
     xlab ="OTU reads (%) Unamended root DNA", 
     ylab ="OTU reads (%) spiked root DNA")

abline(lm.fit2)


## Bacteria E03-E00

data_e03_bac <- dplyr::select(bac, OTUID, E03_bac, E00_bac)

x3 <- data_e03_bac$E03_bac
y3 <- data_e03_bac$E00_bac

cor.test(x3,y3,method="pearson")
lm.fit3 <-lm(y3~0+x3)
coef(lm.fit3)

plot(data_e03_bac$E03_bac,data_e03_bac$E00_bac,
     xlab ="OTU reads (%) Unamended root DNA", 
     ylab ="OTU reads (%) spiked root DNA")

abline(lm.fit3)


## Fungi E05-E00

data_e05_fun <- dplyr::select(fun, OTUID, E05_fun, E00_fun)

x4 <- data_e05_fun$E05_fun
y4 <- data_e05_fun$E00_fun

cor.test(x4,y4,method="pearson")
lm.fit4 <- lm(y4~0+x4)
coef(lm.fit4)

plot(data_e05_fun$E05_fun,data_e05_fun$E00_fun,pch=2,
     xlab ="OTU reads (%) Unamended root DNA", 
     ylab ="OTU reads (%) spiked root DNA")

abline(lm.fit4)


## Fungi E04-E00

data_e04_fun <- dplyr::select(fun, OTUID, E04_fun, E00_fun)

x5 <- data_e04_fun$E04_fun
y5 <- data_e04_fun$E00_fun

cor.test(x5,y5,method="pearson")
lm.fit5 <- lm(y5~0+x5)
coef(lm.fit5)

plot(data_e04_fun$E04_fun,data_e04_fun$E00_fun,pch=2,
     xlab ="OTU reads (%) Unamended root DNA", 
     ylab ="OTU reads (%) spiked root DNA")

abline(lm.fit5)


## Fungi E03-E00

data_e03_fun <- dplyr::select(fun, OTUID, E03_fun, E00_fun)

x6 <- data_e03_fun$E03_fun
y6 <- data_e03_fun$E00_fun

cor.test(x6,y6,method="pearson")
lm.fit6 <- lm(y6~0+x6)
coef(lm.fit6)

plot(data_e03_fun$E03_fun,data_e03_fun$E00_fun,pch=2,
     xlab ="OTU reads (%) Unamended root DNA", 
     ylab ="OTU reads (%) spiked root DNA")

abline(lm.fit6)

coef <- c(coef(lm.fit1), coef(lm.fit2), coef(lm.fit3), coef(lm.fit4), coef(lm.fit5), coef(lm.fit6))
id <- c("E05_bac", "E04_bac", "E03_bac", "E05_fun", "E04_fun", "E03_fun")
pearson_value <- c("0.9992054", "0.9991369", "0.9985348", "0.9735125", "0.1569087", "0.992301")
tab <- rbind(id, coef, pearson_value)

write.table(tab,file = "./result/Fig4a_pearson_value.txt",sep = '\t',row.names = T)

