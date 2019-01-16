setwd("/mnt/bai/qinyuan/xiaoxuan/bac_New/")
print(paste("Your working directory is in",getwd()))
library(dplyr)
library(ggplot2)
library(tidyr)

bac <- read.csv("table/20180808-fig4c-bac-all.csv")
fun <- read.csv("table/20180808-fig4c-fun-all.csv")
fun2<- read.csv("table/20180808-fig4c-fun-piont.csv")

## E05-BAC-RA

t1 <- select(bac, Bac, E05_RA_ONE, E05_RA_TWO)

bac05ra <- gather(t1, key = variable , value = count, `E05_RA_ONE`:`E05_RA_TWO`)

p1 = ggplot(bac05ra, aes(x=variable, y = count, fill = Bac )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_fill_manual(values = c("#ff8364", "#fbd685", 
                               "#7fa99b", "#f6f6e9", "#1d97c1"))+
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+
  ylab("Relative Abundance (%)")+ 
  theme_classic()
p1

ggsave(paste("figure/20180808-fig4c-bac-ra-e05.pdf", sep=""), p1, width = 6.5, height = 6)


## E04-BAC-RA

t2 <- select(bac, Bac, E04_RA_ONE, E04_RA_TWO)

bac04ra <- gather(t2, key = variable , value = count, `E04_RA_ONE`:`E04_RA_TWO`)

p2 = ggplot(bac04ra, aes(x=variable, y = count, fill = Bac )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_fill_manual(values = c("#ff8364", "#fbd685", 
                               "#7fa99b", "#f6f6e9", "#1d97c1"))+
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+
  ylab("Relative Abundance (%)")+ 
  theme_classic()
p2

ggsave(paste("figure/20180808-fig4c-bac-ra-e04.pdf", sep=""), p2 , width = 6.5, height = 6)


## E03-BAC-RA

t3 <- select(bac, Bac, E03_RA_ONE, E03_RA_TWO)

bac03ra <- gather(t3, key = variable , value = count, `E03_RA_ONE`:`E03_RA_TWO`)

p3 = ggplot(bac03ra, aes(x=variable, y = count, fill = Bac )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_fill_manual(values = c("#ff8364", "#fbd685", 
                               "#7fa99b", "#f6f6e9", "#1d97c1"))+
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+
  ylab("Relative Abundance (%)")+ 
  theme_classic()
p3

ggsave(paste("figure/20180808-fig4c-bac-ra-e03.pdf", sep=""), p3 , width = 6.5, height = 6)


## E05-BAC-QA

t4 <- select(bac, Bac, E05_QA_ONE, E05_QA_TWO)

bac05qa <- gather(t4, key = variable , value = count, `E05_QA_ONE`:`E05_QA_TWO`)

p4 = ggplot(bac05qa, aes(x=variable, y = count, fill = Bac )) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_manual(values = c("#ff8364", "#fbd685", 
                               "#7fa99b", "#f6f6e9", "#1d97c1"))+
  xlab("Groups")+
  ylab("Quantitative Abundance (relative to host plant)")+ 
  theme_classic()
p4

ggsave(paste("figure/20180808-fig4c-bac-qa-e05.pdf", sep=""), p4, width = 6.5, height = 6)


## E04-BAC-QA

t5 <- select(bac, Bac, E04_QA_ONE, E04_QA_TWO)

bac04qa <- gather(t5, key = variable , value = count, `E04_QA_ONE`:`E04_QA_TWO`)

p5 = ggplot(bac04qa, aes(x=variable, y = count, fill = Bac )) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_manual(values = c("#ff8364", "#fbd685", 
                               "#7fa99b", "#f6f6e9", "#1d97c1"))+
  xlab("Groups")+
  ylab("Quantitative Abundance (relative to host plant)")+ 
  theme_classic()
p5

ggsave(paste("figure/20180808-fig4c-bac-qa-e04.pdf", sep=""), p5, width = 6.5, height = 6)


## E03-BAC-QA

t6 <- select(bac, Bac, E03_QA_ONE, E03_QA_TWO)

bac03qa <- gather(t6, key = variable , value = count, `E03_QA_ONE`:`E03_QA_TWO`)

p6 = ggplot(bac03qa, aes(x=variable, y = count, fill = Bac )) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_manual(values = c("#ff8364", "#fbd685", 
                               "#7fa99b", "#f6f6e9", "#1d97c1"))+
  xlab("Groups")+
  ylab("Quantitative Abundance (relative to host plant)")+ 
  theme_classic()
p6

ggsave(paste("figure/20180808-fig4c-bac-qa-e03.pdf", sep=""), p6, width = 6.5, height = 6)



## E05-FUN-RA

f1 <- select(fun, FUN, E05_RA_ONE, E05_RA_TWO)

fun05ra <- gather(f1, key = variable , value = count, `E05_RA_ONE`:`E05_RA_TWO`)

fp1 = ggplot(fun05ra, aes(x=variable, y = count, fill = FUN )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_fill_manual(values = c("#8ea6b4", "#f8d5f0","#89a4c7"))+
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+
  ylab("Relative Abundance (%)")+ 
  theme_classic()
fp1

ggsave(paste("figure/20180808-fig4c-fun-ra-e05.pdf", sep=""), fp1, width = 6.5, height = 6)


## E03-FUN-RA

f2 <- select(fun, FUN, E03_RA_ONE, E03_RA_TWO)

fun03ra <- gather(f2, key = variable , value = count, `E03_RA_ONE`:`E03_RA_TWO`)

fp2 = ggplot(fun03ra, aes(x=variable, y = count, fill = FUN )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_fill_manual(values = c("#8ea6b4", "#f8d5f0","#89a4c7"))+
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+
  ylab("Relative Abundance (%)")+ 
  theme_classic()
fp2
ggsave(paste("figure/20180808-fig4c-fun-ra-e03.pdf", sep=""), fp2, width = 6.5, height = 6)



## E05-FUN-QA

f3 <- select(fun, FUN, E05_QA_ONE, E05_QA_TWO)

fun05qa <- gather(f3, key = variable , value = count, `E05_QA_ONE`:`E05_QA_TWO`)

fp3 = ggplot(fun05qa, aes(x=variable, y = count, fill = FUN )) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_manual(values = c("#8ea6b4", "#f8d5f0","#89a4c7"))+
  xlab("Groups")+
  ylab("Quantitative Abundance (relative to host plant)")+ 
  theme_classic()
fp3

ggsave(paste("figure/20180808-fig4c-fun-qa-e05.pdf", sep=""), fp3, width = 6.5, height = 6)



## E03-FUN-QA

f4 <- select(fun, FUN, E03_QA_ONE, E03_QA_TWO)

fun03qa <- gather(f4, key = variable , value = count, `E03_QA_ONE`:`E03_QA_TWO`)

fp4 = ggplot(fun03qa, aes(x=variable, y = count, fill = FUN )) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_manual(values = c("#8ea6b4", "#f8d5f0","#89a4c7"))+
  xlab("Groups")+
  ylab("Quantitative Abundance (relative to host plant)")+ 
  theme_classic()
fp4

ggsave(paste("figure/20180808-fig4c-fun-qa-e03.pdf", sep=""), fp4, width = 6.5, height = 6)


## E03-FUN-QA-POINT

f5 <- select(fun2, FUN, E03_QA_ONE, E03_QA_TWO)

fun03qa2 <- gather(f5, key = variable , value = count, `E03_QA_ONE`:`E03_QA_TWO`)

fp5 = ggplot(fun03qa2, aes(x=variable, y = count, fill = FUN )) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_manual(values = c("#8ea6b4", "#f8d5f0","#89a4c7"))+
  xlab("Groups")+
  ylab("Quantitative Abundance (relative to host plant)")+ 
  theme_classic()
fp5

ggsave(paste("figure/20180808-fig4c-fun-qa-e03-point.pdf", sep=""), fp5, width = 6.5, height = 6)





