# 合并每个级别分类学表

# 0. 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list <- c("dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 1. 读取OTU表
otutab = read.table("~/xiaoxuan/180528/180627_AQ/bac/spikein/otutab_norm_AA.txt", header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
dim(otutab)

# 2. 读取物种注释
tax = read.table("../result/taxonomy_8.txt", header=T, row.names= 1, sep="\t",comment.char = "", stringsAsFactors = F) 
head(tax)

# 标准化，并筛选高丰度菌均值百万分之一0.0001%
# norm = t(t(otutab)/colSums(otutab,na=T))*100

# 筛选OTU
# used_otutab = read.table("AQ1/result/compare/BacHnMH63dry-BacHnMH63wet_all.txt", header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
# otutab = norm[rownames(used_otutab),]

# colSums(norm)
# idx = rowMeans(norm) > 0.0001
# HA = norm[idx,]
HA = otutab
colSums(HA)

# 数据筛选并排序，要求每个OTU必须的注释，可以为空
tax = tax[rownames(HA),]

# 转换为等级|连接格式
tax$Phylum=paste(tax$Kingdom,tax$Phylum,sep = "|")
tax$Class=paste(tax$Phylum,tax$Class,sep = "|")
tax$Order=paste(tax$Class,tax$Order,sep = "|")
tax$Family=paste(tax$Order,tax$Family,sep = "|")
tax$Genus=paste(tax$Family,tax$Genus,sep = "|")
tax$Species=paste(tax$Genus,tax$Species,sep = "|")
# head(tax)

# 按Kingdom合并
grp <- tax[rownames(tax),"Kindom", drop=F]
merge=cbind(HA, grp)
HA_Kingdom = merge %>% group_by(Kindom) %>% summarise_all(sum)
colnames(HA_Kingdom)[1]="Class"
write.table(HA_Kingdom, file="~/xiaoxuan/180528/180627_AQ/bac/spikein/k.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)

# 按Phylum合并
grp <- tax[rownames(tax), "Phylum", drop=F]
merge=cbind(HA, grp)
HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
colnames(HA_Phylum)[1]="Class"
write.table(HA_Phylum, file="~/xiaoxuan/180528/180627_AQ/bac/spikein/p.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)

# 按Class合并
grp <- tax[rownames(tax), "Class", drop=F]
merge=cbind(HA, grp)
HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
colnames(HA_Class)[1]="Class"
write.table(HA_Class, file="~/xiaoxuan/180528/180627_AQ/bac/spikein/c.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)

# 按Order合并
grp <- tax[rownames(tax), "Order", drop=F]
merge=cbind(HA, grp)
HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
colnames(HA_Order)[1]="Class"
write.table(HA_Order, file="~/xiaoxuan/180528/180627_AQ/bac/spikein/o.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)

# 按Family合并
grp <- tax[rownames(tax), "Family", drop=F]
merge=cbind(HA, grp)
HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
colnames(HA_Family)[1]="Class"
write.table(HA_Family, file="~/xiaoxuan/180528/180627_AQ/bac/spikein/f.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)

# 按Genus合并
grp <- tax[rownames(tax), "Genus", drop=F]
merge=cbind(HA, grp)
HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
colnames(HA_Genus)[1]="Class"
write.table(HA_Genus, file="~/xiaoxuan/180528/180627_AQ/bac/spikein/g.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)

# 按Species合并
grp <- tax[rownames(tax), "Species", drop=F]
merge=cbind(HA, grp)
HA_Species = merge %>% group_by(Species) %>% summarise_all(sum)
colnames(HA_Species)[1]="Class"
write.table(HA_Species, file="~/xiaoxuan/180528/180627_AQ/bac/spikein/s.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)

# 合并7个分类级
all = rbind(HA_Kingdom, HA_Phylum, HA_Class, HA_Order, HA_Family, HA_Genus, HA_Species)

# 修改样品名为组名：删除结尾的数字
colnames(all) = gsub("\\d+$","",colnames(all),perl=TRUE)

write.table(all, file="~/xiaoxuan/180528/180627_AQ/bac/spikein/a.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
