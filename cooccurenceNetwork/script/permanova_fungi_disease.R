#### Run PERMANOVA of bac and fungi healthy condition 

rm(list = ls())

## Health
design = read.table("../doc/design.txt", header=T, row.names= 1, sep="\t") 
design$SampleID <- row.names(design)
dist<- bray_curtis <- read.table("../RA/beta_div/bray_curtis_otutab_norm.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
design = design[rownames(design) %in% rownames(bray_curtis),]
# id <- match(row.names(design[design$Site %in% "HaiNan",]),rownames(bray_curtis))
# 
# id
# 
# dist <- bray_curtis[id,id]

# dist<- bray_curtis <- read.table("~/xiaoxuan/180528/180627_AQ/bac/RA_rem_sample_720/beta/bray_curtis_otutab.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# bray_curtis_1 <- read.table("../bac_New/table/beta_div/bray_curtis_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# eulian <- read.table("../bac_New/table/beta_div/euclidean_otu_table_css.txt", sep="\t", header=T, row.names = 1,check.names=F,stringsAsFactors = F)
# 
# design = read.table("../../180627_AQ/bac/doc/design.txt", header=T, row.names= 1, sep="\t") 
# design$SampleID <- row.names(design)

# design = design[rownames(design) %in% rownames(dist),]
map =design


pmanova <- adonis(as.dist(dist) ~  Genotype + Other + Genotype:Other ,  data = map)

map.nbs <- filter(map, Genotype != "Bulksoil")
dist.nbs <- dist[match(map.nbs$SampleID, rownames(dist)), match(map.nbs$SampleID, colnames(dist))]
pmanova.nbs <- adonis(as.dist(dist.nbs) ~ Genotype  ,  data =  map.nbs)

print("permanova result of fungi of HaiNan being as follows:");
pmanova;pmanova.nbs


