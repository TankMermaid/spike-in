
#### Step 1 - Installation and biom format conversion

# library(devtools)
# install_github("zdk123/SpiecEasi")
library(SpiecEasi)

source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)

# install.packages(c('geigen'))
# install.packages("stinepack")
# install_github("hallucigenia-sparsa/seqtime") 
library(seqtime)


#### Step 2 - Data import and preprocessing
input.path="/mnt/bai/xiaoning/xiaoxuan/180213/fun_New/script/conet/"
output.path="/mnt/bai/xiaoning/xiaoxuan/180213/fun_New/script/conet/"
biom.path=file.path(input.path,"arctic_soils_json.biom")


phyloseqobj=import_biom(biom.path)
otus=otu_table(phyloseqobj)
taxa=tax_table(phyloseqobj)


filterobj=filterTaxonMatrix(otus,minocc=20,keepSum = TRUE, return.filtered.indices = TRUE)
otus.f=filterobj$mat
taxa.f=taxa[setdiff(1:nrow(taxa),filterobj$filtered.indices),]
dummyTaxonomy=c("k__dummy","p__","c__","o__","f__","g__","s__")
taxa.f=rbind(taxa.f,dummyTaxonomy)
rownames(taxa.f)[nrow(taxa.f)]="0"
rownames(otus.f)[nrow(otus.f)]="0"


updatedotus=otu_table(otus.f, taxa_are_rows = TRUE)
updatedtaxa=tax_table(taxa.f)
phyloseqobj.f=phyloseq(updatedotus, updatedtaxa)



#### Step 3 - SPIEC-EASI Run

spiec.out=spiec.easi(phyloseqobj.f, method="mb",icov.select.params=list(rep.num=20))


spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(phyloseqobj.f)))
plot_network(spiec.graph, phyloseqobj.f, type='taxa', color="Rank3", label=NULL)



#### Step 4 - SPIEC-EASI Analysis

betaMat=as.matrix(symBeta(getOptBeta(spiec.out)))


clusters=cluster_fast_greedy(spiec.graph)
clusterOneIndices=which(clusters$membership==1)
clusterOneOtus=clusters$names[clusterOneIndices]
clusterTwoIndices=which(clusters$membership==2)
clusterTwoOtus=clusters$names[clusterTwoIndices]

sort(table(getTaxonomy(clusterOneOtus,taxa.f,useRownames = TRUE)),decreasing = TRUE)
sort(table(getTaxonomy(clusterTwoOtus,taxa.f,useRownames = TRUE)),decreasing = TRUE)

### Step 5 - Importing into Cytoscape


write.table(taxa.f,file=file.path(output.path,"arctic_soil_lineages.txt"),sep="\t", quote=FALSE)



