rm(list = ls())

# The following files need to be accessable in the same directory in order
# to allow correct execution of the script:
# 1) OTU_Table_wTax_DilutionExperiment.txt
# 2) Mapping_DilutionExperiment.txt
# 3) designConc_DilutionExperiment.csv
# 4) Total16S_DilutionExperiment.csv
# 5) OTU_Table_wTax_human.txt
# 6) Mapping_human.txt
#

library(hexbin)
library(extrafont)
library(ggplot2)
library(grid)
library(reshape)
library(gridExtra)
library(scales)
library(data.table)
library(plyr)
## ----libraries,message=FALSE,warning=FALSE,results="hide"----------------
libs<-c("ggplot2",
        "grid",
        "gridExtra",
        "reshape",
        "scales",
        "hexbin",
        "data.table",
        "extrafont")
lapply(libs,require,character.only=TRUE)

# install.packages("hexbin")

# install.packages("extrafont")

## ----ImportOTUTable,tidy=F-----------------------------------------------
#Import text formatted OTU Table as derived by biom convert:
scml.OTUs_tax<-read.table("./OTU_Table_wTax_DilutionExperiment.txt",sep="\t",
                          stringsAsFactors=F,header=TRUE,row.names=1)
#Delete taxonomy column
scml.OTUs<-scml.OTUs_tax[,-37] 

## ----CreateTaxaMap-------------------------------------------------------
TaxaMap<-data.frame(OTU=rownames(scml.OTUs),taxonomy=scml.OTUs_tax$taxonomy)

## ----ImportMetaData------------------------------------------------------
#Import MetaData Mapping File
scml.meta<-read.table("./Mapping_DilutionExperiment.txt",sep="\t",header=TRUE)
#Bring OTU table in same order like meta data in respect to sample ID
ord<-match(scml.meta$SampleID,colnames(scml.OTUs))
scml.OTUs<-scml.OTUs[,ord]

## ----DesignConcentration-------------------------------------------------
designConc<-read.table("./designConc_DilutionExperiment.csv",sep="\t",header=TRUE,
                       row.names=1)
designConc.names<-colnames(designConc)
designConc<-designConc[,-c(which(colnames(designConc)=="Sample.87"),
                           which(colnames(designConc)=="Sample.88"))]

## ----Import16SConc-------------------------------------------------------
total16S<-read.table("./Total16S_DilutionExperiment.csv",sep="\t",stringsAsFactors=F,
                     header=TRUE)
ord2<-match(scml.meta$SampleID,total16S$SampleID)
total16S<-total16S[ord2,]

## ----ExcludeMIDs---------------------------------------------------------
midToExclude<-c("MID39","MID40")
midIDX<-which(colnames(scml.OTUs)%in%midToExclude)

#Exclude MIDs from OTU table and convert to matrix:
scml.OTUs.trim<-as.matrix(scml.OTUs[,-midIDX])

#Exclude MIDs from the mapping table
scml.meta.trim<-scml.meta[-midIDX,]

#Exclude MIDs from the total 16S table
midIDX<-which(total16S$SampleID%in%midToExclude)
total16S.trim<-total16S[-midIDX,]

## ----SpikePositions------------------------------------------------------
#Search unique ID for S. Ruber:
SalPos<-grep("AF323500XXXX",rownames(scml.OTUs.trim))
#Search unique ID for R. radiobacter:
RhizPos<-grep("AB247615XXXX",rownames(scml.OTUs.trim))
#Search unique ID for A. acidiphilus:
AliPos<-grep("AB076660XXXX",rownames(scml.OTUs.trim))

SpikeCols<-c("#e41a1c","#4daf4a","#377eb8")

## ----Pseudocount---------------------------------------------------------
scml.OTUs.trim<-scml.OTUs.trim+1

## ----RelativeAbundances--------------------------------------------------
#Divide each sample vector by its sum of read counts:
relAbundance<-sweep(scml.OTUs.trim,2,colSums(scml.OTUs.trim),'/') 

## ----SCML_Salini---------------------------------------------------------
##SCML (S. Ruber)
#Extract S. Ruber reads
Salini<-scml.OTUs.trim[SalPos,]

#Calculate size factor:
sizeFactor_Salini<-Salini/mean(Salini)

#Apply size factor to each sample
scml.counts<-scml.OTUs.trim
for(i in seq(ncol(scml.counts))){
  scml.counts[,i]<-scml.counts[,i]/sizeFactor_Salini[i]
}

## ----AdjustingAAandRR----------------------------------------------------
#read out only spike OTU reads:
spikeReads.Raw<-scml.OTUs.trim[c(AliPos,SalPos,RhizPos),] 
rownames(spikeReads.Raw)<-c("A. acidiphilus","S. ruber","R. radiobacter")
#Create copy of spikeReads.Raw for inserting adjusted spike reads:
spikeReads.Adjusted<-spikeReads.Raw
#Calculate size factors depending on the ratio towards 
#S. ruber (on row 2 of designConc):
sfA<-designConc[1,]/designConc[2,]
sfR<-designConc[3,]/designConc[2,]
#For each spike-in apply the corresponding size factor
spikeReads.Adjusted[1,]<-as.numeric(spikeReads.Adjusted[1,])/as.numeric(sfA)
spikeReads.Adjusted[3,]<-as.numeric(spikeReads.Adjusted[3,])/as.numeric(sfR)

## ----SCMLCombined--------------------------------------------------------
#Sum the raw reads counts of S. ruber and the adjusted read counts of 
#A. acidiphilus and R. radiobacter for each sample:
summedSpikes<-colSums(spikeReads.Adjusted)

#Calculate size factor:
sizeFactor_SCML.combined<-summedSpikes/mean(summedSpikes)

#Apply size factor to each sample:
scml.counts.combined<-scml.OTUs.trim
for(i in seq(ncol(scml.counts.combined))){
  scml.counts.combined[,i]<-scml.counts.combined[,i]/sizeFactor_SCML.combined[i]
}

## ----qRT-PCR_normalization-----------------------------------------------
#Extract SYBR Green pGEM total 16S measurement from the total16S object:
total16S_SYBR<-t(total16S.trim[,1])
#Calculate size factor:
sizeFactor_qRT.PCR<-total16S_SYBR/mean(total16S_SYBR)

#Apply size factor to each sample:
qRTpcr.counts<-scml.OTUs.trim
for(i in seq(ncol(qRTpcr.counts))){
  qRTpcr.counts[,i]<-qRTpcr.counts[,i]*sizeFactor_qRT.PCR[i]
}

## ----ConstructLongTableForSpikeReads-------------------------------------
#Extract scml read counts for the spike bacteria:
spikeReads.SCML<-scml.counts[c(AliPos,SalPos,RhizPos),]
rownames(spikeReads.SCML)<-c("A. acidiphilus","S. ruber","R. radiobacter")
#Melt design concentration to later on add to the long tables:
mConc<-melt(as.matrix(designConc))
#Melt each read count table to long format and add a status and 
#design concentration column, describing the origin of the reads:
mSpikes.raw<-cbind(melt(as.matrix(spikeReads.Raw)),
                   Status="raw",
                   DesignConc=mConc$value)
mSpikes.adjusted<-cbind(melt(as.matrix(spikeReads.Adjusted)),
                        Status="adjusted",
                        DesignConc=mConc$value)
mSpikes.scml<-cbind(melt(as.matrix(spikeReads.SCML)),
                    Status="normalized",
                    DesignConc=mConc$value)
#Unify column names for merging:
colnames(mSpikes.raw)<-
  colnames(mSpikes.adjusted)<-
  colnames(mSpikes.scml)<-
  c("SpikeOTU","SampleID","Reads","Status","DesignConc")
#Combine all three long tables:
mSpikes<-rbind(mSpikes.raw,mSpikes.adjusted,mSpikes.scml)
#Add information on the dilution of each sample:
x<-match(mSpikes$SampleID,scml.meta.trim$SampleID)
#Round Dilution to 2 decimals
mSpikes$Dilution<-round(scml.meta.trim$Dilution,2)[x]
#Reverse the levels to create a ordering for the display of the dilution:
mSpikes$Dilution<-factor(mSpikes$Dilution,levels=rev(levels(factor(mSpikes$Dilution))))

## ----Figure1a_Part1------------------------------------------------------
#Scientific format function for x-axis:
scientific_10 <- function(x) {
    y<-2^x
    parse(text=gsub("e", " %*% 10^", scientific_format()(y)))
}
#Subset for raw reads:
rawData<-droplevels(subset(mSpikes,Status=="raw"))
#Subset for A. acidiphilus:
acidiSpike<-droplevels(subset(rawData,SpikeOTU=="A. acidiphilus"))


#Plot Figure 1a part 1:
Fig1a_p1<-ggplot(acidiSpike,aes(x=log2(DesignConc),y=log2(Reads),
                             colour=Dilution))+
      geom_line()+geom_point()+facet_wrap(~SpikeOTU)+
      theme_bw()+theme(text=element_text(size=13),
                  panel.grid.minor = element_blank(),
                  legend.position = "bottom", 
                  legend.box = "horizontal",
                  axis.text.x = element_text(angle = 90,
                                             hjust = 1,
                                             vjust=0.5),
                  plot.title = element_text(hjust = 0,size=18),
                  axis.title=element_text(size=14))+
      xlab(expression(paste("spiked-in 16S copies by design",sep="")))+
      ylab(expression(paste(log[2]," read counts",sep="")))+
      scale_x_continuous(breaks = pretty_breaks(n=5),labels = scientific_10)+
      guides(colour=guide_legend(title="stool dilution"))

plot(Fig1a_p1+ggtitle("a"))

## ----Figure1a_Part2------------------------------------------------------
#Subset for R. radiobacter:
radioSpike<-droplevels(subset(rawData,SpikeOTU=="R. radiobacter"))
#Plot Figure 1a part 2:
Fig1a_p2<-ggplot(radioSpike,aes(x=log2(DesignConc),y=log2(Reads),
                             colour=Dilution))+
      geom_line()+geom_point()+facet_wrap(~SpikeOTU)+
      theme_bw()+theme(text=element_text(size=13),
                  panel.grid.minor = element_blank(),
                  legend.position = "bottom",
                  legend.box = "horizontal",
                  axis.text.x = element_text(angle = 90,
                                             hjust = 1,
                                             vjust=0.5),#,size=12
                  plot.title = element_text(hjust = 0,size=18),
                  axis.title=element_text(size=14))+
      xlab(expression(paste("spiked-in 16S copies by design",sep="")))+
      ylab(expression(paste(log[2]," read counts",sep="")))+
      scale_x_continuous(breaks = pretty_breaks(n=5),labels = scientific_10)+
      guides(colour=guide_legend(title="stool dilution"))

plot(Fig1a_p2+ggtitle("a"))

## ----Figure1c------------------------------------------------------------
#Plot the adjusted read counts for all spike ins as function of the dilution:
Fig1b<-ggplot(droplevels(subset(mSpikes,Status=="adjusted")),
      aes(x=factor(Dilution,levels=rev(levels(Dilution))),y=log2(Reads)))+
      geom_boxplot(aes(fill=SpikeOTU))+
      theme_bw()+
      theme(text=element_text(size=13),
            panel.grid.major.x = element_blank(),
            legend.title = element_blank(),
            axis.ticks.x=element_blank(),
            plot.title = element_text(hjust = 0),
            legend.position = "bottom",
            legend.box = "horizontal")+
      xlab("stool dilution")+
      ylab(expression(paste(log[2]," read counts",sep="")))+
      scale_fill_manual(values=SpikeCols)+
      geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5))

print(Fig1b+ggtitle("b"))

## ----CorrelationAnalysis-------------------------------------------------
adjusted.Spikes<-droplevels(subset(mSpikes,Status=="adjusted"))

AACor<-droplevels(subset(adjusted.Spikes,SpikeOTU=="A. acidiphilus"))
SRCor<-droplevels(subset(adjusted.Spikes,SpikeOTU=="S. ruber"))
RRCor<-droplevels(subset(adjusted.Spikes,SpikeOTU=="R. radiobacter"))

AAcorrTest<-cor.test(1/as.numeric(as.character(AACor$Dilution)),log2(AACor$Reads))
SRcorrTest<-cor.test(1/as.numeric(as.character(SRCor$Dilution)),log2(SRCor$Reads))
RRcorrTest<-cor.test(1/as.numeric(as.character(RRCor$Dilution)),log2(RRCor$Reads))

## ----CreateReadCountLongTableAllAproaches,eval=FALSE,include=FALSE-------
## #SCML:
## scml.m<-melt(as.matrix(scml.counts))
## scml.m$Approach<-"SCML"
## #Relative Abundance:
## relAb.m<-melt(as.matrix(relAbundance))
## relAb.m$Approach<-"relative abundance"
## #SCML combined:
## scml.combined.m<-melt(as.matrix(scml.counts))
## scml.combined.m$Approach<-"SCML combined"
## #QPCR:
## qRTpcr.m<-melt(as.matrix(qRTpcr.counts))
## qRTpcr.m$Approach<-"total16S"
## 
## mReads<-rbind(scml.m,relAb.m,scml.combined.m,qRTpcr.m)
## colnames(mReads)<-c("OTU_ID","SampleID","Counts","Approach")

## ----CreateAllPossibleCombinations---------------------------------------
#Retrieve all Combinations of Samples:
combinations<-combn(1:ncol(scml.OTUs.trim),2)
#Name Combinations:
comb.Names<-character()
for(i in seq(ncol(combinations))){
  comb.Names[i]<-paste(colnames(scml.OTUs.trim)[combinations[1,i]],
                                   "_vs_",colnames(scml.OTUs.trim)[combinations[2,i]],
                                   sep="")
}
colnames(combinations)<-comb.Names

## ----CalculateExpectedLog2Ratios-----------------------------------------
#Create dummies
expRatio_ali<-expRatio_rhiz<-numeric()
expRatio_Background<-numeric()
#Extract amount of background sample:
background<-scml.meta.trim$Background
#Extract design concentration of A. acidiphilus and R. radiobacter:
spikesExp<-designConc[c(1,3),]
rownames(spikesExp)<-c("A. acidiphilus","R. radiobacter")
#Calculate expected log 2 ratios for A. acidiphilus and R. radiobacter, 
#aswell as the expected ratio for the sample:
for(i in seq(ncol(combinations))){
  expRatio_ali[i]<-spikesExp[1,combinations[1,i]]/spikesExp[1,combinations[2,i]]
  expRatio_rhiz[i]<-spikesExp[2,combinations[1,i]]/spikesExp[2,combinations[2,i]]
  expRatio_Background[i]<-background[combinations[1,i]]/background[combinations[2,i]]
}

## ----CalculateObservedLog2Ratios-----------------------------------------
#Dummies for A. acidiphilus
obsRatio_Ali_relAb<-obsRatio_Ali_SCML<-numeric()
#Dummies for R. radiobacter
obsRatio_Rhiz_relAb<-obsRatio_Rhiz_SCML<-numeric()

#Relative abundance subsets for A. acidiphilus and R. radiobacter:
spikes_relAb<-relAbundance[c(AliPos,RhizPos),]
#Read subsets for A. acidiphilus and R. radiobacter:
spikes_SCML<-scml.counts[c(AliPos,RhizPos),]
rownames(spikes_relAb)<-rownames(spikes_SCML)<-c("Alicyclobacillus","Rhizobium")
unique(expRatio_ali)
#Calculate ratios for A. acidiphilus and R. radiobacter over all possible
#sample comparisons:
for(i in seq(ncol(combinations))){
  #Ratios based on relative abundances:
  obsRatio_Ali_relAb[i]<-
    spikes_relAb[1,combinations[1,i]]/spikes_relAb[1,combinations[2,i]]
  obsRatio_Rhiz_relAb[i]<-
    spikes_relAb[2,combinations[1,i]]/spikes_relAb[2,combinations[2,i]]
  #Ratios based on SCML:
  obsRatio_Ali_SCML[i]<-
    spikes_SCML[1,combinations[1,i]]/spikes_SCML[1,combinations[2,i]]
  obsRatio_Rhiz_SCML[i]<-
    spikes_SCML[2,combinations[1,i]]/spikes_SCML[2,combinations[2,i]]
}

#Build data.frame with all calculated ratios for A. acidiphilus:
ratiosAli<-data.frame(expRatio=expRatio_ali,
                      obsRatio_relAb=obsRatio_Ali_relAb,
                      obsRatio_SCML=obsRatio_Ali_SCML,
                      comparison=colnames(combinations))

# to learn auto seqecence 
ratiosAli$Spike<-"A. acidiphilus"

#Build data.frame with all calculated ratios for R. radiobacter:
ratiosRhiz<-data.frame(expRatio=expRatio_rhiz,
                       obsRatio_relAb=obsRatio_Rhiz_relAb,
                       obsRatio_SCML=obsRatio_Rhiz_SCML,
                       comparison=colnames(combinations))
ratiosRhiz$Spike<-"R. radiobacter"

#Combine the tables of A. acidiphilus and R. radiobacter:
spikeRatios<-rbind(ratiosAli,ratiosRhiz)

## ----Figure2a------------------------------------------------------------
Fig2a<-ggplot(spikeRatios,aes(x=log2(expRatio),y=log2(obsRatio_relAb),
                        colour=Spike,fill=Spike))+
          geom_point(position=position_jitterdodge(dodge.width=0.4,
                                                   jitter.width=0))+
          geom_abline(intercept = 0, slope = 1,col="purple")+
          theme_bw()+annotation_logticks(base=2)+
          ylab(expression(paste("observed ",log[2]," ratio",sep="")))+
          xlab(expression(paste("expected ",log[2]," ratio",sep="")))+
          xlim(c(-11,11))+ylim(c(-11,11))+
          theme(text=element_text(size=13),
                legend.position="bottom",
                plot.title = element_text(hjust = 0))+
          ggtitle("a")+scale_colour_manual(values=SpikeCols[1:2])+
          guides(colour=FALSE,fill=FALSE)+geom_smooth()
plot(Fig2a)

## ----Figure2b------------------------------------------------------------
Fig2b<-ggplot(spikeRatios,aes(x=log2(expRatio),y=log2(obsRatio_SCML),
                        colour=Spike,fill=Spike))+
          geom_point(position=position_jitterdodge(dodge.width=0.4,
                                                   jitter.width=0))+
          geom_abline(intercept = 0, slope = 1,col="purple")+
          theme_bw()+annotation_logticks(base=2)+
          ylab(expression(paste("observed ",log[2]," ratio",sep="")))+
          xlab(expression(paste("expected ",log[2]," ratio",sep="")))+
          xlim(c(-11,11))+ylim(c(-11,11))+
          theme(legend.position="bottom",
                text=element_text(size=13),
                plot.title = element_text(hjust = 0))+
          ggtitle("b")+scale_colour_manual(values=SpikeCols[1:2])+
          guides(colour=FALSE,fill=FALSE)+geom_smooth()
plot(Fig2b)

## ----ErrorSpikes---------------------------------------------------------
#calculate the error (exp vs obs):
errorSCML<-log2(spikeRatios$expRatio)-log2(spikeRatios$obsRatio_SCML)
errorRelAb<-log2(spikeRatios$expRatio)-log2(spikeRatios$obsRatio_relAb)
#Combine Data:
error.Data<-data.frame(SCML=errorSCML,relAb=errorRelAb)

#Convert to long format for ggplot usage:
meD<-melt(error.Data)
meD$Spike<-rep(as.character(spikeRatios$Spike),2)
colnames(meD)<-c("Method","Ratio","Spike")

## ----Figure2c_Part1------------------------------------------------------
#Subset for part 1 of Figure 2c:
meD.relAb<-droplevels(subset(meD,Method=="relAb"))
meD.relAb$Method<-factor(gsub("relAb","relative abundances",meD.relAb$Method))

Fig2c_1<-ggplot(meD.relAb,aes(x=factor(Method),y=Ratio))+
  geom_boxplot(aes(fill=Spike))+xlab("approach")+
  ylab(expression(paste(atop(paste("error",sep=""), "expected vs observed"))))+
  scale_fill_manual(labels=c("A.acidiphilus","R.radiobacter"),
                    values=SpikeCols[1:2])+
  theme_bw()+theme(text=element_text(size=13),
                   legend.position="bottom",
                   axis.text.x  = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text=element_text(size=8),
                   plot.title = element_text(hjust = 0))+
  ggtitle("c")+ylim(c(-11,11))+facet_grid(Method~.,scales="free_x")
print(Fig2c_1)

## ----Figure2c_Part2------------------------------------------------------
#Subset for part 2 of Figure 2c:
meD.SCML<-droplevels(subset(meD,Method=="SCML"))

Fig2c_2<-ggplot(meD.SCML,aes(x=factor(Method),y=Ratio))+
  geom_boxplot(aes(fill=Spike))+xlab("approach")+
  ylab(expression(paste(atop(paste("error",sep=""), "expected vs observed"))))+
  scale_fill_manual(labels=c("A.acidiphilus","R.radiobacter"),
                    values=SpikeCols[1:2])+
  theme_bw()+theme(text=element_text(size=13),
                   legend.position="bottom",
                   axis.text.x  = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   legend.title = element_blank(),
                   legend.text=element_text(size=8),
                   plot.title = element_text(hjust = 0))+
  ggtitle(" ")+ylim(c(-11,11))+facet_grid(Method~.,scales="free_x")+
  guides(fill=guide_legend(ncol=2,byrow=TRUE))
print(Fig2c_2)

## ----ObsRatios_Background------------------------------------------------
#SCML scaled reads (by Salini):
reads_SCML.salini<-scml.counts[-c(SalPos,RhizPos,AliPos),]
#Relative abundances:
reads_relAb<-relAbundance[-c(SalPos,RhizPos,AliPos),]
#SCML scaled reads (by all spike-ins):
reads_SCML.combined<-scml.counts.combined[-c(SalPos,RhizPos,AliPos),]
#Total 16S-rDNA normalized counts:
reads_Qpcr<-qRTpcr.counts[-c(SalPos,RhizPos,AliPos),]
#Number of OTUs
xx<-nrow(reads_SCML.salini)
#Dummies
back_ratio<-numeric()
comp<-character()
obs_SCML<-obs_RelAb<-obs_SCML.combined<-obs_QPCR<-list()
#Calculate ratios for each approach:
for(i in seq(ncol(combinations))){
  comp<-c(comp,rep(colnames(combinations)[i],xx))
  obs_SCML[[i]]<-reads_SCML.salini[,combinations[1,i]]/reads_SCML.salini[,combinations[2,i]]
  obs_RelAb[[i]]<-reads_relAb[,combinations[1,i]]/reads_relAb[,combinations[2,i]]
  obs_SCML.combined[[i]]<-reads_SCML.combined[,combinations[1,i]]/reads_SCML.combined[,combinations[2,i]]
  obs_QPCR[[i]]<-reads_Qpcr[,combinations[1,i]]/reads_Qpcr[,combinations[2,i]]
  back_ratio<-c(back_ratio,rep(background[combinations[1,i]]/background[combinations[2,i]],xx))
}
#Create name vector for the OTUs:
OTU_Names<-rep(rownames(reads_SCML.salini),ncol(combinations))
#Unlist observed ratios:
obs_SCML<-as.character(unlist(obs_SCML))
obs_RelAb<-as.character(unlist(obs_RelAb))
obs_SCML.combined<-as.character(unlist(obs_SCML.combined))
obs_QPCR<-as.character(unlist(obs_QPCR))

#Create data.frame containing all data:
df.BG<-data.frame(OTU=OTU_Names,
                              expRatio=back_ratio,
                              comparison=comp,
                              SCML=as.numeric(obs_SCML),
                              RelAb=as.numeric(obs_RelAb),
                              SCML.combined=as.numeric(obs_SCML.combined),
                              QPCR=as.numeric(obs_QPCR))

#Check correlation between log2 expected ratio and the observed ratios calculated by the 4 different methods
RelAb.corr<-cor.test(log2(df.BG$expRatio),log2(df.BG$RelAb)) #cor=0.3590464
QPCR.corr<-cor.test(log2(df.BG$expRatio),log2(df.BG$QPCR)) #cor=0.7166334 
SCML.corr<-cor.test(log2(df.BG$expRatio),log2(df.BG$SCML)) #cor=0.8334987
SCML.combined.corr<-cor.test(log2(df.BG$expRatio),log2(df.BG$SCML.combined)) #cor=0.8452189

## ----Modified_Binning_Function-------------------------------------------
hexbinIt<-function(x,y){
  ###################################
  #' @title bin the data points to hexagons 
  #'        and recalculate percentages instead of counts for each bin
  #' @param x The value to plot on the X-axis
  #' @param y The value to plot on the Y-axis
  #' @return data.frame containing the percentages, coordinates and cellIDs 
  #'         for each hexagonal bin (ggplot compatible input)
  require(hexbin)
  require(data.table)

  #Count levels of X:
  lx<-levels(factor(x))
  totalComps<-table(x)
  #Bin data into hexgons:
  hexdata<-hexbin(x,y,IDs=TRUE)
  #Create a data.table containing the mapping between the count and the cells:
  idMap<-data.table(count=hexdata@count,cellID=hexdata@cell)
  #Create a data.table containing the mapping between ratios and the cells:
  dt<-data.table(expRatio=x,obsRatio=y,cID=hexdata@cID)
  
  for (i in seq(length(idMap$cellID))){
    #actual cellID:
    cID<-idMap$cellID[i]
    
    #look for position of cellID in dt:
    idx<-which(dt$cID==cID)
    idMap$expR[i]<-expR<-dt$expRatio[idx][1]
    
    #Get total number of values for the level of expected ratio:
    cache<-totalComps[as.character(expR)]
    names(cache)<-NULL
    idMap$total[i]<-cache
    
    #Convert counts to percentage:
    hexdata@count[i]<-round(idMap$count[i]/cache*100,digits=2)
  }
  #Create final data.frame as a mapping between cellID, 
  #percentage for factor of x and the coordinates to plot the hexagons:
  hexdf <- data.frame (cID = hexdata@cell,percentage=hexdata@count,hcell2xy (hexdata)  )
  #Insert NA's for fields with percentage 0 to omit those:
  hexdf$percentage[hexdf$percentage == 0] <- NA
  #Return final ggplot compatible data.frame
  return(hexdf)
}

## ----Figure3a------------------------------------------------------------
#Bin the log2 ratios based on relative abundances:
hexFig3a<-hexbinIt(log2(df.BG$expRatio),log2(df.BG$RelAb))
#Create Figure 3a and omit bins with percentages below 0.05:
Fig3a<-ggplot(subset(hexFig3a,percentage>0.05),aes(x=x,y=y,fill=percentage))+
  geom_hex(stat="identity",colour="grey50")+
  theme_bw()+geom_abline(intercept = 0, slope = 1,col="purple")+
  annotation_logticks(base=2)+
  ylab(expression(paste("observed ",log[2]," ratio",sep="")))+
  xlab(expression(paste("expected ",log[2]," ratio",sep="")))+
  xlim(c(-5,5))+ylim(c(-10,10))+
  scale_fill_gradientn(colours=c("white","darkred"),
                       name="percentage of values\nper exp. ratio",
                       guide="colourbar",na.value="white",
                       limits=c(.05,55))+
  ggtitle("a")+theme(text=element_text(size=13),
                     legend.position="bottom",
                     axis.text.x  = element_blank(),
                     axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     legend.text=element_text(size=8),
                     legend.key.width = unit(3, "line"),
                     legend.key.height = unit(1, "line"),
                     legend.title=element_text(size=8),
                     plot.title = element_text(hjust = 0))+geom_smooth()

print(Fig3a)

## ----Figure3b------------------------------------------------------------
#Bin the log2 ratios based on SCML:
hexFig3b<-hexbinIt(log2(df.BG$expRatio),log2(df.BG$SCML))
#Create Figure 3b and omit bins with percentages below 0.05:
Fig3b<-ggplot(subset(hexFig3b,percentage>0.05),aes(x=x,y=y,fill=percentage))+
  geom_hex(stat="identity",colour="grey50")+
  theme_bw()+geom_abline(intercept = 0, slope = 1,col="purple")+
  annotation_logticks(base=2)+
  ylab(expression(paste("observed ",log[2]," ratio",sep="")))+
  xlab(expression(paste("expected ",log[2]," ratio",sep="")))+
  xlim(c(-5,5))+ylim(c(-10,10))+
  scale_fill_gradientn(colours=c("white","darkred"),
                       name="percentage of values\nper exp. ratio",
                       guide="colourbar",
                       space = "Lab",
                       na.value="white",
                       limits=c(.05,55))+
  ggtitle("b")+theme(legend.position="bottom")+
                theme(text=element_text(size=13),
                      legend.text=element_text(size=8),
                      legend.key.width = unit(3, "line"),
                      legend.key.height = unit(1, "line"),
                      legend.title=element_text(size=8),
                      plot.title = element_text(hjust = 0))+geom_smooth()

print(Fig3b)

## ----ErrorBackground-----------------------------------------------------
#Calculate error between expected and observed log2 ratios:
error_SCML<-log2(df.BG$expRatio)-log2(df.BG$SCML)
error_RelAb<-log2(df.BG$expRatio)-log2(df.BG$RelAb)
error_SCML.combined<-log2(df.BG$expRatio)-log2(df.BG$SCML.combined)
error_QPCR<-log2(df.BG$expRatio)-log2(df.BG$QPCR)
#Create data.frame containing errors of all approaches:
errorData.BG<-data.frame(SCML.solo=error_SCML,
                         RelAb=error_RelAb,
                         SCML.combined=error_SCML.combined,
                         QPCR=error_QPCR)

#Melt data.frame to long format:
m.BG<-melt(errorData.BG)
#Rename for better plotting:
namesCache<-gsub("SCML.combined","SCML\n(combined)",m.BG$variable)
namesCache<-gsub("SCML.solo","SCML",namesCache)
namesCache<-gsub("QPCR","total 16S copies",namesCache)
namesCache<-gsub("RelAb","relative abundances",namesCache)
#use new nomenclature and sort levels:
m.BG$variable<-factor(namesCache,
                        levels=c("relative abundances",
                                 "total 16S copies",
                                 "SCML",
                                 "SCML\n(combined)"))

## ----Figure3c_Part1------------------------------------------------------
#Subset for relative abundances:
m.BG.c1<-droplevels(subset(m.BG,variable=="relative abundances"))

#Create Figure 3c part1:
Fig3c_1<-ggplot(m.BG.c1,aes(x=factor(variable),y=value))+
  geom_boxplot()+xlab("approach")+
  ylab(expression(paste(atop(paste("error",sep=""), "expected vs observed"))))+
  theme_bw()+ggtitle("c")+
  theme(text=element_text(size=13),
        legend.position="bottom",
        axis.text.x  = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0))+
  facet_grid(variable~.,scales="free_x")+ylim(c(-10,10))

print(Fig3c_1)

## ----Figure3c_Part2------------------------------------------------------
#Subset for SCML:
m.BG.c2<-droplevels(subset(m.BG,variable=="SCML"))

#Create Figure 3c part2:
Fig3c_2<-ggplot(m.BG.c2,aes(x=factor(variable),y=value))+
  geom_boxplot()+xlab("approach")+
  ylab(expression(paste(atop(paste("error",sep=""), "expected vs observed"))))+
  theme_bw()+theme(text=element_text(size=13),
                   legend.position="bottom",
                   axis.text.x  = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   legend.title = element_blank(),
                   plot.title = element_text(hjust = 0))+
  facet_grid(variable~.,scales="free_x")+ggtitle(" ")+ylim(c(-10,10))

print(Fig3c_2)

## ----Figure4a------------------------------------------------------------
Fig4a<-Fig3b

print(Fig4a)

## ----Figure4b------------------------------------------------------------
#Bin the log2 ratios based on counts normalized by total 16S copies:
hexFig4b<-hexbinIt(log2(df.BG$expRatio),log2(df.BG$QPCR))
#Create Figure 3b and omit bins with percentages below 0.05:
Fig4b<-ggplot(subset(hexFig4b,percentage>0.05),aes(x=x,y=y,fill=percentage))+
  geom_hex(stat="identity",colour="grey50")+
  theme_bw()+geom_abline(intercept = 0, slope = 1,col="purple")+
  annotation_logticks(base=2)+
  ylab(expression(paste("observed ",log[2]," ratio",sep="")))+
  xlab(expression(paste("expected ",log[2]," ratio",sep="")))+
  xlim(c(-5,5))+ylim(c(-10,10))+
  scale_fill_gradientn(colours=c("white","darkred"),
                       name="percentage of values per exp. ratio",
                       guide="colourbar",
                       space = "Lab",
                       na.value="white",
                       limits=c(.05,55))+
  ggtitle("b")+theme(text=element_text(size=13),
                     legend.position="bottom",
                     legend.text=element_text(size=8),
                     legend.key.width = unit(3, "line"),
                     legend.key.height = unit(1, "line"),
                     legend.title=element_text(size=8),
                     plot.title = element_text(hjust = 0))+geom_smooth()

print(Fig4b)

## ----Figure4c-------------------------------------------------------------
Fig4c<-ggplot(m.BG,aes(x=factor(variable),y=value))+
  geom_boxplot()+xlab("")+
  ylab(expression(paste(atop(paste("error",sep=""), "expected vs observed"))))+
  theme_bw()+theme(text=element_text(size=13),
                   legend.position="bottom",
                   axis.ticks.x = element_blank(),
                   legend.title = element_blank(),
                   plot.title = element_text(hjust = 0))

print(Fig4c)

## Calculate some statistics and summaries on the log2 differences between expected and observed for the different approaches
bgSummary<-ddply(m.BG,.(variable),function(x) summary(m.BG$value))
bgStats<-ddply(m.BG,.(variable),summarize,var=var(value),sd=sd(value),mean=mean(value),max=max(value))

## ----LoadASCTData--------------------------------------------------------
#read count table from QIIME output (converted to classic table from biom)
asct.wTax<-read.table("./OTU_Table_wTax_human.txt",sep="\t",
                      stringsAsFactors=F,header=TRUE,row.names=1)

#remove taxonomy
asct.OTUs<-asct.wTax[,-ncol(asct.wTax)]


#Import MetaData
asct.meta<-read.table("./Mapping_human.txt",sep="\t",header=TRUE)
#Bring OTU table in same order like meta data in respect to sample ID
ord<-match(asct.meta$SampleID,colnames(asct.OTUs))
asct.OTUs<-as.matrix(asct.OTUs[,ord])

## ----SpikePositions_asct-------------------------------------------------
#Search unique ID for S. Ruber:
SalPos_h<-grep("AF323500XXXX",rownames(asct.OTUs))
#Search unique ID for R. radiobacter:
RhizPos_h<-grep("AB247615XXXX",rownames(asct.OTUs))
#Search unique ID for A. acidiphilus:
AliPos_h<-grep("AB076660XXXX",rownames(asct.OTUs))

## ----SCML_Salini_asct----------------------------------------------------
##SCML (S. Ruber)
#Extract S. Ruber reads
Salini_asct<-asct.OTUs[SalPos_h,]

#Calculate size factor for later use:
sizeFactor_Salini_asct<-Salini_asct/mean(Salini_asct)

## ----CreateTaxonomyMapping_asct------------------------------------------
#Use taxonomy column from the OTU table in combination with the OTU ids:
y<-match(rownames(asct.OTUs),rownames(asct.wTax))
TaxaMap_asct<-data.frame(OTU=rownames(asct.OTUs),taxonomy=asct.wTax$taxonomy[y])

## ----BreakTaxonomyStringFunction-----------------------------------------
require(plyr)
#create column names for final table:
taxons<-c(paste("L",seq(7),sep=""))
#split the string by taxonomic levels:
tmp<-strsplit(as.character(TaxaMap_asct$taxonomy),split="; __")
#check the maximum column length and create a data frame:
tmpDF<-t(sapply(tmp, '[', seq(max(sapply(tmp, length)))))
#Subset only for the first 7 taxonomic levels:
tmpDF<-tmpDF[,seq(7)]
#Set row and column names
rownames(tmpDF)<-TaxaMap_asct$OTU
colnames(tmpDF)<-taxons
taxonomy_asct<-as.data.frame(tmpDF)

## ----EnterococcusIDX_asct------------------------------------------------
#Build a subset containing no spike reads:
OTUs_noSpikes<-asct.OTUs[-c(SalPos_h,RhizPos_h,AliPos_h),]
#Identify position of Enterococcus OTU:
ecIDX<-grep("JF896440",rownames(OTUs_noSpikes))
#Save the reads of this OTU for later use:
ecReads<-OTUs_noSpikes[ecIDX,]

## ----ConstructionPhylumReadsMatrix---------------------------------------
#Read Counts contributing to the chosen Enterococcus OTU 
#are excluded from the merge:
cache<-OTUs_noSpikes[-ecIDX,]
#Exclude all 4 excluded OTUs (3 spikes, 1 Enterococcus) 
#from the taxonomy:
taxonomy_Sub<-taxonomy_asct[rownames(cache),]
#Count and save different phyla:
Phylum<-taxonomy_Sub$L2
lvlPhyl<-levels(Phylum)
#Create dummy matrix:
PhylReads<-matrix(ncol=ncol(cache),nrow=length(lvlPhyl),0)
#Fill dummy matrix with read counts collapsed
#to phylum level:
for(i in seq(length(lvlPhyl))){
  idx<-which(Phylum==lvlPhyl[i])
  actOTUS<-rownames(taxonomy_Sub)[idx]
  readCache<-cache[which(rownames(cache)%in%actOTUS),]
  #Handle case were only one OTU contributes to current phylum:
  if(is.null(ncol(readCache))){
    PhylReads[i,]<-colSums(t(as.matrix(readCache)))
  }else{
    PhylReads[i,]<-colSums(readCache)
  }
}
#Set row and column names:
colnames(PhylReads)<-colnames(cache)
rownames(PhylReads)<-lvlPhyl
#Add an extra row for reads of Enterococcus:
PhylReads<-rbind(PhylReads,Enterococcus=ecReads)
#Only include Phyla which have more than 50 reads over 
#all samples:
PhylReads<-PhylReads[rowSums(PhylReads)>50,]

#Rename Firmicutes and Enterococcus to show that Enterococcus is shown separately:
rownames(PhylReads)<-gsub("Enterococcus","Genus: Enterococcus",rownames(PhylReads))
rownames(PhylReads)<-gsub("Firmicutes","Firmicutes w/o Enterococcus",rownames(PhylReads))

## ----ApplyBothApproaches-------------------------------------------------
#Calculate relative abundances:
standard_data<-sweep(PhylReads,2,colSums(PhylReads),'/') 
standard_data<-standard_data*100

#Apply SCML Salini only
SCML_data<-PhylReads
for(i in seq(ncol(PhylReads))){
  SCML_data[,i]<-SCML_data[,i]/sizeFactor_Salini_asct[i]
}

## ----LongTableFormat-----------------------------------------------------
require(reshape)

##Relative Abundance
m.asct_rel<-melt(standard_data)
colnames(m.asct_rel)<-c("Phylum","SampleID","RelativeAbundance")
x<-match(m.asct_rel$SampleID,asct.meta$SampleID)
#Add patient ID:
m.asct_rel$Patient<-asct.meta$Treatment[x]
#Add time point:
m.asct_rel$Time<-factor(asct.meta$Time[x],levels=c("preASCT","d0","d7","d14"))

##SCML
m.asct_scml<-melt(SCML_data)
colnames(m.asct_scml)<-c("Phylum","SampleID","CalibratedCounts")
x<-match(m.asct_scml$SampleID,asct.meta$SampleID)
#Add patient ID:
m.asct_scml$Patient<-asct.meta$Treatment[x]
#Add time point:
m.asct_scml$Time<-factor(asct.meta$Time[x],levels=c("preASCT","d0","d7","d14"))

## ----Figure_5a-----------------------------------------------------------
#Set colors for the bar plots:
asct_colors<-c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628")
#Create Figure 5a:
Figure_5a<-ggplot(data=m.asct_rel,aes(x=Time,y=RelativeAbundance,fill=Phylum))+
  geom_bar(stat="identity")+theme_bw()+
  theme(text=element_text(size=11),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust=.5),
        plot.title = element_text(hjust = 0),
        legend.position = "bottom",
        legend.text=element_text(size=10))+
  facet_wrap(~Patient,ncol = 6,scales = "free_x")+
  ylab("relative abundance")+
  scale_fill_manual(values=asct_colors)+
  ggtitle("a")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE,title=""))

print(Figure_5a)

## ----Figure_5b-----------------------------------------------------------
#Create Figure 5b:
Figure_5b<-ggplot(data=m.asct_scml,aes(x=Time,y=CalibratedCounts,fill=Phylum))+
  geom_bar(stat="identity")+theme_bw()+
  theme(text=element_text(size=11),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust=.5),
        plot.title = element_text(hjust = 0),
        legend.position = "bottom",
        legend.text=element_text(size=10),
        axis.text.y=element_blank())+
  facet_wrap(~Patient,ncol = 6,scales = "free_x")+
  ylab("microbial abundance (scaled)")+
  scale_fill_manual(values=asct_colors)+
  ggtitle("b")+guides(fill=guide_legend(nrow=2,byrow=TRUE,title=""))+
  scale_y_continuous(breaks=NULL)

print(Figure_5b)

## ----SubsettingEnterococcusData------------------------------------------
patientsECpositive<-c("Patient2","Patient4","Patient5")
#Subset the standard data for relative abundances of Enterococcus:
mEC<-subset(m.asct_rel,Phylum=="Genus: Enterococcus")
#Add the calibrated read counts for Enterococcus from the SCML data:
mEC$CalibratedCounts<-subset(m.asct_scml,
                             Phylum=="Genus: Enterococcus")$CalibratedCounts
#Subset only for Enterococcus positive Patients:
mEC<-subset(mEC,Patient%in%patientsECpositive)

## ----CalculateLog2Ratios_Enterococcus------------------------------------
#Calculate log2ratios between the last and first captured time point
#for Patient 2, 4 and 5 based on standard data:
log2ratio_relative<-c(log2(mEC$RelativeAbundance[4]/mEC$RelativeAbundance[1]),
                      log2(mEC$RelativeAbundance[7]/mEC$RelativeAbundance[5]),
                      log2(mEC$RelativeAbundance[10]/mEC$RelativeAbundance[8]))
#Calculate log2ratios between the last and first captured time point
#for Patient 2, 4 and 5 based on SCML data:
log2ratio_scml<-c(log2(mEC$CalibratedCounts[4]/mEC$CalibratedCounts[1]),
                      log2(mEC$CalibratedCounts[7]/mEC$CalibratedCounts[5]),
                      log2(mEC$CalibratedCounts[10]/mEC$CalibratedCounts[8]))

#Create patient vector:
patientVec<-mEC$Patient[c(1,7,10)]
#Create comparison vector:
comparison<-c(paste0(mEC$Time[4],"_vs_",mEC$Time[1]),
              paste0(mEC$Time[7],"_vs_",mEC$Time[5]),
              paste0(mEC$Time[10],"_vs_",mEC$Time[8]))

#Combine to one Table:
log2.ratio<-data.frame(Patient=patientVec,
                       Comparison=comparison,
                       RelativeAbundance=log2ratio_relative,
                       SCML=log2ratio_scml)

m.log2.ratio<-melt(log2.ratio)
colnames(m.log2.ratio)[3:4]<-c("Method","log2ratio")
methodVec<-gsub("RelativeAbundance","relative abundance",m.log2.ratio$Method)
m.log2.ratio$Method<-factor(methodVec,levels=c("relative abundance","SCML"))

## ----Figure_5c-----------------------------------------------------------
Figure_5c<-ggplot(m.log2.ratio,aes(x=Patient,y=log2ratio,fill=Method))+
  geom_bar(stat="identity",position=position_dodge())+
  theme_bw()+xlab("")+
  ylab(expression(paste("Enterococcus ", 
            log[2]," ratio between last and first captured timepoint",
            sep="")))+
  theme(text=element_text(size=13),
        axis.text.x=element_text(angle=90, vjust=0.5, size=10),
        plot.title = element_text(hjust = 0),
        legend.position="bottom",
        legend.text=element_text(size=8))+
  ggtitle("c")+scale_fill_manual(values=c("#bdbdbd", "#636363"),
                                 name="",breaks=c("relative abundance","SCML"))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

print(Figure_5c)

## ----SessionInfo---------------------------------------------------------
sessionInfo()
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 16.04.3 LTS
# 
# Matrix products: default
# BLAS: /usr/lib/openblas-base/libblas.so.3
# LAPACK: /usr/lib/libopenblasp-r0.2.18.so
# 
# locale:
# #   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] plyr_1.8.4          data.table_1.10.4-3 scales_0.4.1        gridExtra_2.3       reshape_0.8.7       ggplot2_2.2.1       extrafont_0.17      hexbin_1.27.1      
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.14         Rttf2pt1_1.3.5       magrittr_1. 5         munsell_0.4.3        colorspace_1.3-2     lattice_0.20-35      rlang_0.1.6          stringr_1.2.0        tools_3.4.1          gtable_0.2.0         lambda.r_1.2        
# [12] futile.logger_1.4.3  extrafontdb_1.0      digest_0.6.13        lazyeval_0.2.1       tibble_1.4.1         Matrix_1.2-12        reshape2_1.4.2       futile.options_1.0.0 labeling_0.3         stringi_1.1.6        compiler_3.4.1      
# [23] pillar_1.0.1        
