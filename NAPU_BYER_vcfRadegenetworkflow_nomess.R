#Namazonurus pustulatus vcfR -> adegenet
#By Jonathan DeBoer, Code sources from Nathan Byer
#Start: 3-31-2021

#######################################################
###    Read in Library - some are older versions    ###
#######################################################

#######packages and functions
library(SNPRelate) #use R 4.0 work around to install (below)
  #if (!requireNamespace("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
  #BiocManager::install("SNPRelate")
library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)
library(RColorBrewer)
library(vegan)
library(hierfstat)
install.packages('hierfstat')
library(apex)
library(mmod)
library(poppr)
library(HWxtest)#removed from repository
#install.packages('HWxtest')
library(plyr)
library(crunch)
library(dartR)
library("devtools")
  #if (!requireNamespace("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
  #BiocManager::install("qvalue")
library(qvalue)
install_github("zhengxwen/gdsfmt") #CHECK
install_github("zhengxwen/SNPRelate")  #CHECK

###################################################
###    set working directory and read in vcf    ###
###################################################

setwd("~/Desktop/cordylid_ddRADSeq/Napu_Floragenex_vcfs")
vcfPR<-read.vcfR("201908201427napu0000_aligned_genotypes_standard.vcf")
pop.data<-read.csv("PopData_400_020620.csv",header=T)
all(colnames(vcfPR@gt)[-1] == pop.data$Sample)

#convert vcf to genlight

aa.genlight<-vcfR2genlight(vcfPR, n.cores = 4)
ploidy(aa.genlight)<-2
pop(aa.genlight)<-pop.data$Pop
aa.genlight
library(dartR)

#keep the vcf if you want - I remove it to save RAM

rm(vcfPR)

##now you can do all sorts of fun stuff with dartR! What follows is my own attempt to replicate
##the workflow of O'Leary et al. (2018).

#some basic summary information of the genlight
nInd(aa.genlight)
nLoc(aa.genlight)

#create a loc.metrics subslot in @other - to my knowledge, this is critical - and required - for all downstream analyses.
#I have not provided this additional code, but if you run into issues with a filter.monomorphs flag - let me know, and 
#I will send you the utils.reset.monomorphs code.

aa.genlight@other$loc.metrics<-data.frame(aa.genlight@loc.names)

#put in the x-y data because why not?
aa.genlight@other$longlat<-data.frame(pop.data$UTMW,pop.data$UTMN)

#filter monomorphs (SNPs that do not vary at all)
gl<-gl.filter.monomorphs(aa.genlight)
gl@other$loc.metrics<-data.frame(gl@loc.names)
gl<-gl.recalc.metrics(gl)

#filter for SNPs < 0.05 maf
gl2<-gl.filter.maf(gl,threshold=0.05)
gl2@other$loc.metrics<-data.frame(gl2@loc.names)

#Obviously a relic of previous workflows - I used to have a subsequent step with some additional filtering here, but I removed
#it.

gl3<-gl2

#iterative, alternating filtering of locus and individual-level missing data

gl4<-gl.filter.callrate(gl3,method="loc",recalc=T,threshold=0.5)
gl5<-gl.filter.callrate(gl4,method="ind",threshold=0.1)
gl6<-gl.filter.callrate(gl5,method="loc",recalc=T,threshold=0.6)
gl7<-gl.filter.callrate(gl6,method="ind",threshold=0.3)
gl8<-gl.filter.callrate(gl7,method="loc",recalc=T,threshold=0.7)
gl9<-gl.filter.callrate(gl8,method="ind",threshold=0.5)
gl10<-gl.filter.callrate(gl9,method="ind",recalc=T,threshold=0.25)
gl11<-gl.filter.callrate(gl10,method="loc",threshold=0.05,mono.rm = T,recalc = T)
gl11<-gl.recalc.metrics(gl11)

#now, I had to come up with some code to trick dartR to actually allow for filtering of hwe. This involved
#mirroring the traditional structure of metadata provided with DART datasets.

locnames<-strsplit(as.character(gl11@loc.names),"_")
locusids<-sapply(locnames,function(x) x[1])
SNPids<-sapply(locnames,function(x) x[2])
gl11@other$loc.metrics$AlleleID<-as.factor(paste(locusids,"|F|",SNPids,sep =""))
gl11@other$loc.metrics$RepAvg<-1

#41 - make sure seed is set for reproducibility
#filter one SNP/locus - the one with the least missing data

set.seed(41)
gl12<-gl.filter.secondaries(gl11,method ="best")

#does not always work as expected - at present, this just filters for global deviations from HWE, 
#but fails when I try to run this for SNPs that deviate from HWE in at least, say, two populations.
#I use some manual methods to remove SNPs out of HWE.
#I also find looking at Fis before/after filtering to be instructive.

gl12_nohwe<-gl.filter.secondaries(gl11,method ="best")
gl12_nohwe<-gl.filter.hwe(gl12_nohwe)
gl12_nohwe_stats<-gl.basic.stats(gl12_nohwe)
gl12hwe<-gl.report.hwe(gl12,subset="all",plot=TR)
gl12hwefail<-subset(gl12hwe,gl12hwe$BonSig=="*")
gl12hwe_2<-gl.report.hwe(gl12,subset = "each")
gl12hwefail_2<-subset(gl12hwe_2,gl12hwe_2$BonSig=="*")
gl12<-gl.drop.loc(gl12,loc.list = gl12hwefail$Locus)
barplot(table(pop(gl12)), las=2)
gl12_stats<-gl.basic.stats(gl12)

#I also filter out SNPs with Heterozygote excess - they might 
#also represent SNPs with sequencing errors or other aberrations

locusexcess<-subset(gl12_stats$perloc,gl12_stats$perloc$Ho>=0.5)
locusexcesslist<-rownames(locusexcess)
gl13<-gl.drop.loc(gl12,loc.list = locusexcesslist)


#some basic summary statistics for the (hopefully final!) dataset!
gl13.stat<-gl.basic.stats(gl13)
summary(gl13@other$loc.metrics$maf)


#Fst calculations - both bootstrapped (with CIs) and single
glFst<-gl.fst.pop(gl13,nboots=100,percent=95)
glFstsingle<-gl.fst.pop(gl13,nboots=1)

#can also run glPca - automated way of running principal components analysis for genlights.
#it takes FOREVER though.

glPca(gl13)



