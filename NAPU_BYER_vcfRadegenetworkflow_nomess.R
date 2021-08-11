#Namazonurus pustulatus vcfR -> adegenet = Post-assemly filtering
#By Jonathan DeBoer, Code sources from Nathan Byer
#Start: 3-31-2021

#######################################################
###    Read in Library - some are older versions    ###
#######################################################

#######packages and functions
library(vcfR)
library(adegenet)
library(adegraphics)
library(dartR)
utils.reset.flags <- function(x, set=FALSE, value=2, verbose=NULL) {
  
  # TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
  # SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
  # FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
  # STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data\n")}
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }
  
  # SCRIPT SPECIFIC ERROR TESTING
  
  if (value < 0 | value > 5){
    cat("  Warning: Parameter 'value' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    value <- 2
  }
  
  # DO THE JOB
  ##################################   
  if (data.type=="SNP"){
    if(verbose >= 2){
      cat("  Resetting flags for AvgPIC, OneRatioRef, OneRatioSnp, PICRef, PICSnp, CallRate, maf, FreqHets, FreqHomRef, FreqHomSnp, monomorphs, OneRatio, PIC to",set,"\n")
      cat("  Resetting SilicoDArT flags for OneRatio, PIC to FALSE\n")
    }
    
    # Check if the x@other$loc.metrics slot exists, if not, create as a dataframe
    if (is.null(x@other$loc.metrics)) {
      x@other$loc.metrics <- as.data.frame(array(NA,nLoc(x)))
    }  
    # Check if the x@other$loc.metrics.flags slot exists, if not, create as a dataframe
    if (is.null(x@other$loc.metrics.flags)) {
      x@other$loc.metrics.flags <- as.data.frame(array(NA,1))
    }
    
    #AvgPIC
    if (is.null(x@other$loc.metrics$AvgPIC)) {
      x@other$loc.metrics$AvgPIC <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric AvgPIC does not exist, creating slot @other$loc.metrics$AvgPIC\n")
      }
    }
    x@other$loc.metrics.flags$AvgPIC <- set
    
    #OneRatioRef  
    if (is.null(x@other$loc.metrics$OneRatioRef)) {
      x@other$loc.metrics$OneRatioRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatioRef does not exist, creating slot @other$loc.metrics$OneRatioRef\n")
      }
    }
    x@other$loc.metrics.flags$OneRatioRef <- set
    
    #OneRatioSnp        
    if (is.null(x@other$loc.metrics$OneRatioSnp)) {
      x@other$loc.metrics$OneRatioSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatioSnp does not exist, creating slot @other$loc.metrics$OneRatioSnp\n")
      }
    }
    x@other$loc.metrics.flags$OneRatioSnp <- set
    
    #PICRef    
    if (is.null(x@other$loc.metrics$PICRef)) {
      x@other$loc.metrics$PICRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PICRef does not exist, creating slot @other$loc.metrics$PICRef\n")
      }
    }
    x@other$loc.metrics.flags$PICRef <- set
    
    #PICSnp
    if (is.null(x@other$loc.metrics$PICSnp)) {
      x@other$loc.metrics$PICSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PICSnp does not exist, creating slot @other$loc.metrics$PICSnp\n")
      }
    }
    x@other$loc.metrics.flags$PICSnp <- set
    
    #CallRate
    if (is.null(x@other$loc.metrics$CallRate)) {
      x@other$loc.metrics$CallRate <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric CallRate does not exist, creating slot @other$loc.metrics$CallRate\n")
      }
    }
    x@other$loc.metrics.flags$CallRate <- set
    
    #FreqHomRef
    if (is.null(x@other$loc.metrics$FreqHomRef)) {
      x@other$loc.metrics$FreqHomRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHomRef does not exist, creating slot @other$loc.metrics$FreqHomRef\n")
      }
    }
    x@other$loc.metrics.flags$FreqHomRef <- set
    
    #FreqHomSnp
    if (is.null(x@other$loc.metrics$FreqHomSnp)) {
      x@other$loc.metrics$FreqHomSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHomSnp does not exist, creating slot @other$loc.metrics$FreqHomSnp\n")
      }
    }
    x@other$loc.metrics.flags$FreqHomSnp <- set
    
    #FreqHets
    if (is.null(x@other$loc.metrics$FreqHets)) {
      x@other$loc.metrics$FreqHets <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHets does not exist, creating slot @other$loc.metrics$FreqHets\n")
      }
    }
    x@other$loc.metrics.flags$FreqHets <- set
    
    #monomorphs
    if (is.null(x@other$loc.metrics$monomorphs)) {
      x@other$loc.metrics$monomorphs <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric monomorphs does not exist, creating slot @other$loc.metrics$monomorphs\n")
      }
    }
    x@other$loc.metrics.flags$monomorphs <- set
    
    #maf
    if (is.null(x@other$loc.metrics$maf)) {
      x@other$loc.metrics$maf <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric maf does not exist, creating slot @other$loc.metrics$maf\n")
      }
    }
    x@other$loc.metrics.flags$maf <- set
    
    #OneRatio
    if (is.null(x@other$loc.metrics$OneRatio)) {
      x@other$loc.metrics$OneRatio <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatio does not exist, creating slot @other$loc.metrics$OneRatio\n")
      }
    }
    x@other$loc.metrics.flags$OneRatio <- FALSE
    
    #PIC
    if (is.null(x@other$loc.metrics$PIC)) {
      x@other$loc.metrics$PIC <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PIC does not exist, creating slot @other$loc.metrics$PIC\n")
      }
    }
    x@other$loc.metrics.flags$PIC <- FALSE
    
    #monomorphs
    # if (is.null(x@other$loc.metrics$monomorphs)) {
    #   x@other$loc.metrics$monomorphs <- array(NA,nLoc(x))
    #   if (verbose >= 3){
    #     cat("  Locus metric monomorphs does not exist, creating slot @other$loc.metrics$monomorphs\n")
    #   }
    # }
    x@other$loc.metrics.flags$monomorphs <- set
    
    #verbosity
    if (is.null(x@other$verbose)) {
      x@other$verbose <- 2
      if (verbose >= 3){
        cat("  Locus metric 'verbose' does not exist, creating slot @other$verbose, setting to default [2]\n")
      }
    }
  }
  ################################## 
  if (data.type=="SilicoDArT"){
    if(verbose >= 2){
      cat("  Resetting flags for CallRate, PIC, OneRatio, monomorphs to",set,"\n")
      cat("  Setting SNP flags for AvgPIC, OneRatioRef, OneRatioSnp, PICRef, PICSnp, maf, FreqHets, FreqHomRef, FreqHomSnp to FALSE\n")
    }
    
    # Check if the x@other$loc.metrics slot exists, if not, create as a dataframe
    if (is.null(x@other$loc.metrics)) {
      x@other$loc.metrics <- as.data.frame(array(NA,nLoc(x)))
    }  
    # Check if the x@other$loc.metrics.flags slot exists, if not, create as a dataframe
    if (is.null(x@other$loc.metrics.flags)) {
      x@other$loc.metrics.flags <- as.data.frame(array(NA,1))
    }
    
    #AvgPIC
    if (is.null(x@other$loc.metrics$AvgPIC)) {
      x@other$loc.metrics$AvgPIC <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric AvgPIC does not exist, creating slot @other$loc.metrics$AvgPIC\n")
      }
    }
    x@other$loc.metrics.flags$AvgPIC <- FALSE
    
    #OneRatioRef  
    if (is.null(x@other$loc.metrics$OneRatioRef)) {
      x@other$loc.metrics$OneRatioRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatioRef does not exist, creating slot @other$loc.metrics$OneRatioRef\n")
      }
    }
    x@other$loc.metrics.flags$OneRatioRef <- FALSE
    
    #OneRatioSnp        
    if (is.null(x@other$loc.metrics$OneRatioSnp)) {
      x@other$loc.metrics$OneRatioSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatioSnp does not exist, creating slot @other$loc.metrics$OneRatioSnp\n")
      }
    }
    x@other$loc.metrics.flags$OneRatioSnp <- FALSE
    
    #PICRef    
    if (is.null(x@other$loc.metrics$PICRef)) {
      x@other$loc.metrics$PICRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PICRef does not exist, creating slot @other$loc.metrics$PICRef\n")
      }
    }
    x@other$loc.metrics.flags$PICRef <- FALSE
    
    #PICSnp
    if (is.null(x@other$loc.metrics$PICSnp)) {
      x@other$loc.metrics$PICSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PICSnp does not exist, creating slot @other$loc.metrics$PICSnp\n")
      }
    }
    x@other$loc.metrics.flags$PICSnp <- FALSE
    
    #CallRate
    if (is.null(x@other$loc.metrics$CallRate)) {
      x@other$loc.metrics$CallRate <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric CallRate does not exist, creating slot @other$loc.metrics$CallRate\n")
      }
    }
    x@other$loc.metrics.flags$CallRate <- set
    
    #FreqHomRef
    if (is.null(x@other$loc.metrics$FreqHomRef)) {
      x@other$loc.metrics$FreqHomRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHomRef does not exist, creating slot @other$loc.metrics$FreqHomRef\n")
      }
    }
    x@other$loc.metrics.flags$FreqHomRef <- FALSE
    
    #FreqHomSnp
    if (is.null(x@other$loc.metrics$FreqHomSnp)) {
      x@other$loc.metrics$FreqHomSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHomSnp does not exist, creating slot @other$loc.metrics$FreqHomSnp\n")
      }
    }
    x@other$loc.metrics.flags$FreqHomSnp <- FALSE
    
    #FreqHets
    if (is.null(x@other$loc.metrics$FreqHets)) {
      x@other$loc.metrics$FreqHets <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHets does not exist, creating slot @other$loc.metrics$FreqHets\n")
      }
    }
    x@other$loc.metrics.flags$FreqHets <- FALSE
    
    #monomorphs
    # if (is.null(x@other$loc.metrics$monomorphs)) {
    #   x@other$loc.metrics$monomorphs <- array(NA,nLoc(x))
    #   if (verbose >= 3){
    #     cat("  Locus metric monomorphs does not exist, creating slot @other$loc.metrics$monomorphs\n")
    #   }
    # }
    x@other$loc.metrics.flags$monomorphs <- set
    
    #maf
    if (is.null(x@other$loc.metrics$maf)) {
      x@other$loc.metrics$maf <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric maf does not exist, creating slot @other$loc.metrics$maf\n")
      }
    }
    x@other$loc.metrics.flags$maf <- FALSE
    
    #OneRatio
    if (is.null(x@other$loc.metrics$OneRatio)) {
      x@other$loc.metrics$OneRatio <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatio does not exist, creating slot @other$loc.metrics$OneRatio\n")
      }
    }
    x@other$loc.metrics.flags$OneRatio <- set
    
    #PIC
    if (is.null(x@other$loc.metrics$PIC)) {
      x@other$loc.metrics$PIC <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PIC does not exist, creating slot @other$loc.metrics$PIC\n")
      }
    }
    x@other$loc.metrics.flags$PIC <- set
    
    #verbosity
    if (is.null(x@other$verbose)) {
      x@other$verbose <- 2
      if (verbose >= 3){
        cat("  Locus metric 'verbose' does not exist, creating slot @other$verbose, setting to default [2]\n")
      }
    }
  }
  
  # ADD TO HISTORY not in utils functions
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat("Completed:", funname, "\n")
  }
  
  return(x)
  
}

###################################################
###    set working directory and read in vcf    ###
###################################################

setwd("~/Desktop/cordylid_ddRADSeq/Napu_Floragenex_vcfs")
#Read in the VCF file - using the relaxed version to heavily filter
vcfPR<-read.vcfR("201908201427napu0000_aligned_genotypes_relaxed.vcf")

#Read in the cleaned metadata, all localities should be cleaned, and added morphological data.
pop.data<-read.csv("napu_popdata_final.csv",header=T)
all(colnames(vcfPR@gt)[-1] == pop.data$Sample)

#convert vcf to genlight
aa.genlight<-vcfR2genlight(vcfPR, n.cores = 4)
ploidy(aa.genlight)<-2
pop(aa.genlight)<-pop.data$Pop
aa.genlight

#Remove vcf to save RAM
rm(vcfPR)

#make sure names line up
pop.data$ID<-as.character(pop.data$ID)
for (i in 1:length(pop.data$ID)){
  pop.data$ID[i]<-paste(pop.data$ID[i],"_1",sep="")
}

pop.data2<-pop.data
for (i in 1:length(pop.data$ID)){
  position<-which(aa.genlight@ind.names==pop.data$ID[i])
  pop.data2[position,]<-pop.data[i,]
}
pop.data<-pop.data2

#some basic summary information of the genlight
nInd(aa.genlight)
nLoc(aa.genlight)

#create a loc.metrics subslot in @other - to my knowledge, this is critical - and required - for all downstream analyses.
#I will send you the utils.reset.monomorphs code.
#I had flag probelms, used utils.reset -> depreciated function
aa.genlight<-utils.reset.flags(aa.genlight)
head(aa.genlight@other$loc.metrics)
aa.genlight@other$loc.metrics<-data.frame(aa.genlight@loc.names)
head(aa.genlight@other$loc.metrics)
aa.genlight<-gl.recalc.metrics(aa.genlight)
head(aa.genlight@other$loc.metrics)

####################################
###   Add Metadata to genlight   ###
####################################
#Add the x-y data -> comes with pop data
aa.genlight@other$latlong<-data.frame(pop.data$Longitude,pop.data$Latitude)
#Add the locality information
aa.genlight@other$locality<-data.frame(pop.data$Locality)
#Add the elevation data
aa.genlight@other$elev<-data.frame(pop.data$Elev)
#Add the sex data 
aa.genlight@other$sex<-data.frame(pop.data$Sex)
#Add the age data
aa.genlight@other$age<-data.frame(pop.data$age)
#Add morphological traits that significantly vary across elevation
#SVL
aa.genlight@other$svl<-data.frame(pop.data$SVL)
#Head width
aa.genlight@other$head_w<-data.frame(pop.data$Head_w)
#Head length
aa.genlight@other$head_l<-data.frame(pop.data$Head_l)
#Head height
aa.genlight@other$head_h<-data.frame(pop.data$Head_h)
#Tail length
aa.genlight@other$tail<-data.frame(pop.data$Tail)

####################################
###           FILTERING          ###
####################################
##now you can do all sorts of fun stuff with dartR! What follows is my own attempt to replicate
##the workflow of O'Leary et al. (2018).

############################################################## GOOD
###   Filter - MONOMORPHS (SNPs that do not vary at all)   ###
##############################################################
gl_no_mono<-gl.filter.monomorphs(aa.genlight)
head(gl_no_mono@other$loc.metrics)

############################################################### GOOD
###      Filter - MAFS (minor allele frequencies)           ###
###############################################################
gl_no_mono_01<-gl.filter.maf(gl_no_mono,threshold=0.01)
gl_no_mono_02<-gl.filter.maf(gl_no_mono,threshold=0.02)
gl_no_mono_03<-gl.filter.maf(gl_no_mono,threshold=0.03)
gl_no_mono_04<-gl.filter.maf(gl_no_mono,threshold=0.04)
gl_no_mono_05<-gl.filter.maf(gl_no_mono,threshold=0.05) 
gl_no_mono_06<-gl.filter.maf(gl_no_mono,threshold=0.06) 

#looking at mafs on charts
#gl.report.maf(gl_no_mono)
#gl.report.maf(gl_no_mono_01)
#gl.report.maf(gl_no_mono_02)
#gl.report.maf(gl_no_mono_03)
#gl.report.maf(gl_no_mono_04)
#gl.report.maf(gl_no_mono_05)
#gl.report.maf(gl_no_mono_06)

#compare the different maf genlights
#gl_no_mono      #SNPs = 48,978
#gl_no_mono_01   #SNPs = 23,857
#gl_no_mono_02   #SNPs = 18,221
#gl_no_mono_03   #SNPs = 12,636
#gl_no_mono_04   #SNPs = 9,438
#gl_no_mono_05   #SNPs = 8,575
#gl_no_mono_06   #SNPs = 7,897

## WHICH IS THE BEST?
#Final MAF filtered genlight is ...
gl_no_mono_MAF02<-gl_no_mono_02
#gl_no_mono_MAF02<-gl_no_mono_01

############################################################### GOOD
### Filter - Missing data (locus and indv-lvl missing data) ###
###############################################################
#Iterative, alternating filtering of locus and individual-level missing data
gl3<-gl.filter.callrate(gl_no_mono_MAF02,method="loc",recalc=T,threshold=0.5)
gl4<-gl.filter.callrate(gl3,method="ind",threshold=0.1)
gl5<-gl.filter.callrate(gl4,method="loc",recalc=T,threshold=0.6)
gl6<-gl.filter.callrate(gl5,method="ind",threshold=0.3)
gl7<-gl.filter.callrate(gl6,method="loc",recalc=T,threshold=0.7)
gl8<-gl.filter.callrate(gl7,method="ind",threshold=0.5)
gl9<-gl.filter.callrate(gl8,method="ind",recalc=T,threshold=0.25)
gl10<-gl.filter.callrate(gl9,method="loc",threshold=0.05,mono.rm = T,recalc = T)
gl_no_mono_MAF02_md10<-gl.recalc.metrics(gl10)

#Plot to compare - missing data
#glPlot(gl_no_mono_MAF02,posi="topleft")      # Pre-Filtering
#glPlot(gl_no_mono_MAF02_md10,posi="topleft") # Post-Filtering

#Plot missing data over entire alignment
#Plot indicates that the missing data is at the beginning of the alignment
#temp <- density(glNA(gl_no_mono_MAF02), bw=10)
#plot(temp, type="n", xlab="Position in the alignment", main="Location of the missing values (NAs)",
#     xlim=c(0,150))
#polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3))
#points(glNA(gl_no_mono_MAF02), rep(0, nLoc(gl_no_mono_MAF02)), pch="|", col="red")
#Plot again after Post-Filtering
#temp <- density(glNA(gl_no_mono_MAF02_md10), bw=10)
#plot(temp, type="n", xlab="Position in the alignment", main="Location of the missing values (NAs)",
#     xlim=c(0,150))
#polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3))
#points(glNA(gl_no_mono_MAF02_md10), rep(0, nLoc(gl_no_mono_MAF02_md10)), pch="|", col="blue")

#dev.off()


################################
###     additional setup     ###
################################
#adding allele IDs
loc_string<-strsplit(as.character(gl_no_mono_MAF02_md10@loc.names),"_")
locus_ids<-sapply(loc_string,function(x) x[3]) #Loci
SNPids<-sapply(loc_string,function(x) x[4])   #SNP 
gl_no_mono_MAF02_md10@other$loc.metrics$AlleleID<-as.factor(paste(locus_ids,"|F|",SNPids,sep =""))
gl_no_mono_MAF02_md10@other$loc.metrics$RepAvg<-1

#assign pops
gl_no_mono_MAF02_md10_names<-sapply(strsplit(indNames(gl_no_mono_MAF02_md10),"_"), `[`, 1)
#names dont match up to reverting.
for (i in 1:length(pop.data$ID)){
  pop.data$ID[i]<-sub("_1","",pop.data$ID[i])
}
gl_localities<-pop.data$Locality[as.character(pop.data$ID) %in% gl_no_mono_MAF02_md10_names]
gl_no_mono_MAF02_md10@pop<-as.factor(gl_localities) #had to make as factor in newer version

################################
###   Filter - secondaries   ###
################################
#secondaries that likely need to be filtered
gl_no_secondaries<-gl.filter.secondaries(gl_no_mono_MAF02_md10,method ="best")

###############################################################
###   Filter - HW (Check for Hardy-Weingberg equilibrium)   ###
###############################################################
#make sure seed is set for reproducibility
set.seed(41)

#Nathan had to come up with some code to trick dartR to actually allow for filtering of hwe. 
#This involved mirroring the traditional structure of metadata provided with DART datasets.

#The gl.filter.hwe function does not always work as expected 
#at present, this just filters for global deviations from HWE, 
#but fails when I try to run this for SNPs that deviate from HWE in at least.

#Use  manual methods to remove SNPs out of HWE by checking for loci out of HWE by EACH population.
#goal is to filter for loci that are out of HWE in at least 2 pops. 
#Looking at Fis before/after filtering to be instructive.
gl_nohwe_stats<-gl.basic.stats(gl_no_secondaries)
#Fis = 0.1386 MAF = 0.2
#find out if loci are out of HWE 
#gl_hwe<-gl.filter.hwe(gl_no_secondaries) #doesnt work.
glhwe_each<-gl.report.hwe(gl_no_secondaries,subset="each",plot=TR)
#subset out loci that signigicantly differ from HWE.
glhwefail_each<-subset(glhwe_each,glhwe_each$BonSig=="*")
#drop said loci
gl_no_hwe_each<-gl.drop.loc(gl_no_secondaries,loc.list = glhwefail_each$Locus)

#check Fis again, previously 0.2986
gl_no_hwe_each_stats<-gl.basic.stats(gl_no_hwe_each)
#Fis = 0.1577


############################################################### NOT COMPLETE
###   Filter - Heterozygotes (SNPs that are heterozygous)   ###
###############################################################
# proportion of samples which score as heterozygous
# plot(gl_no_mono_MAF02_md10$other$loc.metrics$FreqHets)
# which(gl_no_mono_MAF02_md10$other$loc.metrics$FreqHets>0.8) #80% or more
# length(which(gl_no_mono_MAF02_md10$other$loc.metrics$FreqHets>0.8)) #173 individuals greater than 80%

#Visualize with a PCoA - 0 for homo reference, 1 for heterozygous, 2 for alternative homo.
#pc<-gl.pcoa(gl_no_mono_MAF02_md10,nfactor=5)

#NATHAN's CODE
#I also filter out SNPs with Heterozygote excess - they might 
#also represent SNPs with sequencing errors or other aberrations

#gl12_stats<-gl.basic.stats(gl_no_mono_MAF02_md10)

locusexcess<-subset(gl_no_hwe_each_stats$perloc,gl_no_hwe_each_stats$perloc$Ho>=0.5)
locusexcesslist<-rownames(locusexcess)
gl_final<-gl.drop.loc(gl_no_hwe_each,loc.list = locusexcesslist)

############################################################### 
###     FINAL - Summary statististics for final dataset     ###
###############################################################
gl_final_stats<-gl.basic.stats(gl_final)

#at MAF 0.02
#individuals = 177
#SNPs = 8,613
#Missing data = 3.23%


#at MAF 0.01
#individuals = 177
#SNPs = 10,758
#Missing data = 3.21%

#######################################
###         Calculating Fst         ###
#######################################
#Make sure populations are assigned!
#first check the localities
unique(pop.data$Locality)

#running FST with khw (n=2) and elis (n=1)
gl_all_fst<-gl.fst.pop(gl_final, nboots=100,percent = 95, nclusters=1)
gl_all_fst


#dropping khw and elis because very few individuals
gl_final_reduced<-gl.drop.pop(gl_final,c("khw","elis"))
gl_final_reduced@pop

gl_final_reduced_fst<-gl.fst.pop(gl_final_reduced, nboots=100,percent = 95, nclusters=4)
gl_final_reduced_fst
gl_final_reduced_stat<-gl.basic.stats(gl_final_reduced)
gl_final_reduced_stat
gl_final_reduced

#?spca(gl13)

#######################################
###         Calculating IBD         ###
#######################################
#gl.ibd doesnt seem to work 
gl.ibd(gl_final_reduced)
#I changed one part to not subset coordinates using brackets.
gl.ibd_fixed<- function(gl = NULL, Dgen = NULL, Dgeo = NULL, projected = FALSE,
                       permutations = 1000, plot = TRUE)
{
  if (!(requireNamespace("dismo", quietly = TRUE))) {
    stop("Package dismo needed for this function to work. Please install it.")
  }
  else {
    if (!is.null(Dgen) & !is.null(Dgeo))
      cat("Analysis performed on provided genetic and Euclidean distance matrices.")
    if (class(gl) == "genlight") {
      cat("Standard analysis performed on the genlight object. Mantel test and plot will be Fst/1-Fst versus log(distance)\n")
      if (nrow(gl@other$latlong) != nInd(gl))
        stop("Cannot find coordinates for each individual in slot @other$latlong")
      if (sum(match(names(gl@other$latlong), "long"), na.rm = T) == 1)
        gl@other$latlong$lon <- gl@other$latlong$long
      if (!projected) {
        xy <- dismo::Mercator(gl@other$latlong)#[, c("lon", "lat")])
        cat("Coordinates transformed to Mercator (google) projection to calculate distances in meters.\n")
      }
      else {
        xy = gl@other$latlong[, c("lon", "lat")]
        cat("Coordinates not transformed. Distances calculated on the provided coordinates.")
      }
      pop.xy <- apply(xy, 2, function(a) tapply(a, pop(gl),
                                                mean, na.rm = T))
      Dgeo <- dist(pop.xy)
      Dgeo <- log(Dgeo)
      Dgen <- as.dist(StAMPP::stamppFst(gl, nboots = 1))
      Dgen <- Dgen/(1 - Dgen)
      ordering <- levels(pop(gl))
      Dgen <- as.dist(as.matrix(Dgen)[ordering, ordering])
      Dgeo <- as.dist(as.matrix(Dgeo)[ordering, ordering])
    }
    miss = FALSE
    if (sum(is.na(Dgen)) > 0 | sum(is.infinite(Dgen)) > 0 ) {
      miss = TRUE
      cat("There are missing values in the genetic distance matrix. No kernel distance plot is possible.\n")
    }
    if (sum(is.na(Dgeo)) > 0 | sum(is.infinite(Dgeo)) > 0 ) {
      miss = TRUE
      cat("There are missing values in the geographic distance matrix. No kernel distance plot is possible.\n")
    }
    manteltest <- vegan::mantel(Dgen, Dgeo, na.rm = TRUE, permutations = 999)
    print(manteltest)
    if (plot) {
      if (!miss) {
        dens <- MASS::kde2d(as.numeric(Dgeo), as.numeric(Dgen), 
                            n = 300)
        myPal <- colorRampPalette(c("white", "blue", 
                                    "gold", "orange", "red"))
        plot(Dgeo, Dgen, pch = 20, cex = 0.8)
        image(dens, col = transp(myPal(300), 0.7), add = TRUE)
        points(Dgeo, Dgen, pch = 20, cex = 0.8)
        abline(lm(as.numeric(Dgen) ~ as.numeric(Dgeo)))
        title("Isolation by distance")
      }
      else {
        plot(Dgeo, Dgen)
        abline(lm(as.numeric(Dgen) ~ as.numeric(Dgeo)))
        title("Isolation by distance")
      }
    }
    out <- list(Dgen = Dgen, Dgeo = Dgeo, mantel = manteltest)
    return(out)
  }
}

gl_ibd<-gl.ibd_fixed(gl=gl_final_reduced)
gl_ibd

#####################################################
###  Calculating Genetic distance vs Phenotype    ###
#####################################################
#need to examine SVL, HW, and HL by locality because those are the ones that were significant 
#between high and low elevations

#remaining localities
#gh, mkb, lichw, gb, okh, dus, fin
which(gl_final_reduced@other$locality=="gh")
which(gl_final_reduced@other$locality=="mkb")
which(gl_final_reduced@other$locality=="lichw")
which(gl_final_reduced@other$locality=="gb")
which(gl_final_reduced@other$locality=="okh")
which(gl_final_reduced@other$locality=="dus")
which(gl_final_reduced@other$locality=="fin")


#traits
#this is just using samples of individuals sequenced
elev_svl<-gl_final_reduced@other$svl
elev_head_w<-gl_final_reduced@other$head_w
elev_head_l<-gl_final_reduced@other$head_l

#gh
avg_gh_SVL<-mean(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="gh"&gl_final_reduced@other$svl>65)],na.rm=T) #113 samples
max(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="gh"&gl_final_reduced@other$svl>65)])
avg_gh_head_w<-mean(elev_head_w$pop.data.Head_w[which(gl_final_reduced@other$locality=="gh"&gl_final_reduced@other$svl>65)],na.rm=T) #113 samples
avg_gh_head_l<-mean(elev_head_l$pop.data.Head_l[which(gl_final_reduced@other$locality=="gh"&gl_final_reduced@other$svl>65)],na.rm=T) #113 samples

#mkb
avg_mkb_SVL<-mean(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="mkb"&gl_final_reduced@other$svl>65)],na.rm=T) #12 samples
max(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="mkb"&gl_final_reduced@other$svl>65)],na.rm=T) #12 samples
avg_mkb_head_w<-mean(elev_head_w$pop.data.Head_w[which(gl_final_reduced@other$locality=="mkb"&gl_final_reduced@other$svl>65)],na.rm=T) #12 samples
avg_mkb_head_l<-mean(elev_head_l$pop.data.Head_l[which(gl_final_reduced@other$locality=="mkb"&gl_final_reduced@other$svl>65)],na.rm=T) #12 samples

#lichw
avg_lichw_SVL<-mean(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="lichw"&gl_final_reduced@other$svl>65)],na.rm=T) #10 samples
max(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="lichw"&gl_final_reduced@other$svl>65)]) #10 samples
avg_lichw_head_w<-mean(elev_head_w$pop.data.Head_w[which(gl_final_reduced@other$locality=="lichw"&gl_final_reduced@other$svl>65)],na.rm=T) #10 samples
avg_lichw_head_l<-mean(elev_head_l$pop.data.Head_l[which(gl_final_reduced@other$locality=="lichw"&gl_final_reduced@other$svl>65)],na.rm=T) #10 samples

#gb
avg_gb_SVL<-mean(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="gb"&gl_final_reduced@other$svl>65)],na.rm=T) #9 samples
max(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="gb"&gl_final_reduced@other$svl>65)]) #9 samples
avg_gb_head_w<-mean(elev_head_w$pop.data.Head_w[which(gl_final_reduced@other$locality=="gb"&gl_final_reduced@other$svl>65)],na.rm=T) #9 samples
avg_gb_head_l<-mean(elev_head_l$pop.data.Head_l[which(gl_final_reduced@other$locality=="gb"&gl_final_reduced@other$svl>65)],na.rm=T) #9 samples

#okh
avg_okh_SVL<-mean(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="okh"&gl_final_reduced@other$svl>65)],na.rm=T) #9 samples
max(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="okh"&gl_final_reduced@other$svl>65)]) #9 samples
avg_okh_head_w<-mean(elev_head_w$pop.data.Head_w[which(gl_final_reduced@other$locality=="okh"&gl_final_reduced@other$svl>65)],na.rm=T) #9 samples
avg_okh_head_l<-mean(elev_head_l$pop.data.Head_l[which(gl_final_reduced@other$locality=="okh"&gl_final_reduced@other$svl>65)],na.rm=T) #9 samples


#dus
avg_dus_SVL<-mean(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="dus"&gl_final_reduced@other$svl>65)],na.rm=T) #10 samples
max(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="dus"&gl_final_reduced@other$svl>65)]) #10 samples
avg_dus_head_w<-mean(elev_head_w$pop.data.Head_w[which(gl_final_reduced@other$locality=="dus"&gl_final_reduced@other$svl>65)],na.rm=T) #10 samples
avg_dus_head_l<-mean(elev_head_l$pop.data.Head_l[which(gl_final_reduced@other$locality=="dus"&gl_final_reduced@other$svl>65)],na.rm=T) #10 samples

#fin
avg_fin_SVL<-mean(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="fin"&gl_final_reduced@other$svl>65)],na.rm=T) #11 samples
max(elev_svl$pop.data.SVL[which(gl_final_reduced@other$locality=="fin"&gl_final_reduced@other$svl>65)]) #11 samples
avg_fin_head_w<-mean(elev_head_w$pop.data.Head_w[which(gl_final_reduced@other$locality=="fin"&gl_final_reduced@other$svl>65)],na.rm=T) #11 samples
avg_fin_head_l<-mean(elev_head_l$pop.data.Head_l[which(gl_final_reduced@other$locality=="fin"&gl_final_reduced@other$svl>65)],na.rm=T) #11 samples

####using all samples###
elev_svl<-pop.data$SVL
elev_head_w<-pop.data$Head_l
elev_head_l<-pop.data$Head_w

#gh
avg_gh_SVL<-mean(elev_svl[which(pop.data$Locality=="gh"&pop.data$SVL>65)],na.rm=T) #113 samples
avg_gh_head_w<-mean(elev_head_w[which(pop.data$Locality=="gh"&pop.data$SVL>65)],na.rm=T) #113 samples
avg_gh_head_l<-mean(elev_head_l[which(pop.data$Locality=="gh"&pop.data$SVL>65)],na.rm=T) #113 samples

#mkb
avg_mkb_SVL<-mean(elev_svl[which(pop.data$Locality=="mkb"&pop.data$SVL>65)],na.rm=T) #12 samples
avg_mkb_head_w<-mean(elev_head_w[which(pop.data$Locality=="mkb"&pop.data$SVL>65)],na.rm=T) #12 samples
avg_mkb_head_l<-mean(elev_head_l[which(pop.data$Locality=="mkb"&pop.data$SVL>65)],na.rm=T) #12 samples

#lichw
avg_lichw_SVL<-mean(elev_svl[which(pop.data$Locality=="lichw"&pop.data$SVL>65)],na.rm=T) #10 samples
avg_lichw_head_w<-mean(elev_head_w[which(pop.data$Locality=="lichw"&pop.data$SVL>65)],na.rm=T) #10 samples
avg_lichw_head_l<-mean(elev_head_l[which(pop.data$Locality=="lichw"&pop.data$SVL>65)],na.rm=T) #10 samples

#gb
avg_gb_SVL<-mean(elev_svl[which(pop.data$Locality=="gb"&pop.data$SVL>65)],na.rm=T) #9 samples
avg_gb_head_w<-mean(elev_head_w[which(pop.data$Locality=="gb"&pop.data$SVL>65)],na.rm=T) #9 samples
avg_gb_head_l<-mean(elev_head_l[which(pop.data$Locality=="gb"&pop.data$SVL>65)],na.rm=T) #9 samples

#okh
avg_okh_SVL<-mean(elev_svl[which(pop.data$Locality=="okh"&pop.data$SVL>65)],na.rm=T) #9 samples
avg_okh_head_w<-mean(elev_head_w[which(pop.data$Locality=="okh"&pop.data$SVL>65)],na.rm=T) #9 samples
avg_okh_head_l<-mean(elev_head_l[which(pop.data$Locality=="okh"&pop.data$SVL>65)],na.rm=T) #9 samples


#dus
avg_dus_SVL<-mean(elev_svl[which(pop.data$Locality=="dus"&pop.data$SVL>65)],na.rm=T) #10 samples
avg_dus_head_w<-mean(elev_head_w[which(pop.data$Locality=="dus"&pop.data$SVL>65)],na.rm=T) #10 samples
avg_dus_head_l<-mean(elev_head_l[which(pop.data$Locality=="dus"&pop.data$SVL>65)],na.rm=T) #10 samples

#fin
avg_fin_SVL<-mean(elev_svl[which(pop.data$Locality=="fin"&pop.data$SVL>65)],na.rm=T) #11 samples
avg_fin_head_w<-mean(elev_head_w[which(pop.data$Locality=="fin"&pop.data$SVL>65)],na.rm=T) #11 samples
avg_fin_head_l<-mean(elev_head_l[which(pop.data$Locality=="fin"&pop.data$SVL>65)],na.rm=T) #11 samples



#order must be in the following to be compared directly with the Fst values:
#gh, mkb, lichw, okh, gb, dus, fin.
#missing are kwh and elis because too few samples and they were filtered out.

#dist matrix of average measurements by locality
D_svl<-as.matrix(c(avg_gh_SVL, avg_mkb_SVL, avg_lichw_SVL, avg_okh_SVL, avg_gb_SVL, avg_dus_SVL, avg_fin_SVL))
row.names(D_svl)<-c("gh", "mkb", "lichw", "okh", "gb", "dus", "fin")
D_svl<-dist(D_svl)

D_head_w<-as.matrix(c(avg_gh_head_w, avg_mkb_head_w, avg_lichw_head_w, avg_okh_head_w, avg_gb_head_w, avg_dus_head_w, avg_fin_head_w))
row.names(D_head_w)<-c("gh", "mkb", "lichw", "okh", "gb", "dus", "fin")
D_head_w<-dist(D_head_w)

D_head_l<-as.matrix(c(avg_gh_head_l, avg_mkb_head_l, avg_lichw_head_l, avg_okh_head_l, avg_gb_head_l, avg_dus_head_l, avg_fin_head_l))
row.names(D_head_l)<-c("gh", "mkb", "lichw", "okh", "gb", "dus", "fin")
D_head_l<-dist(D_head_l)

Fst_matrix<-as.dist(gl_final_reduced_fst$Fsts)

library("ade4")
set.seed(1000)
mantel.rtest(Fst_matrix,D_svl)
set.seed(1000)
mantel.rtest(D_svl,Fst_matrix)
plot(c(Fst_matrix),c(D_svl))
plot(r1 <- mantel.rtest(D_svl,Fst_matrix), main = "Mantel's test")

mantel.rtest(Fst_matrix,D_head_l)
mantel.rtest(Fst_matrix,D_head_w)

##########################################################################
### Comparing allele frequencies vs morphological traits by individual ###
##########################################################################

#SVL
svl<-as.data.frame(gl_final_reduced@other$svl[c(which(gl_final_reduced@other$svl>65)),])
rownames(svl)<-c(gl_final_reduced@ind.names[c(which(gl_final_reduced@other$svl>65))])
colnames(svl)<-"svl"
complete_svl_names<-rownames(svl)[c(which(complete.cases(svl)==T))]
svl_data<-as.matrix(svl$svl[c(which(complete.cases(svl)==T))])
row.names(svl_data)<-c(complete_svl_names)
length(svl_data)#127 indviduals with SNPs and SVL data
D_svl<-dist(svl_data)

#Head_w
head_w<-as.data.frame(c((gl_final_reduced@other$head_w[which(gl_final_reduced@other$svl>65),])))
rownames(head_w)<-c(gl_final_reduced@ind.names[which(gl_final_reduced@other$svl>65)])
colnames(head_w)<-"head_w"
complete_head_w_names<-rownames(head_w)[c(which(complete.cases(head_w)==T))]
head_w_data<-as.matrix(head_w$head_w[c(which(complete.cases(head_w)==T))])
row.names(head_w_data)<-c(complete_head_w_names)
length(head_w_data)#124 indviduals with SNPs and head width data
D_head_w<-dist(head_w_data)

#Head_l
head_l<-as.data.frame(c(gl_final_reduced@other$head_l[which(gl_final_reduced@other$svl>65),]))
rownames(head_l)<-c(gl_final_reduced@ind.names[which(gl_final_reduced@other$svl>65)])
colnames(head_l)<-"head_l"
complete_head_l_names<-rownames(head_l)[c(which(complete.cases(head_l)==T))]
head_l_data<-as.matrix(head_l$head_l[c(which(complete.cases(head_l)==T))])
row.names(head_l_data)<-c(complete_head_l_names)
length(head_l_data)#125 indviduals with SNPs and SVL data
D_head_l<-dist(head_l_data)

#making allel frequency  matrix
gl_final_reduced@gen
#gl.write.csv(gl_final_reduced,"alelle_freq.csv",outpath = getwd())
allele_freq<-read.csv("alelle_freq.csv")
#subset columns we want to use
allele_freq<-allele_freq[,c(17:ncol(allele_freq))]
corrected_names<-substr(colnames(allele_freq),2,25)
#manually fix errors
corrected_names[78]<-"napu1_1"
corrected_names[154]<-"napu2_1"
corrected_names[173]<-"napu3_1"
corrected_names[174]<-"napu4_1"
#reassign the names
colnames(allele_freq)<-corrected_names
#all genotypes with corrected names
head(allele_freq)
#remove locality information
allele_freq<-allele_freq[-1,]
AF2<-t(allele_freq)
rownames(AF2)<-corrected_names
#This is the entire allele frequency matrix that we will use dist on. but first... 
#subset out only the ones with the correct names per phenotype
AF2_svl<-AF2[which(rownames(AF2)%in%rownames(svl_data)),]
AF2_head_w<-AF2[which(rownames(AF2)%in%rownames(head_w_data)),]
AF2_head_l<-AF2[which(rownames(AF2)%in%rownames(head_l_data)),]
#verify that they are the same
nrow(AF2_svl)==nrow(svl_data) #127
nrow(AF2_head_w)==nrow(head_w_data) #124
nrow(AF2_head_l)==nrow(head_l_data) #125
#make distance matrices
D_allele_4_svl<-dist(AF2_svl)
D_allele_4_head_w<-dist(AF2_head_w)
D_allele_4_head_l<-dist(AF2_head_l)

#Run mantel's test on variables
set.seed(999)
man_svl<-mantel.rtest(D_allele_4_svl,D_svl)
man_head_w<-mantel.rtest(D_allele_4_head_w,D_head_w)
man_head_l<-mantel.rtest(D_allele_4_head_l,D_head_l)

#Notes
#Observation = Mantel's R
man_svl
man_head_w
man_head_l

mantel.rtest(gl_ibd$Dgeo,gl_ibd$Dgen)
Dgeo<-as.matrix(gl_ibd$Dgeo)
Dgen<-as.matrix(gl_ibd$Dgen)
plot(gl_ibd$Dgeo,gl_ibd$Dgen)
abline(lm(gl_ibd$Dgen~gl_ibd$Dgeo))
gl_ibd


#######################################
###     Check Head H and Tail       ###
#######################################

#Head_h
head_h<-as.data.frame(c(gl_final_reduced@other$head_h[which(gl_final_reduced@other$svl>65),]))
rownames(head_h)<-c(gl_final_reduced@ind.names[which(gl_final_reduced@other$svl>65)])
colnames(head_h)<-"head_h"
complete_head_h_names<-rownames(head_h)[c(which(complete.cases(head_h)==T))]
head_h_data<-as.matrix(head_h$head_h[c(which(complete.cases(head_h)==T))])
row.names(head_h_data)<-c(complete_head_h_names)
length(head_h_data)#131 indviduals with SNPs and HH data
D_head_h<-dist(head_h_data)

#Tail
tail<-as.data.frame(c(gl_final_reduced@other$tail[which(gl_final_reduced@other$svl>65),]))
rownames(tail)<-c(gl_final_reduced@ind.names[which(gl_final_reduced@other$svl>65)])
colnames(tail)<-"tail"
complete_tail_names<-rownames(tail)[c(which(complete.cases(tail)==T))]
tail_data<-as.matrix(tail$tail[c(which(complete.cases(tail)==T))])
row.names(tail_data)<-c(complete_tail_names)
length(tail_data)#59 indviduals with SNPs and Tail data
D_tail<-dist(tail_data)

#This is the entire allele frequency matrix that we will use dist on. but first... 
#subset out only the ones with the correct names per phenotype
AF2_head_h<-AF2[which(rownames(AF2)%in%rownames(head_h_data)),]
AF2_tail<-AF2[which(rownames(AF2)%in%rownames(tail_data)),]
#verify that they are the same
nrow(AF2_head_h)==nrow(head_h_data) #124
nrow(AF2_tail)==nrow(tail_data) #125
#make distance matrices
D_allele_4_head_h<-dist(AF2_head_h)
D_allele_4_tail<-dist(AF2_tail)

#Run mantel's test on variables
set.seed(999)
man_head_h<-mantel.rtest(D_allele_4_head_h,D_head_h)
man_tail<-mantel.rtest(D_allele_4_tail,D_tail)

#Notes
#Observation = Mantel's R
man_head_h
man_tail

#######################################
###           Outliers              ###
#######################################
library("devtools")
library("BiocInstaller")
library("pcadapt") #install.packages("pcadapt")
library("qvalue") #BiocManager::install("whitlock/OutFLANK")
library("OutFLANK") #BiocManager::install("whitlock/OutFLANK")
library("ggplot2")

#unable to install bioclite - depreciated because installs bad stuff..
source("http://bioconductor.org/biocLite.R")
BiocManager::valid("biocLite.R")
BiocManager::install("biocLite.R")

#SECTION 1: LOAD THE DATA
sel <- read.table("data/SNPselection1.txt", head = TRUE)
dim(sel)

#######################################
###              PRDA              ###
#######################################


devtools::install_github("bcm-uga/TESS3_encho_sen")


  
#######################################
###              AMOVA              ###
#######################################
PrDA

gl_final_reduced_AMOVA<-gl.amova(gl_final_reduced,permutations=100)
gl_final_reduced_AMOVA


###############################################
### Exporting files for FastStructure file  ###
###############################################
#Make a population file based on the Genlight object being exported. 
#This file is used for the visualizations of the fastStructure analysis using distruct.py
write.table(gl_final_reduced@pop, file = "napu_pop_reduced.txt")

#write out structure file. 
gl2faststructure(x = gl_final_reduced, outfile = "napu_reduced.str", outpath = getwd())

#######################################################
###  Using DAPC for alternative clustering method  ###
#######################################################

species.genlight<-gl_final_reduced

#Cluster
clustersoverall<-find.clusters(species.genlight,max.n.clust=20,n.iter=1e6,n.start=1000,n.pca=100,choose.n.clust=F)
plot(clustersoverall$Kstat,type="o",xlab="number of clusters (K)", ylab="BIC",col="blue",main="Detection based on BIC")

#How many Clusters?
order(clustersoverall$Kstat)

#Discriminant Analyses of Principal Compentent Analyses (DAPC)
dapc.napu<-dapc(species.genlight,clustersoverall$grp)
dapc.napu2<-dapc(species.genlight,clustersoverall$grp,n.pca = 40, n.da = 40,var.contrib = T,scale = F,pop = species.genlight$pop)

scatter(dapc.napu,cell=0, pch=18:23,cstar=0,mstree=T,lwd =2, lty=2)
scatter(dapc.napu2,cell=0, pch=18:23,cstar=0,mstree=T,lwd =2, lty=2)

library(varhandle)
compoplot(dapc.napu,legend=F,lab=unfactor(species.genlight@pop),show.lab=T)
compoplot(dapc.napu2,legend=T,lab=indNames(species.genlight))
dev.off()

####################################
###### RELATEDNESS ANALYSIS ########
####################################
#Installation of SNPRelate - https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#installation-of-the-package-snprelate
library(gdsfmt)
library(SNPRelate)
library(SeqArray)

#Export data to .gds file (uses the faststructure pkg but may need to use snpgdsVCF2GDS() function)
setwd("/Users/jonathandeboer727/Desktop/cordylid_ddRADSeq/Napu_gds")
gl2gds(x = gl_final_reduced, outfile = "napu_rel.gds", outpath = "/Users/jonathandeboer727/Desktop/cordylid_ddRADSeq/Napu_gds")
snpgdsSummary("napu_rel.gds")

#For relatedness analysis, identity-by-descent (IBD) estimation 
#in SNPRelate using moments (MoM) (Purcell et al., 2007) 
napufile <- snpgdsOpen("napu_rel.gds",allow.duplicate = T)
napu_snp.id<-read.gdsn(index.gdsn(napufile,"snp.id"))
napu_sample.id<-read.gdsn(index.gdsn(napufile,"sample.id"))

#use the method of moments to find identity by descent.
rel_napu <- snpgdsIBDMoM(napufile, sample.id=napu_sample.id, snp.id=napu_snp.id,
                          num.thread=2,autosome.only = F,remove.monosnp = F)
#WITHOUT cutoff
ibd.coeff <- snpgdsIBDSelection(rel_napu)
head(ibd.coeff)
#how many kinship relationships?
length(which(ibd.coeff$kinship>0))

#plot relatedness
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="Napu samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

##### WITH CUTOFF
ibd.coeff2 <- snpgdsIBDSelection(rel_napu,kinship.cutoff=1/32)
head(ibd.coeff2)
#how many kinship relationships?
length(which(ibd.coeff2$kinship>0.25))
ibd.coeff2[which(ibd.coeff2$kinship>0.25),]

#plot relatedness
plot(ibd.coeff2$k0, ibd.coeff2$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="Napu samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)




#gl2gds(gl13,outfile = "allsites_allyear.glds",outpath="G://Postdoc//pygmy rabbits//3RAD_output//3RAD_all_pools_min400_outfiles//",)
#gdsallyear<-snpgdsOpen("allsites_allyear.glds",allow.duplicate = T)
#snp.id<-read.gdsn(index.gdsn(gdsallyear,"snp.id"))
#sample.id<-read.gdsn(index.gdsn(gdsallyear,"sample.id"))
#relallyear <- snpgdsIBDMoM(gdsallyear, sample.id=sample.id, snp.id=snp.id,
#                           maf=0.01, missing.rate=0.20, num.thread=2,autosome.only = F,remove.monosnp = F)
#relallyear.coeff <- snpgdsIBDSelection(relallyear)

# gh samples only?
sample.id <- read.gdsn(index.gdsn(napufile, "sample.id"))
gh.id <- sample.id[napu_pop_code == "gh"]



#############################
###### snpgdsIBDKING ########
#############################

######NATHAN, what about the family ID thing?

##### WITH CUT OFF
napu_robust<-snpgdsIBDKING(napufile,sample.id=sample.id,autosome.only = F)
names(napu_robust)
dat<-snpgdsIBDSelection(napu_robust,kinship.cutoff=1/32)
head(dat)
length(which(dat$kinship>0.25))
plot(dat$IBS0, dat$kinship, xlab="Proportion of Zero IBS",
     ylab="Estimated Kinship Coefficient (KING-robust)")

##### WITHOUT CUTOFF
king_napu<-snpgdsIBDKING(napufile,sample.id=sample.id,snp.id=snp.id,
                         num.thread=2,autosome.only = F,remove.monosnp = F)
king.coeff <- snpgdsIBDSelection(king_napu)

head(king.coeff)
#how many kinship relationships?
length(which(king.coeff$kinship>0.25))
length(which(king.coeff$kinship<0.25&king.coeff$kinship>0.20))
length(which(king.coeff$kinship==0))
length(which(king.coeff$kinship<0))

#plot relatedness
plot(king.coeff$IBS0, king.coeff$kinship, xlab="Proportion of Zero IBS",
     ylab="Estimated Kinship Coefficient (KING-robust)")

##### AUTOSOME.ONLY = T
{
  king_napu2<-snpgdsIBDKING(napufile,sample.id=sample.id,snp.id=snp.id,
                          num.thread=2,autosome.only = T,remove.monosnp = F)
king.coeff2 <- snpgdsIBDSelection(king_napu2)
head(king.coeff2)
#how many kinship relationships?
length(which(king.coeff2$kinship>0))
length(which(king.coeff2$kinship==0))
length(which(king.coeff2$kinship<0))

#plot relatedness
plot(king.coeff2$IBS0, king.coeff2$kinship, xlab="Proportion of Zero IBS",
     ylab="Estimated Kinship Coefficient (KING-robust)")
}


######################
###  ResistanceGA  ###
######################
library(ResistanceGA)

 devtools::install_github("wpeterman/ResistanceGA", build_vignettes = TRUE)


#something about Julia...

#making directories
if("ResistanceGA_Examples"%in%dir("C:/")==FALSE)
  dir.create(file.path("~/ResistanceGA_Examples"))
# Create a subdirectory for the first example
setwd("~/ResistanceGA_Examples")
dir.create(file.path("~/ResistanceGA_Examples/","SingleSurface"))
# Directory to write .asc files and results
write.dir <- "~/ResistanceGA_Examples/SingleSurface/"
# Give path to CIRCUITSCAPE .exe file
# Default = '"C:/Program Files/Circuitscape/cs_run.exe"'
CS.program <- paste('"~/Program Files/Circuitscape/cs_run.exe"')

#Load resistance surfaces and export as .asc file for use with CIRCUITSCAPE. These surfaces were made
#using the RandomFields package
data(resistance_surfaces)
continuous <- resistance_surfaces[[2]]
writeRaster(continuous,
            filename = paste0(write.dir,"cont.asc"),
            overwrite = TRUE)
#Load the example sample location data and export as .txt file. 
#This is formatted for input into CIRCUITSCAPE
data(samples)
write.table(samples,file=paste0(write.dir,"samples.txt"),sep="\t",col.names=F,row.names=F)
# Create a spatial points object for plotting
sample.locales <- SpatialPoints(samples[,c(2,3)])

#Plot surface and overlay the sample points
plot(continuous)
plot(sample.locales, pch = 16, col = "blue", add = TRUE) # Add points

#####Prepare data for optimization#####
#Run the GA.prep and CS.prep functions
# Set the random number seed to reproduce the results presented
GA.inputs <- GA.prep(ASCII.dir = write.dir,
                     max.cat = 500,
                     max.cont = 500,
                     select.trans = "M",
                     method = "LL",
                     seed = 555)
CS.inputs <- CS.prep(n.Pops = length(sample.locales),
                     CS_Point.File = paste0(write.dir,"samples.txt"),
                     CS.program = CS.program)


###################################################
###  Using Pophelper to visulize the Structure files ###
###################################################
#install_github('royfrancis/pophelper')
library(pophelper)

setwd("~/proj/fastStructure/napu_faststructure_final")

#Take only the mean Q's. I used K = 1:10. It was not working with only K = 2:10. 
sfiles <- list.files(path="~/proj/fastStructure/napu_faststructure_final/meanQs", full.names=T)
sfiles

slist <- readQ(sfiles)
slist[[3]][1,]
slist1<- alignK(slist[c(3,6,10)])

#This is K=2, 5, and 9
p1 <- plotQ(slist[c(3,6,10)],imgoutput="join",returnplot=T,exportplot=F,basesize=11)
#need to force to use gridExtra package
gridExtra::grid.arrange(p1$plot[[1]])


p1 <- plotQ(slist1,imgoutput="join",returnplot=T,exportplot=F,basesize=11)
gridExtra::grid.arrange(p1$plot[[1]])


#######################################################
### Using tess3r for alternative clustering method  ###
#######################################################
install.packages("tess3")
devtools::install_github("bcm-uga/TESS3_encho_sen")

library(tess3r)

##############################
### Calculating Gene Flow  ###
##############################

##### Nathan said individual based approaches. mantel correlegrams, characterize potenital space of dispersal. 

install.packages("admixr")
library(admixr) #qp3Pop, qpDstat, qpF4ratio, qpAdm, qpWave
library(qp3pop)
library(tidyverse)

#test it
prefix <- download_data(dirname = tempdir())
list.files(path = dirname(prefix), pattern = basename(prefix), full.names = TRUE)
snps<-eigenstrat(prefix)
pops <- c("French", "Sardinian", "Han", "Papuan", "Khomani_San", "Mbuti", "Dinka")
result <- d(W = pops, X = "Yoruba", Y = "Vindija", Z = "Chimp", data = snps)

#my turn
#1 get files into eigenstrat

#make pops list
pops<-c("dus","fin","gb","gh","lichw","mkb","okh")

#calculate D-statistic
dstats <- d(data=snps, W=pops)

#############EXCESS CODE
#can also run glPca - automated way of running principal components analysis for genlights.
#it takes FOREVER though.
#### THIS IS ON SECOND SCRIPT
pca1<-glPca(gl4)
pca1
scatter(pca1,posi="topleft")

#NJ - check with Kajo for crazy outliers?
library(ape)
tre <- njs(dist(as.matrix(gl13)))
plot(tre,typ="fan")

gl.tree.nj(gl13,type = "phylogram")


#PC using color
myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
#add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)
plot(tre, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=myCol, cex=4)


gl.ibd_fixed2<-function(gl = NULL, Dgen = NULL, Dgeo = NULL, projected = FALSE,
                     permutations = 1000, plot = TRUE)
{
  if (!(requireNamespace("dismo", quietly = TRUE))) {
    stop("Package dismo needed for this function to work. Please install it.")
  }
  else {
    if (!is.null(Dgen) & !is.null(Dgeo))
      cat("Analysis performed on provided genetic and Euclidean distance matrices.")
    if (class(gl) == "genlight") {
      cat("Standard analysis performed on the genlight object. Mantel test and plot will be Fst/1-Fst versus log(distance)\n")
      if (nrow(gl@other$latlong) != nInd(gl))
        stop("Cannot find coordinates for each individual in slot @other$latlong")
      if (sum(match(names(gl@other$latlong), "long"), na.rm = T) == 1)
        gl@other$latlong$lon <- gl@other$latlong$long
      if (!projected) {
        xy <- dismo::Mercator(gl@other$latlong)#[, c("lon", "lat")])
        cat("Coordinates transformed to Mercator (google) projection to calculate distances in meters.\n")
      }
      else {
        xy = gl@other$latlong[, c("lon", "lat")]
        cat("Coordinates not transformed. Distances calculated on the provided coordinates.")
      }
      pop.xy <- apply(xy, 2, function(a) tapply(a, pop(gl),
                                                mean, na.rm = T))
      Dgeo <- dist(pop.xy)
      Dgeo <- log(Dgeo)
      Dgen <- as.dist(StAMPP::stamppFst(gl, nboots = 1))
      Dgen <- Dgen/(1 - Dgen)
      ordering <- levels(pop(gl))
      Dgen <- as.dist(as.matrix(Dgen)[ordering, ordering])
      Dgeo <- as.dist(as.matrix(Dgeo)[ordering, ordering])
    }
    miss = FALSE
    if (sum(is.na(Dgen)) > 0 | sum(is.infinite(Dgen)) > 0 ) {
      miss = TRUE
      cat("There are missing values in the genetic distance matrix. No kernel distance plot is possible.\n")
    }
    if (sum(is.na(Dgeo)) > 0 | sum(is.infinite(Dgeo)) > 0 ) {
      miss = TRUE
      cat("There are missing values in the geographic distance matrix. No kernel distance plot is possible.\n")
    }
    manteltest <- vegan::mantel(Dgen, Dgeo, na.rm = TRUE, permutations = 999)
    print(manteltest)
    if (plot) {
      if (!miss) {
        dens <- MASS::kde2d(as.numeric(Dgeo), as.numeric(Dgen), 
                            n = 300)
        myPal <- colorRampPalette(c("white", "blue", 
                                    "gold", "orange", "red"))
        plot(Dgeo, Dgen, pch = 20, cex = 0.8)
        #image(dens, col = transp(myPal(300), 0.7), add = TRUE)
        points(Dgeo, Dgen, pch = 20, cex = 2.0)
        abline(lm(as.numeric(Dgen) ~ as.numeric(Dgeo)))
        title("Isolation by distance")
      }
      else {
        plot(Dgeo, Dgen)
        abline(lm(as.numeric(Dgen) ~ as.numeric(Dgeo)))
        title("Isolation by distance")
      }
    }
    out <- list(Dgen = Dgen, Dgeo = Dgeo, mantel = manteltest)
    return(out)
  }
}
gl_ibd2<-gl.ibd_fixed2(gl=gl_final_reduced)
gl_ibd
