library(methylKit)

####Male DMR Analyses####
###50 ug/kg BPA vs. Corn Oil Analyses###
##Preparing the methylation dataset for 50 ug/kg BPA vs. Corn Oil comparisons
file.list=list(file.path("./data/methylation/33M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/14M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/49M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/50M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/52M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/5M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/69M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/73M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/84M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/93M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/9M1.sorted.dedup_resorted.bismark.cov"))


myobj_50BPA_M = methRead(file.list,
                         sample.id=list("33M1","14M1","49M1","50M1","52M1","5M1","69M1","73M1","84M1","93M1","9M1" ),
                         assembly="rn7",
                         treatment=c(0,1,0,0,0,1,1,1,0,1,0), #0 is Corn Oil, 1 is prenatal BP
                         context="CpG",
                         mincov = 1,
                         pipeline = "bismarkCoverage")

#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_50BPA_M,win.size=500,step.size=500,cov.bases = 10)
meth=unite(tiles,min.per.group=4L)

##Calculating DMRs without licking/grooming covariate
covariates=data.frame(batch=c(3,4,2,3,4,1,5,3,3,5,2))

myDiff=calculateDiffMeth(meth,
                         overdispersion = "MN",
                         covariates=covariates,
                         adjust="SLIM",
                         test="Chisq")

#Extract DMRs with q-value < 0.05 and >5% methylation difference
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs

#Save output file for assessing overlap in DMRs and gene annotation
write.csv(myDiff5p, file="DMR_50BPA_M.csv", quote=FALSE)

##Calculating DMRs with licking/grooming covariate
covariates=data.frame(LG=c(-0.0793,1.3241,0.6706,1.5728,-0.1859,0.5683,-1.1157,-0.0222,-0.4216,0.3745,-1.8635),batch=c(3,4,2,3,4,1,5,3,3,5,2))

myDiff=calculateDiffMeth(meth,
                         overdispersion = "MN",
                         covariates=covariates,
                         adjust="SLIM",
                         test="Chisq")

#Extract DMRs with q-value < 0.05 and >5% methylation difference
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs

#Save output file for assessing overlap in DMRs and gene annotation
write.csv(myDiff5p, file="DMR_50BPA_M_LG-Covariate.csv", quote=FALSE)

###50 ug/kg Mixed BP vs. Corn Oil Analyses###
##Preparing the methylation dataset for 50 ug/kg Mixed BP vs. Corn Oil comparisons
file.list=list(file.path("./data/methylation/33M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/49M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/50M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/103M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/38M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/84M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/39M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/53M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/52M1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/9M1.sorted.dedup_resorted.bismark.cov"))


myobj_50BP_M = methRead(file.list,
                        sample.id=list("33M1","49M1","50M1","103M1","38M1","84M1","39M1","53M1","52M1","9M1" ),
                        assembly="rn7",
                        treatment=c(0,0,0,1,1,0,1,1,0,0), #0 is Corn Oil, 1 is prenatal BP
                        context="CpG",
                        mincov = 1,
                        pipeline = "bismarkCoverage")

#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_50BP_M,win.size=500,step.size=500,cov.bases = 10)
meth=unite(tiles,min.per.group=4L)

##Calculating DMRs without licking/grooming covariate
covariates=data.frame(batch=c(3,2,3,5,3,3,2,3,4,2))

myDiff=calculateDiffMeth(meth,
                         overdispersion = "MN",
                         covariates=covariates,
                         adjust="SLIM",
                         test="Chisq")

#Extract DMRs with q-value < 0.05 and >5% methylation difference
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs

#Save output file for assessing overlap in DMRs and gene annotation
write.csv(myDiff5p, file="DMR_50BP_M.csv", quote=FALSE)

##Calculating DMRs with licking/grooming covariate
covariates=data.frame(LG=c(-0.0793,0.6706,1.5728,-1.1313,0.4293,-0.4216,0.5953,-1.9328,-0.1859,-1.8635),batch=c(3,2,3,5,3,3,2,3,4,2))

myDiff=calculateDiffMeth(meth,
                         overdispersion = "MN",
                         covariates=covariates,
                         adjust="SLIM",
                         test="Chisq")

#Extract DMRs with q-value < 0.05 and >5% methylation difference
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs

#Save output file for assessing overlap in DMRs and gene annotation
write.csv(myDiff5p, file="DMR_50BP_M_LG-Covariate.csv", quote=FALSE)

###150 ug/kg Mixed BP vs. Corn Oil Analyses###
##Preparing the methylation dataset for 150 ug/kg Mixed BP vs. Corn Oil comparisons
 file.list=list(file.path("./data/methylation/33M1.sorted.dedup_resorted.bismark.cov"),
                file.path("./data/methylation/46M1.sorted.dedup_resorted.bismark.cov"),
                file.path("./data/methylation/49M1.sorted.dedup_resorted.bismark.cov"),
                file.path("./data/methylation/50M1.sorted.dedup_resorted.bismark.cov"),
                file.path("./data/methylation/52M1.sorted.dedup_resorted.bismark.cov"),
                file.path("./data/methylation/61M1.sorted.dedup_resorted.bismark.cov"),
                file.path("./data/methylation/62M1.sorted.dedup_resorted.bismark.cov"),
                file.path("./data/methylation/63M1.sorted.dedup_resorted.bismark.cov"),
                file.path("./data/methylation/84M1.sorted.dedup_resorted.bismark.cov"),
                file.path("./data/methylation/87M1.sorted.dedup_resorted.bismark.cov"),
                file.path("./data/methylation/91M1.sorted.dedup_resorted.bismark.cov"),
                file.path("./data/methylation/9M1.sorted.dedup_resorted.bismark.cov"))
 
 
 myobj_150BP_M = methRead(file.list,
                          sample.id=list("33M1","46M1","49M1","50M1","52M1","61M1","62M1","63M1","84M1","87M1","91M1","9M1" ),
                          assembly="rn7",
                          treatment=c(0,1,0,0,0,1,1,1,0,1,1,0), #0 is Corn Oil, 1 is prenatal BP
                          context="CpG",
                          mincov = 1,
                          pipeline = "bismarkCoverage")
 
 #Creating regions with 500 bp length
 tiles = tileMethylCounts(myobj_150BP_M,win.size=500,step.size=500,cov.bases = 10)
 meth=unite(tiles,min.per.group=4L)
 
 ##Calculating DMRs without licking/grooming covariate
 covariates=data.frame(batch=c(3,1,2,3,4,1,3,4,3,2,5,2))
 
 myDiff=calculateDiffMeth(meth,
                          overdispersion = "MN",
                          covariates=covariates,
                          adjust="SLIM",
                          test="Chisq")
 
 #Extract DMRs with q-value < 0.05 and >5% methylation difference
 myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
 getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
 getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs
 
 #Save output file for assessing overlap in DMRs and gene annotation
 write.csv(myDiff5p, file="DMR_150BP_M.csv", quote=FALSE)
 
 ##Calculating DMRs with licking/grooming covariate
 covariates=data.frame(LG=c(-0.0793,0.0820,0.6706,1.5728,-0.1859,0.6268,0.1139,-2.0646,-0.4216,-0.88954,0.1461,-1.8635),batch=c(3,1,2,3,4,1,3,4,3,2,5,2))
 
 myDiff=calculateDiffMeth(meth,
                          overdispersion = "MN",
                          covariates=covariates,
                          adjust="SLIM",
                          test="Chisq")
 
 #Extract DMRs with q-value < 0.05 and >5% methylation difference
 myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
 getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
 getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs
 
 #Save output file for assessing overlap in DMRs and gene annotation
 write.csv(myDiff5p, file="DMR_150BP_M_LG-Covariate.csv", quote=FALSE)
 