library(methylKit)

####Female DMR Analyses####
###50 ug/kg BPA vs. Corn Oil Analyses###
##Preparing the methylation dataset for 50 ug/kg BPA vs. Corn Oil comparisons
file.list=list(file.path("./data/methylation/2F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/33F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/49F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/50F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/14F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/5F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/73F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/84F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/93F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/9F1.sorted.dedup_resorted.bismark.cov"))


myobj_50BPA_F = methRead(file.list,
                         sample.id=list("2F1","33F1","49F1","50F1","14F1","5F1","73F1","84F1","93F1","9F1" ),
                         assembly="rn7",
                         treatment=c(0,0,0,0,1,1,1,0,1,0), #0 is Corn Oil, 1 is prenatal BP
                         context="CpG",
                         mincov = 1,
                         pipeline = "bismarkCoverage")

#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_50BPA_F,win.size=500,step.size=500,cov.bases = 10)
meth=unite(tiles,min.per.group=4L)

##Calculating DMRs without licking/grooming covariate
covariates=data.frame(batch=c(1,3,2,1,2,4,5,5,5,2))

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
write.csv(myDiff5p, file="DMR_50BPA_F.csv", quote=FALSE)

##Calculating DMRs with licking/grooming covariate
covariates=data.frame(LG=c(0.0777,0.401,1.338,2.466,2.155,1.210,0.472,-0.0271,0.968,-1.829),batch=c(1,3,2,1,2,4,5,5,5,2))

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
write.csv(myDiff5p, file="DMR_50BPA_F_LG-Covariate.csv", quote=FALSE)

###50 ug/kg Mixed BP vs. Corn Oil Analyses###
##Preparing the methylation dataset for 50 ug/kg Mixed BP vs. Corn Oil comparisons

file.list=list(file.path("./data/methylation/2F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/33F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/49F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/50F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/100F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/103F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/38F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/84F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/39F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/53F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/55F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/9F1.sorted.dedup_resorted.bismark.cov"))


myobj_50BP_F = methRead(file.list,
                        sample.id=list("2F1","33F1","49F1","50F1","100F1","103F1","38F1","84F1","39F1","53F1","55F1","9F1" ),
                        assembly="rn7",
                        treatment=c(0,0,0,0,1,1,1,0,1,1,1,0), #0 is Corn Oil, 1 is prenatal BP
                        context="CpG",
                        mincov = 1,
                        pipeline = "bismarkCoverage")

#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_50BP_F,win.size=500,step.size=500,cov.bases = 10)
meth=unite(tiles,min.per.group=4L)

##Calculating DMRs without licking/grooming covariate
covariates=data.frame(batch=c(1,3,2,1,5,5,4,5,2,3,2,2))

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
write.csv(myDiff5p, file="DMR_50BP_F.csv", quote=FALSE)

##Calculating DMRs with licking/grooming covariate
covariates=data.frame(LG=c(0.0777,0.401,1.338,2.466,0.276,-0.914,1.037,-0.0271,1.244,-1.916,0.0741,-1.829),batch=c(1,3,2,1,5,5,4,5,2,3,2,2))

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
write.csv(myDiff5p, file="DMR_50BP_F_LG-Covariate.csv", quote=FALSE)

###150 ug/kg Mixed BP vs. Corn Oil Analyses###
##Preparing the methylation dataset for 150 ug/kg Mixed BP vs. Corn Oil comparisons
file.list=list(file.path("./data/methylation/2F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/33F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/49F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/50F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/61F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/62F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/77F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/84F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/87F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/88F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/89F1.sorted.dedup_resorted.bismark.cov"),
               file.path("./data/methylation/9F1.sorted.dedup_resorted.bismark.cov"))

myobj_150BP_F = methRead(file.list,
                         sample.id=list("2F1","33F1","49F1","50F1","61F1","62F1","77F1","84F1","87F1","88F1","89F1","9F1" ),
                         assembly="rn7",
                         treatment=c(0,0,0,0,1,1,1,0,1,1,1,0), #0 is Corn Oil, 1 is prenatal BP
                         context="CpG",
                         mincov = 1,
                         pipeline = "bismarkCoverage")

#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_150BP_F,win.size=500,step.size=500,cov.bases = 10)
meth=unite(tiles,min.per.group=4L)

##Calculating DMRs without licking/grooming covariate
covariates=data.frame(batch=c(1,3,2,1,1,3,5,5,2,2,5,2))

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
write.csv(myDiff5p, file="DMR_150BP_F.csv", quote=FALSE)

##Calculating DMRs with licking/grooming covariate
covariates=data.frame(LG=c(0.0777,0.401,1.338,2.466,1.283,0.642,0.610,-0.0271,-0.619,-0.510,-1.145,-1.829),batch=c(1,3,2,1,1,3,5,5,2,2,5,2))

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
write.csv(myDiff5p, file="DMR_150BP_F_LG-Covariate.csv", quote=FALSE)
