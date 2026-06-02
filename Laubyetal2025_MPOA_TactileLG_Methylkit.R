library(methylKit)

####Female Tactile LG DMR Analyses####
####Bisphenol Treatment vs. Corn Oil within Supplemental Tactile Stimulation Condition Analyses####

###150 ug/kg BPA vs. Corn Oil within nonstimulated group###
##Preparing the methylation dataset for nonstimulated 150 ug/kg BPA vs. Corn Oil comparisons
file.list=list(file.path("./124F12_S4_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./125F12_S1_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./126F12_S11_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./128F12_S15_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./150F12_S19_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./156F12_S29_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./111F12_S1_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./112F12_S3_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./113F12_S7_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./115F12_S13_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./133F12_S21_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./140F12_S25_L002.sorted.dedup_resorted.bismark.cov"))

myobj_150BPA_Nonstim = methRead(file.list,
                                sample.id=list("124F12","125F12","126F12","128F12","150F12","156F12","111F12","112F12","113F12","115F12","133F12","140F12"),
                                assembly="rn7",
                                treatment=c(0,0,0,0,0,0,1,1,1,1,1,1), #0 is Corn Oil, 1 is prenatal BP
                                context="CpG",
                                mincov = 1,
                                pipeline = "bismarkCoverage")

#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_150BPA_Nonstim,win.size=500,step.size=500,cov.bases = 10)
meth=methylKit::unite(tiles, min.per.group = 4L)

##Calculating DMRs
covariates=data.frame(batch=c(1,1,2,2,3,3,1,1,2,2,3,3))

myDiff=calculateDiffMeth(meth,
                         overdispersion = "MN",
                         covariates=covariates,
                         adjust="SLIM",
                         test="Chisq")

#Extract DMRs with q-value < 0.05 and >5% methylation difference
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs

#Save output file for gene annotation
write.csv(myDiff5p, file="DMR_150BPA_Control_Nonstim.csv", quote=FALSE)

###150 ug/kg BPA vs. Corn Oil within stimulated group###
##Preparing the methylation dataset for stimulated 150 ug/kg BPA vs. Corn Oil comparisons
file.list=list(file.path("./124F22_S5_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./125F22_S2_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./126F22_S12_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./128F22_S16_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./150F22_S20_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./156F22_S30_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./111F22_S2_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./112F22_S4_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./113F22_S8_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./115F22_S14_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./133F23_S22_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./140F22_S26_L002.sorted.dedup_resorted.bismark.cov"))

myobj_150BPA_Stim = methRead(file.list,
                             sample.id=list("124F22","125F22","126F22","128F22","150F22","156F22","111F22","112F22","113F22","115F22","133F22","140F22"),
                             assembly="rn7",
                             treatment=c(0,0,0,0,0,0,1,1,1,1,1,1), #0 is Corn Oil, 1 is prenatal BP
                             context="CpG",
                             mincov = 1,
                             pipeline = "bismarkCoverage")
                             
#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_150BPA_Stim,win.size=500,step.size=500,cov.bases = 10)
meth=methylKit::unite(tiles, min.per.group = 4L)

##Calculating DMRs
covariates=data.frame(batch=c(1,1,2,2,3,3,1,1,2,2,3,3))

myDiff=calculateDiffMeth(meth,
                         overdispersion = "MN",
                         covariates=covariates,
                         adjust="SLIM",
                         test="Chisq")

#Extract DMRs with q-value < 0.05 and >5% methylation difference
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs

#Save output file for gene annotation
write.csv(myDiff5p, file="DMR_150BPA_Control_Stim.csv", quote=FALSE)

###150 ug/kg Mixed BP vs. Corn Oil within nonstimulated group###
##Preparing the methylation dataset for nonstimulated 150 ug/kg Mixed BP vs. Corn Oil comparisons
file.list=list(file.path("./124F12_S4_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./125F12_S1_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./126F12_S11_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./128F12_S15_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./150F12_S19_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./156F12_S29_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./118F12_S3_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./119F12_S9_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./120F12_S5_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./141F12_S17_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./145F12_S23_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./146F12_S27_L002.sorted.dedup_resorted.bismark.cov"))

myobj_150BP_Nonstim = methRead(file.list,
                               sample.id=list("124F12","125F12","126F12","128F12","150F12","156F12","118F12","119F12","120F12","141F12","145F12","146F12"),
                               assembly="rn7",
                               treatment=c(0,0,0,0,0,0,1,1,1,1,1,1), #0 is Corn Oil, 1 is prenatal BP
                               context="CpG",
                               mincov = 1,
                               pipeline = "bismarkCoverage")
                               
#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_150BP_Nonstim,win.size=500,step.size=500,cov.bases = 10)
meth=methylKit::unite(tiles, min.per.group = 4L)

##Calculating DMRs
covariates=data.frame(batch=c(1,1,2,2,3,3,1,2,1,2,3,3))

myDiff=calculateDiffMeth(meth,
                         overdispersion = "MN",
                         covariates=covariates,
                         adjust="SLIM",
                         test="Chisq")

#Extract DMRs with q-value < 0.05 and >5% methylation difference
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs

#Save output file for gene annotation
write.csv(myDiff5p, file="DMR_150BP_Control_Nonstim.csv", quote=FALSE)

###150 ug/kg Mixed BP vs. Corn Oil within stimulated group###
##Preparing the methylation dataset for stimulated 150 ug/kg Mixed BP vs. Corn Oil comparisons
file.list=list(file.path("./124F22_S5_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./125F22_S2_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./126F22_S12_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./128F22_S16_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./150F22_S20_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./156F22_S30_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./119F22_S10_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./120F22_S6_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./141F22_S18_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./145F22_S24_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./146F22_S28_L002.sorted.dedup_resorted.bismark.cov"))

myobj_150BP_Stim = methRead(file.list,
                            sample.id=list("124F22","125F22","126F22","128F22","150F22","156F22","119F22","120F22","141F22","145F22","146F22"),
                            assembly="rn7",
                            treatment=c(0,0,0,0,0,0,1,1,1,1,1), #0 is Corn Oil, 1 is prenatal BP
                            context="CpG",
                            mincov = 1,
                            pipeline = "bismarkCoverage")
                            
#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_150BP_Stim,win.size=500,step.size=500,cov.bases = 10)
meth=methylKit::unite(tiles, min.per.group = 4L)

##Calculating DMRs
covariates=data.frame(batch=c(1,1,2,2,3,3,2,1,2,3,3))

myDiff=calculateDiffMeth(meth,
                         overdispersion = "MN",
                         covariates=covariates,
                         adjust="SLIM",
                         test="Chisq")

#Extract DMRs with q-value < 0.05 and >5% methylation difference
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs

#Save output file for gene annotation
write.csv(myDiff5p, file="DMR_150BP_Control_Stim.csv", quote=FALSE)

####Nonstimulated vs. Stimulated within Prenatal Treatment Group Analyses####
###Nonstimulated vs. Stimulated within Corn Oil group###
##Preparing the methylation dataset for Corn Oil nonstimulated vs. stimulated comparisons
file.list=list(file.path("./124F12_S4_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./124F22_S5_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./125F12_S1_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./125F22_S2_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./126F12_S11_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./126F22_S12_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./128F12_S15_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./128F22_S16_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./150F12_S19_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./150F22_S20_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./156F12_S29_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./156F22_S30_L002.sorted.dedup_resorted.bismark.cov"))

myobj_Stim_Control = methRead(file.list,
                              sample.id=list("124F12","124F22","125F12","125F22","126F12","126F22","128F12","128F22","150F12","150F22","156F12","156F22" ),
                              assembly="rn7",
                              treatment=c(0,1,0,1,0,1,0,1,0,1,0,1), #0 is nonstimulated, 1 is stimulated
                              context="CpG",
                              mincov = 1,
                              pipeline = "bismarkCoverage")
                              
#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_Stim_Control,win.size=500,step.size=500,cov.bases = 10)
meth=methylKit::unite(tiles, min.per.group = 4L)

##Calculating DMRs
covariates=data.frame(batch=c(1,1,1,1,2,2,2,2,3,3,3,3))

myDiff=calculateDiffMeth(meth,
                         overdispersion = "MN",
                         covariates=covariates,
                         adjust="SLIM",
                         test="Chisq")

#Extract DMRs with q-value < 0.05 and >5% methylation difference
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs

#Save output file for gene annotation
write.csv(myDiff5p, file="DMR_Control_Stim.csv", quote=FALSE)

###Nonstimulated vs. Stimulated within 150 ug/kg BPA group###
##Preparing the methylation dataset for 150 ug/kg BPA nonstimulated vs. stimulated comparisons
file.list=list(file.path("./111F12_S1_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./111F22_S2_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./112F12_S3_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./112F22_S4_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./113F12_S7_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./113F22_S8_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./115F12_S13_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./115F22_S14_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./133F12_S21_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./133F23_S22_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./140F12_S25_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./140F22_S26_L002.sorted.dedup_resorted.bismark.cov"))

myobj_Stim_150BPA = methRead(file.list,
                             sample.id=list("111F12","111F22","112F12","112F22","113F12","113F22","115F12","115F22","133F12","133F22","140F12","140F22" ),
                             assembly="rn7",
                             treatment=c(0,1,0,1,0,1,0,1,0,1,0,1), #0 is nonstimulated, 1 is stimulated
                             context="CpG",
                             mincov = 1,
                             pipeline = "bismarkCoverage")
                             
#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_Stim_150BPA,win.size=500,step.size=500,cov.bases = 10)
meth=methylKit::unite(tiles, min.per.group = 4L)

##Calculating DMRs
covariates=data.frame(batch=c(1,1,1,1,2,2,2,2,3,3,3,3))

myDiff=calculateDiffMeth(meth,
                         overdispersion = "MN",
                         covariates=covariates,
                         adjust="SLIM",
                         test="Chisq")

#Extract DMRs with q-value < 0.05 and >5% methylation difference
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs

#Save output file for gene annotation
write.csv(myDiff5p, file="DMR_150BPA_Stim.csv", quote=FALSE)

###Nonstimulated vs. Stimulated within 150 ug/kg Mixed BP group###
##Preparing the methylation dataset for 150 ug/kg Mixed BP nonstimulated vs. stimulated comparisons
file.list=list(file.path("./118F12_S3_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./119F12_S9_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./119F22_S10_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./120F12_S5_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./120F22_S6_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./141F12_S17_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./141F22_S18_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./145F12_S23_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./145F22_S24_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./146F12_S27_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./146F22_S28_L002.sorted.dedup_resorted.bismark.cov"))

myobj_Stim_150BP = methRead(file.list,
                            sample.id=list("118F12","119F12","119F22","120F12","120F22","141F12","141F22","145F12","145F22","146F12","146F22" ),
                            assembly="rn7",
                            treatment=c(0,0,1,0,1,0,1,0,1,0,1), #0 is nonstimulated, 1 is stimulated
                            context="CpG",
                            mincov = 1,
                            pipeline = "bismarkCoverage")
                            
#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_Stim_150BP,win.size=500,step.size=500,cov.bases = 10)
meth=methylKit::unite(tiles, min.per.group = 4L)

##Calculating DMRs
covariates=data.frame(batch=c(1,2,2,1,1,2,2,3,3,3,3))

myDiff=calculateDiffMeth(meth,
                         overdispersion = "MN",
                         covariates=covariates,
                         adjust="SLIM",
                         test="Chisq")

#Extract DMRs with q-value < 0.05 and >5% methylation difference
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.05)
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper") #View hypermethylated DMRs
getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo") #View hypomethylated DMRs

#Save output file for gene annotation
write.csv(myDiff5p, file="DMR_150BP_Stim.csv", quote=FALSE)


####Covariate/interaction analysis####
###150 ug/kg BPA vs. Corn Oil Without Supplemental Tactile Stimulation Covariate###
##Preparing the methylation dataset for 150 ug/kg BPA vs. Corn Oil comparisons
file.list=list(file.path("./124F12_S4_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./124F22_S5_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./125F12_S1_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./125F22_S2_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./126F12_S11_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./126F22_S12_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./128F12_S15_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./128F22_S16_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./150F12_S19_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./150F22_S20_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./156F12_S29_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./156F22_S30_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./111F12_S1_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./111F22_S2_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./112F12_S3_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./112F22_S4_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./113F12_S7_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./113F22_S8_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./115F12_S13_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./115F22_S14_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./133F12_S21_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./133F23_S22_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./140F12_S25_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./140F22_S26_L002.sorted.dedup_resorted.bismark.cov"))

myobj_150BPA = methRead(file.list,
                        sample.id=list("124F12","124F22","125F12","125F22","126F12","126F22","128F12","128F22","150F12","150F22","156F12","156F22","111F12","111F22","112F12","112F22","113F12","113F22","115F12","115F22","133F12","133F23","140F12","140F22"),
                        assembly="rn7",
                        treatment=c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1), #0 is Corn Oil, 1 is prenatal BP
                        context="CpG",
                        mincov = 1,
                        pipeline = "bismarkCoverage")

#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_150BPA,win.size=500,step.size=500,cov.bases = 10)
meth=methylKit::unite(tiles, min.per.group = 10L)

##Calculating DMRs without supplemental tactile stimulation covariate
covariates=data.frame(batch=c(1,1,1,1,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,3,3,3,3))

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
write.csv(myDiff5p, file="DMR_150BPA_Control_No-Stim-Covariate.csv", quote=FALSE)

###150 ug/kg BPA vs. Corn Oil With Supplemental Tactile Stimulation Covariate###
##Calculating DMRs with supplemental tactile stimulation covariate
covariates=data.frame(batch=c(1,1,1,1,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,3,3,3,3), 
                      stim=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1)) #0 is nonstimulated, 1 is stimulated

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
write.csv(myDiff5p, file="DMR_150BPA_Control_Stim-Covariate.csv", quote=FALSE)

###150 ug/kg Mixed BP vs. Corn Oil Without Supplemental Tactile Stimulation Covariate###
##Preparing the methylation dataset for 150 ug/kg Mixed BP vs. Corn Oil comparisons
file.list=list(file.path("./124F12_S4_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./124F22_S5_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./125F12_S1_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./125F22_S2_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./126F12_S11_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./126F22_S12_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./128F12_S15_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./128F22_S16_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./150F12_S19_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./150F22_S20_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./156F12_S29_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./156F22_S30_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./118F12_S3_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./119F12_S9_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./119F22_S10_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./120F12_S5_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./120F22_S6_L001.sorted.dedup_resorted.bismark.cov"),
               file.path("./141F12_S17_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./141F22_S18_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./145F12_S23_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./145F22_S24_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./146F12_S27_L002.sorted.dedup_resorted.bismark.cov"),
               file.path("./146F22_S28_L002.sorted.dedup_resorted.bismark.cov"))

myobj_150BP = methRead(file.list,
                       sample.id=list("124F12","124F22","125F12","125F22","126F12","126F22","128F12","128F22","150F12","150F22","156F12","156F22","118F12","119F12","119F22","120F12","120F22","141F12","141F22","145F12","145F22","146F12","146F22"),
                       assembly="rn7",
                       treatment=c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1), #0 is Corn Oil, 1 is prenatal BP
                       context="CpG",
                       mincov = 1,
                       pipeline = "bismarkCoverage")

#Creating regions with 500 bp length
tiles = tileMethylCounts(myobj_150BP,win.size=500,step.size=500,cov.bases = 10)
meth=methylKit::unite(tiles, min.per.group = 10L)

##Calculating DMRs without supplemental tactile stimulation covariate
covariates=data.frame(batch=c(1,1,1,1,2,2,2,2,3,3,3,3,1,2,2,1,1,2,2,3,3,3,3))

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
write.csv(myDiff5p, file="DMR_150BP_Control_No-Stim-Covariate.csv", quote=FALSE)

###150 ug/kg Mixed BP vs. Corn Oil Without Supplemental Tactile Stimulation Covariate###
##Calculating DMRs with supplemental tactile stimulation covariate
covariates=data.frame(batch=c(1,1,1,1,2,2,2,2,3,3,3,3,1,2,2,1,1,2,2,3,3,3,3), 
                      stim=c(0,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,0,1)) #0 is nonstimulated, 1 is stimulated

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
write.csv(myDiff5p, file="DMR_150BP_Control_Stim-Covariate.csv", quote=FALSE)

