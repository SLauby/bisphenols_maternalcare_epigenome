library(WGCNA)
library(flashClust)
library(readr)
library(datasets)
library(ggplot2)
library(tidyverse)

#####WGCNA Analysis#####
#Female MPOA
##Import normalized counts file that was created from DESeq2
datExpr_F <- read_csv("./DESeq2-MPOA-F-fullmodel.csv")
#An 'X' needs to be added to the top row of the first column of the file before reformatting the file
colnames(datExpr_F)[1] <- "X"
View(datExpr_F)

#Reformatting file for use in WGCNA
datExpr_F_matrix <- dat_Expr_F[ , -1]
row.names(datExpr_F_matrix) <- dat_Expr_F$X
datExpr_F_matrix <- as.matrix(datExpr_F_matrix)

datExpr_F_log <- log2(datExpr_F_matrix+1) #Log transform the count data

datExpr_F_log = as.data.frame(t(datExpr_F_log)) #Samples are rows and genes are columns
dim(datExpr_F_log)

#Run this to check if there are gene outliers (there always are)
gsg = goodSamplesGenes(datExpr_F_log, verbose = 3)
gsg $allOK

#Remove gene outliers (low counts and/or variability)

if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
  printFlush(paste("Removing genes:", paste(names(datExpr_F_log)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr_F_log)[!gsg$goodSamples], collapse=", ")))
  datExpr_F_log= datExpr_F_log[gsg$goodSamples, gsg$goodGenes]
}

#Import file that contains the study variables
#Prenatal treatment group is dummy coded, all study variables need to be in numeric format otherwise the WGCNA code will not work
datTraits_F <- read_csv("./data/datTraits_FM.csv", 
                        col_types = cols(ID = col_character(), 
                                         Group_1 = col_number(), Group_2 = col_number(), 
                                         Group_3 = col_number(), Nest_attendance_P1_5_scaled = col_number(), 
                                         LG_P1_5_scaled = col_number()))
View(datTraits_F)

#Reformatting file
datTraits_F_matrix <- datTraits_F[ , -1]
row.names(datTraits_F_matrix) <- datTraits_F$ID

head(datTraits_F)

#form a data frame analogous to the count data
dim(datTraits_F)
table(rownames(datTraits_F_matrix)==rownames(datExpr_F_log)) #should return TRUE if datasets align correctly

#####Cluster count data#####

A = adjacency(t(datExpr_F_log),type="signed") #this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 #standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")

#Convert traits to a color representation where red indicates high values
traitColors = data.frame(numbers2colors(datTraits_F_matrix,signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits_F_matrix))
datColors = data.frame(outlier = outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram and Trait Heatmap") #no sample outliers were removed in this study

#Soft power was determined by a prior analysis with male mPFC data
enableWGCNAThreads()
softPower = 8
adjacency = adjacency(datExpr_F_log, power = softPower, type = "signed")

##Construct Networks

#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed")
dissTOM = 1-TOM


##Generate co-expressed eigengene modules

#Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")

plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#Minimum number of genes to cluster into a module was set to 30
minModuleSize = 30 
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)

#Create module eigengene lists
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr_F_log, colors= dynamicColors,softPower = 8)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

#Plot dendrogram with module colors below
plotDendroAndColors(geneTree, cbind(dynamicColors), c("Eigengene Modules"), rowText=names(moduleLabels), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

#####Statistical Analyses#####

##Import and reformat study traits file for statistical analysis
MPOA_design_F <- read_excel("./data/MPOA_tagseq_variables_F_outliersremoved.xlsx")
View(MPOA_design_F)

#Factor prenatal BP groups, control group first
MPOA_design_F$Group <- factor(MPOA_design_F$Group, levels = c('Corn Oil', '50 μg/kg BPA', '50 μg/kg Mixed BP', '150 μg/kg Mixed BP'))

#Scaling maternal care measures so 0 is approximately the mean and 1 is approximately the stdev
MPOA_design_F$Nest_attendance_P1_5_scaled <- (MPOA_design_F$Nest_attendance_P1_5-1600)/400
MPOA_design_F$LG_P1_5_scaled <- (MPOA_design_F$LG_P1_5-300)/100

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr_F_log, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
colnames(MEs) #ME color names

#Incorporate study variables into MEs data
MEs$Group <- MPOA_design_F$Group
MEs$LG_P1_5_scaled <- MPOA_design_F$LG_P1_5_scaled
MEs$Nest_attendance_P1_5_scaled <- MPOA_design_F$Nest_attendance_P1_5_scaled

#Factor prenatal BP groups, control group first
MEs$Group <- factor(MEs$Group, levels = c('Corn Oil', '50 μg/kg BPA', '50 μg/kg Mixed BP', '150 μg/kg Mixed BP'))

##Significance testing for prenatal BP treatment group, maternal care, and their interactions

moduleTraitLm = lm(MEblue ~ Group + LG_P1_5_scaled + Group:LG_P1_5_scaled + Nest_attendance_P1_5_scaled + Group:Nest_attendance_P1_5_scaled, data = MEs)
moduleTraitLm = lm(MEbrown ~ Group + LG_P1_5_scaled + Group:LG_P1_5_scaled + Nest_attendance_P1_5_scaled + Group:Nest_attendance_P1_5_scaled, data = MEs)
moduleTraitLm = lm(MEmagenta ~ Group + LG_P1_5_scaled + Group:LG_P1_5_scaled + Nest_attendance_P1_5_scaled + Group:Nest_attendance_P1_5_scaled, data = MEs)
moduleTraitLm = lm(MEturquoise ~ Group + LG_P1_5_scaled + Group:LG_P1_5_scaled + Nest_attendance_P1_5_scaled + Group:Nest_attendance_P1_5_scaled, data = MEs)
moduleTraitLm = lm(MEyellow ~ Group + LG_P1_5_scaled + Group:LG_P1_5_scaled + Nest_attendance_P1_5_scaled + Group:Nest_attendance_P1_5_scaled, data = MEs)

summary(moduleTraitLm)

##Calculating kME connectivity with estrogen receptors for each module

datKME = signedKME(datExpr_F_log, MEs, outputColumnName = "kME")

datKME %>%
  filter(row.names(datKME) %in% c('Esr1','Esr2','Esrrg','Esrra','Esrrb')) -> Esr_Hub_F_MPOA

write.csv(Esr_Hub_F_MPOA, "WGCNA-Female-MPOA-Esr-kME.csv")

#Exporting gene names within modules as a txt file
bluemodule<-as.data.frame(names(datExpr_F_log)[moduleColors=="blue"])
brownmodule<-as.data.frame(names(datExpr_F_log)[moduleColors=="brown"])
magentamodule<-as.data.frame(names(datExpr_F_log)[moduleColors=="magenta"])
turquoisemodule<-as.data.frame(names(datExpr_F_log)[moduleColors=="turquoise"])
yellowmodule<-as.data.frame(names(datExpr_F_log)[moduleColors=="yellow"])

write_delim(bluemodule, file = "F-MPOA-bluemodule.txt", delim = "\n")
write_delim(brownmodule, file = "F-MPOA-brownmodule.txt", delim = "\n")
write_delim(magentamodule, file = "F-MPOA-magentamodule.txt", delim = "\n")
write_delim(turquoisemodule, file = "F-MPOA-turquoisemodule.txt", delim = "\n")
write_delim(yellowmodule, file = "F-MPOA-yellowmodule.txt", delim = "\n")

allgenes<-as.data.frame(names(datExpr_F_log)) #all genes analyzed for WGCNA
write_delim(allgenes, file = "F-MPOA-allgenes.txt", delim = "\n")

####Graphing MEs Data#####

#Subset the data by prenatal BP group so individual regression lines can be displayed
MEs_F_Control <- subset(MEs, Group == "Corn Oil")
MEs_F_BPA <- subset(MEs, Group == "50 μg/kg BPA")
MEs_F_BP_Low <- subset(MEs, Group == "50 μg/kg Mixed BP")
MEs_F_BP_High <- subset(MEs, Group == "150 μg/kg Mixed BP")

##ggplot2 scatterplot with points and regression lines colored by treatment group
WGCNA_MPOA_F <- ggplot(MEs, aes(x = LG_P1_5_scaled, y = MEblue, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#2c98b0',"#fb945c",'#967be4'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Blue Eigengene Value") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Licking/Grooming received per day from\nPND1-5 (scaled)")+
  geom_smooth(data = MEs_F_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BPA, method = 'lm', se= FALSE, colour = "#2c98b0", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BP_Low, method = 'lm', se= FALSE, colour = "#fb945c", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BP_High, method = 'lm', se= FALSE, colour = "#967be4", fullrange = TRUE)

WGCNA_MPOA_F

WGCNA_MPOA_F <- ggplot(MEs, aes(x = LG_P1_5_scaled, y = MEbrown, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#2c98b0',"#fb945c",'#967be4'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Brown Eigengene Value") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Licking/Grooming received per day from\nPND1-5 (scaled)")+
  geom_smooth(data = MEs_F_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BPA, method = 'lm', se= FALSE, colour = "#2c98b0", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BP_Low, method = 'lm', se= FALSE, colour = "#fb945c", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BP_High, method = 'lm', se= FALSE, colour = "#967be4", fullrange = TRUE)

WGCNA_MPOA_F

WGCNA_MPOA_F <- ggplot(MEs, aes(x = LG_P1_5_scaled, y = MEmagenta, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#2c98b0',"#fb945c",'#967be4'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Magenta Eigengene Value") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Licking/Grooming received per day from\nPND1-5 (scaled)")+
  geom_smooth(data = MEs_F_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BPA, method = 'lm', se= FALSE, colour = "#2c98b0", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BP_Low, method = 'lm', se= FALSE, colour = "#fb945c", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BP_High, method = 'lm', se= FALSE, colour = "#967be4", fullrange = TRUE)

WGCNA_MPOA_F

WGCNA_MPOA_F <- ggplot(MEs, aes(x = LG_P1_5_scaled, y = MEturquoise, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#2c98b0',"#fb945c",'#967be4'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Turquoise Eigengene Value") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Licking/Grooming received per day from\nPND1-5 (scaled)")+
  geom_smooth(data = MEs_F_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BPA, method = 'lm', se= FALSE, colour = "#2c98b0", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BP_Low, method = 'lm', se= FALSE, colour = "#fb945c", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BP_High, method = 'lm', se= FALSE, colour = "#967be4", fullrange = TRUE)

WGCNA_MPOA_F

WGCNA_MPOA_F <- ggplot(MEs, aes(x = LG_P1_5_scaled, y = MEyellow, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#2c98b0',"#fb945c",'#967be4'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Yellow Eigengene Value") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Licking/Grooming received per day from\nPND1-5 (scaled)")+
  geom_smooth(data = MEs_F_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BPA, method = 'lm', se= FALSE, colour = "#2c98b0", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BP_Low, method = 'lm', se= FALSE, colour = "#fb945c", fullrange = TRUE)+
  geom_smooth(data = MEs_F_BP_High, method = 'lm', se= FALSE, colour = "#967be4", fullrange = TRUE)

WGCNA_MPOA_F
