library(WGCNA)
library(flashClust)
library(readr)
library(datasets)
library(ggplot2)
library(tidyverse)

#####WGCNA Analysis#####
#Male MPOA
##Import normalized counts file that was created from DESeq2
datExpr_M <- read_csv("./DESeq2-MPOA-M-fullmodel.csv")
#An 'X' needs to be added to the top row of the first column of the file before reformatting the file
colnames(datExpr_M)[1] <- "X"
View(datExpr_M)

#Reformatting file for use in WGCNA
datExpr_M_matrix <- dat_Expr_M[ , -1]
row.names(datExpr_M_matrix) <- dat_Expr_M$X
datExpr_M_matrix <- as.matrix(datExpr_M_matrix)

datExpr_M_log <- log2(datExpr_M_matrix+1) #Log transform the count data

datExpr_M_log = as.data.frame(t(datExpr_M_log)) #Samples are rows and genes are columns
dim(datExpr_M_log)

#Run this to check if there are gene outliers (there always are)
gsg = goodSamplesGenes(datExpr_M_log, verbose = 3)
gsg $allOK

#Remove gene outliers (low counts and/or variability)

if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
  printFlush(paste("Removing genes:", paste(names(datExpr_M_log)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr_M_log)[!gsg$goodSamples], collapse=", ")))
  datExpr_M_log= datExpr_M_log[gsg$goodSamples, gsg$goodGenes]
}

#Import file that contains the study variables
#Prenatal treatment group is dummy coded, all study variables need to be in numeric format otherwise the WGCNA code will not work
datTraits_M <- read_csv("./data/datTraits_MPOA_M.csv", 
                        col_types = cols(ID = col_character(), 
                                         Group_1 = col_number(), Group_2 = col_number(), 
                                         Group_3 = col_number(), Nest_attendance_P1_5_scaled = col_number(), 
                                         LG_P1_5_scaled = col_number()))
View(datTraits_M)

#Reformatting file
datTraits_M_matrix <- datTraits_M[ , -1]
row.names(datTraits_M_matrix) <- datTraits_M$ID

head(datTraits_M)

#form a data frame analogous to the count data
dim(datTraits_M)
table(rownames(datTraits_M_matrix)==rownames(datExpr_M_log)) #should return TRUE if datasets align correctly

#####Cluster count data#####

A = adjacency(t(datExpr_M_log),type="signed") #this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 #standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")

#Convert traits to a color representation where red indicates high values
traitColors = data.frame(numbers2colors(datTraits_M_matrix,signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits_M_matrix))
datColors = data.frame(outlier = outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram and Trait Heatmap") #no sample outliers were removed in this study

#Soft power was determined by a prior analysis with male mPFC data
enableWGCNAThreads()
softPower = 8
adjacency = adjacency(datExpr_M_log, power = softPower, type = "signed")

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
MEList= moduleEigengenes(datExpr_M_log, colors= dynamicColors,softPower = 8)
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
MPOA_design_M <- read_csv("./data/datTraits_MM.csv")
View(MPOA_design_M)

#Factor prenatal BP groups, control group first
MPOA_design_M$Group <- factor(MPOA_design_M$Group, levels = c('Corn Oil', '50 μg/kg BPA', '50 μg/kg Mixed BP', '150 μg/kg Mixed BP'))

#Scaling maternal care measures so 0 is approximately the mean and 1 is approximately the stdev
MPOA_design_M$Nest_attendance_P1_5_scaled <- (MPOA_design_M$Nest_attendance_P1_5-1600)/400
MPOA_design_M$LG_P1_5_scaled <- (MPOA_design_M$LG_P1_5-300)/100

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr_M_log, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
colnames(MEs) #ME color names

#Incorporate study variables into MEs data
MEs$Group <- MPOA_design_M$Group
MEs$LG_P1_5_scaled <- MPOA_design_M$LG_P1_5_scaled
MEs$Nest_attendance_P1_5_scaled <- MPOA_design_M$Nest_attendance_P1_5_scaled

#Factor prenatal BP groups, control group first
MEs$Group <- factor(MEs$Group, levels = c('Corn Oil', '50 μg/kg BPA', '50 μg/kg Mixed BP', '150 μg/kg Mixed BP'))

##Significance testing for prenatal BP treatment group, maternal care, and their interactions
#Outlier/insignificant modules were excluded
moduleTraitLm = lm(MEsaddlebrown ~ Group, data = MEs)

summary(moduleTraitLm)

#Exporting gene names within modules as a txt file
saddlebrownmodule<-as.data.frame(names(datExpr_M_log)[moduleColors=="saddlebrown"])

write_delim(saddlebrownmodule, file = "M-MPOA-saddlebrownmodule.txt", delim = "\n")

allgenes<-as.data.frame(names(datExpr_M_log)) #all genes analyzed for WGCNA
write_delim(allgenes, file = "M-MPOA-allgenes.txt", delim = "\n")

####Graphing MEs Data#####

#Subset the data by prenatal BP group so individual regression lines can be displayed
MEs_M_Control <- subset(MEs, Group == "Corn Oil")
MEs_M_BPA <- subset(MEs, Group == "50 μg/kg BPA")
MEs_M_BP_Low <- subset(MEs, Group == "50 μg/kg Mixed BP")
MEs_M_BP_High <- subset(MEs, Group == "150 μg/kg Mixed BP")

##ggplot2 scatterplot with points and regression lines colored by treatment group
WGCNA_MPOA_M <- ggplot(MEs, aes(x = LG_P1_5_scaled, y = MEsaddlebrown, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#2c98b0',"#fb945c",'#967be4'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Saddlebrown Eigengene Value") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Licking/Grooming received per day from\nPND1-5 (scaled)")+
  geom_smooth(data = MEs_M_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = MEs_M_BPA, method = 'lm', se= FALSE, colour = "#2c98b0", fullrange = TRUE)+
  geom_smooth(data = MEs_M_BP_Low, method = 'lm', se= FALSE, colour = "#fb945c", fullrange = TRUE)+
  geom_smooth(data = MEs_M_BP_High, method = 'lm', se= FALSE, colour = "#967be4", fullrange = TRUE)

WGCNA_MPOA_M
