library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

#####Multiomic Factor Analysis#####
##Male MPOA
#Transcriptome data preparation
##Import normalized counts file from DESeq2
mRNA_M <- read_csv("./MOFA-tagseq-M.csv")

#Filter out top 5000 most variable genes
mRNA_M_filtered <- mRNA_M %>%
  rowwise() %>%
  mutate(variance = c_across(2:19) %>% var(na.rm=TRUE)) %>%
  ungroup() %>%
  slice_max(n = 5000, order_by = variance)

mRNA_M_matrix <- mRNA_M_filtered[ , -1]
mRNA_M_matrix <- mRNA_M_matrix[ , 1:18]
row.names(mRNA_M_matrix) <- mRNA_M_filtered$X
mRNA_M_matrix <- as.matrix(mRNA_M_matrix)

#DNA Methylome data preparation
##Import normalized DNA methylation file from Methylkit
Methylation_M <- read_csv("./MOFA-methylation-M-transformed.csv")

#Filter out top 10000 most variable regions
Methylation_M_filtered <- Methylation_M %>%
  rowwise() %>%
  mutate(variance = c_across(2:19) %>% var(na.rm=TRUE)) %>%
  ungroup() %>% 
  slice_max(n = 10000, order_by = variance)

Methylation_M_matrix <- Methylation_M_filtered[ , -1]
Methylation_M_matrix <- Methylation_M_matrix[ , 1:18]
row.names(Methylation_M_matrix) <- Methylation_M_filtered$X
Methylation_M_matrix <- as.matrix(Methylation_M_matrix)

##Integrate Transcriptome and DNA methylome datasets
data <- list(mRNA_M_matrix,Methylation_M_matrix)
lapply(data,dim)

MOFAobject <- create_mofa_from_matrix(data, groups = NULL)
views_names(MOFAobject) <- c("mRNA", "Methylation")

plot_data_overview(MOFAobject)

#View and change MOFA parameters
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE
data_opts

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5
model_opts

train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)

##Run MOFA modeling, basilisk does not work properly for MOFA2 so it is advisable to use reticulate and provide your own python directory with numpy and MOFA included
reticulate::use_python("/usr/bin/python3", required = TRUE)  ##use own directory where python is saved
outfile = file.path(getwd(),"model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile)
model <- load_model(file.path("./model.hdf5"))

##Import and reformat study metadata file for statistical analysis
sample_metadata_M <- read_excel("./data/MOFA-variables-M.xlsx")
View(sample_metadata_M)

#Scaling maternal care measures so 0 is approximately the mean and 1 is approximately the stdev
sample_metadata_M$Nest_attendance_P1_5 <- (sample_metadata_M$Nest_attendance_P1_5-1600)/400
sample_metadata_M$LG_P1_5 <- (sample_metadata_M$LG_P1_5-300)/100

#Examining how much variance is explained by each factor and view
head(model@cache$variance_explained$r2_total[[1]])

plot_variance_explained(model, x="view", y="factor")

plot_factor_cor(model)
plot_variance_explained(model, plot_total = T)[[2]]

#####Statistical Analyses#####

samples_metadata(model) <- sample_metadata_M

##Extract and format factor data for analysis
data_M <- get_factors(model, as.data.frame = T)

data_M_wide <- data_M %>% pivot_wider(id_cols = c("sample"), 
                                      names_from = "factor", 
                                      values_from = "value")

#Incorporate study variables into factor data
data_M_wide$Group <- sample_metadata_M$Group
data_M_wide$LG <- sample_metadata_M$LG_P1_5
data_M_wide$Nest <- sample_metadata_M$Nest_attendance_P1_5

#Factor prenatal BP groups, control group first
data_M_wide$Group <- factor(data_M_wide$Group, levels = c('Corn Oil', '50 μg/kg BPA', '50 μg/kg Mixed BP', '150 μg/kg Mixed BP'))

#Significance testing for prenatal BP treatment group, maternal care, and their interactions
feature_lm = lm(Factor1 ~ Group + LG + Group:LG + Nest + Group:Nest, data = data_M_wide)
feature_lm = lm(Factor2 ~ Group + LG + Group:LG + Nest + Group:Nest, data = data_M_wide)
feature_lm = lm(Factor3 ~ Group + LG + Group:LG + Nest + Group:Nest, data = data_M_wide)
feature_lm = lm(Factor4 ~ Group + LG + Group:LG + Nest + Group:Nest, data = data_M_wide)
feature_lm = lm(Factor5 ~ Group + LG + Group:LG + Nest + Group:Nest, data = data_M_wide)

summary(feature_lm)

#Save all model results
write.csv(data_M, file="MOFA-results-M.csv", quote=FALSE)

####Graphing Factor Data#####

#Subset the data by prenatal BP group so individual regression lines can be displayed
Factor_M_Control <- subset(data_M_wide, Group == "Corn Oil")
Factor_M_BPA <- subset(data_M_wide, Group == "50 μg/kg BPA")
Factor_M_BP_Low <- subset(data_M_wide, Group == "50 μg/kg Mixed BP")
Factor_M_BP_High <- subset(data_M_wide, Group == "150 μg/kg Mixed BP")

##ggplot2 scatterplot with points and regression lines colored by treatment group
MOFA_Factor_M <- ggplot(data_M_wide, aes(x = LG, y = Factor3, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#2c98b0',"#fb945c",'#967be4'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Factor 3 Value") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Licking/Grooming received per day from\nPND1-5 (scaled)")+
  geom_smooth(data = Factor_M_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = Factor_M_BPA, method = 'lm', se= FALSE, colour = "#2c98b0", fullrange = TRUE)+
  geom_smooth(data = Factor_M_BP_Low, method = 'lm', se= FALSE, colour = "#fb945c", fullrange = TRUE)+
  geom_smooth(data = Factor_M_BP_High, method = 'lm', se= FALSE, colour = "#967be4", fullrange = TRUE)

MOFA_Factor_M

MOFA_Factor_M <- ggplot(data_M_wide, aes(x = LG, y = Factor4, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#2c98b0',"#fb945c",'#967be4'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Factor 4 Value") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Licking/Grooming received per day from\nPND1-5 (scaled)")+
  geom_smooth(data = Factor_M_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = Factor_M_BPA, method = 'lm', se= FALSE, colour = "#2c98b0", fullrange = TRUE)+
  geom_smooth(data = Factor_M_BP_Low, method = 'lm', se= FALSE, colour = "#fb945c", fullrange = TRUE)+
  geom_smooth(data = Factor_M_BP_High, method = 'lm', se= FALSE, colour = "#967be4", fullrange = TRUE)

MOFA_Factor_M

MOFA_Factor_M <- ggplot(data_M_wide, aes(x = Nest, y = Factor4, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#2c98b0',"#fb945c",'#967be4'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Factor 4 Value") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Time spent on nest per day from\nPND1-5 (scaled)")+
  geom_smooth(data = Factor_M_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = Factor_M_BPA, method = 'lm', se= FALSE, colour = "#2c98b0", fullrange = TRUE)+
  geom_smooth(data = Factor_M_BP_Low, method = 'lm', se= FALSE, colour = "#fb945c", fullrange = TRUE)+
  geom_smooth(data = Factor_M_BP_High, method = 'lm', se= FALSE, colour = "#967be4", fullrange = TRUE)

MOFA_Factor_M
