setwd("/Users/clairegreen/Downloads/")

ALS <- read.csv("targets_associated_with_amyotrophic_lateral_sclerosis.csv")
FTLD <- read.csv("targets_associated_with_Frontotemporal_dementia.csv")
AD <- read.csv("targets_associated_with_Alzheimers_disease.csv")
LBD <- read.csv("targets_associated_with_Lewy_body_dementia.csv")
IBM <- read.csv("targets_associated_with_Inclusion_body_myopathy_with_Paget_disease_of_bone_and_frontotemporal_dementia.csv")

ALS_list <- ALS$target.gene_info.symbol
FTLD_list <- FTLD$target.gene_info.symbol
AD_list <- AD$target.gene_info.symbol
LBD_list <- LBD$target.gene_info.symbol
IBM_list <- IBM$target.gene_info.symbol

consensus <- Reduce(intersect, list(ALS_list, FTLD_list, AD_list, LBD_list, IBM_list))
consensus
