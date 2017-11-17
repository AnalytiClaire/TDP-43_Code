setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/")
C9_cor.5 <- read.csv("C9_cor.5_cytoscape_combo.csv")
sALS_cor.5 <-read.csv("sALS_cor.5_cytoscape_combo.csv")
FTLD_cor.5 <-read.csv("FTLD_cor.5_cytoscape_combo.csv")
VCP_cor.5 <-read.csv("VCP_cor.5_cytoscape_combo.csv")
PET_cor.5 <-read.csv("PET_cor.5_cytoscape_combo.csv")
RAV_cor.5 <-read.csv("RAV_cor.5_cytoscape_combo.csv")

Commonedge <- Reduce(intersect, list(C9_cor.5$combo, sALS_cor.5$combo, FTLD_cor.5$combo,
                                     VCP_cor.5$combo, PET_cor.5$combo, RAV_cor.5$combo))


#Subset each dataset with these common names so they are all the same size
C9_CE <- subset(C9_cor.5, C9_cor.5$combo %in% Commonedge)
C9_CE <- C9_CE[order(C9_CE$combo),]
sALS_CE <- subset(sALS_cor.5, sALS_cor.5$combo %in% Commonedge)
sALS_CE <- sALS_CE[order(sALS_CE$combo),]
FTLD_CE <- subset(FTLD_cor.5, FTLD_cor.5$combo %in% Commonedge)
FTLD_CE <- FTLD_CE[order(FTLD_CE$combo),]
VCP_CE <- subset(VCP_cor.5, VCP_cor.5$combo %in% Commonedge)
VCP_CE <- VCP_CE[order(VCP_CE$combo),]
PET_CE <- subset(PET_cor.5, PET_cor.5$combo %in% Commonedge)
PET_CE <- PET_CE[order(PET_CE$combo),]
RAV_CE <- subset(RAV_cor.5, RAV_cor.5$combo %in% Commonedge)
RAV_CE <- RAV_CE[order(RAV_CE$combo),]



CommonGroup <- data.frame(row.names = C9_CE$combo,
                         C9 = C9_CE$reg.mat,
                         sALS = sALS_CE$reg.mat,
                         FTLD = FTLD_CE$reg.mat,
                         VCP = VCP_CE$reg.mat,
                         PET = PET_CE$reg.mat,
                         RAV = RAV_CE$reg.mat)


CG_conserved_up <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x > 0)), ]
CG_conserved_down <- CommonGroup[apply(CommonGroup, MARGIN = 1, function(x) all(x < 0)), ]

CG_samedir <- rbind(CG_conserved_up, CG_conserved_down)
CG_samedir$corMean <- rowMeans(CG_samedir, na.rm = FALSE, dims = 1)
write.csv(CG_samedir, "RMETHOD_samedir_mean.csv", quote = F)
