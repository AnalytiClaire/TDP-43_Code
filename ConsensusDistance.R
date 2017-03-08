##ConsensusDistance##

#Run consensus on pathprint
Consensus <- consensusFingerprint(VCP_pathprint, threshold = 0.9)

#Calculate distance from the consensus
Distance<-consensusDistance(Consensus, GEO.fingerprint.matrix)

#Plot histogram
# plot histograms
# par(mfcol = c(2,1), mar = c(0, 4, 4, 2))
# C9Distance.hist<-hist(C9orf72Distance[,"distance"],
#                                    nclass = 50, xlim = c(0,1), main = "Distance from C9orf72 consensus")
# par(mar = c(7, 4, 4, 2))
# hist(C9orf72Distance[C9.LCM_pathprint, "distance"],
#      breaks = C9Distance.hist$breaks, xlim = c(0,1), 
#      main = "", xlab = "above: all GEO, below: curated C9orf72 samples")
