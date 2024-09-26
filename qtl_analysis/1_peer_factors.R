###PEER factor analysis

library(limma)
library(edgeR)
library(peer)
#load expression matrix
MODELED_EX_MAT="gaba/DAC_AnalysisGABA.Rdata"
load(MODELED_EX_MAT)


dir.create(paste0(celltype,"/Peer_results/"), showWarnings = FALSE)
dir.create(paste0(celltype,"/Peer_plots/"), showWarnings = FALSE)

## DesignCovars is provided within the loaded objects as designFinalCovars
DesignCovars = designFinalCovars

#try(if(any(DesignCovars %in% colnames(allInfo) == FALSE)) stop("variables not present"))

## need to create the design matrix:
design_peer <- model.matrix(as.formula(paste(" ~ 0 + ", paste(DesignCovars, collapse = "+"))), allInfo)

svdval = svd(design_peer)
svdval$d

if(length(which(svdval$d<10^-10))>0){
  designsvd = svdval$u[, -which(svdval$d<10^-10)]
  cat("design colinearity","\n")
}

StartingTime = Sys.time()
print(paste("Peer update intialized at:", StartingTime))

model = PEER()
# Set the expression and uncertainty values
PEER_setPhenoMean(model, t(modeledVoomObj$E))
# dim(PEER_getPhenoMean(model)) # [1]    122 18812

PEER_setPhenoVar(model,t(1/modeledVoomObj$weights))
PEER_setCovariates( model, design_peer)
# Run PEER with 1 latent variables for the test:

#nSVs = i = as.numeric(as.character(args[[1]]))
i=10
nSVs=10
PEER_setNk(model, 10)
PEER_update(model)

print(paste("Peer update completed at:", Sys.time()))

k = ncol(design_peer) + 1
SV = PEER_getX(model)[, k:(k + nSVs - 1)]

CompletionTime = Sys.time()
print(paste("Peer update completed at:", CompletionTime))
print(CompletionTime - StartingTime)

peer_results <- ls()
peer_results$factors <- SV  # samples x PEER factors
peer_results$precision <- PEER_getAlpha(model)  # PEER factors x 1
peer_results$residuals <- t(PEER_getResiduals(model))  # peaks x samples
peer_results$run_time <- paste0("peer with ", i, " latent variables took ", CompletionTime - StartingTime, "to complete")

saveRDS(peer_results, paste0(celltype, "/Peer_results/Results_peer_", i,"_LatVar.RDs",sep=""))

# save.image(paste0(WorkingFolder, "Peer_results/Workspace_peer_for_eQTL", i,"_LatVar.Rda",sep=""))

pdf(paste0(celltype,"/Peer_plots/Plot_peer_for_eQTL_", i,"_LatVar.pdf"))
PEER_plotModel(model)
dev.off()

rm(list = ls())
# file.remove(".RData")
