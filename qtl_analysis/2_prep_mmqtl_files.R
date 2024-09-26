#files needed for MMQTL
# 1. normalized residualized matrix 
# 2. GRM using gemma (done)
#3. plink file with imputed genotypes, filtered for 5% MAF and 95% completion. (done)
#4. ene annotation in BED format: chr, starting pos, ending position, name 
'

## 1. expression matrix
library(mgsub)
library(limma)
library(biomaRt)
library(edgeR)
# library(data.table)
celltype="gaba"
LV <- 10
## get the estimated peer factors:
PEER_results<-readRDS(paste0(celltype,"/Peer_results/Results_peer_10_LatVar.RDs"))
colnames(PEER_results$factors) <- paste0("peerLV_", seq(1,ncol(PEER_results$factors),1))
dim(PEER_results$factors)


load("gaba/DAC_Analysis.Rdata")
##add sample names back to peer
#should be same as initialDegobj
rownames(PEER_results$factors)=colnames(initialDgeObj)


## keep only the 848 individuals:
ids<-read.delim("ind_for_eqtl.txt") #848 inds

allInfo = allInfo[allInfo$Sample_RNA_ID%in%ids$Sample_RNA_ID,]
dim(allInfo)

initialDgeObj = initialDgeObj[, colnames(initialDgeObj) %in% ids$Sample_RNA_ID ]
dim(initialDgeObj) 

PEER_results$factors<-PEER_results$factors[rownames(PEER_results$factors)%in%ids$Sample_RNA_ID,]
dim(PEER_results$factors)


#####covariates
Model_withPEER <- cbind(allInfo[, designFinalCovars], PEER_results$factors)

designFinalCovars = colnames(Model_withPEER)
designFinalString = paste("~ 0 +", paste(designFinalCovars, collapse = " + "))
designFinalMatrix = model.matrix(as.formula(designFinalString), Model_withPEER)
modeledVoomObj = voomWithQualityWeights(initialDgeObj , designFinalMatrix, plot=T, normalize.method = "none") 
dim(modeledVoomObj)

# get the calcResiduals function (from Jaro's scripts):
calcResiduals <- function(geneBySampleValues, samplesByCovariates, varsToAddBackIn=NULL, sampleWeights=NULL) {
  # Use the calcResiduals() code after this section in a for loop:
  if (is.matrix(sampleWeights)) {
    residualizedMat = matrix(NA, nrow=nrow(geneBySampleValues), ncol=ncol(geneBySampleValues), dimnames=dimnames(geneBySampleValues))
    for (gInd in 1:nrow(geneBySampleValues)) {
      gRow = calcResiduals(geneBySampleValues[gInd, , drop=FALSE], samplesByCovariates, varsToAddBackIn, sampleWeights[gInd, ])
      residualizedMat[gInd, ] = gRow
    }
    return(residualizedMat)
  }
  result.lm = lm(t(geneBySampleValues) ~ 0 + samplesByCovariates, weights=sampleWeights)
  covarNames = colnames(samplesByCovariates)
  
  coef = result.lm$coefficients
  isMatrixForm = is.matrix(coef)
  if (isMatrixForm) {    rownames(coef) = covarNames  }
  else {    names(coef) = covarNames  }
  
  allVarsToAddBack = '(Intercept)'
  if (!is.null(varsToAddBackIn)) {    allVarsToAddBack = c(allVarsToAddBack, varsToAddBackIn)  }
  allVarsToAddBack = intersect(allVarsToAddBack, covarNames)
  
  residualizedMat = result.lm$residuals
  for (v in allVarsToAddBack) {
    if (isMatrixForm) {      multCoef = coef[v, , drop=FALSE]    }
    else {      multCoef = coef[v]    }
    residualizedMat = residualizedMat + samplesByCovariates[, v, drop=FALSE] %*% multCoef
  }
  residualizedMat = t(residualizedMat)
  rownames(residualizedMat) = rownames(geneBySampleValues)
  colnames(residualizedMat) = colnames(geneBySampleValues)
  
  return(residualizedMat)
}
############

residualizedVoomExp = calcResiduals(
  modeledVoomObj$E,
  designFinalMatrix,
  sampleWeights = modeledVoomObj$weights
)
dim(residualizedVoomExp) # [1] 18898   135


## and now get the needed gene info:
# library("biomaRt")
ensembl = useDataset("hsapiens_gene_ensembl", mart = useMart("ensembl") )
genelist <- rownames(residualizedVoomExp)
table(duplicated(genelist))
#'FALSE 
#18898'

GeneInfoTable <- getBM(attributes = c(
  'chromosome_name', 'start_position', 'end_position', 'ensembl_gene_id','strand', 'external_gene_name', 'hgnc_symbol'),
  filters = 'ensembl_gene_id', 
  values = genelist, 
  mart = ensembl)
colnames(GeneInfoTable) <- c("Chr", "start", "end", "EnsampleID", "strand", "gene_Symbol", "gene_HGNC")
dim(GeneInfoTable)
# [1] 18857     7

# GeneInfoTable[1:5,]

table(duplicated( GeneInfoTable$EnsampleID))
#'FALSE  TRUE 
#18856     1 '

## remove the dups
GeneInfoTable <- GeneInfoTable[-which(duplicated( GeneInfoTable$EnsampleID)),]
rownames(GeneInfoTable) <- GeneInfoTable$EnsampleID

GeneInfoTable$strand <- mgsub(GeneInfoTable$strand, pattern = c(1, -1), replacement = c("+","-"))
GeneInfoTable$Chr <- paste0("chr", GeneInfoTable$Chr)

table(GeneInfoTable$EnsampleID %in%  rownames(residualizedVoomExp))
#' TRUE 
#18856'
table(rownames(residualizedVoomExp) %in% GeneInfoTable$EnsampleID)


## will sort by chromosome and positions:
GeneInfoTable <- GeneInfoTable[order(GeneInfoTable$Chr, GeneInfoTable$start),]

residualizedVoomExp <- residualizedVoomExp[GeneInfoTable$EnsampleID,]
dim(residualizedVoomExp) # [1] 18856   101

## make sure to use the correct sample names to be compatible with the genetic data:
NewSampleNames = allInfo$Genotype_ID[match(colnames(residualizedVoomExp), allInfo$ID)]
NewSampleNames = ids$snpArray[match(colnames(residualizedVoomExp), ids$Sample_RNA_ID)]
# paste(colnames(residualizedVoomExp), NewSampleNames, sep = "_____")
colnames(residualizedVoomExp) = NewSampleNames

## save the objects:
OutDir <- paste0("qtl", "/gaba_exp/")
dir.create(OutDir, showWarnings = F, recursive = T)

saveRDS(residualizedVoomExp, file = paste0(OutDir,"Residualized_Peer10_gaba.RDs"))


saveRDS(GeneInfoTable[, c("Chr", "start", "end", "EnsampleID")], file = paste0(OutDir,"GeneInfo.RDs"))



write.table(GeneInfoTable, file = paste0(OutDir,"GeneInfo.txt"), col.names = TRUE,row.names=F,quote=F,sep="\t")


write.table(x, "Residualized_Peer10_gaba.txt", col.names = TRUE,row.names=T,quote=F,sep="\t")
