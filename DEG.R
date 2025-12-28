library(ggplot2)
library(edgeR)
library(Rtsne)
library(variancePartition)

########################################################################################
##### CONFIG & PROLOGUE ################################################################

{
  ROOT = "~/Desktop/molecular_profiling/molecular_profiling/" # !!! FIXME: SET TO YOUR CUSTOM DIRECTORY !!!
  CPU_CORES = 8                                               # !!! FIXME: SET THE DESIRED NUMBER OF CORES !!!
  
  # Load metadata (+ custom fixes)
  METADATA = file.path(ROOT, "inputs", "qc_all_rna.csv")      # Pre-calculated QC metrics for RNA-seq samples from processing computational pipeline
  allInfo = read.csv(METADATA)
  rownames(allInfo) = allInfo$ID
  allInfo$Dx_asFactor = ordered(allInfo$Dx, levels=c("Control", "SCZ"))
  allInfo$cell_subtype_asFactor = ordered(allInfo$cell_subtype_abbreviation, levels=c("GABA", "GLU", "OLIG", "MGAS"))
  allInfo$Gender_asFactor = ordered(allInfo$Gender, levels=c("Male", "Female"))
  allInfo$Person_ID_asFactor = as.factor(allInfo$Person_ID)
  
  # Load raw read count matrix
  RNASEQ_COUNT_MATRIX_RAW = file.path(ROOT, "inputs", "rnaseq_count_matrix_raw.RDS")
  countMatrixRaw = readRDS(RNASEQ_COUNT_MATRIX_RAW)

  # Load Ensembl
  gtf = readRDS(file.path(ROOT, "inputs", "genes_annot.RDS"))
  rownames(gtf) = gtf$PeakID
  
  # Default output dir
  outDir = file.path(ROOT, "outputs")
  
  # BIC-related settings for covariate selection: the selection procedure is repeated until the best-performing covariates from selection pool improved at least (my_covarFracCutOff e.g. 5%)% of genes/peaks by (my_deltaBicCutOff, e.g. 2)
  my_covarFracCutOff = 0.05
  my_deltaBicCutOff = 2
  
  # Load helper scripts
  source(file.path(ROOT, "helper_functions.R"))
  options(mc.cores=CPU_CORES)
}

########################################################################################
##### SELECTION OF COVARIATES USING BAYESIAN INFORMATION CRITERION (BIC) APPROACH ######

{
  ### Initial set of covariates for exploration (RNA-seq)
  covariateInfo = list(
    isBiologicalNumeric = c("AOD_scaled", "geno_PC1", "geno_PC2", "geno_PC3", "deconvolution_GABA", "deconvolution_GLU", "deconvolution_AST", "deconvolution_MG", "deconvolution_ODC"),
    isBiologicalFactor = c("cell_subtype_asFactor", "Dx_asFactor", "Gender_asFactor"),
    isTechnicalNumeric = c("insertMetrics_MEDIAN_INSERT_SIZE_cScaled",  "insertMetrics_MEDIAN_ABSOLUTE_DEVIATION_cScaled", "insertMetrics_MEAN_INSERT_SIZE_cScaled", "insertMetrics_STANDARD_DEVIATION_cScaled",
                           "insertMetrics_WIDTH_OF_10_PERCENT_cScaled", "insertMetrics_WIDTH_OF_20_PERCENT_cScaled", "insertMetrics_WIDTH_OF_30_PERCENT_cScaled", "insertMetrics_WIDTH_OF_40_PERCENT_cScaled", 
                           "insertMetrics_WIDTH_OF_50_PERCENT_cScaled", "insertMetrics_WIDTH_OF_60_PERCENT_cScaled", "insertMetrics_WIDTH_OF_70_PERCENT_cScaled", "insertMetrics_WIDTH_OF_80_PERCENT_cScaled",
                           "insertMetrics_WIDTH_OF_90_PERCENT_cScaled", "picard_meanGcContent_cScaled", "picard_AT_DROPOUT_cScaled", "picard_GC_NC_0_19_cScaled", "picard_GC_NC_20_39_cScaled",
                           "picard_GC_NC_40_59_cScaled", "picard_GC_NC_60_79_cScaled", "picard_GC_NC_80_100_cScaled", "star_pct_of_reads_mapped_to_multiple_loci_scaled", "star_pct_of_reads_mapped_to_too_many_loci_scaled",
                           "star_pct_of_reads_unmapped_other_scaled", "star_pct_of_reads_unmapped_too_short_scaled", "star_Uniquely_mapped_reads_pct_scaled", "finalReadCount_scaled",
                           "finalReadCountFrac_scaled", "picard_PERCENT_DUPLICATION_scaled", "star_Average_mapped_length_scaled", "star_Number_of_input_reads_scaled", 
                           "star_Number_of_reads_mapped_to_multiple_loci_scaled", "star_Number_of_reads_mapped_to_too_many_loci_scaled", "star_Uniquely_mapped_reads_number_scaled", "PMI..in.hours.",
                           "rnaseqc_PCT_R1_TRANSCRIPT_STRAND_READS", "rnaseqc_PCT_R2_TRANSCRIPT_STRAND_READS", "rnaseqc_PCT_RIBOSOMAL_BASES", "rnaseqc_PCT_CODING_BASES", "rnaseqc_PCT_UTR_BASES", 
                           "rnaseqc_PCT_INTRONIC_BASES", "rnaseqc_PCT_INTERGENIC_BASES", "rnaseqc_PCT_MRNA_BASES", "rnaseqc_PCT_USABLE_BASES"),
    isTechnicalFactor = c())
}

########################################################################################
##### READ COUNT NORMALIZATION (TMM METHOD) & MDS & tSNE & ANALYSIS OF COVARIATES ######

{
  # Make edgeR object:
  expObjAll = DGEList(counts=countMatrixRaw, genes=gtf[match(rownames(countMatrixRaw),rownames(gtf)),])
  
  # Keep genes with at least 1 count-per-million reads (cpm) in at least 20% of the samples:
  fracSamplesWithMinCPM = rowMeans(cpm(expObjAll) >= 1)
  isNonLowExpr = fracSamplesWithMinCPM >= 0.2
  expObjNonLow = expObjAll[isNonLowExpr, , keep.lib.sizes=F]
  
  geneNormObj = calcNormFactors(expObjNonLow, method="TMM")
  
  # Voom normalization
  initialDgeObj = geneNormObj
  
  mpdf("misc_DEG_voomFirstPlot", outDir=file.path(ROOT, "outputs"))
  initialVoomObj = voom(initialDgeObj, design=NULL, plot=T)
  dev.off()
  
  ######
  ### Plotting MDS
  {
    myDist = as.dist(sqrt(1-cor(initialVoomObj$E)^2)) # # Squared distance correlation
    mdsResults = cmdscale(myDist, k=2, eig=T)
    colnames(mdsResults$points) = c("Coordinate_1", "Coordinate_2")
    mdsResults$points = cbind(mdsResults$points, allInfo)
    
    mdsPreCovsPlot = ggplot(mdsResults$points, aes(x=Coordinate_1, y=Coordinate_2, shape=Dx, color=cell_subtype_abbreviation)) + geom_point(size=2) + ggtitle("MDS") +
      coord_fixed() + xlab("Coordinate 1") + ylab("Coordinate 2") + theme_classic() + theme(axis.text=element_text(colour="black"))
    
    mpdf("misc_DEG_MDS_preCovs", outDir=file.path(ROOT, "outputs")); print(mdsPreCovsPlot); dev.off()
  }
  
  ######
  ### Plotting tSNE
  {
    set.seed(sample(1:100, 1))
    tsne = Rtsne(t(initialVoomObj$E), perplexity=floor(ncol(initialVoomObj$E)/3), verbose=T)
    allInfo$tsne_dim_1 = tsne$Y[,1]
    allInfo$tsne_dim_2 = tsne$Y[,2]
    
    myPlot = ggplot(allInfo, aes(x=tsne_dim_1, y=tsne_dim_2, shape=Dx, color=cell_subtype_abbreviation)) +
      geom_point(size=2) + ggtitle("tSNE") + coord_fixed() + xlab("Dimension 1") + ylab("Dimension 2") + theme_classic() +
      theme(axis.text=element_text(colour="black"))
    mpdf("misc_DEG_tSNE_preCovs", outDir=file.path(ROOT, "outputs")); print(myPlot); dev.off()
  }
  
  ######
  ### Calculating which covariates are significantly associated (at FDR<0.20) with at least one PC of gene expression variance that explain at least 1pct of variance
  {
    covariatesForExploration = unname(unlist(covariateInfo[c("isTechnicalNumeric","isTechnicalFactor")]))
    
    myName = "RNA-seq covariate exploration"
    myNameSuffix="aggPca"
    
    preCovPca = list()
    preCovPca$pca_all = calcAndPlotPCAs(unlist(unname(covariateInfo)), initialVoomObj$E, allInfo, "all", myName, covariateColsToAccountForInMultipleTesting=covariatesForExploration)
    preCovPca$aggPcaCorrDat=preCovPca$pca_all$pcaCorr
    preCovPca$aggPcaCorrDat$origPC=preCovPca$aggPcaCorrDat$PC
    preCovPca$aggPcaCorrDat$PC=paste(preCovPca$aggPcaCorrDat$subset,preCovPca$aggPcaCorrDat$origPC)
    preCovPca$aggPcaCorrDat$origPC_with_pct=preCovPca$aggPcaCorrDat$PC_with_pct
    preCovPca$aggPcaCorrDat$PC_with_pct=paste(preCovPca$aggPcaCorrDat$subset,preCovPca$aggPcaCorrDat$origPC_with_pct)
    
    mpdf("misc_DEG_PCA_correl", outDir=file.path(ROOT, "outputs"), width=10, height=7)
    plotPcaCovariateCorr(preCovPca$aggPcaCorrDat,paste0(myName,"_",myNameSuffix,"_unfiltered"))
    plotPcaCovariateCorr(preCovPca$aggPcaCorrDat[preCovPca$aggPcaCorrDat$minBH_AdjP<0.05,],paste0(myName,"_",myNameSuffix,"_FDR_0p05"))
    plotPcaCovariateCorr(preCovPca$aggPcaCorrDat[preCovPca$aggPcaCorrDat$minBH_AdjP<0.1,],paste0(myName,"_",myNameSuffix,"_FDR_0p1"))
    plotPcaCovariateCorr(preCovPca$aggPcaCorrDat[preCovPca$aggPcaCorrDat$minBH_AdjP<0.2,],paste0(myName,"_",myNameSuffix,"_FDR_0p2"))
    dev.off()
    
    signifCovs = unique(as.character(preCovPca$aggPcaCorrDat$covar[preCovPca$aggPcaCorrDat$targetedBH_AdjP < 0.05 & !is.na(preCovPca$aggPcaCorrDat$targetedBH_AdjP)]))
    
    z=unname(unlist(covariateInfo))
    z=data.frame(
      covar=z,
      dof=sapply(z,function(x){if(is.factor(allInfo[,x])){length(unique(allInfo[,x]))-1}else{1}},USE.NAMES=F),
      isBiologicalNumeric=z %in% covariateInfo$isBiologicalNumeric,
      isBiologicalFactor=z %in% covariateInfo$isBiologicalFactor,
      isTechnicalNumeric=z %in% covariateInfo$isTechnicalNumeric,
      isTechnicalFactor=z %in% covariateInfo$isTechnicalFactor,
      stringsAsFactors=F
    )
    
    myCats=c("isBiologicalNumeric","isBiologicalFactor","isTechnicalNumeric","isTechnicalFactor")
    for(i in c("covar",myCats))
      z=z[order(z[,i]),]
    if(any(!rowSums(z[,myCats])==1)) stop("Something is wrong with the classification of covariates")
    rm(i,myCats)
    
    preCovPca$signifCovs = z[z$covar %in% signifCovs,]
  }
}

########################################################################################
##### SELECTION OF COVARIATES USING BAYESIAN INFORMATION CRITERION (BIC) APPROACH ######

{
  # Initial model is "gene expression = Cell_type + Disease + Cell_type:Disease + Sex" ... note that "Cell_type + Disease + Cell_type + Cell_type:Disease" is encoded as Groups variable 
  baseModel = c("Groups", "Gender_asFactor")
  
  ######
  ### Variables of baseline model
  firstBaseModel = baseModel[1]
  remainingBaseModel = baseModel[-1]
  bicBaseModel = genericBicCalc(initialVoomObj$E, allInfo, remainingBaseModel, "baseModelEval", firstBaseModel, iterative=T, covarFracCutOff=-Inf)
  
  ######
  ### Single BIC model using numeric technical covariates
  techNumBaseModel = baseModel
  techNumModels = setdiff(preCovPca$signifCovs$covar[preCovPca$signifCovs$isTechnicalNumeric],techNumBaseModel)
  bicSingleTechNumModels = genericBicCalc(initialVoomObj$E,allInfo, techNumModels, "misc_BIC_techNumModels", techNumBaseModel, iterative=T, covarFracCutOff=my_covarFracCutOff, deltaBicCutOff=my_deltaBicCutOff)
  selTechNumModel = unname(unlist(bicSingleTechNumModels$selectedCovars))
  
  ######
  ### Single BIC model using categorical technical covariates
  techFactBaseModel = c(baseModel,selTechNumModel)
  techFactModels = setdiff(preCovPca$signifCovs$covar[preCovPca$signifCovs$isTechnicalFactor],techFactBaseModel)
  bicSingletechFactModels = genericBicCalc(initialVoomObj$E,allInfo,techFactModels,"misc_BIC_techFactModels",techFactBaseModel,iterative=T,covarFracCutOff=my_covarFracCutOff,deltaBicCutOff=my_deltaBicCutOff)
  selTechFactModel = unname(unlist(bicSingletechFactModels$selectedCovars))
  
  ######
  ### Put it all together
  finalBicModel = c(baseModel, selTechNumModel, selTechFactModel)
  finalBicModelCorrelation = generalizedColumnCorr(allInfo[,finalBicModel,drop=F],plotTitle="finalBicModel")
}

########################################################################################
##### PERFORM VARIANCE PARTITIONING (BEFORE ADDING COVARIATES) #########################

{
  varPartPreCovsModel = paste("~", paste(c(sapply(setdiff(c("Dx_asFactor", "cell_subtype_asFactor", finalBicModel), "Groups"), function(covar) ifelse(is.factor(allInfo[,covar]), paste0("(1|", covar,")"), covar)), "(1|Person_ID_asFactor)"), collapse=" + "))
  preCovVarPart = fitExtractVarPartModel(initialVoomObj, varPartPreCovsModel, allInfo, showWarnings=F)
  mpdf("misc_varPart_RNAseq_PRE_COVS",width=10,height=5); plotVarPart(preCovVarPart); dev.off()
}

########################################################################################
##### DIFFERENTIAL ANALYSIS ############################################################

{
  # Create Dx contrasts (per cell type as well as merged) and formula for differential analysis
  designFinalCovars = union("Groups", finalBicModel)
  designFinalString = paste("~ 0 +", paste(designFinalCovars,collapse=" + "))
  designFinalMatrix = model.matrix(as.formula(designFinalString),allInfo)
  contrastStringDf = genDxContrastStrings(allInfo, dxColumn="Dx")
  contrastMatrix = makeContrasts(contrasts=contrastStringDf$eq,levels=designFinalMatrix)
  colnames(contrastMatrix) = rownames(contrastStringDf)
  dreamFinalString = paste0(designFinalString, " + (1|Person_ID)")   # final formula also correct for inter-individual variance
  
  ######
  ### Voom normalizatino (with covariates)
  mpdf("misc_DEG_voomFinal"); modeledVoomObj = voomWithDreamWeights(initialDgeObj, dreamFinalString, allInfo, plot=T); dev.off()
  
  ######
  ### Differential analysis via dream
  fitDream = variancePartition::dream(modeledVoomObj, dreamFinalString, allInfo, contrastMatrix)
  degResults = analyzeAndPlotFit(fitDream, contrastStringDf, "misc_DEG_analyseDreamResults", allInfo, pAdjustVal=0.05, forceNoPlots=T)
  saveRDS(degResults, file=file.path(ROOT, "outputs", "misc_DEG_full_results.RDS"))
}

########################################################################################
##### PERFORM VARIANCE PARTITIONING (AFTER REGRESSING OUT COVARIATES) ##################

{
  # Helper function
  eval_residuals = function(form, form_full, vobj, METADATA, CPU_CORES=2) {
    f = function(fit) {
      residuals(fit) + variancePartition::get_prediction(fit, form)
    }
    i = match(colnames(vobj), rownames(METADATA))
    info = METADATA[i,]
    resid.lst = fitVarPartModel(vobj, form_full, info, showWarnings=T, fxn = f)
    do.call(rbind, resid.lst)
  }
  
  # Regress out the effect of covariates but (1) keep the effect Dx+Cell_type and(2) Cell_type
  residualizedModel_withGroups = paste("~", c(paste(c(designFinalCovars, "(1|Person_ID_asFactor)"), collapse=" + ")))
  count_matrix_residualized_Dx_CellType_kept = eval_residuals(~ Groups + (1|Person_ID_asFactor), residualizedModel_withGroups, modeledVoomObj, allInfo, CPU_CORES=5)
  saveRDS(count_matrix_residualized_Dx_CellType_kept, file=file.path(ROOT, "outputs", "misc_RNAseq_residualized_Dx_CellType_kept.RDS"))
  
  # Variance partition (using residualized matrix)
  postCovVarPart = fitExtractVarPartModel(count_matrix_residualized_Dx_CellType_kept, varPartPostCovsModel, allInfo, showWarnings=F)
  mpdf("misc_varPart_RNAseq_POST_COVS",width=10,height=5); plotVarPart(postCovVarPart); dev.off()
}

########################################################################################
##### ADDITIONAL PLOTTING (USING RESIDUALIZED COUNT MATRICES) ##########################

{
  ######
  ### Plotting MDS
  {
    myDist = as.dist(sqrt(1-cor(count_matrix_residualized_Dx_CellType_kept)^2)) # # Squared distance correlation
    mdsResults = cmdscale(myDist, k=2, eig=T)
    colnames(mdsResults$points) = c("Coordinate_1", "Coordinate_2")
    mdsResults$points = cbind(mdsResults$points, allInfo)
    
    mdsPostCovsPlot = ggplot(mdsResults$points, aes(x=Coordinate_1, y=Coordinate_2, shape=Dx, color=cell_subtype_abbreviation)) + geom_point(size=2) + ggtitle("MDS") +
      coord_fixed() + xlab("Coordinate 1") + ylab("Coordinate 2") + theme_classic() + theme(axis.text=element_text(colour="black"))
    
    mpdf("misc_DEG_MDS_postCovs", outDir=file.path(ROOT, "outputs")); print(mdsPostCovsPlot); dev.off()
  }
  
  ######
  ### Plotting tSNE
  {
    set.seed(sample(1:100, 1))
    tsne = Rtsne(t(count_matrix_residualized_Dx_CellType_kept), perplexity=floor(ncol(count_matrix_residualized_Dx_CellType_kept)/6), verbose=T)
    allInfo$tsne_dim_1 = tsne$Y[,1]
    allInfo$tsne_dim_2 = tsne$Y[,2]
    
    postCovTsne = ggplot(allInfo, aes(x=tsne_dim_1, y=tsne_dim_2, shape=Dx, color=cell_subtype_abbreviation)) +
      geom_point(size=2) + ggtitle("tSNE") + coord_fixed() + xlab("Dimension 1") + ylab("Dimension 2") + theme_classic() +
      theme(axis.text=element_text(colour="black"))
    mpdf("misc_DEG_tSNE_postCovs", outDir=file.path(ROOT, "outputs")); print(postCovTsne); dev.off()
  }
}
