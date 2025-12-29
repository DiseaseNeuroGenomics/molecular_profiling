library(ggplot2)
library(edgeR)
library(Rtsne)
library(variancePartition)

########################################################################################
##### CONFIG & PROLOGUE ################################################################

{
  ROOT = "~/Desktop/molecular_profiling/molecular_profiling/" # !!! FIXME: SET TO YOUR CUSTOM DIRECTORY !!!
  CPU_CORES = 8                                               # !!! FIXME: SET THE DESIRED NUMBER OF CORES !!!
  
  # Load metadata
  METADATA = file.path(ROOT, "inputs", "qc_all_atac.csv")
  allInfo = read.csv(METADATA)
  rownames(allInfo) = allInfo$ID
  allInfo$Dx_asFactor = ordered(allInfo$Dx, levels=c("Control", "SCZ"))
  allInfo$cell_subtype_asFactor = ordered(allInfo$cell_subtype_abbreviation, levels=c("GABA", "GLU", "OLIG", "MGAS"))
  allInfo$Gender_asFactor = ordered(allInfo$Gender, levels=c("Male", "Female"))
  allInfo$Person_ID_asFactor = as.factor(allInfo$Person_ID)
  
  # Load read count matrix
  ATACSEQ_COUNT_MATRIX_RAW = file.path(ROOT, "inputs", "atacseq_count_matrix_raw.RDS")     # Raw read count matrix for ATAC-seq data
  countMatrixRaw = readRDS(ATACSEQ_COUNT_MATRIX_RAW)
  
  # Load peaks
  peaks = readRDS(file.path(ROOT, "inputs", "atacseq_peaks.RDS"))
  rownames(peaks) = peaks$PeakID
  
  # Default output dir
  outDir=file.path(ROOT, "outputs")
  
  # BIC-related settings for covariate selection: the selection procedure is repeated until the best-performing covariates from selection pool improved at least (my_covarFracCutOff e.g. 5%)% of genes/peaks by (my_deltaBicCutOff, e.g. 2)
  my_covarFracCutOff = 0.05
  my_deltaBicCutOff = 2
  
  ### Initial set of covariates for exploration (ATAC-seq)
  covariateInfo = list(
    isBiologicalNumeric = c("ageOfDeath_scaled", "geno_PC1", "geno_PC2", "geno_PC3", "deconvolution_GABA", "deconvolution_GLU", "deconvolution_OLIG", "deconvolution_MGAS"),
    isBiologicalFactor = c("cell_subtype_asFactor", "Dx_asFactor", "Gender_asFactor"),
    isTechnicalNumeric = c("Yield_of_nuclei_cScaled", "chrMFrac_cScaled", "chrMReads_cScaled", "fracReadsInNonBlacklistedPeaks_cScaled", "fracReadsInOnlyBlacklistedPeaks_cScaled",
                           "fracReadsInOnlyBlacklistedPeaksVsInAllPeaks_cScaled", "peakGappedCount_cScaled", "peakNarrowFDR1pctCount_cScaled", "ctcf_fos_cScaled", "spp_NSC_cScaled", "spp_RSC_cScaled",
                           "insertMetrics_MEDIAN_INSERT_SIZE_cScaled",  "insertMetrics_MEDIAN_ABSOLUTE_DEVIATION_cScaled", "insertMetrics_MEAN_INSERT_SIZE_cScaled", "insertMetrics_STANDARD_DEVIATION_cScaled",
                           "insertMetrics_WIDTH_OF_10_PERCENT_cScaled", "insertMetrics_WIDTH_OF_20_PERCENT_cScaled", "insertMetrics_WIDTH_OF_30_PERCENT_cScaled", "insertMetrics_WIDTH_OF_40_PERCENT_cScaled", 
                           "insertMetrics_WIDTH_OF_50_PERCENT_cScaled", "insertMetrics_WIDTH_OF_60_PERCENT_cScaled", "insertMetrics_WIDTH_OF_70_PERCENT_cScaled", "insertMetrics_WIDTH_OF_80_PERCENT_cScaled",
                           "insertMetrics_WIDTH_OF_90_PERCENT_cScaled", "picard_meanGcContent_cScaled", "picard_AT_DROPOUT_cScaled", "picard_GC_NC_0_19_cScaled", "picard_GC_NC_20_39_cScaled",
                           "picard_GC_NC_40_59_cScaled", "picard_GC_NC_60_79_cScaled", "picard_GC_NC_80_100_cScaled", "star_pct_of_reads_mapped_to_multiple_loci_scaled", "star_pct_of_reads_mapped_to_too_many_loci_scaled",
                           "star_pct_of_reads_unmapped_other_scaled", "star_pct_of_reads_unmapped_too_short_scaled", "star_Uniquely_mapped_reads_pct_scaled", "Total_PCR_cycles_scaled", "finalReadCount_scaled",
                           "finalReadCountFrac_scaled", "pbc_scaled", "picard_PERCENT_DUPLICATION_scaled", "star_Average_mapped_length_scaled", "star_Number_of_input_reads_scaled", 
                           "star_Number_of_reads_mapped_to_multiple_loci_scaled", "star_Number_of_reads_mapped_to_too_many_loci_scaled", "star_Uniquely_mapped_reads_number_scaled", "PMI_scaled"),
    isTechnicalFactor = c("Barcode_combination_asFactor", "Illumina_index_1_asFactor", "Illumina_index_2_asFactor", "PCR_date_asFactor", "Pool_asFactor"))
  
  
  # Load helper scripts
  source(file.path(ROOT, "helper_functions.R"))
  options(mc.cores=CPU_CORES)
}


########################################################################################
##### READ COUNT NORMALIZATION (TMM METHOD) & tSNE & ANALYSIS OF COVARIATES ############

{
  # Make edgeR object:
  expObjAll = DGEList(counts=countMatrixRaw, genes=peaks[match(rownames(countMatrixRaw),rownames(peaks)),])
  
  # Keep peaks with at least 1 count-per-million reads (cpm) in at least 20% of the samples
  fracSamplesWithMinCPM = rowMeans(cpm(expObjAll) >= 1)
  isNonLowExpr = fracSamplesWithMinCPM >= 0.2
  expObjNonLow = expObjAll[isNonLowExpr, , keep.lib.sizes=F]
  
  geneNormObj = calcNormFactors(expObjNonLow, method="TMM")
  
  # Voom normalization
  initialDgeObj = geneNormObj
  
  mpdf("misc_DAC_voomFirstPlot", outDir=file.path(ROOT, "outputs"))
  initialVoomObj = voom(initialDgeObj, design=NULL, plot=T)
  dev.off()
  
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
    mpdf("misc_DAC_tSNE_preCovs", outDir=file.path(ROOT, "outputs")); print(myPlot); dev.off()
  }
  
  ######
  ### Calculating which covariates are significantly associated (at FDR<0.20) with at least one PC of peak expression variance that explain at least 1pct of variance
  {
    covariatesForExploration = unname(unlist(covariateInfo[c("isTechnicalNumeric","isTechnicalFactor")]))
    
    myName = "ATAC-seq covariate exploration"
    myNameSuffix="aggPca"
    
    preCovPca = list()
    preCovPca$pca_all = calcAndPlotPCAs(unlist(unname(covariateInfo)), initialVoomObj$E, allInfo, "all", myName, covariateColsToAccountForInMultipleTesting=covariatesForExploration)
    preCovPca$aggPcaCorrDat=preCovPca$pca_all$pcaCorr
    preCovPca$aggPcaCorrDat$origPC=preCovPca$aggPcaCorrDat$PC
    preCovPca$aggPcaCorrDat$PC=paste(preCovPca$aggPcaCorrDat$subset,preCovPca$aggPcaCorrDat$origPC)
    preCovPca$aggPcaCorrDat$origPC_with_pct=preCovPca$aggPcaCorrDat$PC_with_pct
    preCovPca$aggPcaCorrDat$PC_with_pct=paste(preCovPca$aggPcaCorrDat$subset,preCovPca$aggPcaCorrDat$origPC_with_pct)
    
    mpdf("misc_DAC_PCA_correl", outDir=file.path(ROOT, "outputs"), width=10, height=7)
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
  # Initial model is "chromatin accessibility = Cell_type + Disease + Cell_type:Disease + Sex" ... note that "Cell_type + Disease + Cell_type + Cell_type:Disease" is encoded as Groups variable 
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
  modelForVarpart = setdiff(c(finalBicModel, "Dx"), "Groups")
  
  varPartPreCovsModel = paste("~", paste(c(sapply(setdiff(c("Dx_asFactor", "cell_subtype_asFactor", finalBicModel), "Groups"), function(covar) ifelse(is.factor(allInfo[,covar]), paste0("(1|", covar,")"), covar)), "(1|Person_ID_asFactor)"), collapse=" + "))
  preCovVarPart = fitExtractVarPartModel(initialVoomObj, varPartPreCovsModel, allInfo, showWarnings=F)
  mpdf("misc_varPart_ATACseq_PRE_COVS",width=10,height=5); plotVarPart(preCovVarPart); dev.off()
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
  mpdf("misc_DAC_voomFinal"); modeledVoomObj = voomWithDreamWeights(initialDgeObj, dreamFinalString, allInfo, plot=T); dev.off()
  
  ######
  ### Differential analysis via dream
  fitDream = variancePartition::dream(modeledVoomObj, dreamFinalString, allInfo, contrastMatrix)
  dacResults = analyzeAndPlotFit(fitDream, contrastStringDf, "misc_DAC_analyseDreamResults", allInfo, pAdjustVal=0.05, forceNoPlots=T)
  saveRDS(dacResults, file=file.path(ROOT, "outputs", "misc_DAC_full_results.RDS"))
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
  saveRDS(count_matrix_residualized_Dx_CellType_kept, file=file.path(ROOT, "outputs", "misc_ATACseq_residualized_Dx_CellType_kept.RDS"))
  
  # Variance partition (using residualized matrix)
  postCovVarPart = fitExtractVarPartModel(count_matrix_residualized_Dx_CellType_kept, varPartPreCovsModel, allInfo, showWarnings=F)
  mpdf("misc_varPart_ATACseq_POST_COVS",width=10,height=5); plotVarPart(postCovVarPart); dev.off()
}

########################################################################################
##### ADDITIONAL PLOTTING (USING RESIDUALIZED COUNT MATRICES) ##########################

{
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
    mpdf("misc_DAC_tSNE_postCovs", outDir=file.path(ROOT, "outputs")); print(postCovTsne); dev.off()
  }
}
