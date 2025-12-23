library(ggplot2)

########################################################################################
##### CONFIG & PROLOGUE ################################################################


{
  ROOT = "~/Desktop/molecular_profiling/molecular_profiling/" # !!! FIXME: SET TO YOUR CUSTOM DIRECTORY !!!
  
  QC_ATACSEQ = file.path(ROOT, "inputs", "qc_all_atac.tsv")  # Pre-calculated QC metrics for ATAC-seq samples from processing computational pipeline
  QC_RNASEQ = file.path(ROOT, "inputs", "qc_all_rna.tsv")    # Pre-calculated QC metrics for RNA-seq samples from processing computational pipeline
  
  RNASEQ_COUNT_MATRIX_RAW = file.path(ROOT, "inputs", "rnaseq_count_matrix_raw.RDS")     # Raw read count matrix for RNA-seq data
  
  # Load metadata and count matrix
  rnaseq_countMatrixRaw = readRDS(RNASEQ_COUNT_MATRIX_RAW)
  allInfo_rnaseq = read.csv(QC_RNASEQ, sep="\t")
  allInfo_atacseq = read.csv(QC_ATACSEQ, sep="\t")
}


########################################################################################
##### HELPER FUNCTIONS #################################################################

{
  
}

########################################################################################
##### SELECTION OF COVARIATES USING BAYESIAN INFORMATION CRITERION (BIC) APPROACH ######

{
  covariateInfoRnaseq = list(
    ibBiologicalNumeric = c("AOD_scaled", "geno_PC1", "geno_PC2", "geno_PC3", "deconvolution_GABA", "deconvolution_GLU", "deconvolution_AST", "deconvolution_MG", "deconvolution_ODC"),
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
                           "rnaseqc_PCT_INTRONIC_BASES", "rnaseqc_PCT_INTERGENIC_BASES", "rnaseqc_PCT_MRNA_BASES", "rnaseqc_PCT_USABLE_BASES", "rnaseqc_MEDIAN_CV_COVERAGE", "rnaseqc_MEDIAN_5PRIME_BIAS", 
                           "rnaseqc_MEDIAN_3PRIME_BIAS", "rnaseqc_MEDIAN_5PRIME_TO_3PRIME_BIAS"),
    isTechnicalFactor = c())
  
  table(unlist(covariateInfoRnaseq) %in% colnames(allInfo_rnaseq))
  unlist(covariateInfoRnaseq)[!unlist(covariateInfoRnaseq) %in% colnames(allInfo_rnaseq)]
  table(allInfo_rnaseq$RIN)
  
  covariateInfoAtacseq = list(
    ibBiologicalNumeric = c("ageOfDeath_scaled", "geno_PC1", "geno_PC2", "geno_PC3", "deconvolution_GABA", "deconvolution_GLU", "deconvolution_OLIG", "deconvolution_MGAS"),
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
                           "star_Number_of_reads_mapped_to_multiple_loci_scaled", "star_Number_of_reads_mapped_to_too_many_loci_scaled", "star_Uniquely_mapped_reads_number_scaled", "PMI_scaled",),
    isTechnicalFactor = c("Barcode_combination_asFactor", "Illumina_index_1_asFactor", "Illumina_index_2_asFactor", "PCR_date_asFactor", "Pool_asFactor"))
  
  table(unlist(covariateInfo) %in% colnames(allInfo_atacseq))
  unlist(covariateInfo)[!unlist(covariateInfo) %in% colnames(allInfo_atacseq)]
  
  
}

########################################################################################
##### SELECTION OF COVARIATES USING BAYESIAN INFORMATION CRITERION (BIC) APPROACH ######

{
  # Initial model is "gene expression = Cell_type + Disease + Cell_type:Disease + Sex" ... note that "Cell_type + Disease + Cell_type + Cell_type:Disease" is encoded as Groups variable 
  baseModel = c("Groups", "Gender_asFactor")
  
  ############################################
  ## Variables of baseline model
  firstBaseModel = baseModel[1]
  remainingBaseModel = baseModel[-1]
  bicBaseModel = genericBicCalc(initialVoomObj$E,allInfo,remainingBaseModel,"baseModelEval",firstBaseModel,iterative=T,covarFracCutOff=-Inf)
  
  ############################################
  ## Single BIC model using numeric technical covariates
  techNumBaseModel=baseModel
  techNumModels=setdiff(preCovPca$signifCovs$covar[preCovPca$signifCovs$isTechnicalNumeric],techNumBaseModel)
  bicSingleTechNumModels=genericBicCalc(initialVoomObj$E,allInfo,techNumModels,"techNumModels",techNumBaseModel,iterative=T,covarFracCutOff=my_covarFracCutOff,deltaBicCutOff=my_deltaBicCutOff)
  selTechNumModel=unname(unlist(bicSingleTechNumModels$selectedCovars))
  
  ############################################
  ## Single BIC model using categorical technical covariates
  techFactBaseModel=c(baseModel,selTechNumModel)
  techFactModdels=setdiff(preCovPca$signifCovs$covar[preCovPca$signifCovs$isTechnicalFactor],techFactBaseModel)
  bicSingletechFactModels=genericBicCalc(initialVoomObj$E,allInfo,techFactModdels,"techFactModdels",techFactBaseModel,iterative=T,covarFracCutOff=my_covarFracCutOff,deltaBicCutOff=my_deltaBicCutOff)
  selTechFactModel=unname(unlist(bicSingletechFactModels$selectedCovars))
  
  ############################################
  ## Single BIC model using numeric technical covariates squared
  techNumSqBaseModel=c(baseModel,selTechNumModel,selTechFactModel)
  allInfo=squareModelCols(allInfo,selTechNumModel)$myDf #add squared cols to table
  techNumSqModels=squareModelCols(allInfo,selTechNumModel)$newCol #get the names of the squared cols
  bicSingletechNumSqModels=genericBicCalc(initialVoomObj$E,allInfo,techNumSqModels,"techNumSqModels",techNumSqBaseModel,iterative=T,covarFracCutOff=my_covarFracCutOff,deltaBicCutOff=my_deltaBicCutOff)
  selTechNumSqModel=unname(unlist(bicSingletechNumSqModels$selectedCovars))
  
  ############################################
  ## Put it all together
  finalBicModel=c(baseModel,selTechNumModel,selTechFactModel,selTechNumSqModel)
  finalBicModelCorrelation=generalizedColumnCorr(allInfo[,finalBicModel,drop=F],plotTitle="finalBicModel")
  
  
}

