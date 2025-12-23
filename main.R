library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggsci)
library(MASS)
library(viridis)
library(reshape)
library(readxl)
library(RColorBrewer)
library(patchwork)
library(tidyr)

########################################################################################
##### CONFIG ###########################################################################

{
  ROOT = "~/Desktop/molecular_profiling/molecular_profiling/" # !!! FIXME: SET TO YOUR CUSTOM DIRECTORY !!!
  
  QC_ATACSEQ = file.path(ROOT, "inputs", "qc_all_atac.csv")  # Pre-calculated QC metrics for ATAC-seq samples from processing computational pipeline
  QC_RNASEQ = file.path(ROOT, "inputs", "qc_all_rna.csv")    # Pre-calculated QC metrics for RNA-seq samples from processing computational pipeline
  KINSHIP_ATACSEQ_SNPPARRAY = file.path(ROOT, "inputs", "kinship_atacseq_snparray.csv")  # Pre-calculated comparison between SNPs called from ATAC-seq reads and SNParrays 
  KINSHIP_RNASEQ_SNPPARRAY = file.path(ROOT, "inputs", "kinship_rnaseq_snparray.csv")    # Pre-calculated comparison between SNPs called from RNA-seq reads and SNParrays 
  
  ATACSEQ_PEAKS =  file.path(ROOT, "inputs", "atacseq_peaks.RDS")                        # Peaks called from ATAC-seq data
  ATACSEQ_COUNT_MATRIX_RAW = file.path(ROOT, "inputs", "atacseq_count_matrix_raw.RDS")   # Raw read count matrix for ATAC-seq data
  ATACSEQ_COUNT_MATRIX_ADJ = file.path(ROOT, "inputs", "atacseq_count_matrix_adj.RDS")   # Covariate-adjusted read count matrix for ATAC-seq data
  ATACSEQ_COUNT_MATRIX_RESIDUALIZED_DX_CELLTYPE_KEPT = file.path(ROOT, "inputs", "atacseq_count_matrix_residualized_CellType_kept.RDS")   # Count matrix from which the effect of technical covariates were regressed out, but Dx & Cell type effect kept
  RNASEQ_COUNT_MATRIX_RAW = file.path(ROOT, "inputs", "rnaseq_count_matrix_raw.RDS")     # Raw read count matrix for RNA-seq data
  RNASEQ_COUNT_MATRIX_ADJ = file.path(ROOT, "inputs", "rnaseq_count_matrix_adj.RDS")     # Covariate-adjusted read count matrix for RNA-seq data
  RNASEQ_COUNT_MATRIX_RESIDUALIZED_DX_CELLTYPE_KEPT = file.path(ROOT, "inputs", "rnaseq_count_matrix_residualized_CellType_kept.RDS")   # Count matrix from which the effect of technical covariates were regressed out, but Dx & Cell type effect kept
  
  DAC_ANALYSIS = file.path(ROOT, "inputs", "DAC_Analysis.Rdata")  # Pre-calculated results for analysis of differential chromatin accessibility
  DEG_ANALYSIS = file.path(ROOT, "inputs", "DEG_Analysis.Rdata")  # Pre-calculated results for analysis of differential gene expression 
  DET_ANALYSIS = file.path(ROOT, "inputs", "DET_Analysis.Rdata")  # Pre-calculated results for analysis of differential transcript expression 
  REMACOR_ANALYSIS = file.path(ROOT, "inputs", "REMACOR_ANALYSIS.xlsx")
  
  DEG_ANALYSIS_PSYCHAD_c07x = file.path(ROOT, "inputs", "DEG_Analysis_PsychAD_c07x.Rdata")  # Pre-calculated results for differential SCZ case-control analysis from PsychAD paper (Lee et atl 2025); contrast c07x
  
  METADATA_HAUBERG_2020 = file.path(ROOT, "inputs", "hauberg_2020_metadata.csv")  # Metadata for samples from Hauberg et al 2020 (dataset used for comparison)
  GEXPR_HAUBERG_2020 = file.path(ROOT, "inputs", "hauberg_2020_gExpr.RDS")        # Covariate-adjusted read count matrix for FANS ATAC-seq data from Hauberg et al 2020
  PEAKS_HAUBERG_2020 = file.path(ROOT, "inputs", "hauberg_2020_peaks.RDS")        # Peaks called from ATAC-seq data from Hauberg et al 2020
  
  METADATA_COLEMAN_2023 = file.path(ROOT, "inputs", "coleman_2023_metadata.csv")  # Metadata for samples from Coleman et al 2023 (dataset used for comparison)
  GEXPR_COLEMAN_2023 = file.path(ROOT, "inputs", "coleman_2023_gExpr.RDS")        # Covariate-adjuste read count matrix from FANS RNA-seq data from Coleman et al 2023
  
  ENSEMBL_INFO = file.path(ROOT, "inputs", "muchEnsemblInfo_hg38.tsv.gz")    # Ensembl (canonical genes, not transcripts)
  GENCODE = file.path(ROOT, "inputs", "gencode.v30.annotation.gtf.gz")       # GENCODE (version 30)
  
  npgList = list("NEURON"="#B2182B", "GLIA"="#2166AC",             # (a) ATAC-seq cell type
                 "green"="#67A61A", "yellow"="#E4AB00", "pink"="#E4288A", "gray"="727272",
                 "DLPFC"  ="#E7298A", "ACC"="#66A61E",             # (b) ATAC-seq brain region
                 "DLPFC_NEURON"="#F39B7F", "ACC_NEURON"="#C77B85", # (c) ATAC-seq cell type & brain region
                 "DLPFC_GLIA"  ="#4DBBD5", "ACC_GLIA"="#9DB9D4",   # (c) ATAC-seq cell type & brain region
                 "FP"="#E6AB02", "ACC"="#E7298A",                  # (d) RNA-seq brain region
                 "DLPFC"="#66A61E", "IFG"="#7570B3",               # (d) RNA-seq brain region
                 "Promoter"="#E6AB02", "Intron"="#66A61E",
                 "Distal_Intergenic"="#E7298A", "Exon"="#7570B3",
                 "unique"="#946317",         "all"="#666666",      # (e) Differentially accessible peaks
                 "SCZ_Control"="#7570B3", "BP_Control"="#E7298A",  # (f) Phenotypes
                 "GLU" = "#E6AB02", "GLU2", "#FFF4D6",
                 "GABA" = "#66A61E", "GABA2", "#DBEFC4",
                 "OLIG" = "#E7298A", "OLIG2", "#FFC9E5",
                 "MGAS" = "#7570B3", "MGAS2", "#C1BEE2")
  
  myPalette = colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"), space = "Lab")
  
  dir.create(file.path(ROOT, "outputs"))
}

########################################################################################
##### HELPER FUNCTIONS #################################################################

{
  mpdf = function(x, width=7,height=7, outDir=outDir, onefile=T) eval.parent(substitute({ pdf(paste0(outDir, "/", make.names(x),".pdf"), useDingbats=F, width=width, height=height, onefile=onefile) }))
  mtsv = function(x, myHeader=T, filename=NULL, outDir=outDir, myRownames=F){ if(is.null(filename)) { filename = make.names(deparse(substitute(x)))}; write.table(x, file=paste0(outDir, "/", filename,".tsv"), na="", sep="\t", quote=F, row.names=myRownames, col.names=myHeader) }
  
  # Convenient way to handle errors
  myStop=function(...) eval.parent(substitute({
    debugEnv <<- as.environment(as.list(environment(), all.names=T))
    stop(paste0(...,". Script aborted. Variables saved to the 'debugEnv' environment for debugging purposes. To access the variables use the normal variable name preceeded by this and a dollar sign. For instance you can have a look at the big table with much information using 'debugEnv$allInfo'."),call.=F)
  }))
  
  # Function to create DGE object with filtering, normalization, and plotting
  getGeneFilteredGeneExprMatrix=function(
      rawReadCountsNoBlacklist,
      allInfo,
      qcPeakAnno,
      MIN_GENE_CPM=1,
      MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.1,
      calcNormFactors.method="codingPromoterTMM", #alternatives are standard edgeR norm methods
      geneTssPeakMapping=NULL,
      housekeepingPeakInfo=NULL, #currently just used for plotting
      plotName=NULL #if provided make a plot with this name
  ){
    w=list()
    
    #allInfo must be sorted
    if(any(order(allInfo$ID)!=seq(nrow(allInfo)))) myStop("allInfo must be sorted by ID column")
    
    #align metadata and expression mat
    if(any(!allInfo$ID %in% colnames(rawReadCountsNoBlacklist))) myStop("metadata and expression data does not line up")
    datExpr=data.frame(rawReadCountsNoBlacklist[,allInfo$ID])
    
    #check that we have metadata
    if(any(!rownames(rawReadCountsNoBlacklist) %in% qcPeakAnno$PeakID)) myStop("one or more peak wasn't found in the provided qcPeakAnno")
    
    #Make edgeR object:
    expObjAll=DGEList(
      counts=datExpr,
      genes=qcPeakAnno[match(rownames(datExpr),qcPeakAnno$PeakID),]
    )
    
    #Keep genes with at least MIN_GENE_CPM count-per-million reads (cpm) in at least (MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM)% of the samples:
    fracSamplesWithMinCPM=rowMeans(cpm(expObjAll) >= MIN_GENE_CPM)
    isNonLowExpr=fracSamplesWithMinCPM >= MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM
    expObjNonLow=expObjAll[isNonLowExpr, , keep.lib.sizes=F]
    
    message(paste0("\nWill normalize expression counts for ", sum(isNonLowExpr), " out of ", length(isNonLowExpr), " OCRs (", sum(!isNonLowExpr), " OCRs discarded, which is ", sprintf("%.2f", 100 *  sum(!isNonLowExpr)/length(isNonLowExpr)), "%)"))
    message(paste0("The OCRs that we keep are those with a minimum of ", MIN_GENE_CPM, " CPM in at least ", sprintf("%.2f", 100 * MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM), "% of the ", ncol(expObjNonLow), " samples."))
    
    if(calcNormFactors.method=="codingPromoterTMM"){
      if(is.null(geneTssPeakMapping)) stop("For codingPromoterTMM normalization, a geneTssPeakMapping must be provided")
      w$protCodingPeaks=geneTssPeakMapping$genePeakMapping$PeakID[geneTssPeakMapping$genePeakMapping$Gene.type=="protein_coding"]
      w=c(w,targetedNormalization(expObjNonLow,w$protCodingPeaks))
    }else{
      w$dgeObj=calcNormFactors(expObjNonLow, method=calcNormFactors.method)
    }
    
    #normalization plots
    if(!is.null(plotName) & !is.null(housekeepingPeakInfo))
      normalizationPlots(w$dgeObj,allInfo,plotName,paste0(outDir,"/normPlots"),peaksToHighlight=housekeepingPeakInfo,highlightingName="housekeeping")
    
    #cpm plots
    gRes=ggplot(data.frame(fracSamplesWithMinCPM=as.numeric(fracSamplesWithMinCPM),stringsAsFactors=F), aes(x=fracSamplesWithMinCPM)) +
      geom_vline(xintercept=MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM, linetype="solid", col="red") +
      geom_histogram(color="black", fill="white", binwidth=0.02) +
      xlab(paste0("Fraction of samples with at least ", MIN_GENE_CPM, " CPM")) + ylab("# of OCRs")
    
    #Optionally make that plot
    if(!is.null(plotName)){
      mpdf(paste0("OCR_CPM_HIST_",plotName), width=8, height=8); print(gRes); dev.off()
    }
    
    return(w)
  }
  
  # Calculate density for density plots
  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  shrinkAtacToptables = function(x) {
    sapply(x,function(x)x[,!colnames(x) %in% c("annotation", "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand", "transcriptId", "distanceToTSS", "Gene.type", "Description")],simplify=F)
  }
  
  addDacExclStats = function(df,dac,setSize=NA,alpha=0.05){
    df$signifCount=sapply(df$Set,function(x)sum(dac[[x]]$adj.P.Val < alpha),USE.NAMES=F)
    df$allCount=setSize # sapply(df$Set,function(x)nrow(dac[[x]]),USE.NAMES=F)
    df$frac=df$signifCount/df$allCount
    df
  }
  
  # Ensembl Gene name and indentifier converter
  symbol_to_ensembl <- function(gene_symbols, anno_df) {
    ids <- anno_df$PeakID[match(gene_symbols, anno_df$gene_name)]
    # drop NAs and keep unique
    ids <- unique(ids[!is.na(ids)])
    ids
  }
  
  # Load Ensembl (and custom adjustment)
  ensemblInfo = read.csv(ENSEMBL_INFO, sep="\t", header=T, stringsAsFactors=F)
  ensemblInfo = ensemblInfo[,c("Ensembl.Gene.ID", "Transcript.length..including.UTRs.and.CDS.", "Gene.type", "Gene...GC.content")]
  colnames(ensemblInfo) = c("gene_id", "Length", "transcript_biotype", "gcContent")
  
  # Load GENCODE (and custom adjustment)
  gtf = data.frame(rtracklayer::import(GENCODE))
  gtf = gtf[which(gtf$type == "gene"),]
  gtf = gtf[,c("gene_id", "seqnames", "start", "end", "strand", "gene_name", "transcript_name", "width", "gene_type")]
  gtf$exonLength = ensemblInfo[match(sapply(strsplit(as.character(gtf[,"gene_id"]), "\\."), "[[", 1), ensemblInfo[,"gene_id"]),]$Length
  gtf$biotype = ensemblInfo[match(sapply(strsplit(as.character(gtf[,"gene_id"]), "\\."), "[[", 1), ensemblInfo[,"gene_id"]),]$transcript_biotype
  gtf$gcContent = ensemblInfo[match(sapply(strsplit(as.character(gtf[,"gene_id"]), "\\."), "[[", 1), ensemblInfo[,"gene_id"]),]$gcContent
  gtf = gtf[!is.na(match(sapply(strsplit(as.character(gtf[,"gene_id"]), "\\."), "[[", 1), ensemblInfo[,"gene_id"])),]
  gtf$PeakID = sapply(strsplit(as.character(gtf$PeakID), "\\."), "[[", 1)
  gtf = gtf[!duplicated(gtf$PeakID),]
  rownames(gtf) = gtf$PeakID
}

####################################################################################################
##### FIG. S1 :: DEMOGRAPHIC AND CLINICAL CHARACTERISTICS OF SCZ CASES AND CONTROLS ################

{
  tmpRna = read.csv(QC_RNASEQ, sep="\t")
  tmpAtac = read.csv(QC_ATACSEQ, sep="\t")
  
  ###################################
  #TODO: custom fixes: correct CSV files later
  deconvRnaseq = read.csv("/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/rnaseq/analysis/gene_cell_all__norm_none__BIC_2_in_0.05__CPM_1__in_0.2/deconv_atacseq.csv")
  genoPcRnaseq = read.csv("/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/rnaseq/analysis/gene_cell_all__norm_none__BIC_2_in_0.05__CPM_1__in_0.2/files/reviewer_ancestryPC/allInfo.csv")
  qcRna$deconvolution_GABA = deconvRnaseq[match(qcRna$ID, deconvRnaseq$ID),"deconvolution_GABA"]
  qcRna$deconvolution_GLU = deconvRnaseq[match(qcRna$ID, deconvRnaseq$ID),"deconvolution_GLU"]
  qcRna$deconvolution_AST = deconvRnaseq[match(qcRna$ID, deconvRnaseq$ID),"deconvolution_AST"]
  qcRna$deconvolution_MG = deconvRnaseq[match(qcRna$ID, deconvRnaseq$ID),"deconvolution_MG"]
  qcRna$deconvolution_ODC = deconvRnaseq[match(qcRna$ID, deconvRnaseq$ID),"deconvolution_ODC"]
  qcRna$geno_PC1 = genoPcRnaseq[match(qcRna$ID, genoPcRnaseq$ID),"geno_PC1"]
  qcRna$geno_PC2 = genoPcRnaseq[match(qcRna$ID, genoPcRnaseq$ID),"geno_PC2"]
  qcRna$geno_PC3 = genoPcRnaseq[match(qcRna$ID, genoPcRnaseq$ID),"geno_PC3"]
  qcRna$cell_subtype_abbreviation = toupper(qcRna$cell_subtype_abbreviation)
  write.csv(qcRna, file="/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/atacseq/analysis/all/deconv_rnaseq.csv", row.names=F)
  write.csv(qcRna, file="~/Desktop/molecular_profiling/molecular_profiling//inputs/qc_all_rna.csv", row.names=F)
  
  deconvAtacseq = read.csv("/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/atacseq/analysis/all/deconv_rnaseq.csv")
  genoPcAtacseq = read.csv("/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/atacseq/analysis/all/reviewer_ancestryPC/allInfo.csv")
  qcAtac$deconvolution_GABA = deconvAtacseq[match(qcAtac$ID, deconvAtacseq$ID),"deconvolution_GABA"]
  qcAtac$deconvolution_GLU = deconvAtacseq[match(qcAtac$ID, deconvAtacseq$ID),"deconvolution_GLU"]
  qcAtac$deconvolution_OLIG = deconvAtacseq[match(qcAtac$ID, deconvAtacseq$ID),"deconvolution_OLIG"]
  qcAtac$deconvolution_MGAS = deconvAtacseq[match(qcAtac$ID, deconvAtacseq$ID),"deconvolution_MGAS"]
  qcAtac$geno_PC1 = genoPcAtacseq[match(qcAtac$ID, gsub("Sample_", "", genoPcAtacseq$ID)),"geno_PC1"]
  qcAtac$geno_PC2 = genoPcAtacseq[match(qcAtac$ID, gsub("Sample_", "", genoPcAtacseq$ID)),"geno_PC2"]
  qcAtac$geno_PC3 = genoPcAtacseq[match(qcAtac$ID, gsub("Sample_", "", genoPcAtacseq$ID)),"geno_PC3"]
  qcAtac$cell_subtype_abbreviation = gsub("GABAergic", "GABA", gsub("glutamatergic", "GLU", gsub("oligodendrocytes", "OLIG", gsub("microgliaAndAstrocytes", "MGAS", qcAtac$cell_subtype_abbreviation))))
  write.csv(qcAtac, file="/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/rnaseq/analysis/gene_cell_all__norm_none__BIC_2_in_0.05__CPM_1__in_0.2/deconv_atacseq.csv", row.names=F)
  write.csv(qcAtac, file="~/Desktop/molecular_profiling/molecular_profiling//inputs/qc_all_atac.csv", row.names=F)
  ###################################
  
  # Custom fixes
  tmpRna$RIN = ifelse(tmpRna$RIN>10, NA, tmpRna$RIN)        # there was an outlier value "105" which was clearly a typo
  cols = intersect(colnames(tmpRna), colnames(tmpAtac))
  tmpAll = rbind.data.frame(tmpRna[,cols], tmpAtac[,cols])
  tmpAll = tmpAll[!duplicated(tmpAll$Person_ID),]
  tmpAll$Ancestry = gsub("Hispanic", "AMR", gsub("African-American", "AFR", gsub("Asian", "AS", gsub("Caucasian", "EUR", tmpAll$Ethnicity))))

  # Plot (density): Sex-by-Age
  tmpFiltered <- tmpAll %>% filter(Sex %in% c("XX", "XY")) %>% mutate(SexLabel = ifelse(Sex == "XX", "Female", "Male"))
  sexByAgePlot = ggplot(tmpFiltered, aes(x = ageOfDeath, fill = SexLabel, color = SexLabel)) + geom_density(alpha = 0.4, adjust = 1.2) + 
    geom_vline(data = tmpFiltered %>% group_by(SexLabel) %>% summarise(m = mean(ageOfDeath)),
               aes(xintercept = m, color = SexLabel), linetype = "dashed", size = 1) + scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
    scale_color_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) + labs(x = "Age of death", y = "Density", fill = NULL, color = NULL) +
    theme_minimal(base_size = 14) + theme(legend.position = "bottom", legend.title = element_blank(), panel.grid.minor = element_blank())

  # Plot (pie): Ancestry distribution
  ancestry_counts <- tmpAll %>% filter(Ancestry %in% c("AFR", "AMR", "EUR", "AS")) %>% count(Ancestry) %>% mutate(pct = round(n / sum(n) * 100), label = paste0(pct, "%"))
  ancestry_colors <- c("AFR" = "#E69F00", "AMR" = "#56B4E9", "EUR" = "#009E73", "AS" = "#F0E442")
  ancestryPiePlot <- ggplot(ancestry_counts, aes(x = "", y = n, fill = Ancestry)) + geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar(theta = "y") + geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5) + scale_fill_manual(values = ancestry_colors) +
    theme_void(base_size = 14) + theme(legend.position = "bottom", legend.title = element_blank())
  
  # Plot (pie): Sex distribution
  sex_counts <- tmpAll %>% filter(Sex %in% c("XX", "XY")) %>% mutate(SexLabel = ifelse(Sex == "XX", "Female", "Male")) %>% count(SexLabel) %>%
    mutate(pct = round(n / sum(n) * 100), label = paste0(pct, "%"))
  sex_colors <- c("Female" = "#E69F00", "Male" = "#56B4E9") # Define colors
  sexPiePlot <- ggplot(sex_counts, aes(x = "", y = n, fill = SexLabel)) + geom_bar(stat = "identity", width = 1, color = "black") + coord_polar(theta = "y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5) + scale_fill_manual(values = sex_colors) + theme_void(base_size = 14) +
    theme(legend.position = "bottom", legend.title = element_blank())

  # Plot (density): Dx-by-age plot
  dx_filtered <- tmpAll %>% filter(Dx %in% c("Control", "SCZ")) %>% mutate(DxLabel = ifelse(Dx == "SCZ", "SCZ", "Control"))
  dx_colors <- c("SCZ" = "#56B4E9", "Control" = "#E69F00") # Define colors matching the figure
  dxByAgePlot <- ggplot(dx_filtered, aes(x = ageOfDeath, fill = DxLabel, color = DxLabel)) + geom_density(alpha = 0.4, adjust = 1.2) +
    geom_vline(data = dx_filtered %>% group_by(DxLabel) %>% summarise(m = mean(ageOfDeath)), aes(xintercept = m, color = DxLabel), linetype = "dashed", size = 1) +
    scale_fill_manual(values = dx_colors) + scale_color_manual(values = dx_colors) + labs(x = "Age of death", y = "Density", fill = NULL, color = NULL) +
    theme_minimal(base_size = 14) + theme(legend.position = "bottom", legend.title = element_blank(), panel.grid.minor = element_blank())

  # Plot (density): Dx-by-pH plot
  ph_filtered <- tmpAll %>% filter(Dx %in% c("Control", "SCZ"), !is.na(pH)) %>% mutate(DxLabel = ifelse(Dx == "SCZ", "SCZ", "Control"))
  dx_colors <- c("SCZ" = "#56B4E9", "Control" = "#E69F00") # Define colors
  dxByPhPlot <- ggplot(ph_filtered, aes(x = pH, fill = DxLabel, color = DxLabel)) + geom_density(alpha = 0.4, adjust = 1.2) +
    geom_vline(data = ph_filtered %>% group_by(DxLabel) %>% summarise(m = mean(pH)), aes(xintercept = m, color = DxLabel), linetype = "dashed", size = 1) +
    scale_fill_manual(values = dx_colors) + scale_color_manual(values = dx_colors) + labs(x = "pH", y = "Density", fill = NULL, color = NULL) +
    theme_minimal(base_size = 14) + theme(legend.position = "bottom", legend.title = element_blank(), panel.grid.minor = element_blank())

  # Plot (pie): Antipsychotics use for SCZ samples
  ap_pie_data <- tmpAll %>% filter(Dx == "SCZ") %>% mutate(AP_Category = case_when( AntipsychAtyp & AntipsychTyp ~ "Both",
                                                                                    AntipsychAtyp & !AntipsychTyp ~ "Atyp only", !AntipsychAtyp & AntipsychTyp ~ "Typ only", TRUE ~ "None" )) %>%
  count(AP_Category) %>% mutate(pct = round(n / sum(n) * 100), label = paste0(pct, "%"))
  ap_colors <- c("Both" = "#D55E00", "Atyp only" = "#E69F00", "Typ only" = "#56B4E9", "None" = "#009E73")  # Define colors for each category
  antipsychPlot <- ggplot(ap_pie_data, aes(x = "", y = n, fill = AP_Category)) + geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar(theta = "y") + geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5) + scale_fill_manual(values = ap_colors) +
    theme_void(base_size = 14) + theme(legend.position = "bottom", legend.title = element_blank())
  
  # Fig. S1
  fig_S1 = ancestryPiePlot + antipsychPlot + sexPiePlot + sexByAgePlot + dxByAgePlot + plot_layout(nrow = 1)
  mpdf("Fig_S1", outDir=file.path(ROOT, "outputs"), width=12, height=8); print(fig_SX1); dev.off()
}

####################################################################################################
##### FIG. S2 :: QUALITY CONTROL FOR RNA-SEQ AND ATAC-SEQ DATA #####################################

{
  qcRna = read.csv(QC_RNASEQ, sep="\t")
  qcRna$RIN = ifelse(qcRna$RIN > 10, NA, qcRna$RIN) # Manual fix of the value 105 which was clearly a nonsense.
  qcRna$mergingDesigns = paste0(qcRna$cell_subtype, "_", qcRna$Dx)
  
  qcRnaSum = ddply(qcRna, "mergingDesigns", summarize,
                     `Cell type` = unique(cell_subtype),
                     `Diagnosis` = unique(Dx),
                     `Sample count` = length(ID),
                     `RIN` = mean(RIN),
                     `Aligned reads` = mean(finalReadCount),
                     `Fraction of duplicated reads` = mean(DuplicateReadFrac),
                     `Insert size` = mean(insertMetrics_MEDIAN_INSERT_SIZE),
                     `GC content in consensus peaks` = mean(picard_meanGcContent))
  
  qcRnaSum = qcRnaSum[,2:ncol(qcRnaSum)] %>% mutate_if(is.numeric, round, 3)
  mtsv(qcRnaSum, filename="Fig_S2_background_rnaseq", outDir=file.path(ROOT, "outputs"), myHeader=T)
  
  selectedCols = c("RIN"="RNA integrity number",
                   "finalReadCount"="Number of uniquely mapped reads",
                   "star_Uniquely_mapped_reads_pct"="Fraction of uniquely mapped reads",
                   "picard_PERCENT_DUPLICATION"="Fraction of duplicated reads", 
                   "rnaseqc_PCT_INTERGENIC_BASES"="Intergenic rate",
                   "rnaseqc_PCT_INTRONIC_BASES"="Intronic rate",
                   "picard_meanGcContent"="GC content in consensus peaks",
                   "insertMetrics_MEDIAN_INSERT_SIZE"="Median insert size")
  
  allPlots = list()
  for(name in names(selectedCols)) {
    myPlot = ggplot(qcRna, aes_string(x = "mergingDesigns", y = name)) +
      geom_boxplot(aes(fill = factor(mergingDesigns)), outlier.size=-1) + scale_fill_brewer(palette="Set3") + labs(x="", y="") +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1,
                         axis.text.y=element_text(colour = "black"), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none") + ylab(selectedCols[name])
    allPlots[[name]] = myPlot
  }
  
  plotQc2 = ggdraw() +
    draw_plot(allPlots[["RIN"]],   .00, .66, .33, .33) +
    draw_plot(allPlots[["finalReadCount"]],                   .33, .66, .33, .33) +
    draw_plot(allPlots[["star_Uniquely_mapped_reads_pct"]],   .66, .66, .33, .33) +
    draw_plot(allPlots[["picard_meanGcContent"]],       .00, .33, .33, .33) +
    draw_plot(allPlots[["picard_PERCENT_DUPLICATION"]],     .33, .33, .33, .33) +
    draw_plot(allPlots[["rnaseqc_PCT_INTERGENIC_BASES"]],       .66, .33, .33, .33) +
    draw_plot_label(c("d", "e", "f", "g", "h", "i"),
                    c(.00, .33, .66, .00, .33, .66),
                    c(.99, .99, .99, .66, .66, .66),
                    size = 15)
  
  mpdf("Fig_S2_d_e_f_g_h_i", outDir=file.path(ROOT, "outputs"), width=11, height=11); print(plotQc2); dev.off()

  #####
  # Fig. S2a :: Median read insert size distribution
  GABA = cbind.data.frame(unlist(sapply(unique(qcRna[qcRna$cell_subtype=="GABA","insertMetrics_MEDIAN_INSERT_SIZE"]), function(x) rep(x, sum(qcRna[qcRna$cell_subtype=="GABA","insertMetrics_MEDIAN_INSERT_SIZE"]==x)))), "GABA")
  colnames(GABA) = c("insertSize", "type")
  GLU = cbind.data.frame(unlist(sapply(unique(qcRna[qcRna$cell_subtype=="GLU","insertMetrics_MEDIAN_INSERT_SIZE"]), function(x) rep(x, sum(qcRna[qcRna$cell_subtype=="GLU","insertMetrics_MEDIAN_INSERT_SIZE"]==x)))), "GLU")
  colnames(GLU) = c("insertSize", "type")
  Olig = cbind.data.frame(unlist(sapply(unique(qcRna[qcRna$cell_subtype=="Olig","insertMetrics_MEDIAN_INSERT_SIZE"]), function(x) rep(x, sum(qcRna[qcRna$cell_subtype=="Olig","insertMetrics_MEDIAN_INSERT_SIZE"]==x)))), "Olig")
  colnames(Olig) = c("insertSize", "type")
  MgAs = cbind.data.frame(unlist(sapply(unique(qcRna[qcRna$cell_subtype=="MgAs","insertMetrics_MEDIAN_INSERT_SIZE"]), function(x) rep(x, sum(qcRna[qcRna$cell_subtype=="MgAs","insertMetrics_MEDIAN_INSERT_SIZE"]==x)))), "MgAs")
  colnames(MgAs) = c("insertSize", "type")
  histMedianInsertSizeDf = rbind.data.frame(GABA, GLU, Olig, MgAs)
  histMedianInsertSizeDf$type = ordered(histMedianInsertSizeDf$type, levels=c("GABA", "GLU", "Olig", "MgAs"))
  
  histMedianInsertSize = ggplot(histMedianInsertSizeDf, aes(insertSize, colour = type)) + geom_density(size=1) + coord_equal()  + 
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = c(0.5, 0.85), 
                       axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) +
    scale_colour_manual(labels=c("GABA neurons", "GLU neurons", "Olig", "Microglia & Astrocytes"), values=c(npgList$GABA, npgList$GLU, npgList$OLIG, npgList$MGAS))# + xlab("Median insert size [bp]") + ylab("Density")
  mpdf("Fig_S2_a", outDir=file.path(ROOT, "outputs")); print(histMedianInsertSize); dev.off()
  
  #####
  # Fig. S2b :: Sex check based on measuring the number reads mapped on chromosome Y
  chrY_genes = qcRna$qcPeakAnno[(qcRna$qcPeakAnno$seqnames=="chrY") & (qcRna$qcPeakAnno$PeakID %in% qcRna$initialDgeObj$genes$PeakID),]
  chrY_genes = chrY_genes[(chrY_genes$end < 10001 | chrY_genes$start > 2781479) & (chrY_genes$end < 155701383 | chrY_genes$start > 156030895),]
  
  #####  Fig. S2c :: Genotype check based on pair-wise comparison of genotypes called from RNA-seq samples with SNP-arrays
  kinshipRnaSnparray = read.csv(KINSHIP_RNASEQ_SNPPARRAY)
  z = kinshipRnaSnparray
  z$`Same person`= ordered(ifelse(z$samePerson, "yes", "no"), levels=c("yes", "no"))
  kinshipRnaSnparray = ggplot(z, aes(Kinship)) + scale_color_manual(name="Same person", labels=c("yes", "no"), values = c(npgList[["NEURON"]], npgList[["GLIA"]])) + 
    theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    geom_density(alpha=0, aes(color=`Same person`, fill=`Same person`), size=1) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = c(0.4, 0.85),
          axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + coord_equal() + xlab("Kinship score") + ylab("Density")
  kinshipRnaSnparray
  mpdf("Fig_S2_c", outDir=file.path(ROOT, "outputs")); print(kinshipRnaSnparray); dev.off()
}

####################################################################################################
##### FIG. S2 (ATAC-SEQ) / PART II :: QUALITY CONTROL ##############################################

{
  qcAtac = read.csv(QC_ATACSEQ, sep="\t")
  qcAtac$mergingDesigns = paste0(qcAtac$cell_subtype, "_", qcAtac$Dx)
  qcAtac$mergingDesigns = ordered(qcAtac$mergingDesigns, levels=c("GABAergic_Control", "GABAergic_SCZ", "glutamatergic_Control", "glutamatergic_SCZ", "microgliaAndAstrocytes_Control", "microgliaAndAstrocytes_SCZ", "oligodendrocytes_Control", "oligodendrocytes_SCZ"))
  
  qcAtacSum = ddply(qcAtac, "mergingDesigns", summarize,
                    `Cell type` = unique(cell_subtype),
                    `Diagnosis` = unique(Dx),
                    `Sample count` = length(ID),
                    `pH` = mean(na.omit(pH)),
                    `PMI [hours]` = mean(na.omit(PMI)),
                    `Age of death` = mean(na.omit(ageOfDeath)),
                    `Ratio of male samples` = sum(Sex == "Male") / length(Sex), 
                    `Ratio of Caucasian samples` = sum(Ethnicity == "Caucasian") / length(Race),
                    `CDR` = mean(na.omit(CDR)))
  
  qcAtacSum = qcAtacSum[,2:ncol(qcAtacSum)] %>% mutate_if(is.numeric, round, 3)
  mtsv(qcAtacSum, filename="Fig_S2_background_atacseq", outDir=file.path(ROOT, "outputs"), myHeader=T)
  
  selectedCols = c("star_Uniquely_mapped_reads_pct"="Fraction of uniquely mapped reads",
                   "finalReadCount"="Number of uniquely mapped reads",
                   "picard_PERCENT_DUPLICATION"="Fraction of duplicated reads", 
                   "chrMFrac"="Fraction of mitDNA reads",
                   "peakNarrowFDR1pctCount"="Number of narrow peaks",
                   "fracReadsInNonBlacklistedPeaks"="Fraction of reads in peaks (FRiP)",
                   "picard_meanGcContent"="GC content in consensus peaks",
                   "insertMetrics_MEDIAN_INSERT_SIZE"="Median insert size", 
                   "pbc"="PCR Bottleneck Coefficient")
  
  allPlots = list()
  for(name in names(selectedCols)) {
    myPlot = ggplot(qcAtac, aes_string(x = "mergingDesigns", y = name)) +
      geom_boxplot(aes(fill = factor(mergingDesigns)), outlier.size=-1) + scale_fill_brewer(palette="Set3") + labs(x="", y="") +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1,
                         axis.text.y=element_text(colour = "black"), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none") + ylab(selectedCols[name])
    allPlots[[name]] = myPlot
  }
  
  plotQc2 = ggdraw() +
    draw_plot(allPlots[["star_Uniquely_mapped_reads_pct"]],   .00, .66, .33, .33) +
    draw_plot(allPlots[["finalReadCount"]],                   .33, .66, .33, .33) +
    draw_plot(allPlots[["picard_PERCENT_DUPLICATION"]],       .66, .66, .33, .33) +
    draw_plot(allPlots[["peakNarrowFDR1pctCount"]],                         .00, .33, .33, .33) +
    draw_plot(allPlots[["fracReadsInNonBlacklistedPeaks"]],           .33, .33, .33, .33) +
    draw_plot(allPlots[["picard_meanGcContent"]],   .66, .33, .33, .33) +
    draw_plot(allPlots[["chrMFrac"]],             .00, .00, .33, .33) +
    draw_plot(allPlots[["pbc"]], .33, .00, .33, .33) +
    draw_plot_label(c("p", "q", "r", "s", "t", "u", "v", "w"),
                    c(.00, .33, .66, .00, .33, .66, .00, .33),
                    c(.99, .99, .99, .66, .66, .66, .33, .33),
                    size = 15)
  
  mpdf("Fig_S2_p_q_r_s_t_u_v_w", outDir=file.path(ROOT, "outputs"), width=11, height=11); print(plotQc2); dev.off()
  
  #####
  # Fig. S2j :: Median read insert size distribution
  gabaergic = cbind.data.frame(unlist(sapply(unique(qcAtac[qcAtac$cell_subtype=="GABAergic","insertMetrics_MEDIAN_INSERT_SIZE"]), function(x) rep(x, sum(qcAtac[qcAtac$cell_subtype=="GABAergic","insertMetrics_MEDIAN_INSERT_SIZE"]==x)))), "GABAergic")
  colnames(gabaergic) = c("insertSize", "type")
  glutamatergic = cbind.data.frame(unlist(sapply(unique(qcAtac[qcAtac$cell_subtype=="glutamatergic","insertMetrics_MEDIAN_INSERT_SIZE"]), function(x) rep(x, sum(qcAtac[qcAtac$cell_subtype=="glutamatergic","insertMetrics_MEDIAN_INSERT_SIZE"]==x)))), "glutamatergic")
  colnames(glutamatergic) = c("insertSize", "type")
  oligodendrocytes = cbind.data.frame(unlist(sapply(unique(qcAtac[qcAtac$cell_subtype=="oligodendrocytes","insertMetrics_MEDIAN_INSERT_SIZE"]), function(x) rep(x, sum(qcAtac[qcAtac$cell_subtype=="oligodendrocytes","insertMetrics_MEDIAN_INSERT_SIZE"]==x)))), "oligodendrocytes")
  colnames(oligodendrocytes) = c("insertSize", "type")
  microgliaAndAstrocytes = cbind.data.frame(unlist(sapply(unique(qcAtac[qcAtac$cell_subtype=="microgliaAndAstrocytes","insertMetrics_MEDIAN_INSERT_SIZE"]), function(x) rep(x, sum(qcAtac[qcAtac$cell_subtype=="microgliaAndAstrocytes","insertMetrics_MEDIAN_INSERT_SIZE"]==x)))), "microgliaAndAstrocytes")
  colnames(microgliaAndAstrocytes) = c("insertSize", "type")
  histMedianInsertSizeDf = rbind.data.frame(gabaergic, glutamatergic, oligodendrocytes, microgliaAndAstrocytes)
  histMedianInsertSizeDf$type = gsub("GABAergic", "GABA", gsub("glutamatergic", "GLU", gsub("oligodendrocytes", "OLIG", gsub("microgliaAndAstrocytes", "MGAS", histMedianInsertSizeDf$type))))
  histMedianInsertSizeDf$type = ordered(histMedianInsertSizeDf$type, levels=c("GABA", "GLU", "OLIG", "MGAS"))
  
  histMedianInsertSize = ggplot(histMedianInsertSizeDf, aes(insertSize, colour = type)) + geom_density(size=1) + coord_equal()  + 
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = c(0.5, 0.85), 
                       axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) +
    scale_colour_manual(labels=c("GABA", "GLU", "OLIG", "MGAS"), values=c(npgList$GABA, npgList$GLU, npgList$OLIG, npgList$MGAS))# + xlab("Median insert size [bp]") + ylab("Density")
  mpdf("Fig_S2_j", outDir=file.path(ROOT, "outputs")); print(histMedianInsertSize); dev.off()
  
  #####
  # Fig. S2k :: Distance of OCRs from the closest TSS
  maxDistance = 1E5
  breaksVector = seq(-20,20)*(1E5/20)
  gabaPeaks = read.csv(file.path(ROOT, "inputs", "peaks_GABA.bed"), sep="\t", header=F)
  gabaHistogram = hist(gabaPeaks$V14[abs(gabaPeaks$V14) <= maxDistance], breaks=breaksVector, plot=F)
  gabaHistogram$proportions = gabaHistogram$counts / sum(gabaHistogram$counts)
  gluPeaks = read.csv(file.path(ROOT, "inputs", "peaks_GLU.bed"), sep="\t", header=F)
  gluHistogram = hist(gluPeaks$V14[abs(gluPeaks$V14) <= maxDistance], breaks=breaksVector, plot=F)
  gluHistogram$proportions = gluHistogram$counts / sum(gluHistogram$counts)
  oligPeaks = read.csv(file.path(ROOT, "inputs", "peaks_OLIG.bed"), sep="\t", header=F)
  oligHistogram = hist(oligPeaks$V14[abs(oligPeaks$V14) <= maxDistance], breaks=breaksVector, plot=F)
  oligHistogram$proportions = oligHistogram$counts / sum(oligHistogram$counts)
  mgasPeaks = read.csv(file.path(ROOT, "inputs", "peaks_MGAS.bed"), sep="\t", header=F)
  mgasHistogram = hist(mgasPeaks$V14[abs(mgasPeaks$V14) <= maxDistance], breaks=breaksVector, plot=F)
  mgasHistogram$proportions = mgasHistogram$counts / sum(mgasHistogram$counts)
  
  histDf = data.frame(t(rbind(gabaHistogram$proportions, gluHistogram$proportions, oligHistogram$proportions, mgasHistogram$proportions)))
  colnames(histDf) = c("GABA", "GLU", "OLIG", "MGAS")
  histDf$breaks = gabaHistogram$mids
  
  histTssDist = ggplot() +
    geom_line(data=histDf, aes(x=breaks, y=GABA, color=npgList$GABA), linetype="solid", size=0.5) +
    geom_point(data=histDf, aes(x=breaks, y=GABA, color=npgList$GABA), size=1) +
    geom_line(data=histDf, aes(x=breaks, y=GLU, color=npgList$GLU), linetype="solid", size=0.5) +
    geom_point(data=histDf, aes(x=breaks, y=GLU, color=npgList$GLU), size=1) +
    geom_line(data=histDf, aes(x=breaks, y=OLIG, color=npgList$OLIG), linetype="solid", size=0.5) +
    geom_point(data=histDf, aes(x=breaks, y=OLIG, color=npgList$OLIG), size=1) +
    geom_line(data=histDf, aes(x=breaks, y=MGAS, color=npgList$MGAS), linetype="solid", size=0.5) +
    geom_point(data=histDf, aes(x=breaks, y=MGAS, color=npgList$MGAS), size=1) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = c(0.75, 0.85), 
                       axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + 
    xlab("Distance to TSS [bp]") + ylab("Proportion of OCRs") + coord_equal()  + 
    scale_colour_manual(labels=c("GABA", "GLU", "OLIG", "MGAS"), values=c(npgList$GABA, npgList$GLU, npgList$OLIG, npgList$MGAS))# + xlab("Median insert size [bp]") + ylab("Density")
  mpdf("Fig_S2_k", outDir=file.path(ROOT, "outputs")); print(histTssDist); dev.off()
  
  #####
  # Fig. S2n :: Sex check based on measuring the number reads mapped on OCRs located at chromosome Y (pseudoautosomal regions not counted)
  chrYplot = ggplot(qcAtac, aes(fracReadsInNonBlacklistedPeaks, chryCounts, color=Gender)) + 
    geom_point() + scale_color_manual(name="Sex", labels=c("Female", "Male "), values = c(npgList[["NEURON"]], npgList[["GLIA"]])) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = c(0.15, 0.85),
                       axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + 
    xlab("Fraction of OCRs in peaks of open chromatin") + ylab("chrY read count") + coord_equal() 
  mpdf("Fig_S2_n", outDir=file.path(ROOT, "outputs")); print(chrYplot); dev.off()
  
  #####
  # Fig. S2o :: Genotype check based on pair-wise comparison of genotypes called from ATAC-seq samples with SNP-arrays
  kinshipAtacSnparray = read.csv(KINSHIP_ATACSEQ_SNPPARRAY)
  z = kinshipAtacSnparray
  z$ID1x = sapply(z$ID2, function(x) {
    xx = strsplit(x, "_")[[1]]
    paste0(xx[2:length(xx)], collapse="_")
  })
  z = z[(z$ID1x %in% qcAtac$ID),]
  z = z[!(is.na(z$samePerson)),]
  z$`Same person`= ordered(ifelse(z$samePerson, "yes", "no"), levels=c("yes", "no"))
  kinshipAtacSnparray = ggplot(z, aes(Kinship)) + scale_color_manual(name="Same person", labels=c("yes", "no"), values = c(npgList[["NEURON"]], npgList[["GLIA"]])) + 
    theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    geom_density(alpha=0, aes(color=`Same person`, fill=`Same person`), size=1) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = c(0.4, 0.85),
          axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + coord_equal() + xlab("Kinship score") + ylab("Density")
  mpdf("Fig_S2_o", outDir=file.path(ROOT, "outputs")); print(kinshipAtacSnparray); dev.off()
}

####################################################################################################
##### FIG. S4 (ATAC-SEQ) :: COMPARISON WITH HAUBERG ET AL. 2023 AND COLEMAN ET AL. 2023 ############

{
  # Fig. S4a :: Correlation of log2(cpm+1) counts between our ATAC-seq data and external ATAC-seq data from 4 cell types from the prefrontal cortex (Hauberg et al 2020)
  atacseq_countMatrixRaw = readRDS(ATACSEQ_COUNT_MATRIX_RAW)
  atacseq_countMatrixAdj = readRDS(ATACSEQ_COUNT_MATRIX_ADJ)
  rnaseq_countMatrixRaw = readRDS(RNASEQ_COUNT_MATRIX_RAW)
  rnaseq_countMatrixAdj = readRDS(RNASEQ_COUNT_MATRIX_ADJ)
  qcPeakAnno = readRDS(ATACSEQ_PEAKS)
  
  ############################
  # ATAC-seq :: CPM comparison with Hauberg et al 2020
  ggom_allInfo = read.csv(METADATA_HAUBERG_2020)
  ggom_countMatrix = readRDS(GEXPR_HAUBERG_2020)
  ggom_qcPeakAnno = readRDS(PEAKS_HAUBERG_2020)
  ggom_qcPeakAnno$PeakID = qcPeakAnno$PeakID
  rownames(ggom_qcPeakAnno) = qcPeakAnno$PeakID
  
  # Load metadata from Hauberg et al 2020
  rownames(ggom_allInfo) = ggom_allInfo$ID
  ggom_allInfo = ggom_allInfo[order(ggom_allInfo$ID),]
  
  # Load count matrix from Hauberg et al 2020
  rownames(ggom_qcPeakAnno) = qcPeakAnno$PeakID
  ggom_countMatrix = ggom_countMatrix[rownames(atacseq_countMatrixRaw),]
  outDir = file.path(ROOT, "tmp")
  dir.create(outDir)
  geneNormObj = getGeneFilteredGeneExprMatrix(ggom_countMatrix, ggom_allInfo, ggom_qcPeakAnno[rownames(ggom_qcPeakAnno) %in% rownames(ggom_countMatrix),], plotName="PRE_COVS", geneTssPeakMapping=NULL,housekeepingPeakInfo=NULL, MIN_GENE_CPM=0, MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0, calcNormFactors.method="TMM")
  colnames(geneNormObj$dgeObj) = gsub("^X", "", colnames(geneNormObj$dgeObj))
  
  initialDgeObj = geneNormObj$dgeObj
  geneNormObj$dgeObj = NULL
  initialVoomObj = voom(initialDgeObj, design=NULL, plot=F)
  
  currentDf = do.call("cbind.data.frame", lapply(tolower(unique(qcAtac$cell_subtype)), function(ctype) {
    rowMeans(log2(2^atacseq_countMatrixAdj[,qcAtac[tolower(qcAtac$cell_subtype) == ctype,"ID"]]+1))
  }))
  colnames(currentDf) = paste0("current_", tolower(unique(qcAtac$cell_subtype)))
  ggomDf = do.call("cbind.data.frame", lapply(tolower(unique(ggom_allInfo$cell_subtype)), function(ctype) {
    rowMeans(log2(2^initialVoomObj$E[,ggom_allInfo[tolower(ggom_allInfo$cell_subtype)==ctype,"ID"]]+1))
  }))
  colnames(ggomDf) = paste0("ggom_", tolower(unique(qcAtac$cell_subtype)))
  
  for(ctype in tolower(unique(qcAtac$cell_subtype))) {
    for(ctype2 in tolower(unique(qcAtac$cell_subtype))) {
      df = cbind.data.frame(currentDf, ggomDf)
      df$density = get_density(df[,paste0("current_", ctype)], df[,paste0("ggom_", ctype2)], n = 100)
      axisMax = round(max(min(abs(df[,paste0("current_", ctype)])),max(abs(df[,paste0("ggom_", ctype2)]))+0.5))
      pearson = cor.test(df[,paste0("current_", ctype)], df[,paste0("ggom_", ctype2)], method="pearson")
      spearman = cor.test(df[,paste0("current_", ctype)], df[,paste0("ggom_", ctype2)], method="spearman")
      print((paste0("> ", ctype, " / ", ctype2, " :: Pearson / Spearman = ", round(pearson$estimate, 3), " / ", round(spearman$estimate, 3))))
      
      densityScatter_current_ggom = ggplot(df, aes_string(x=paste0("current_", ctype), y=paste0("ggom_", ctype2))) + geom_point(aes_string(x=paste0("current_", ctype), y=paste0("ggom_", ctype2), color="density")) + scale_color_viridis() +
        coord_equal() + theme_classic() + theme(axis.text.y=element_text(colour="black")) + 
        xlab(paste0("log2(cpm+1); this study - ", ctype)) + ylab(paste0("log2(cpm+1); GGOM - ", ctype2)) + xlim(c(0,axisMax)) + ylim(c(0,axisMax)) +
        geom_abline(intercept=0, slope=1, color="gray", linetype="dashed") + geom_hline(yintercept=0, color="gray", linetype="dashed") + geom_smooth(method=lm, se=FALSE) +
        ggtitle(paste0("Pearson / Spearman = ", round(pearson$estimate, 3), " / ", round(spearman$estimate, 3)))
      mpdf(paste0("Fig_S4_hauberg__", ctype, "__", ctype2), outDir=file.path(ROOT, "outputs")); print(densityScatter_current_ggom); dev.off()
    }
  }  
  
  #####
  # Fig. S4b :: Correlation of log2(cpm+1) counts between our RNA-seq data and external RNA-seq dataset of 3 cell types from parahippocampal gyrus (Coleman et al 2023)
  phg_allInfo = read.csv(METADATA_COLEMAN_2023)
  phg_initDgeObj = readRDS("/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/tmp/phg_rnaseqInitialDgeObj.RDS")
  phg_initVoomObj = readRDS(GEXPR_COLEMAN_2023)
  
  currentDf = do.call("cbind.data.frame", lapply(tolower(unique(qcRna$cell_subtype)), function(ctype) {
    rowMeans(log2(2^rnaseq_countMatrixAdj[,qcRna[tolower(qcRna$cell_subtype) == ctype,"ID"]]+1))
  }))
  colnames(currentDf) = paste0("current_", tolower(unique(qcRna$cell_subtype)))
  phgDf = do.call("cbind.data.frame", lapply(tolower(unique(phg_allInfo$cell_subtype)), function(ctype) {
    rowMeans(log2(2^phg_initVoomObj$E[,phg_allInfo[tolower(phg_allInfo$cell_subtype)==ctype,"ID"]]+1))
  }))
  colnames(phgDf) = paste0("phg_", tolower(unique(unique(phg_allInfo$cell_subtype))))
  
  isectGenes = intersect(rownames(currentDf), rownames(phgDf))
  dfx = cbind.data.frame(currentDf[isectGenes,], phgDf[isectGenes,])
  for(ctype in tolower(unique(qcRna$cell_subtype))) {
    for(ctype2 in tolower(unique(qcRna$cell_subtype))) {
      ctype2 = ifelse(ctype2 == "gaba", "neuron", ifelse(ctype2 == "glu", "neuron", ifelse(ctype2 == "olig", "oligodendrocytes", ifelse(ctype2 == "mgas", "astroandmicroglia", NA))))
      df = dfx
      df$density = get_density(df[,paste0("current_", ctype)], df[,paste0("phg_", ctype2)], n = 100)
      axisMax = round(max(min(abs(df[,paste0("current_", ctype)])),max(abs(df[,paste0("phg_", ctype2)]))+0.5))
      pearson = cor.test(df[,paste0("current_", ctype)], df[,paste0("phg_", ctype2)], method="pearson")
      spearman = cor.test(df[,paste0("current_", ctype)], df[,paste0("phg_", ctype2)], method="spearman")
      print((paste0("> ", ctype, " / ", ctype2, " :: Pearson / Spearman = ", round(pearson$estimate, 3), " / ", round(spearman$estimate, 3))))
      
      densityScatter_current_phg = ggplot(df, aes_string(x=paste0("current_", ctype), y=paste0("phg_", ctype2))) + geom_point(aes_string(x=paste0("current_", ctype), y=paste0("phg_", ctype2), color="density")) + scale_color_viridis() +
        coord_equal() + theme_classic() + theme(axis.text.y=element_text(colour="black")) + 
        xlab(paste0("log2(cpm+1); this study - ", ctype2)) + ylab(paste0("log2(cpm+1); PHG - ", ctype2)) + xlim(c(0,axisMax)) + ylim(c(0,axisMax)) +
        geom_abline(intercept=0, slope=1, color="gray", linetype="dashed") + geom_hline(yintercept=0, color="gray", linetype="dashed") + geom_smooth(method=lm, se=FALSE) +
        ggtitle(paste0("Pearson / Spearman = ", round(pearson$estimate, 3), " / ", round(spearman$estimate, 3)))
      mpdf(paste0("Fig_S4_coleman__", ctype, "__", ctype2), outDir=file.path(ROOT, "outputs")); print(densityScatter_current_phg); dev.off()
    }
  }
}

####################################################################################################
##### FIG. S5 :: COMPARISON OF SIGNAL FOR MARKER GENES IN RNA-SEQ AND PROMOTER OCRS ################

{
  # Define marker genes and corresponding OCRs
  CELL_TYPES = c("GABA", "GLU", "OLIG", "MGAS")
  CELL_MARKER_GENES = list(
    "GABA" = c("PVALB", "SST", "SLC32A1", "DLX6"), 
    "GLU" = c("NEUROD6", "BDNF", "TBR1", "SLC17A6"),
    "OLIG" = c("OLIG1", "OPALIN", "MBP", "PDGFRA"),
    "MGAS" = c("TYROBP", "TREM2", "CX3CR1", "ALDH1L1")
  )
  CELL_MARKER_PROMOTERS = list(
    "GABA" = c("Peak_346651", "Peak_400760", "Peak_330236", "Peak_608799"), 
    "GLU" = c("Peak_590189", "Peak_86350", "Peak_301041", "Peak_83899"), 
    "OLIG" = c("Peak_342301", "Peak_72780", "Peak_259883", "Peak_419664"), 
    "MGAS" = c("Peak_262781", "Peak_536771", "Peak_360255", "Peak_382793")
  )
  
  atacseq_countMatrixResi = readRDS(ATACSEQ_COUNT_MATRIX_RESIDUALIZED_DX_CELLTYPE_KEPT)
  atacseq_expr = log2(2^atacseq_countMatrixResi + 1)
  
  rnaseq_countMatrixResi = readRDS(RNASEQ_COUNT_MATRIX_RESIDUALIZED_DX_CELLTYPE_KEPT)
  rnaseq_expr = log2(2^rnaseq_countMatrixResi + 1)
  
  ###
  # Plot gene expression markers
  qcRna$cell_subtype = toupper(qcRna$cell_subtype)
  
  CELL_MARKER_GENES_ENSEMBL = lapply(CELL_MARKER_GENES, symbol_to_ensembl, anno_df = gtf) # Convert marker gene names to Ensembl gene identifiers
  
  genes_to_plot = unique(unlist(CELL_MARKER_GENES_ENSEMBL)) |> intersect(rownames(rnaseq_expr))
  samples_all = qcRna %>% filter(cell_subtype %in% CELL_TYPES) %>% pull(ID)
  expr_sub = rnaseq_expr[genes_to_plot, samples_all, drop = FALSE]
  
  plot_df = as.data.frame(expr_sub) %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
    left_join(qcRna[, c("ID","Dx","cell_subtype")], by = c("sample" = "ID")) %>%
    mutate(
      cell_subtype = factor(toupper(cell_subtype), levels = CELL_TYPES),
      Dx = factor(Dx, levels = c("Control","SCZ")),
      group = factor(paste0(cell_subtype, "_", Dx),
                     levels = c("GABA_Control","GABA_SCZ",
                                "GLU_Control","GLU_SCZ",
                                "OLIG_Control","OLIG_SCZ",
                                "MGAS_Control","MGAS_SCZ")),
      gene_symbol = unlist(CELL_MARKER_GENES)[match(gene, unlist(CELL_MARKER_GENES_ENSEMBL))],
      gene_symbol = ordered(gene_symbol, levels = unlist(CELL_MARKER_GENES)))
  
  final_plot_facet = ggplot(plot_df, aes(x = group, y = expr, fill = cell_subtype)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", size = 0.3) +
    scale_fill_manual(values = c("GABA" = npgList$GABA, "GLU"  = npgList$GLU, "OLIG" = npgList$OLIG, "MGAS" = npgList$MGAS)) +
    facet_wrap(~ gene_symbol, ncol = 4, scales = "free_y") +
    theme_minimal(base_size = 9) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 9),
      axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 8),
      legend.position = "none",
      strip.text = element_text(size = 9, face = "bold")
    ) +
    labs(y = "log2(expr+1)")
  
  mpdf(paste0("Fig_S5a_rnaseq"), outDir=file.path(ROOT, "outputs")); print(final_plot_facet); dev.off()
  
  ###
  # Plot chromatin accessibility (peak) markers
  qcAtac$cell_subtype = gsub("GABAergic", "GABA", gsub("glutamatergic", "GLU", gsub("oligodendrocytes", "OLIG", gsub("microgliaAndAstrocytes", "MGAS", qcAtac$cell_subtype))))
  
  peaks_to_plot = unique(unlist(CELL_MARKER_PROMOTERS)) |> intersect(rownames(atacseq_expr))
  samples_all = qcAtac %>% filter(cell_subtype %in% CELL_TYPES) %>% pull(ID)
  expr_sub = atacseq_expr[unlist(CELL_MARKER_PROMOTERS), samples_all, drop = FALSE]
  
  plot_df = as.data.frame(expr_sub) %>%
    tibble::rownames_to_column("peak") %>%
    pivot_longer(-peak, names_to = "sample", values_to = "expr") %>%
    left_join(qcAtac[, c("ID","Dx","cell_subtype")], by = c("sample" = "ID")) %>%
    mutate(
      cell_subtype = factor(toupper(cell_subtype), levels = CELL_TYPES),
      Dx = factor(Dx, levels = c("Control","SCZ")),
      group = factor(paste0(cell_subtype, "_", Dx),
                     levels = c("GABA_Control","GABA_SCZ",
                                "GLU_Control","GLU_SCZ",
                                "OLIG_Control","OLIG_SCZ",
                                "MGAS_Control","MGAS_SCZ")),
      gene_symbol = unlist(CELL_MARKER_GENES)[match(peak, unlist(CELL_MARKER_PROMOTERS))],
      gene_symbol = ordered(gene_symbol, levels = unlist(CELL_MARKER_GENES)))
  
  final_plot_facet = ggplot(plot_df, aes(x = group, y = expr, fill = cell_subtype)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", size = 0.3) +
    scale_fill_manual(values = c("GABA" = npgList$GABA, "GLU"  = npgList$GLU, "OLIG" = npgList$OLIG, "MGAS" = npgList$MGAS)) +
    facet_wrap(~ gene_symbol, ncol = 4, scales = "free_y") +
    theme_minimal(base_size = 9) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 9),
      axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 8),
      legend.position = "none",
      strip.text = element_text(size = 9, face = "bold")
    ) +
    labs(y = "log2(expr+1)")
  
  mpdf(paste0("Fig_S5b_atacseq"), outDir=file.path(ROOT, "outputs")); print(final_plot_facet); dev.off()
}

####################################################################################################
##### FIG. S6 :: ESTIMATED CELL TYPE COMPOSITION FOR SCZ CASE AND CONTROL SAMPLES ##################

{
  #####
  # RNA-seq
  plot_df = qcRna %>%
    dplyr::select(ID, Dx, cell_subtype, deconvolution_GABA, deconvolution_GLU, deconvolution_ODC, deconvolution_AST, deconvolution_MG) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("deconvolution_"),
      names_to = "deconv_target",
      names_prefix = "deconvolution_",
      values_to = "estimate"
    ) %>%
    dplyr::mutate(
      cell_subtype  = case_when(
        cell_subtype == "GABA"  ~ "GABA",
        cell_subtype == "GLU"  ~ "GLU",
        cell_subtype == "ODC"  ~ "OLIG",
        cell_subtype == "OLIG" ~ "OLIG",
        cell_subtype == "MGAS" ~ "MGAS",
        TRUE ~ cell_subtype
      ),
      cell_subtype  = factor(cell_subtype, levels = CELL_TYPES),
      deconv_target = factor(toupper(deconv_target),
                             levels = c("GABA", "GLU", "ODC", "AST", "MG", "OLIG", "MGAS")),
      Dx = factor(Dx, levels = c("Control","SCZ")),
      group = factor(paste0(deconv_target, "_", Dx),
                     levels = c("GABA_Control", "GABA_SCZ", "GLU_Control", "GLU_SCZ", "ODC_Control", "ODC_SCZ", "AST_Control", "AST_SCZ", "MG_Control", "MG_SCZ"))
    )
  
  make_violin_deconv <- function(df, title = NULL) {
    ggplot(df, aes(x = group, y = estimate, fill = deconv_target)) + # geom_violin(trim = TRUE, scale = "width", color = "gray40", size = 0.3) +
      geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", size = 0.3) +
      scale_fill_manual(values = npgList) + coord_cartesian(ylim = c(0, NA)) +
      theme_minimal(base_size = 9) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
      ) +
      labs(y = "Deconvolution estimate", title = title)
  }
  
  plots_rna = list(
    make_violin_deconv(filter(plot_df, cell_subtype == "GABA"), "Sorted nuclei: GABA"),
    make_violin_deconv(filter(plot_df, cell_subtype == "GLU"),  "Sorted nuclei: GLU"),
    make_violin_deconv(filter(plot_df, cell_subtype == "OLIG"), "Sorted nuclei: OLIG"),
    make_violin_deconv(filter(plot_df, cell_subtype == "MGAS"), "Sorted nuclei: MGAS")
  )
  
  final_plot = wrap_plots(plots_rna, ncol = 4)
  mpdf(paste0("Fig_S6a_rnacseq"), outDir=file.path(ROOT, "outputs"), width = 9.5, height = 3); print(final_plot); dev.off()
  
  #####
  # ATAC-seq
  plot_df = qcAtac %>%   # Prepare data for plotting
    dplyr::select(ID, Dx, cell_subtype, dplyr::starts_with("deconvolution_")) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("deconvolution_"), names_to = "deconv_target", names_prefix = "deconvolution_", values_to = "estimate") %>%
    dplyr::mutate(
      cell_subtype = factor(toupper(cell_subtype), levels = CELL_TYPES),
      deconv_target = factor(toupper(deconv_target), levels = CELL_TYPES),
      Dx = factor(Dx, levels = c("Control","SCZ")),
      group = factor(paste0(deconv_target, "_", Dx), levels = c("GABA_Control", "GABA_SCZ", "GLU_Control", "GLU_SCZ", "OLIG_Control", "OLIG_SCZ", "MGAS_Control", "MGAS_SCZ")))
  
  # Violin + boxplot function
  make_violin_deconv = function(df, title = NULL) {
    ggplot(df, aes(x = group, y = estimate, fill = deconv_target)) +
      #geom_violin(trim = TRUE, scale = "width", color = "gray40", size = 0.3) +
      geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", size = 0.3) +
      scale_fill_manual(values = npgList) +
      theme_minimal(base_size = 9) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
      ) +
      labs(y = "Deconvolution estimate", title = title)
  }
  
  # Make 4 panels — one per sorted subtype
  plots_atac = lapply(CELL_TYPES, function(ct) {
    make_violin_deconv(filter(plot_df, cell_subtype == ct), paste0("Sorted nuclei: ", ct))
  })
  
  # Combine into one grid
  final_plot = wrap_plots(plots_atac, ncol = 4) + plot_annotation(title = "ATAC-seq: Deconvolution-based cell type composition per sorted population")
  mpdf(paste0("Fig_S6b_atacseq"), width = 9.5, height = 3, outDir=file.path(ROOT, "outputs")); print(final_plot); dev.off()
}

{ # Deconvolution plot :: RNA-seq
  library(tidyverse)
  library(patchwork)
  
  # Define canonical order and palette
  celltypes <- c("GABA","GLU","OLIG","MGAS")
  npgList <- list(
    "GABA"="#66A61E",
    "GLU" ="#E6AB02",
    "OLIG"="#E7298A",
    "AST" ="#7570B3",
    "MG"  ="#C1BEE2" 
  )
  
  # Standardize subtype labels
  rnaseqAll$allInfo$cell_subtype = toupper(rnaseqAll$allInfo$cell_subtype)
  rnaseqAll$allInfo$deconvolution_best <- NULL
  rnaseqAll$allInfo = allInfo[,!colnames(allInfo) %in% c("deconvolution_Sample")]
  
  # Prepare data (keep all five deconv targets)
  plot_df <- rnaseqAll$allInfo %>%
    dplyr::select(ID, Dx, cell_subtype, dplyr::starts_with("deconvolution_")) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("deconvolution_"),
      names_to = "deconv_target",
      names_prefix = "deconvolution_",
      values_to = "estimate"
    ) %>%
    dplyr::mutate(
      cell_subtype  = case_when(
        cell_subtype == "ODC"  ~ "OLIG",
        cell_subtype == "OLIG" ~ "OLIG",
        cell_subtype == "MGAS" ~ "MGAS",
        TRUE ~ cell_subtype
      ),
      cell_subtype  = factor(cell_subtype, levels = celltypes),
      deconv_target = factor(toupper(deconv_target),
                             levels = c("GABA","GLU","ODC","AST","MG")),
      Dx = factor(Dx, levels = c("Control","SCZ")),
      group = factor(paste0(deconv_target, "_", Dx),
                     levels = c("GABA_Control","GABA_SCZ",
                                "GLU_Control","GLU_SCZ",
                                "ODC_Control","ODC_SCZ",
                                "AST_Control","AST_SCZ",
                                "MG_Control","MG_SCZ"))
    )
  
  # ---- Plot function ----
  make_violin_deconv <- function(df, title = NULL) {
    ggplot(df, aes(x = group, y = estimate, fill = deconv_target)) +
      geom_violin(trim = TRUE, scale = "width", color = "gray40", size = 0.3) +
      geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", size = 0.3) +
      scale_fill_manual(values = npgList) +
      theme_minimal(base_size = 9) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
      ) +
      labs(y = "Deconvolution estimate", title = title)
  }
  
  # ---- Create panels ----
  plots_rna <- list(
    make_violin_deconv(filter(plot_df, cell_subtype == "GABA"), "Sorted nuclei: GABA"),
    make_violin_deconv(filter(plot_df, cell_subtype == "GLU"),  "Sorted nuclei: GLU"),
    make_violin_deconv(filter(plot_df, cell_subtype == "OLIG"), "Sorted nuclei: OLIG"),
    make_violin_deconv(filter(plot_df, cell_subtype == "MGAS" &
                                deconv_target %in% c("ODC","AST","MG")),
                       "Sorted nuclei: MGAS")
  )
  
  # ---- Combine grid ----
  final_plot <- wrap_plots(plots_rna, ncol = 2) +
    plot_annotation(title = "RNA-seq: Deconvolution-based cell type estimates")
  
  pdf("~/Desktop/deconv_RNA_4panels_mergedMGAS.pdf", width = 8, height = 6); print(final_plot); dev.off()
  
  final_plot
}


####################################################################################################
##### FIG. S9 :: CONCORDANCE BETWEEN DIFF GENE EXPRESSION RESULTS FROM THIS STUDY AND PSYCHAD ######

{
  # Load the results of our and PsychAD differential analysis
  degAnalysis = new.env(); load(DEG_ANALYSIS, envir=degAnalysis)
  psychadDf = readRDS(DEG_ANALYSIS_PSYCHAD_c07x)
  
  # Set the corresponding cell types between PsychAD and our study. Note that our OLIG and MGAS correspond to multiple PsychAD cell populations
  PAIRS = list(
    "GLU.SCZ_Control" = c("EN"),
    "GABA.SCZ_Control" = c("IN"),
    "Olig.SCZ_Control" = c("Oligo", "OPC"),
    "MgAs.SCZ_Control" = c("Astro", "Immune")   # Theoretically also "Mural" and "Endo" but we ignore them as those are very small cell populations
  )
  
  # PsychAD gene identifiers is a mix of Ensembl gene names and identifiers -> let's convert it to identifiers
  psychadDf$Ensembl_IDx = gtf[match(psychadDf$ID, gtf$gene_name), "PeakID"]
  psychadDf$Ensembl_ID = ifelse(is.na(psychadDf$Ensembl_IDx), psychadDf$ID, psychadDf$Ensembl_IDx)
  
  # Iterate over all possible pairs (our t-stats -vs- PsychAD t-stats for all combinations of cell populations) and save their comparison
  DENSITY_PLOTS = list()
  for(pairName in names(PAIRS)) {
    pair = PAIRS[[pairName]]
    for(psychAD_ctype in pair) {
      ours = degAnalysis$dacResults$dac[[pairName]][degAnalysis$dacResults$dac[[pairName]]$adj.P.Val < 0.05,]  # Only genes that are FDR significant in our study are being tested
      psychad_subsetDf = psychadDf[(psychadDf$assay == psychAD_ctype),]
      psychad_subsetDf = psychad_subsetDf[!duplicated(psychad_subsetDf$Ensembl_ID),]
      rownames(psychad_subsetDf) = psychad_subsetDf$Ensembl_ID
      isect = intersect(ours$PeakID, psychad_subsetDf$Ensembl_ID)
      colnames(psychad_subsetDf) = paste0("psychad_", colnames(psychad_subsetDf))
      df = cbind.data.frame(ours[isect,c("t", "AveExpr", "logFC")], psychad_subsetDf[isect,])
      if(nrow(df) < 10)
        next
      
      df$density = get_density(df$t, df[,"psychad_t"], n = 100)
      axisMax = round(max((abs(df$t)),max(abs(df[,"psychad_t"]))+0.5))
      
      axisMaxX = round(max(abs(df$t)) + 0.5)
      axisMaxY = round(max(abs(df[,"psychad_t"])) + 0.5)
      densityScatter = ggplot(df, aes_string(x="t", y="psychad_t")) + geom_point(aes_string(x="t", y="psychad_t", color="density")) + scale_color_viridis() +
        coord_cartesian(xlim = c(-axisMaxX, axisMaxX), ylim = c(-axisMaxY, axisMaxY)) + theme_classic() + theme(aspect.ratio = 1, axis.text.y=element_text(colour="black")) + 
        xlab(paste0("t-stats; this study - ", pairName)) + ylab(paste0("t-stats; PsychAD ", psychAD_ctype))  +
        geom_abline(intercept=0, slope=1, color="gray", linetype="dashed") + geom_hline(yintercept=0, color="gray", linetype="dashed") + geom_smooth(method=lm, se=FALSE) +
        ggtitle(paste0("R=", round(cor.test(df$t, df[,"psychad_t"])$estimate, 3)))
      print(densityScatter)
      
      DENSITY_PLOTS[[paste0("densityScatter_PsychAD_", pairName, "__", psychAD_ctype)]] = densityScatter
    }
  }
  
  combined_plot = wrap_plots(DENSITY_PLOTS, ncol = 3)
  mpdf(paste0("Fig_S9"), outDir=file.path(ROOT, "outputs"), width=9, height=6); print(combined_plot); dev.off()
}

####################################################################################################
##### FIG. 2A :: NUMBERS OF DIFFERENTIALLY ACCESSIBLE PEAKS ########################################

{
  dacAnalysis = new.env(); load(DAC_ANALYSIS, envir=dacAnalysis)
  
  atacMerged = list(
    "GABAergic neurons" = dacAnalysis$dacResults$dac$GABAergic.SCZ_Control[dacAnalysis$dacResults$dac$GABAergic.SCZ_Control$adj.P.Val<0.05,],
    "Glutamatergic neurons" = dacAnalysis$dacResults$dac$glutamatergic.SCZ_Control[dacAnalysis$dacResults$dac$glutamatergic.SCZ_Control$adj.P.Val<0.05,],
    "Oligodendrocytes" = dacAnalysis$dacResults$dac$oligodendrocytes.SCZ_Control[dacAnalysis$dacResults$dac$oligodendrocytes.SCZ_Control$adj.P.Val<0.05,],
    "Microglia+Astrocytes" = dacAnalysis$dacResults$dac$microgliaAndAstrocytes.SCZ_Control[dacAnalysis$dacResults$dac$microgliaAndAstrocytes.SCZ_Control$adj.P.Val<0.05,],
    "Merged" = dacAnalysis$dacResults$dac$SCZ_Control[dacAnalysis$dacResults$dac$SCZ_Control$adj.P.Val<0.05,]
  )
  atacMerged = shrinkAtacToptables(atacMerged)
  z=names(atacMerged)
  atacMerged_metadata = data.frame(Set=z, overallAssay="ATAC", cell="Merged", assay="ATAC", contrast=gsub("(atacGABA__|_All$|_DLPFC|_ACC$)","",z), brainRegion=gsub(".+_","",z), stringsAsFactors=F)
  atacMerged_metadata$SetFullName = with(atacMerged_metadata, paste(overallAssay, cell, contrast, brainRegion))
  atacMerged_metadata = addDacExclStats(atacMerged_metadata, atacMerged, nrow(dacAnalysis$modeledVoomObj$E))
  rm(z)
  
  atacUp = list(
    "GABAergic neurons" = dacAnalysis$dacResults$dacUp$GABAergic.SCZ_Control[dacAnalysis$dacResults$dacUp$GABAergic.SCZ_Control$adj.P.Val<0.05,],
    "Glutamatergic neurons" = dacAnalysis$dacResults$dacUp$glutamatergic.SCZ_Control[dacAnalysis$dacResults$dacUp$glutamatergic.SCZ_Control$adj.P.Val<0.05,],
    "Oligodendrocytes" = dacAnalysis$dacResults$dacUp$oligodendrocytes.SCZ_Control[dacAnalysis$dacResults$dacUp$oligodendrocytes.SCZ_Control$adj.P.Val<0.05,],
    "Microglia+Astrocytes" = dacAnalysis$dacResults$dacUp$microgliaAndAstrocytes.SCZ_Control[dacAnalysis$dacResults$dacUp$microgliaAndAstrocytes.SCZ_Control$adj.P.Val<0.05,],
    "Merged" = dacAnalysis$dacResults$dacUp$SCZ_Control[dacAnalysis$dacResults$dacUp$SCZ_Control$adj.P.Val<0.05,]
  )
  atacUp = shrinkAtacToptables(atacUp)
  z=names(atacUp)
  atacUp_metadata = data.frame(Set=z, overallAssay="ATAC", cell="Up", assay="ATAC", contrast=gsub("(atacGABA__|_All$|_DLPFC|_ACC$)","",z), brainRegion=gsub(".+_","",z), stringsAsFactors=F)
  atacUp_metadata$SetFullName = with(atacUp_metadata, paste(overallAssay, cell, contrast, brainRegion))
  atacUp_metadata = addDacExclStats(atacUp_metadata, atacUp, nrow(dacAnalysis$modeledVoomObj$E))
  rm(z)
  
  atacDown = list(
    "GABAergic neurons" = dacAnalysis$dacResults$dacUp$GABAergic.Control_SCZ[dacAnalysis$dacResults$dacUp$GABAergic.Control_SCZ$adj.P.Val<0.05,],
    "Glutamatergic neurons" = dacAnalysis$dacResults$dacUp$glutamatergic.Control_SCZ[dacAnalysis$dacResults$dacUp$glutamatergic.Control_SCZ$adj.P.Val<0.05,],
    "Oligodendrocytes" = dacAnalysis$dacResults$dacUp$oligodendrocytes.Control_SCZ[dacAnalysis$dacResults$dacUp$oligodendrocytes.Control_SCZ$adj.P.Val<0.05,],
    "Microglia+Astrocytes" = dacAnalysis$dacResults$dacUp$microgliaAndAstrocytes.Control_SCZ[dacAnalysis$dacResults$dacUp$microgliaAndAstrocytes.Control_SCZ$adj.P.Val<0.05,],
    "Merged" = dacAnalysis$dacResults$dacUp$Control_SCZ[dacAnalysis$dacResults$dacUp$Control_SCZ$adj.P.Val<0.05,]
  )
  atacDown = shrinkAtacToptables(atacDown)
  z=names(atacDown)
  atacDown_metadata = data.frame(Set=z, overallAssay="ATAC", cell="Down", assay="ATAC", contrast=gsub("(atacGABA__|_All$|_DLPFC|_ACC$)","",z), brainRegion=gsub(".+_","",z), stringsAsFactors=F)
  atacDown_metadata$SetFullName = with(atacDown_metadata, paste(overallAssay, cell, contrast, brainRegion))
  atacDown_metadata = addDacExclStats(atacDown_metadata, atacDown, nrow(dacAnalysis$modeledVoomObj$E))
  rm(z)
  
  z = rbind.fill(atacUp_metadata, atacDown_metadata, atacMerged_metadata)
  #z = z[z$contrast != "Glutamatergic neurons",]
  z$cell = ordered(z$cell, levels=c("Up", "Down", "Merged"))
  z$contrast = ordered(z$contrast, levels=c("GABAergic neurons", "Glutamatergic neurons", "Oligodendrocytes", "Microglia+Astrocytes", "Merged"))
  unique(z$contrast)
  z$assayAndRegion = ordered(paste(z$assay, z$cell), levels=rev(c("ATAC Up", "ATAC Down", "ATAC Merged")))
  z$assayAndContrast = paste(z$assay,z$contrast)
  
  #####
  # Fig. 2a :: Numbers of differentially accessible OCRs (BH-adjusted P-value < 0.05) stratified by cell type and direction of change
  mpdf("Fig_2_a", outDir=file.path(ROOT, "outputs"));
  print(ggplot(z,aes(contrast,assayAndRegion),height=40) +
          ggtitle("signifHeatMap") +
          scale_y_discrete(expand=c(0, 0)) +
          scale_x_discrete(expand=c(0, 0)) +
          theme_classic() +
          theme(axis.text=element_text(colour="black")) +
          coord_fixed() +
          theme(axis.text.x=element_text(angle=45, hjust=1)) +
          theme(axis.title=element_blank()) +
          theme(legend.title=element_text(size=12, face="bold")) +
          geom_tile(aes(fill=frac)) + scale_fill_gradientn(colours=myPalette(100),name="frac significant") +
          geom_text(aes(label=signifCount))
  )
  dev.off()
}

####################################################################################################
##### FIG. 2B-C :: HERITABILITY ENRICHMENTS OF SCZ AND OTHER TRAITS IN DIFF OCR SETS ###############

{
  # Read results of LDsc run
  ldscScores = read.csv(file.path(ROOT, "inputs", "ldsc_results.tsv"), sep="\t", stringsAsFactors=F)
  
  # Custom fixes for more interpretable labels in plots
  ldscScores$analysisType = as.factor(sapply(ldscScores$annoID, function(x) ifelse(grepl(x=x, "P_01"), "P_01", ifelse(grepl(x=x, "P_05"), "P_05", "FDR"))))
  ldscScores$annoID = gsub("\\.1000bp.all", "", ldscScores$annoID)
  
  # Keeping only version with padding = 1000b on either side (other versions of padding: 0bp, 500bp) and only FDR / P-value<0.05 sets (we previously also tried P<0.01); GLU excluded entirely due to low coverage
  ldscScores = ldscScores[(ldscScores$ldscPadding == 1000) & (ldscScores$analysisType %in% c("FDR", "P_05")) & (!startsWith(ldscScores$annoName, prefix="GLU")),]
  ldscScores$direction = sapply(ldscScores$annoID, function(x) x=strsplit(x, "_")[[1]][length(strsplit(x, "_")[[1]])])
  
  # Preparation for plotting
  colorScheme = c(Specific="#9C3A1F",All="#666666") #1B9E77 #A6761D #752D19
  colScale = scale_colour_manual(name = "Open\nchromatin",values = colorScheme)
  ldsc = ldscScores[ldscScores$gwasAcronym=="sz3",] 
  ldsc$P.value = 10^-ldsc$minus_log10_p_regression
  ldsc$adj.P.value = p.adjust(ldsc$P.value, method="BH")
  ldsc$myLabel = ""
  ldsc$myLabel[ldsc$P.value<0.05] = "·"
  ldsc$myLabel[ldsc$adj.P.value<0.05] = "#"
  ldsc$myCoefficient=ldsc$Coefficient
  ldsc$error_left = ldsc$Coefficient - ldsc$Coefficient_std_error
  ldsc$error_right = ldsc$Coefficient + ldsc$Coefficient_std_error
  
  # Plot Fig. 2b :: Heritability coefficients for SCZ risk variants across various sets of differentially accessible OCRs, stratified by direction of regulation (upregulated, downregulated, or both) and cell type. 
  fig2b_plot = ggplot(data = ldsc,aes(x = myCoefficient,y = annoName, color=analysisType)) + 
    geom_vline(xintercept =0, alpha = 0.5, linetype = "dotted") +
    geom_point(show.legend=T) +
    labs(size="-logP") +
    colScale +
    geom_errorbarh(height=0,size=1,aes(xmin = error_left,xmax = error_right,color=analysisType)) +
    geom_text(size=5,aes(label=myLabel)) +
    theme_classic() +
    theme(axis.text=element_text(colour="black")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.title = element_text(size = 12, face = "bold")) +
    xlab("Heritability Coefficient") + ylab("Annotation") +
    ggtitle(paste(unique(ldsc$sumstatName),collapse=" ")) +
    facet_grid(direction ~ .,scales="free_y",space="free_y") +
    theme(strip.text=element_text(face="bold",size=9.2)) +
    theme(axis.title.y = element_blank())+
    theme(plot.title = element_text(face="bold", size=13)) + scale_color_manual(values = c("FDR" = "#00441A", "P_05" = "#58B567"))
  mpdf("Fig_2_b", outDir=file.path(ROOT, "outputs"), width=7, height=5); print(fig2b_plot); dev.off(); 
  
  #####
  # Fig. 2c :: Enrichment of SCZ-upregulated OCRs (FDR < 0.05) per cell type for common variants associated with various brain-related disorders
  # Select traits for plotting
  SELECTED_TRAITS = c("pd_without_23andMe", "als2", "alzBellenguez", "bip2", "sz3", "mdd_without_23andMe", "eduAttainment", "cd", "bmi", "uc")
  
  # Prepare df for plotting (calc adj.p-val, adjust labels etc)
  ldsc = ldscScores[(ldscScores$gwasAcronym %in% SELECTED_TRAITS) & (ldscScores$analysisType %in% c("FDR", "P_05")) & (ldscScores$direction %in% c("up")),]
  ldsc$gwasAcronym = ordered(ldsc$gwasAcronym, levels=SELECTED_TRAITS)
  ldsc = ldsc[order(ldsc$gwasAcronym),]
  ldsc$sumstatName = ordered(ldsc$sumstatName, levels=unique(rev(ldsc$sumstatName)))
  ldsc$minus_log10_p_regression = -log10(ldsc$p_regression)
  ldsc$plotLabel = ""
  ldsc$plotLabel[ldsc$p_regression < 0.05] = "·"
  ldsc$plotLabel[p.adjust(ldsc$p_regression, method="BH") < 0.05] = "#"
  plotTextSize = 9
  
  # Plot Fig. 2c :: Enrichment of SCZ-upregulated OCRs (FDR < 0.05) per cell type for common variants associated with various brain-related disorders
  fig2c_plot = ggplot(ldsc, aes(sumstatName, annoName, fill = minus_log10_p_regression)) + geom_tile() + scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + ylab("Trait") + 
    xlab("Annotation") + 
    theme_classic(base_size = plotTextSize) + 
    theme(axis.text = element_text(colour = "black")) + 
    coord_fixed() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.title = element_text(size = 10, face = "bold")) + 
    geom_tile(aes(fill = minus_log10_p_regression)) + 
    scale_fill_gradientn(colours = myPalette(100), name = "-logP") + 
    geom_text(aes(label = plotLabel), size = plotTextSize * 0.55)
  mpdf("Fig_2_c", outDir=file.path(ROOT, "outputs"), width=7, height=5); print(fig2c_plot); dev.off();
  
  #####
  # Table SX2 :: Heritability enrichment for sets of SCZ-associated OCRs.
  # Custom fixes of GWAS naming and re-formatting columns for printing
  fullsumstatConvertor = list("Bipolar Disorder 2"="Bipolar disorder", "BMI"="Body mass index", "Alzheimer's Disease Bellenguez"="Alzheimer's disease",
                              "Crohn's Disease"="Crohn's disease", "Schizophrenia 3"="Schizophrenia",
                              "Depression (without 23andMe)"="Depression", "Parkinson Disease (without 23andMe)"="Parkinson disease",
                              "Ulcerative Colitis"="Ulcerative colitis")
  table_ldsc = ldscScores[ldscScores$gwasAcronym %in% SELECTED_TRAITS,]
  table_ldsc$sumstatName = sapply(table_ldsc$sumstatName, function(x) ifelse(x %in% names(fullsumstatConvertor), fullsumstatConvertor[[x]], x))
  table_ldsc$celltype = sapply(table_ldsc$annoID, function(x) strsplit(x,"_")[[1]][1])
  selColumns = c("annoID", "celltype", "analysisType", "direction", "sumstatName", "gwasPubmed", "Prop._h2", "Prop._h2_std_error", "Enrichment", "Enrichment_std_error", "Enrichment_p", "Coefficient", "Coefficient_std_error", "p_regression")
  table_ldsc = table_ldsc[,selColumns]
  colnames(table_ldsc) = c("Category - fullname", "Cell type", "Significance", "Direction", "GWAS Trait", "GWAS Pubmed", "Prop. h2", "Prop. h2 std error", "Enrichment", "Enrichment std error", "Enrichment p", "Coefficient", "Coefficient std error", "Coefficient p")
  
  # Write Table SX2
  mtsv(table_ldsc, filename="Table_SX2_ldsc", outDir=file.path(ROOT, "outputs"), myHeader=T)
}

####################################################################################################
##### FIG. 2F :: PRIORITIZED TF GENES ##############################################################

{
  # Load HOMER & TOBIAS results
  tfList = readRDS(file.path(ROOT, "inputs", "homer_and_tobias.RDS"))
  
  # Definition of "blacklisted motifs that we don't use because there are better alternatives in the results for the same TFs
  blacklistedMotifName = c("ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer", "ETS:E-box(ETS,bHLH)/HPC7-Scl-ChIP-Seq(GSE22178)/Homer", "ETS(ETS)/Promoter/Homer", 
                           "OCT:OCT(POU,Homeobox)/NPC-OCT6-ChIP-Seq(GSE43916)/Homer", "OCT:OCT(POU,Homeobox)/NPC-Brn1-ChIP-Seq(GSE35496)/Homer", "OCT:OCT(POU,Homeobox,IR1)/NPC-Brn2-ChIP-Seq(GSE35496)/Homer", 
                           "OCT:OCT-short(POU,Homeobox)/NPC-OCT6-ChIP-Seq(GSE43916)/Homer", "RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer", "Tcf3(HMG)/mES-Tcf3-ChIP-Seq(GSE11724)/Homer", "E2A(bHLH),near_PU.1/Bcell-PU.1-ChIP-Seq(GSE21512)/Homer",
                           "Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer", "RBPJ:Ebox(?,bHLH)/Panc1-Rbpj1-ChIP-Seq(GSE47459)/Homer", "Stat3+il21(Stat)/CD4-Stat3-ChIP-Seq(GSE19198)/Homer", 
                           "STAT6(Stat)/Macrophage-Stat6-ChIP-Seq(GSE38377)/Homer", "Tcf12(bHLH)/GM12878-Tcf12-ChIP-Seq(GSE32465)/Homer", "THRb(NR)/HepG2-THRb.Flag-ChIP-Seq(Encode)/Homer")
  
  # First result-filtering, i.e. keep only "up" & "down" (not "all" which is "up"+"down") & remove suboptimal motifs that have better alternatives in the results
  tfDf_complete = do.call("rbind", tfList)
  tfDf_complete$pc1_corr_pearsonAbs = abs(tfDf_complete$pc1_corr_pearson)
  tfDf_complete = tfDf_complete[(tfDf_complete$Direction %in% c("up", "down")) & (!tfDf_complete$Motif.Name %in% blacklistedMotifName),]
  
  # Second result-filtering, i.e. keeping only motif that are significantly enriched (FDR<0.05) and correlated with expression of predicted downstream genes
  tfDf = tfDf_complete[(tfDf_complete$adj.P.value < 0.05) & (tfDf_complete$pc1_corr_pearson_pval < 0.05),]
  
  # Perform hierarchical clustering on -log(P-val) scores (to reorder dot-heatmap) 
  myPalette = colorRampPalette(brewer.pal(9, "Greens")[3:9], space="Lab")
  mat = do.call("cbind.data.frame", lapply(unique(tfDf_complete$Gene), function(gene) {
    sapply(unique(tfDf_complete$cat), function(cat) ifelse(length(tfDf_complete[(tfDf_complete$Gene == gene) & (tfDf_complete$cat == cat),"Log.P.value"]) == 0, 0, tfDf_complete[(tfDf_complete$Gene == gene) & (tfDf_complete$cat == cat),"Log.P.value"]))
  }))
  mat = data.frame(t(mat))
  colnames(mat) = unique(tfDf_complete$cat)
  rownames(mat) = unique(tfDf_complete$Gene)
  clust = hclust(dist(mat %>% as.matrix()))
  
  # Plot Fig. 2f :: Prioritized TF genes
  tfDf_complete = tfDf_complete[(tfDf_complete$Gene %in% unique(tfDf$Gene)) & (tfDf_complete$Motif.Name %in% unique(tfDf$Motif.Name)) & (tfDf_complete$pc1_corr_pearson_pval < 0.05), ]
  fig2f_plot = tfDf_complete %>%  
    mutate(pc1_corr_pearsonAbs, Gene = factor(Gene, levels = clust$labels[clust$order]), visible = ifelse(adj.P.value < 0.05, TRUE, FALSE)) %>% 
    ggplot(aes(y=cat, x=Gene, color = pc1_corr_pearsonAbs, size = minus_Log.P.value)) + geom_point(aes(size = -Log.P.value, alpha = visible)) + 
    cowplot::theme_cowplot() + theme(axis.line  = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab('') + theme(axis.ticks = element_blank()) + scale_color_gradientn(colours = myPalette(100),name="abs(TF target corr)") + coord_flip()
  mpdf("Fig_2_f", outDir=file.path(ROOT, "outputs"), width=4, height=12); print(fig2f_plot); dev.off();
  
  #####
  # Table SX3 :: Summary of TF motif enrichment and footprinting analysis results across SCZ-associated OCRs
  selCols = c("Motif", "Consensus", "Gene_ID", "Gene", "ctype", "cat", "Direction", "P.value", "adj.P.value", "Number.Target.Sequences.with.Motif", "Number.Background.Sequences.with.Motif", "Pct.Target.Sequences.with.Motif", "Pct.Background.Sequences.with.Motif", "pc1_corr_pearson", "pc1_corr_pearson_pval", "pc2_corr_pearson", "pc2_corr_pearson_pval", "pc3_corr_pearson", "pc3_corr_pearson_pval")
  table_tf = do.call("rbind.data.frame", lapply(unique(tfDf$Gene), function(gene) {
    tfDf[which(tfDf$Gene == gene),selCols]
  }))
  
  # Write Table SX3
  mtsv(table_tf, filename="Table_SX3_TF", outDir=file.path(ROOT, "outputs"), myHeader=T)
}

####################################################################################################
##### FIG. 4A :: NUMBERS OF DIFFERENTIALLY EXPRESSED GENES #########################################

{
  degAnalysis = new.env(); load(DEG_ANALYSIS, envir=degAnalysis)
  
  rnaMerged = list(
    "GABAergic neurons" = degAnalysis$dacResults$dac$GABA.SCZ_Control[degAnalysis$dacResults$dac$GABA.SCZ_Control$adj.P.Val<0.05,],
    "Glutamatergic neurons" = degAnalysis$dacResults$dac$GLU.SCZ_Control[degAnalysis$dacResults$dac$GLU.SCZ_Control$adj.P.Val<0.05,],
    "Oligodendrocytes" = degAnalysis$dacResults$dac$Olig.SCZ_Control[degAnalysis$dacResults$dac$Olig.SCZ_Control$adj.P.Val<0.05,],
    "Microglia+Astrocytes" = degAnalysis$dacResults$dac$MgAs.SCZ_Control[degAnalysis$dacResults$dac$MgAs.SCZ_Control$adj.P.Val<0.05,],
    "Merged" = degAnalysis$dacResults$dac$SCZ_Control[degAnalysis$dacResults$dac$SCZ_Control$adj.P.Val<0.05,]
  )
  rnaMerged = shrinkAtacToptables(rnaMerged)
  z=names(rnaMerged)
  rnaMerged_metadata = data.frame(Set=z, overallAssay="atac", cell="Merged", assay="rna", contrast=gsub("(rnaGABA__|_All$|_DLPFC|_ACC$)","",z), brainRegion=gsub(".+_","",z), stringsAsFactors=F)
  rnaMerged_metadata$SetFullName = with(rnaMerged_metadata, paste(overallAssay, cell, contrast, brainRegion))
  rnaMerged_metadata = addDacExclStats(rnaMerged_metadata, rnaMerged, nrow(degAnalysis$modeledVoomObj$E))
  rm(z)
  
  rnaUp = list(
    "GABAergic neurons" = degAnalysis$dacResults$dacUp$GABA.SCZ_Control[degAnalysis$dacResults$dacUp$GABA.SCZ_Control$adj.P.Val<0.05,],
    "Glutamatergic neurons" = degAnalysis$dacResults$dacUp$GLU.SCZ_Control[degAnalysis$dacResults$dacUp$GLU.SCZ_Control$adj.P.Val<0.05,],
    "Oligodendrocytes" = degAnalysis$dacResults$dacUp$Olig.SCZ_Control[degAnalysis$dacResults$dacUp$Olig.SCZ_Control$adj.P.Val<0.05,],
    "Microglia+Astrocytes" = degAnalysis$dacResults$dacUp$MgAs.SCZ_Control[degAnalysis$dacResults$dacUp$MgAs.SCZ_Control$adj.P.Val<0.05,],
    "Merged" = degAnalysis$dacResults$dacUp$SCZ_Control[degAnalysis$dacResults$dacUp$SCZ_Control$adj.P.Val<0.05,]
  )
  rnaUp = shrinkAtacToptables(rnaUp)
  z=names(rnaUp)
  rnaUp_metadata = data.frame(Set=z, overallAssay="rna", cell="Up", assay="rna", contrast=gsub("(rnaGABA__|_All$|_DLPFC|_ACC$)","",z), brainRegion=gsub(".+_","",z), stringsAsFactors=F)
  rnaUp_metadata$SetFullName = with(rnaUp_metadata, paste(overallAssay, cell, contrast, brainRegion))
  rnaUp_metadata = addDacExclStats(rnaUp_metadata, rnaUp, nrow(degAnalysis$modeledVoomObj$E))
  rm(z)
  
  rnaDown = list(
    "GABAergic neurons" = degAnalysis$dacResults$dacUp$GABA.Control_SCZ[degAnalysis$dacResults$dacUp$GABA.Control_SCZ$adj.P.Val<0.05,],
    "Glutamatergic neurons" = degAnalysis$dacResults$dacUp$GLU.Control_SCZ[degAnalysis$dacResults$dacUp$GLU.Control_SCZ$adj.P.Val<0.05,],
    "Oligodendrocytes" = degAnalysis$dacResults$dacUp$Olig.Control_SCZ[degAnalysis$dacResults$dacUp$Olig.Control_SCZ$adj.P.Val<0.05,],
    "Microglia+Astrocytes" = degAnalysis$dacResults$dacUp$MgAs.Control_SCZ[degAnalysis$dacResults$dacUp$MgAs.Control_SCZ$adj.P.Val<0.05,],
    "Merged" = degAnalysis$dacResults$dacUp$Control_SCZ[degAnalysis$dacResults$dacUp$Control_SCZ$adj.P.Val<0.05,]
  )
  rnaDown = shrinkAtacToptables(rnaDown)
  z=names(rnaDown)
  rnaDown_metadata = data.frame(Set=z, overallAssay="rna", cell="Down", assay="rna", contrast=gsub("(rnaGABA__|_All$|_DLPFC|_ACC$)","",z), brainRegion=gsub(".+_","",z), stringsAsFactors=F)
  rnaDown_metadata$SetFullName = with(rnaDown_metadata, paste(overallAssay, cell, contrast, brainRegion))
  rnaDown_metadata = addDacExclStats(rnaDown_metadata, rnaDown, nrow(degAnalysis$modeledVoomObj$E))
  rm(z)
  
  z=rbind.fill(rnaUp_metadata, rnaDown_metadata, rnaMerged_metadata)
  z$cell = ordered(z$cell, levels=c("Up", "Down", "Merged"))
  z$contrast = ordered(z$contrast, levels=c("GABAergic neurons", "Glutamatergic neurons", "Oligodendrocytes", "Microglia+Astrocytes", "Merged"))
  z$assayAndRegion = paste(z$assay, z$cell)
  z$assayAndContrast = paste(z$assay,z$contrast)
  z$assayAndRegion = ordered(paste(z$assay, z$cell), levels=rev(c("rna Up", "rna Down", "rna Merged")))
  
  #####
  # Fig. 4a :: Numbers of differentially expressed genes (DEGs) per cell type. b, Overlap of DEGs between cell types
  mpdf("Fig_4_a", outDir=file.path(ROOT, "outputs"));
  print(ggplot(z,aes(contrast,assayAndRegion),height=40) +
          ggtitle("signifHeatMap") +
          scale_y_discrete(expand=c(0, 0)) +
          scale_x_discrete(expand=c(0, 0)) +
          theme_classic() +
          theme(axis.text=element_text(colour="black")) +
          coord_fixed() +
          theme(axis.text.x=element_text(angle=45, hjust=1)) +
          theme(axis.title=element_blank()) +
          theme(legend.title=element_text(size=12, face="bold")) +
          geom_tile(aes(fill=frac)) + scale_fill_gradientn(colours=myPalette(100),name="frac significant") +
          geom_text(aes(label=signifCount))
  )
  dev.off()
}
  
####################################################################################################
##### FIG. 5A :: NUMBERS OF DIFFERENTIALLY EXPRESSED TRANSCRIPTS ###################################
  
{
  rnaseqDET = new.env(); load(DET_ANALYSIS, envir=rnaseqDET)
  sheet_names = excel_sheets(REMACOR_ANALYSIS)
  list_of_data_frames = lapply(sheet_names, read_excel, path="~/drive_lab/ROUSSOS_LAB_SHARED/Manuscripts/Molecular_Profiling/Tables/Table_S8.xlsx")
  table_s7 = lapply(list_of_data_frames, as.data.frame)
  names(table_s7) = sheet_names
  table_s7$Combined = do.call("rbind.data.frame", table_s7[c("GABA_meta", "GLU_meta", "OLIG_meta", "MGAS_meta")])
  table_s7$Combined = table_s7$Combined[order(table_s7$Combined$BH_RE2C),]
  table_s7$Combined = table_s7$Combined[!duplicated(table_s7$Combined$gene),]
  
  colnames(table_s7$GABA_meta) = gsub("BH_RE2C", "adj.P.Val", colnames(table_s7$GABA_meta))
  colnames(table_s7$GLU_meta) = gsub("BH_RE2C", "adj.P.Val", colnames(table_s7$GLU_meta))
  colnames(table_s7$OLIG_meta) = gsub("BH_RE2C", "adj.P.Val", colnames(table_s7$OLIG_meta))
  colnames(table_s7$MGAS_meta) = gsub("BH_RE2C", "adj.P.Val", colnames(table_s7$MGAS_meta))
  colnames(table_s7$Combined) = gsub("BH_RE2C", "adj.P.Val", colnames(table_s7$MGAS_meta))
  
  rnaRemacor = list(
    "GABAergic neurons" = table_s7$GABA_meta,
    "Glutamatergic neurons" = table_s7$GLU_meta,
    "Oligodendrocytes" = table_s7$OLIG_meta,
    "Microglia+Astrocytes" = table_s7$MGAS_meta
  )
  rnaRemacor = shrinkAtacToptables(rnaRemacor)
  z = names(rnaRemacor)
  rnaRemacor_metadata = data.frame(Set=z, overallAssay="rna", cell="Remacor", assay="rna", contrast=gsub("(rnaGABA__|_All$|_DLPFC|_ACC$)","",z), brainRegion=gsub(".+_","",z), stringsAsFactors=F)
  rnaRemacor_metadata$SetFullName = with(rnaRemacor_metadata, paste(overallAssay, cell, contrast, brainRegion))
  rnaRemacor_metadata = addDacExclStats(rnaRemacor_metadata, rnaRemacor, nrow(rnaseq_countMatrixRaw))
  rm(z)
  
  ##########
  sheet_names = excel_sheets("~/drive_lab/ROUSSOS_LAB_SHARED/Manuscripts/Molecular_Profiling/Tables/Table_S8.xlsx")
  list_of_data_frames = lapply(sheet_names, read_excel, path="~/drive_lab/ROUSSOS_LAB_SHARED/Manuscripts/Molecular_Profiling/Tables/Table_S8.xlsx")
  table_s7 = lapply(list_of_data_frames, as.data.frame)
  names(table_s7) = sheet_names
  table_s7$Combined = do.call("rbind.data.frame", table_s7[c("GABA_meta", "GLU_meta", "OLIG_meta", "MGAS_meta")])
  table_s7$Combined = table_s7$Combined[order(table_s7$Combined$BH_RE2C),]
  table_s7$Combined = table_s7$Combined[!duplicated(table_s7$Combined$gene),]
  
  rnaMerged = list(
    "GABAergic neurons" = table_s7$GABA[(table_s7$GABA$adj.P.Val<0.05),],
    "Glutamatergic neurons" = table_s7$GLU[(table_s7$GLU$adj.P.Val<0.05),],
    "Oligodendrocytes" =  table_s7$OLIG[(table_s7$OLIG$adj.P.Val<0.05),],
    "Microglia+Astrocytes" = table_s7$MGAS[(table_s7$MGAS$adj.P.Val<0.05),]
  )
  rnaMerged = shrinkAtacToptables(rnaMerged)
  z=names(rnaMerged)
  rnaMerged_metadata = data.frame(Set=z, overallAssay="atac", cell="Merged", assay="rna", contrast=gsub("(rnaGABA__|_All$|_DLPFC|_ACC$)","",z), brainRegion=gsub(".+_","",z), stringsAsFactors=F)
  rnaMerged_metadata$SetFullName = with(rnaMerged_metadata, paste(overallAssay, cell, contrast, brainRegion))
  rnaMerged_metadata = addDacExclStats(rnaMerged_metadata, rnaMerged, nrow(rnaseqDET$initialDgeObj))
  rm(z)
  
  rnaUp = list(
    "GABAergic neurons" = table_s7$GABA[(table_s7$GABA$adj.P.Val<0.05) & (table_s7$GABA$logFC > 0),],
    "Glutamatergic neurons" = table_s7$GLU[(table_s7$GLU$adj.P.Val<0.05) & (table_s7$GLU$logFC > 0),],
    "Oligodendrocytes" =  table_s7$OLIG[(table_s7$OLIG$adj.P.Val<0.05) & (table_s7$OLIG$logFC > 0),],
    "Microglia+Astrocytes" = table_s7$MGAS[(table_s7$MGAS$adj.P.Val<0.05) & (table_s7$MGAS$logFC > 0),]
  )
  rnaUp = shrinkAtacToptables(rnaUp)
  z=names(rnaUp)
  rnaUp_metadata = data.frame(Set=z, overallAssay="rna", cell="Up", assay="rna", contrast=gsub("(rnaGABA__|_All$|_DLPFC|_ACC$)","",z), brainRegion=gsub(".+_","",z), stringsAsFactors=F)
  rnaUp_metadata$SetFullName = with(rnaUp_metadata, paste(overallAssay, cell, contrast, brainRegion))
  rnaUp_metadata = addDacExclStats(rnaUp_metadata, rnaUp, nrow(rnaseqDET$initialDgeObj))
  rm(z)
  
  rnaDown = list(
    "GABAergic neurons" = table_s7$GABA[(table_s7$GABA$adj.P.Val<0.05) & (table_s7$GABA$logFC < 0),],
    "Glutamatergic neurons" = table_s7$GLU[(table_s7$GLU$adj.P.Val<0.05) & (table_s7$GLU$logFC < 0),],
    "Oligodendrocytes" =  table_s7$OLIG[(table_s7$OLIG$adj.P.Val<0.05) & (table_s7$OLIG$logFC < 0),],
    "Microglia+Astrocytes" = table_s7$MGAS[(table_s7$MGAS$adj.P.Val<0.05) & (table_s7$MGAS$logFC < 0),]
  )
  rnaDown = shrinkAtacToptables(rnaDown)
  z=names(rnaDown)
  rnaDown_metadata = data.frame(Set=z, overallAssay="rna", cell="Down", assay="rna", contrast=gsub("(rnaGABA__|_All$|_DLPFC|_ACC$)","",z), brainRegion=gsub(".+_","",z), stringsAsFactors=F)
  rnaDown_metadata$SetFullName = with(rnaDown_metadata, paste(overallAssay, cell, contrast, brainRegion))
  rnaDown_metadata = addDacExclStats(rnaDown_metadata, rnaDown, nrow(rnaseqDET$initialDgeObj))
  rm(z)
  
  z=rbind.fill(rnaUp_metadata, rnaDown_metadata, rnaMerged_metadata, rnaRemacor_metadata)
  z$cell = ordered(z$cell, levels=c("Up", "Down", "Merged"))
  z$contrast = ordered(z$contrast, levels=c("GABAergic neurons", "Glutamatergic neurons", "Oligodendrocytes", "Microglia+Astrocytes", "Merged"))
  z$assayAndRegion = paste(z$assay, z$cell)
  z$assayAndContrast = paste(z$assay,z$contrast)
  z$assayAndRegion = ordered(paste(z$assay, z$cell), levels=rev(c("rna Up", "rna Down", "rna Merged")))
  
  #####
  # Fig. 5a :: Numbers of dysregulated transcripts as well as genes with at least one differentially expressed transcript (FDR < 0.05) stratified by cell type
  mpdf("Fig_5_a", outDir=file.path(ROOT, "outputs"));
  print(ggplot(z,aes(contrast,assayAndRegion),height=40) +
          ggtitle("signifHeatMap") +
          scale_y_discrete(expand=c(0, 0)) +
          scale_x_discrete(expand=c(0, 0)) +
          theme_classic() +
          theme(axis.text=element_text(colour="black")) +
          coord_fixed() +
          theme(axis.text.x=element_text(angle=45, hjust=1)) +
          theme(axis.title=element_blank()) +
          theme(legend.title=element_text(size=12, face="bold")) +
          geom_tile(aes(fill=frac)) + scale_fill_gradientn(colours=myPalette(100),name="frac significant") +
          geom_text(aes(label=signifCount))
  )
  dev.off()
}

####################################################################################################
##### FIG. 5D,F :: COMPARISON WITH KOZLENKOV ET AL. 2023 ##############################################

{
  # Load pre-processed data from Kozlenkov et al. 2023
  TRANSCRIPT_ANALYSIS = file.path(ROOT, "inputs", "transcript_analysis.Rdata")
  rnaseqTranscriptEnv = new.env(); load(TRANSCRIPT_ANALYSIS, envir=rnaseqTranscriptEnv)
  
  # Get normalized count matrix (effect of technical covariates was regressed out)
  mx = rnaseqTranscriptEnv$residualized_DxBrainRegion_EffectKept
  
  # Get IDs of samples depending on their cell types (ODC / OPC) and age of donor (Adult / Infant)
  sampleGroups = sapply(unique(rnaseqTranscriptEnv$allInfo$Groups), function(group) rnaseqTranscriptEnv$allInfo[(rnaseqTranscriptEnv$allInfo$Groups == group),"ID"])
  names(sampleGroups) = unique(rnaseqTranscriptEnv$allInfo$Groups)
  
  # Define transcripts of interests for CACNA1C and KMT5A; note that not all transcripts were sufficiently expressed (CPM > 1 in at least 20% of samples)
  cacna1c_transcripts = c("ENST00000496818", "ENST00000491104", "ENST00000492150", "ENST00000483136", "ENST00000399655", "ENST00000480911", "ENST00000465278", "ENST00000541871")
  cacna1c_transcripts = cacna1c_transcripts[cacna1c_transcripts %in% rownames(mx)]
  trim2_transcripts = c("ENST00000338700", "ENST00000502281", "ENST00000460908", "ENST00000494872", "ENST00000482578", "ENST00000632856", "ENST00000437508", "ENST00000491446")
  trim2_transcripts = trim2_transcripts[trim2_transcripts %in% rownames(mx)]
  
  #####
  # Fig. 5d :: Comparison of ENST00000465278 and ENST00000483136 expression in the study profiling OPC and mature oligodendrocytes (MO) in infants and adults
  cacna1c_exp = do.call("rbind.data.frame", lapply(cacna1c_transcripts, function(transcriptId) {
    unlist(sapply(sampleGroups, function(sGroup) { 
      mean(mx[transcriptId, sGroup])
    }))
  }))
  colnames(cacna1c_exp) = names(sampleGroups)
  rownames(cacna1c_exp) = cacna1c_transcripts
  
  cacna1c_exp2 = reshape::melt(mx[cacna1c_transcripts,])
  cacna1c_exp2$group = rnaseqTranscriptEnv$allInfo[match(cacna1c_exp2$X2, rnaseqTranscriptEnv$allInfo$ID), "Groups"]
  
  cacna1c_plot = ggplot(data=cacna1c_exp2[cacna1c_exp2$X1 %in% c("ENST00000483136", "ENST00000465278"),], aes(x=X1, y=value, fill=group)) + geom_boxplot() + labs(title="", x="Transcript ID", y="Expression") + 
    theme_classic() + theme(axis.text.x = element_text(angle = 90))
  mpdf("Fig_5_d", outDir=file.path(ROOT, "outputs"), width=8, height=5); print(cacna1c_plot); dev.off();
  
  #####
  # Fig. 5f :: Comparison of ENST00000338700 and ENST00000460908 expression in the study profiling OPC and MO in infants and adults
  trim2_exp = do.call("rbind.data.frame", lapply(trim2_transcripts, function(transcriptId) {
    unlist(sapply(sampleGroups, function(sGroup) { 
      mean(mx[transcriptId, sGroup])
    }))
  }))
  colnames(trim2_exp) = names(sampleGroups)
  rownames(trim2_exp) = trim2_transcripts
  print(trim2_exp)
  
  trim2_exp2 = reshape::melt(mx[trim2_transcripts,])
  trim2_exp2$group = rnaseqTranscriptEnv$allInfo[match(trim2_exp2$X2, rnaseqTranscriptEnv$allInfo$ID), "Groups"]
  
  trim2_plot = ggplot(data=trim2_exp2[trim2_exp2$X1 %in% c("ENST00000338700", "ENST00000460908"),], aes(x=X1, y=value, fill=group)) + geom_boxplot() + labs(title="", x="Transcript ID", y="Expression") + 
    theme_classic() + theme(axis.text.x = element_text(angle = 90))
  trim2_plot
  mpdf("Fig_5_f", outDir=file.path(ROOT, "outputs"), width=8, height=5); print(trim2_plot); dev.off();
}
