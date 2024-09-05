library(ggplot2)
library(data.table)
library(dplyr)
library(plyr)
library(ggpubr)
library(cowplot)
library(ggsci)
library(MASS)
library(viridis)
library(reshape)
library(readxl)

########################################################################################
##### CONFIG ###########################################################################

{
  ROOT = "~/Desktop/molecular_profiling/" # !!! FIXME: SET TO YOUR CUSTOM DIRECTORY !!!
  
  QC_ATACSEQ = file.path(ROOT, "inputs", "qc_all_atac.tsv")  # Pre-calculated QC metrics for ATAC-seq samples from processing computational pipeline
  QC_RNASEQ = file.path(ROOT, "inputs", "qc_all_rna.tsv")    # Pre-calculated QC metrics for RNA-seq samples from processing computational pipeline
  KINSHIP_ATACSEQ_SNPPARRAY = file.path(ROOT, "inputs", "kinship_atacseq_snparray.csv")  # Pre-calculated comparison between SNPs called from ATAC-seq reads and SNParrays 
  KINSHIP_RNASEQ_SNPPARRAY = file.path(ROOT, "inputs", "kinship_rnaseq_snparray.csv")    # Pre-calculated comparison between SNPs called from RNA-seq reads and SNParrays 
  
  ATACSEQ_PEAKS =  file.path(ROOT, "inputs", "atacseq_peaks.RDS")                        # Peaks called from ATAC-se data
  ATACSEQ_COUNT_MATRIX_RAW = file.path(ROOT, "inputs", "atacseq_count_matrix_raw.RDS")   # Raw read count matrix for ATAC-seq data
  ATACSEQ_COUNT_MATRIX_ADJ = file.path(ROOT, "inputs", "atacseq_count_matrix_adj.RDS")   # Covariate-adjusted read count matrix for ATAC-seq data
  RNASEQ_COUNT_MATRIX_RAW = file.path(ROOT, "inputs", "rnaseq_count_matrix_raw.RDS")     # Raw read count matrix for RNA-seq data
  RNASEQ_COUNT_MATRIX_ADJ = file.path(ROOT, "inputs", "rnaseq_count_matrix_adj.RDS")     # Covariate-adjusted read count matrix for RNA-seq data
  
  DAC_ANALYSIS = file.path(ROOT, "inputs", "DAC_Analysis.Rdata")  # Pre-calculated results for analysis of differential chromatin accessibility
  DEG_ANALYSIS = file.path(ROOT, "inputs", "DEG_Analysis.Rdata")  # Pre-calculated results for analysis of differential gene expression 
  DET_ANALYSIS = file.path(ROOT, "inputs", "DET_Analysis.Rdata")  # Pre-calculated results for analysis of differential transcript expression 
  REMACOR_ANALYSIS = file.path(ROOT, "inputs", "REMACOR_ANALYSIS.xlsx")
  
  METADATA_HAUBERG_2020 = file.path(ROOT, "inputs", "hauberg_2020_metadata.csv")  # Metadata for samples from Hauberg et al 2020 (dataset used for comparison)
  GEXPR_HAUBERG_2020 = file.path(ROOT, "inputs", "hauberg_2020_gExpr.RDS")        # Covariate-adjusted read count matrix for FANS ATAC-seq data from Hauberg et al 2020
  PEAKS_HAUBERG_2020 = file.path(ROOT, "inputs", "hauberg_2020_peaks.RDS")        # Peaks called from ATAC-seq data from Hauberg et al 2020
  
  METADATA_COLEMAN_2023 = file.path(ROOT, "inputs", "coleman_2023_metadata.csv")  # Metadata for samples from Coleman et al 2023 (dataset used for comparison)
  GEXPR_COLEMAN_2023 = file.path(ROOT, "inputs", "coleman_2023_gExpr.RDS")        # Covariate-adjuste read count matrix from FANS RNA-seq data from Coleman et al 2023
  
  phg_initDgeObj = readRDS("/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/tmp/phg_rnaseqInitialDgeObj.RDS")
  phg_initVoomObj = readRDS("/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/tmp/phg_rnaseqInitialVoomObj.RDS")
  
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
  
  myPalette = colorRampPalette(c("#F7FCF5", "#E5F5E0", 
                                 "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", 
                                 "#006D2C", "#00441B"), space = "Lab")
  
  dir.create(file.path(ROOT, "outputs"))
}

########################################################################################
##### HELPFR FUNCTIONS #################################################################

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
}

####################################################################################################
##### FIG. S1 (RNA-SEQ):: QUALITY CONTROL FOR ######################################################

{
  qcRna = read.csv(QC_RNASEQ, sep="\t")
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
  mtsv(qcRnaSum, filename="Fig_S1_background_rnaseq", outDir=file.path(ROOT, "outputs"), myHeader=T)
  
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
    draw_plot(allPlots[["picard_PERCENT_DUPLICATION"]],       .00, .33, .33, .33) +
    draw_plot(allPlots[["rnaseqc_PCT_INTERGENIC_BASES"]],     .33, .33, .33, .33) +
    draw_plot(allPlots[["rnaseqc_PCT_INTRONIC_BASES"]],       .66, .33, .33, .33) +
    draw_plot(allPlots[["picard_meanGcContent"]],             .00, .00, .33, .33) +
    draw_plot(allPlots[["insertMetrics_MEDIAN_INSERT_SIZE"]], .33, .00, .33, .33) +
    draw_plot_label(c("A", "B", "C", "D", "E", "F", "G", "H"),
                    c(.00, .33, .66, .00, .33, .66, .00, .33),
                    c(.99, .99, .99, .66, .66, .66, .33, .33),
                    size = 15)
  
  mpdf("Fig_S1_d_e_f_g_h_i", outDir=file.path(ROOT, "outputs"), width=11, height=11); print(plotQc2); dev.off()

  #####
  # Fig. S1a :: Median read insert size distribution
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
  histMedianInsertSize
  mpdf("Fig_S1_a", outDir=file.path(ROOT, "outputs")); print(histMedianInsertSize); dev.off()
  
  #####
  # Fig. S1b :: Sex check based on measuring the number reads mapped on chromosome Y
  chrY_genes = qcRna$qcPeakAnno[(qcRna$qcPeakAnno$seqnames=="chrY") & (qcRna$qcPeakAnno$PeakID %in% qcRna$initialDgeObj$genes$PeakID),]
  chrY_genes = chrY_genes[(chrY_genes$end < 10001 | chrY_genes$start > 2781479) & (chrY_genes$end < 155701383 | chrY_genes$start > 156030895),]
  
  #####  # Fig. S1c :: Genotype check based on pair-wise comparison of genotypes called from RNA-seq samples with SNP-arrays
  kinshipRnaSnparray = read.csv(KINSHIP_RNASEQ_SNPPARRAY)
  z = kinshipRnaSnparray
  z$`Same person`= ordered(ifelse(z$samePerson, "yes", "no"), levels=c("yes", "no"))
  kinshipRnaSnparray = ggplot(z, aes(Kinship)) + scale_color_manual(name="Same person", labels=c("yes", "no"), values = c(npgList[["NEURON"]], npgList[["GLIA"]])) + 
    theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    geom_density(alpha=0, aes(color=`Same person`, fill=`Same person`), size=1) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = c(0.4, 0.85),
          axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + coord_equal() + xlab("Kinship score") + ylab("Density")
  kinshipRnaSnparray
  mpdf("Fig_S1_c", outDir=file.path(ROOT, "outputs")); print(kinshipRnaSnparray); dev.off()
}

####################################################################################################
##### FIG. S1 (ATAC-SEQ) :: QUALITY CONTROL ########################################################

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
  mtsv(qcAtacSum, filename="Fig_S1_background_atacseq", outDir=file.path(ROOT, "outputs"), myHeader=T)
  
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
    draw_plot(allPlots[["chrMFrac"]],                         .00, .33, .33, .33) +
    draw_plot(allPlots[["peakNarrowFDR1pctCount"]],           .33, .33, .33, .33) +
    draw_plot(allPlots[["fracReadsInNonBlacklistedPeaks"]],   .66, .33, .33, .33) +
    draw_plot(allPlots[["pbc"]],             .00, .00, .33, .33) +
    draw_plot(allPlots[["insertMetrics_MEDIAN_INSERT_SIZE"]], .33, .00, .33, .33) +
    draw_plot_label(c("A", "B", "C", "D", "E", "F", "G", "H"),
                    c(.00, .33, .66, .00, .33, .66, .00, .33),
                    c(.99, .99, .99, .66, .66, .66, .33, .33),
                    size = 15)
  
  mpdf("Fig_S1_p_r_s_t_u_v_w", outDir=file.path(ROOT, "outputs"), width=11, height=11); print(plotQc2); dev.off()
  
  #####
  # Fig. S1j :: Median read insert size distribution
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
  mpdf("Fig_S1_j", outDir=file.path(ROOT, "outputs")); print(histMedianInsertSize); dev.off()
  
  #####
  # Fig. S1k :: Distance of OCRs from the closest TSS
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
  mpdf("Fig_S1_k", outDir=file.path(ROOT, "outputs")); print(histTssDist); dev.off()
  
  #####
  # Fig. S1n :: Sex check based on measuring the number reads mapped on OCRs located at chromosome Y (pseudoautosomal regions not counted)
  chrYplot = ggplot(qcAtac, aes(fracReadsInNonBlacklistedPeaks, chryCounts, color=Gender)) + 
    geom_point() + scale_color_manual(name="Sex", labels=c("Female", "Male "), values = c(npgList[["NEURON"]], npgList[["GLIA"]])) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1, legend.position = c(0.15, 0.85),
                       axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black")) + 
    xlab("Fraction of OCRs in peaks of open chromatin") + ylab("chrY read count") + coord_equal() 
  mpdf("Fig_S1_n", outDir=file.path(ROOT, "outputs")); print(chrYplot); dev.off()
  
  #####
  # Fig. S1o :: Genotype check based on pair-wise comparison of genotypes called from ATAC-seq samples with SNP-arrays
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
  mpdf("Fig_S1_o", outDir=file.path(ROOT, "outputs")); print(kinshipAtacSnparray); dev.off()
}

####################################################################################################
##### FIG. S3 (ATAC-SEQ) :: COMPARISON WITH HAUBERG ET AL. 2023 AND COLEMAN ET AL. 2023 ############

{
  # Fig. S3a :: Correlation of log2(cpm+1) counts between our ATAC-seq data and external ATAC-seq data from 4 cell types from the prefrontal cortex (Hauberg et al 2020)
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
      mpdf(paste0("Fig_S3_hauberg__", ctype, "__", ctype2), outDir=file.path(ROOT, "outputs")); print(densityScatter_current_ggom); dev.off()
    }
  }  
  
  #####
  # Fig. S3b :: Correlation of log2(cpm+1) counts between our RNA-seq data and external RNA-seq dataset of 3 cell types from parahippocampal gyrus (Coleman et al 2023)
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
      mpdf(paste0("Fig_S3_coleman__", ctype, "__", ctype2), outDir=file.path(ROOT, "outputs")); print(densityScatter_current_phg); dev.off()
    }
  }
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
##### FIG. 4A :: NUMBERS OF DIFFERENTIALLY EXPRESSED GENES ###################################

{
  DEG_ANALYSIS = file.path(ROOT, "inputs", "DEG_Analysis.Rdata")
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
