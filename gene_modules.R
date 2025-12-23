library(WGCNA)
library(dplyr)
library(ggplot2)
library(data.table)

################################################################################
##### CONFIG ###################################################################

{
  ROOT = "~/Desktop/molecular_profiling/molecular_profiling/" # !!! FIXME: SET TO YOUR CUSTOM DIRECTORY !!
  QC_ATACSEQ = file.path(ROOT, "inputs", "qc_all_atac.csv")   # Pre-calculated QC metrics for ATAC-seq samples from processing computational pipeline
  QC_RNASEQ = file.path(ROOT, "inputs", "qc_all_rna.csv")     # Pre-calculated QC metrics for RNA-seq samples from processing computational pipeline
  ATACSEQ_COUNT_MATRIX_RESIDUALIZED_DX_CELLTYPE_KEPT = file.path(ROOT, "inputs", "atacseq_count_matrix_residualized_Dx_CellType_kept.RDS") # Count matrix from which the effect of technical covariates were regressed out, but Dx & Cell type effect kept
  RNASEQ_COUNT_MATRIX_RESIDUALIZED_DX_CELLTYPE_KEPT = file.path(ROOT, "inputs", "rnaseq_count_matrix_residualized_Dx_CellType_kept.RDS")   # Count matrix from which the effect of technical covariates were regressed out, but Dx & Cell type effect kept
  ABC = file.path(ROOT, "inputs", "abc.RData")                    # Output of activity-by-contact (ABC) analysis containing E-P interactions between OCRs and genes
  ATACSEQ_PEAKS =  file.path(ROOT, "inputs", "atacseq_peaks.RDS") # Peaks called from ATAC-seq data
  DAC_ANALYSIS = file.path(ROOT, "inputs", "DAC_Analysis.Rdata")  # Pre-calculated results for analysis of differential chromatin accessibility
  
  allowWGCNAThreads(nThreads = 5)  # TODO: Check that this is fine for your machine
  
  options(stringsAsFactors = FALSE)
  
  CELL_TYPES = c("GABA", "GLU", "OLIG", "MGAS")
  
  # Load metadata and residualized count matrices for both ATAC-seq and RNA-seq data
  qcAtac = read.csv(QC_ATACSEQ)
  rownames(qcAtac) = qcAtac$ID
  qcRna = read.csv(QC_RNASEQ)
  rownames(qcRna) = qcRna$ID
  atacseq_countMatrixResi = readRDS(ATACSEQ_COUNT_MATRIX_RESIDUALIZED_DX_CELLTYPE_KEPT)
  rnaseq_countMatrixResi = readRDS(RNASEQ_COUNT_MATRIX_RESIDUALIZED_DX_CELLTYPE_KEPT)
  
  # Helper function
  mpdf = function(x, width=7,height=7, outDir=outDir, onefile=T) eval.parent(substitute({ pdf(paste0(outDir, "/", make.names(x),".pdf"), useDingbats=F, width=width, height=height, onefile=onefile) }))

  # Load E-P interactions defined by activity-by-contact (ABC) method
  abcEnv = new.env(); load(ABC, envir=abcEnv)
  
  # Load all OCRs & differential OCRs
  qcPeakAnno = readRDS(ATACSEQ_PEAKS)
  rownames(qcPeakAnno) = qcPeakAnno$PeakID
  dacAnalysis = new.env(); load(DAC_ANALYSIS, envir=dacAnalysis)
  
  # Load pre-calculated MAGMA results ()
  final = read.csv("/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/publication_plots/final.csv")  # Load ABC
  
  # Prepare for GSEA
  customBlacklist = file.path(ROOT, "inputs", "hg38.blacklist_combination_kundaje_new_and_custom.bed")
  geneSetFile = file.path(ROOT, "inputs", "gseafunctions_hg38.Rdata")
  selectedGeneMetaSets = c("msigdbSetsPruned", "msigdbSetsStrictlyPruned", "msigdbSetsVeryStrictlyPruned") # other choices: , "msigdbSets100Pruned", "brainSets.humanSingleCell", "brainSets.humanSingleLake", "brainSets.humanSingleLakeDetailed", "brainSets.humanSingleHabib", "brainSets.humanSingleHabibDetailed", "brainSets.lake", "brainSets.zhangAndZeisel", "brainSets.candidate_gene_sets", "brainSets.mckenzie", "brainSets.syngoToplevel", "brainSets.syngoAll", "brainSets.velmeshev", "brainSets.velmeshevGeneral"
  z = local({load(geneSetFile); environment()})
  standardGeneSets = z$standardFisher$standardGeneSets
  pasteNull=function(...){x=list(...); if(any(sapply(x,is.null))){NULL}else{paste(unlist(x),collapse="")}}; #return NULL if any arg is NULL. otherwise paste0
}

################################################################################
##### DATA PREPARATION #########################################################

{
  # Per-cell type gene coexpression results
  MEs_LIST = list()
  NETs_LIST = list()
  DO_GSEA = T             # Perform GSEA analysis
  GSEA = list()           # for saving per-cell-type GSEA results
  DAC_ENRICH = list()
  GWAS_ENRICH = list()
  MAGMA_ENRICH = list()    # for saving per-module MAGMA GWAS enrichment (enrichment is calculated for all genes belonging to the module and SCZ GWAS)
  EIGENCORREL = list()     # for saving correlations between OCRs associated with each module and module eigenvector
  EIGENCORRELSUM = list()  # for saving per-module summaries of OCRs that are correlated with their module eigenvector
  MODULE_COLOR = list()
  
  # Let's run gene coexpression analysis separately for each cell type
  for(CTYPE in CELL_TYPES) {
    ovrl_subjects = intersect(qcAtac[qcAtac$cell_subtype_abbreviation == CTYPE, "Subject_ID"], qcRna[qcRna$cell_subtype_abbreviation == CTYPE, "Subject_ID"])
    ids_atacseq = qcAtac[(qcAtac$cell_subtype == CTYPE),]
    ids_atacseq = ids_atacseq[match(ovrl_subjects, ids_atacseq$Subject_ID), "ID"]
    ids_rnaseq = qcRna[match(ovrl_subjects, qcRna$Subject_ID), "ID"]
    
    ###
    ### Network construction and module detection (incl. dendrogram plotting)
    {
      expr = t(rnaseq_countMatrixResi[,ids_rnaseq])
      chosen_power = 7
      
      net = blockwiseModules(expr,
                              power = chosen_power,
                              networkType = "signed",
                              TOMType = "signed",
                              corType = "bicor",
                              blockSize = 4000,
                              deepSplit = 2,
                              minModuleSize = 50,
                              mergeCutHeight = 0.25,
                              reassignThreshold = 0,
                              pamRespectsDendro = FALSE,
                              saveTOMs = FALSE,
                              verbose = 3)
      
      moduleColors = labels2colors(net$colors)
      NETs_LIST[[CTYPE]] = net
      
      # Quantify module sizes: You should see: (1) no module with <20 genes, (2) no single module with 90% of genes
      print(table(moduleColors))
      
      # Eigengenes and disease association
      MEs = moduleEigengenes(expr, moduleColors)$eigengenes
      MEs = orderMEs(MEs)
      MEs_LIST[[CTYPE]] = MEs
      
      # Gene dendrogram + module colors: this is the canonical WGCNA visualization
      mpdf(paste0("misc_gene_module_dendrogram_", CTYPE), outDir=file.path(ROOT, "outputs"), width=7, height=4);
      plotDendroAndColors(net$dendrograms[[1]],
                          moduleColors[net$blockGenes[[1]]],
                          "Modules",
                          dendroLabels = FALSE,
                          hang = 0.03,
                          addGuide = TRUE,
                          guideHang = 0.05,
                          main = paste0(CTYPE, " co-expression modules"))
      dev.off()
    }
    
    ###
    ### Eigengenes and disease association
    {
      datME = as.data.frame(MEs)
      datME$Dx = as.integer(ifelse(qcRna[ids_rnaseq, "Dx"] == "Control", 0, 1))
      
      moduleStats = sapply(colnames(MEs), function(m) {
        # fit ME_m ~ Dx using the columns of datME
        fit = lm(datME[, m] ~ datME$Dx)
        co = summary(fit)$coefficients
        beta = co[2, "Estimate"]   # effect of Dx
        p = co[2, "Pr(>|t|)"]      # p-value for Dx
        c(beta = beta, p = p)
      })
      
      moduleStats = as.data.frame(t(moduleStats))
      moduleStats$FDR = p.adjust(moduleStats$p, "BH")
      
      moduleStats_df = data.frame(
        Module = rownames(moduleStats),
        beta   = moduleStats[,"beta"],
        p   = moduleStats[,"p"],
        FDR    = moduleStats[,"FDR"],
        row.names = NULL,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
      eigenPlot = ggplot(moduleStats_df, aes(x = reorder(Module, beta), y = beta, fill = FDR < 0.05)) +
        geom_col() + coord_flip() + theme_minimal() + labs(title = "Diagnosis effect on module eigengenes",
                                                           y = "Beta (Dx effect)", x = "Module") + scale_fill_manual(values = c("grey70", "firebrick"), name = "FDR < 0.05")
      mpdf(paste0("misc_eigen_dx_correl_", CTYPE, ".pdf"), outDir=file.path(ROOT, "outputs")); print(eigenPlot); dev.off()
    }
    
    ###
    ### Pathway enrichment for genesets (genes belonging to each module)
    {
      moduleColors = labels2colors(net$colors)
      names(moduleColors) = colnames(expr)   # ensure gene names are attached
      
      geneModuleTable = data.frame(
        gene = names(moduleColors),
        module = moduleColors
      )
      moduleGeneList = split(geneModuleTable$gene, geneModuleTable$module)
      
      source(file.path(ROOT, "helper_functions.R"))
      if(DO_GSEA) {
        cellGseaAllWithBg = universalGsea(
          testMethod = "fisher",
          inputForTest = moduleGeneList, 
          inputForTestMetadata = NULL,
          GENOME_VERSION = "hg38",
          background = as.character(unique(unlist(moduleGeneList))), 
          myDataType = "genes",
          myGeneMetaSets = standardGeneSets[selectedGeneMetaSets],
          outDir = file.path(ROOT, "outputs", paste0("WGCNA_gsea_", CTYPE)),
          furtherArgs__geneRegDomains = list(blacklistFile=customBlacklist),
          shrinkOutput = T,
          forceNoOutDirCheck = T
        )
        GSEA[[CTYPE]] = cellGseaAllWithBg
      }
    }
    
    ###
    ### Define OCRs associated with genes of each module (either E-P linked OCRs or those in intronic, UTR, exon or promoter regions (only to the distance of 50bp from promoter))
    modulePeakList = lapply(names(moduleGeneList), function(moduleName) {
      abcPeaks = unique(abcEnv$abcResults_cell[[CTYPE]][abcEnv$abcResults_cell[[CTYPE]]$TargetGene %in% moduleGeneList[[moduleName]], "PeakID"])
      
      proxPeaks = qcPeakAnno[(qcPeakAnno$geneId %in% moduleGeneList[[moduleName]]),]
      proxPeaks = proxPeaks[(proxPeaks$annotationSimple != "Distal Intergenic"),]
      proxPeaks = proxPeaks[(proxPeaks$annotationSimple != "Promoter") | ((proxPeaks$annotationSimple == "Promoter") & (proxPeaks$distanceToTSS < 50)),]
      unique(c(proxPeaks$PeakID, abcPeaks))
    })
    names(modulePeakList) = names(moduleGeneList)
    
    atacOrig = t(atacseq_countMatrixResi[,ids_atacseq])
    
    # Eigengene–peak correlations (define "regulatory elements for a module")
    corResults = lapply(names(modulePeakList), function(mod) {
      peaks = modulePeakList[[mod]]
      peaks = intersect(peaks, colnames(atacOrig))
      if (length(peaks) == 0) return(NULL)
      
      atac = atacOrig[, peaks]
      
      MEvec = MEs[ids_rnaseq, paste0("ME", mod)]  # adjust if names differ
      A = atac[, peaks, drop = FALSE]
      
      r = apply(A, 2, function(x) cor(x, MEvec, use = "pairwise"))
      n = sum(!is.na(MEvec) & !is.na(A[,1]))
      p = 2 * pt(-abs(r * sqrt((n - 2)/(1 - r^2))), df = n - 2)
      
      out = data.frame(peak = peaks, module = mod, r = r, p = p)
      out$FDR = p.adjust(out$p, "BH")
      out
    })
    corResults = do.call(rbind, corResults)
    corResults$ctype = CTYPE
    EIGENCORREL[[CTYPE]] = corResults
    write.csv(corResults, file=file.path(ROOT, "outputs", paste0("misc_WGCNA_eigen_allcorrel_", CTYPE, ".csv")), row.names=F)
    
    ###
    # Calculate per-module summaries
    summaryByModule = aggregate(cbind(nSig = FDR < 0.05) ~ module, data = corResults, sum)
    summaryByModule$fracSig = summaryByModule$nSig / table(corResults$module)[summaryByModule$module]
    summaryByModule$all_GENE_count = table(moduleColors)[summaryByModule$module]
    summaryByModule$all_OCR_count = table(corResults$module)[summaryByModule$module]
    summaryByModule$ctype = CTYPE
    EIGENCORRELSUM[[CTYPE]] = summaryByModule
    
    ###
    # Permutation tests to calculate whether we have significantly more modules (in each cell type) that have more than one OCRs significantly associated with module eigenvector
    {
      n_perm = 10
      perm_summary = vector("list", n_perm)
      for (perm in 1:n_perm) {
        print(perm)
        perm_corResults = lapply(names(modulePeakList), function(mod) {
          peaks = modulePeakList[[mod]]
          peaks = intersect(peaks, colnames(atacOrig))
          if (length(peaks) == 0) return(NULL)
          
          atac = atacOrig[, peaks]
          # Permute eigengene
          perm_MEvec = sample(MEs[ids_rnaseq, paste0("ME", mod)])
          
          A = atac[, peaks, drop = FALSE]
          r = apply(A, 2, function(x) cor(x, perm_MEvec, use = "pairwise"))
          n = sum(!is.na(perm_MEvec) & !is.na(A[,1]))
          p = 2 * pt(-abs(r * sqrt((n - 2)/(1 - r^2))), df = n - 2)
          
          df = data.frame(module = mod, p = p)
          df$FDR = p.adjust(df$p, "BH")
          return(df)
        })
        perm_corResults = do.call(rbind, perm_corResults)
        
        perm_summary[[perm]] = aggregate(FDR < 0.05 ~ module, data = perm_corResults, FUN = sum)
      }
      
      # Count how many modules have ≥1 significant OCR in each perm
      null_sig_module_counts = sapply(perm_summary, function(df) sum(df[[2]] > 0))
      
      # Observed count
      observed_sig_module_count = sum(aggregate(FDR < 0.05 ~ module, data = corResults, FUN = sum)[[2]] > 0)
      
      # Empirical p-value
      empirical_p = mean(null_sig_module_counts >= observed_sig_module_count)
      cat("Empirical p-value for OCR–eigengene coupling:\n", empirical_p, "\n")
      
      df = data.frame(null_counts = null_sig_module_counts)
      permutPlot = ggplot(df, aes(x = null_counts)) + geom_histogram(binwidth = 1, fill = "gray70", color = "black") +
        geom_vline(xintercept = observed_sig_module_count, color = "red", linetype = "dashed", size = 1) +
        theme_classic() + labs(title = "Null distribution of regulatory module coupling",
                               x = "# of modules with ≥1 significant OCR (per permutation)",y = "Frequency",
                               subtitle = paste("Observed =", observed_sig_module_count, "; empirical p < ", empirical_p))
      pdf(file=file.path(ROOT, "outputs", paste0("misc_wgcna_permutPlot_", CTYPE, ".pdf")), width=4, height=2); print(permutPlot); dev.off()
    }
    write.csv(summaryByModule, file=file.path(ROOT, "outputs", paste0("misc_WGCNA_eigen_summary_", CTYPE, ".csv")), row.names=F)
    
    ###
    # Calculate enrichment of OCRs in each module and differential (SCZ-associated) OCRs
    modules = unique(corResults$module)
    dacPeaks = dacAnalysis$dacResults$dac[[paste0(CTYPE_ALT, ".SCZ_Control")]][(dacAnalysis$dacResults$dac[[paste0(CTYPE_ALT, ".SCZ_Control")]]$adj.P.Val<0.05),"PeakID"]
    enrichTable = do.call(rbind, lapply(modules, function(mod) {
      modPeaks = subset(corResults, module == mod & FDR < 0.05)$peak
      
      inMod = colnames(atacOrig) %in% modPeaks
      isDAC = colnames(atacOrig) %in% dacPeaks
      
      tbl = table(Module = inMod, DAC = isDAC)
      if (all(dim(tbl) == c(2,2))) {
        ft = fisher.test(tbl)
        
        data.frame(
          Module = mod,
          Peaks_in_Module = sum(inMod),
          Peaks_in_DAC = sum(isDAC),
          Overlap = tbl["TRUE","TRUE"],
          Odds_Ratio = ft$estimate,
          P_value = ft$p.value
        )
      }
    }))
    enrichTable$FDR = p.adjust(enrichTable$P_value, "BH")
    write.csv(enrichTable, file=file.path(ROOT, "outputs", paste0("Table_S11_", CTYPE, "_enrich.csv")), row.names=F)
    enrichTable$ctype = CTYPE
    DAC_ENRICH[[CTYPE]] = enrichTable
    
    ###
    # Calculate enrichment of genes in each module in SCZ GWAS (using MAGMA)
    magma = data.frame(fread(file.path(ROOT, "inputs", paste0("magma_", CTYPE , ".tsv.gz")), sep="\t", stringsAsFactors=F))
    magmaScz = magma[(magma$gwasAcronym == "sz3"),c("VARIABLE", "NGENES", "BETA", "BETA_STD", "SE", "P")]
    magmaScz$FDR = p.adjust(magmaScz$P, method="BH")
    magmaScz$ctype = CTYPE
    MAGMA_ENRICH[[CTYPE]] = magmaScz
    write.csv(magmaScz, file=file.path(ROOT, "outputs", paste0("Table_S11_", CTYPE, "_magma.csv")), row.names=F)
  }
  
  ###
  # Plotting
  {
    df = eigenVectorCorrelSum
    df$ctype = factor(df$ctype, levels = c("GABA","GLU","OLIG","MGAS"))
    df_out = df %>% group_by(ctype) %>%
      mutate(Q1 = quantile(all_GENE_count, 0.25), Q3 = quantile(all_GENE_count, 0.75),
             IQR = Q3 - Q1, outlier = (all_GENE_count < (Q1 - 1.5 * IQR)) | (all_GENE_count > (Q3 + 1.5 * IQR))
      ) %>% ungroup()
    
    # Median bars
    df$all_GENE_count = as.numeric(unlist(df$all_GENE_count, use.names = FALSE))
    df$all_OCR_count = as.numeric(unlist(df$all_OCR_count, use.names = FALSE))
    df$fracSig = as.numeric(unlist(df$fracSig, use.names = FALSE))
    df_stat = df %>% group_by(ctype) %>% summarise(median = median(all_GENE_count), .groups = "drop")
    df_stat_OCRs = df %>% group_by(ctype) %>% summarise(median = median(all_OCR_count), .groups = "drop")
    
    wgcnaPlot = ggplot(df_out, aes(x = ctype, y = all_GENE_count, fill = ctype)) + geom_violin(trim = FALSE, alpha = 0.9, width = 0.9, color = NA) +
      geom_boxplot(width = 0.18, outlier.shape = NA, color = "black", alpha = 0.35) + geom_point(data = subset(df_out, outlier), 
                                                                                                 aes(x = ctype, y = all_GENE_count), shape = 21, fill = "white", color = "black", size = 2, stroke = 0.4, position = position_jitter(width = 0.05),
                                                                                                 inherit.aes = FALSE) + scale_fill_npg() + scale_y_continuous(labels = scales::comma) +
      labs(x = NULL, y = "Number of genes per module") +
      theme_classic(base_size = 12) + theme(legend.position = "none", axis.text.x = element_text(face = "bold"), plot.title = element_text(face = "bold")) + coord_flip()
    wgcnaPlot
    pdf(file=file.path(outDir, "wgcna_modules_by_ctype.pdf"), width=4, height=2); print(wgcnaPlot); dev.off()
  }
  
  # Text for answer, part 1
  {
    print(paste0("> Modules per cell type: ", min(table(df$ctype)), " - ", max(table(df$ctype))))
    print("> Median number of genes by cell type")
    print(df_stat)
    print("> Median number of OCRs by cell type")
    print(df_stat_OCRs)
    print(paste0("> Modules per cell type: ", min(table(df$ctype)), " - ", max(table(df$ctype))))
    print(paste0("> Every module had approximately ", round(sum(df$all_OCR_count)/sum(df$all_GENE_count), 2), "x more OCRs than genes."))
  }
  
  # Correlation between modules
  {
    for(ctype in c("GABA", "GLU", "OLIG", "MGAS")) {
      MEcor = cor(MEs_LIST[[ctype]], use = "pairwise.complete.obs")
      hc = hclust(as.dist(1 - abs(MEcor)), method = "average")
      pdf(file=file.path("~/Desktop/", paste0("wgcna_modules_", ctype, ".pdf")), width=4, height=4) 
      heatmap(
        MEcor,
        Rowv = as.dendrogram(hc),
        Colv = as.dendrogram(hc),
        symm = TRUE,
        scale = "none",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        margins = c(8, 8),
        main = paste0("Eigengene correlation structure: ", CTYPE)
      )
      dev.off()
      
      MEcor = cor(MEs_LIST$MGAS)
      diag(MEcor) = NA
      mean(abs(MEcor), na.rm = TRUE)
      write.csv(file=file.path(outDir, paste0("wgcna_correl_modules_", ctype, ".csv")), MEcor)
      
      me = MEs_LIST[[ctype]]
      names(me) = gsub("^ME", "", names(me))
      
      pdf(file=file.path(outDir, paste0("wgcna_eigen_network_", ctype, ".pdf")), width=8, height=7) 
      WGCNA::plotEigengeneNetworks(
        me,
        setLabels = CTYPE,
        marDendro = c(3,3,2,4),
        marHeatmap = c(3,4,2,2),
        cex.lab = 0.8,
        xLabelsAngle = 90
      )
      dev.off()
    }
  }
  
  # Plot for reviewer, panel B ::: Intramodular coherence across cell types
  {
    hub_thresh = 0.6
    
    for(ctype in c("GABA", "GLU", "OLIG", "MGAS")) {
      kME = cor(t(rnaseqList[[ctype]]$residualized_DxBrainRegion_EffectKept[,rownames(MEs_LIST[[ctype]])]), MEs_LIST[[ctype]], use = "pairwise.complete.obs")
      kME_df = data.frame(
        gene = rownames(rnaseqList[[ctype]]$residualized_DxBrainRegion_EffectKept[,rownames(MEs_LIST[[ctype]])]),
        module = labels2colors(NETs_LIST[[ctype]]$colors)
      )
      
      kME_df$kME = mapply(function(g, m) {
        if(m == "grey") return(NA)
        kME[g, paste0("ME", m)]
      }, g = rownames(rnaseqList[[ctype]]$residualized_DxBrainRegion_EffectKept[,rownames(MEs_LIST[[ctype]])]), m = labels2colors(NETs_LIST[[ctype]]$colors))
      
      kmePlot = ggplot(kME_df, aes(x = kME)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "white") +
        facet_wrap(~ module, scales = "free_y") +
        theme_classic() +
        labs(x = "Module membership (kME)", y = "Gene count",
             title = paste("kME distributions:", CTYPE))
      pdf(file=file.path(outDir, paste0("wgcna_module_kME_full_", ctype, ".pdf")), width=8, height=6); print(kmePlot); dev.off()
      
      hub_summary = kME_df %>%
        filter(!is.na(kME)) %>%
        group_by(module) %>%
        summarise(
          frac_hubs = mean(kME > hub_thresh),
          .groups = "drop"
        )
      
      # Intramodular hub enrichment
      hubSumPlot = ggplot(hub_summary, aes(x = "", y = frac_hubs)) + geom_violin(fill = "#d73027", alpha = 0.6, width = 0.8) +
        geom_point(size = 2, position = position_jitter(width = 0.05)) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        coord_cartesian(ylim = c(0, 1)) + theme_classic() + labs(y = "Fraction of hub genes (kME > 0.6)",x = NULL)
      pdf(file=file.path(outDir, paste0("wgcna_module_kME_threshold_", ctype, ".pdf")), width=2, height=3); print(hubSumPlot); dev.off()
    }
  }
  
  ###
  # Plot for reviewer, panel C ::: Fraction of modules with ≥1 regulatory OCR
  {
    df = eigenVectorCorrelSum
    df = df %>% mutate(has_reg_OCR = nSig > 0)
    prop_nonzero = df %>% group_by(ctype) %>% summarise(
      n_modules = n(), n_regulated = sum(has_reg_OCR), frac_modules_regulated = mean(has_reg_OCR), .groups = "drop")
    
    propNonZeroPlot = ggplot(prop_nonzero, aes(x = ctype, y = frac_modules_regulated, fill = ctype)) +
      geom_col(width = 0.7) +
      scale_fill_npg() +
      scale_y_continuous(labels = percent_format()) +
      theme_classic() +
      labs(
        y = "Fraction of modules with ≥1 regulatory OCR",
        x = NULL,
        title = "Prevalence of chromatin–transcription coupling"
      )
    pdf(file=file.path(outDir, "wgcna_addAtac_propNonZero.pdf"), width=3, height=3); print(propNonZeroPlot); dev.off()
    
    ###
    # Plot for reviewer, panel D ::: Distribution of OCR burden per module
    signOcrPerModPlot = ggplot(df, aes(x = ctype, y = nSig, fill = ctype)) +
      geom_violin(trim = FALSE) + geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6) +
      scale_fill_npg() + scale_y_continuous(trans = "log10", labels = comma) +
      theme_classic() + labs(y = "Number of regulatory OCRs per module (log scale)", x = NULL, title = "Regulatory burden per module")
    pdf(file=file.path(outDir, "wgcna_addAtac_signOcrPerModPlot.pdf"), width=5, height=3); print(signOcrPerModPlot); dev.off()
    
    pdf(file=file.path(outDir, "wgcna_minorRank.pdf"), width=5, height=3)
    df %>% group_by(ctype) %>% arrange(desc(nSig)) %>% mutate(rank = row_number()) %>%
      ggplot(aes(x = rank, y = nSig, color = ctype)) + geom_line() + scale_color_npg() + scale_y_continuous(trans = "log10") +
      theme_classic() + labs(x = "Module rank (by regulatory OCRs)", y = "nSig (log scale)", title = "A minority of modules account for most regulatory signal")
    dev.off()
    
    ###
    # Strength of regulation: fraction per module (currently not used)
    strengthPlot = ggplot(df, aes(x = ctype, y = fracSig, fill = ctype)) + geom_violin(trim = FALSE) +
      geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5) + scale_fill_npg() +
      scale_y_continuous(labels = percent_format()) + theme_classic() +
      labs(y = "Fraction of regulatory OCRs per module", x = NULL, title = "Strength of chromatin–transcription coupling")
    pdf(file=file.path(outDir, "wgcna_addAtac_strength.pdf"), width=5, height=3); print(strengthPlot); dev.off()
    
    strengthPlot2 = ggplot(df, aes(x = ctype, y = fracSig)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
      theme_classic() +
      labs(
        x = "Cell type",
        y = "Fraction of eigengene-correlated OCRs",
        title = "Proportion of module-regulating OCRs by cell type"
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1))
    pdf(file=file.path(outDir, "wgcna_addAtac_strength_2.pdf"), width=5, height=3); print(strengthPlot2); dev.off()
    
    print("> Medians of module-level proportions of significantly associated OCRs:")
    df$fracSig = as.numeric(df$fra)
    medians = df %>%
      group_by(ctype) %>%
      summarize(med = median(fracSig, na.rm = TRUE))
    print(medians)
    
    strengthPlot3 = ggplot(df, aes(x = fracSig, fill = ctype, color = ctype)) +
      geom_density(alpha = 0.35, size = 1) + geom_vline(data = medians, aes(xintercept = med, color = ctype),
                                                        linetype = "solid", size = 1.1, show.legend = FALSE) + scale_fill_npg() + scale_color_npg() +
      theme_classic(base_size = 13) + labs(x = "Fraction of eigengene-correlated OCRs per module", y = "Density", fill = "Cell type",
                                           title = "Chromatin–transcription coupling across cell types", subtitle = "Fraction of regulatory elements whose accessibility tracks module eigengene"
      ) + theme(legend.position = "right", plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 10))
    pdf(file=file.path(outDir, "wgcna_addAtac_strength_3.pdf"), width=5, height=3); print(strengthPlot3); dev.off()
  }
  
  # Plot of FDR-sign enrichment for MAGMA (currently not shown)
  sigM = magmaMerged %>% filter(FDR < 0.05)
  magmaGlobalPlot = ggplot(sigM, aes(x = VARIABLE, y = -log10(FDR), fill = ctype)) +
    geom_col(position = "dodge") + coord_flip() + facet_wrap(~ ctype, scales = "free_y") +
    theme_classic() + labs(title = "Disease-enriched modules by cell type")
  pdf(file=file.path(outDir, "wgcna_magma_global_fdr.pdf"), width=5, height=3); print(magmaGlobalPlot); dev.off()
  
  ###
  # Plot for reviewer, panel F ::: OLIG GSEA
  #df = read.csv("/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/publication_plots/WGCNA_gsea_OLIG/universalBackground/result_textFiles/msigdbSetsVeryStrictlyPruned.tsv.gz", sep="\t")
  df = read.csv("/sc/arion/projects/CommonMind/roussp01a/MOLECULAR_PROFILING/publication_plots/WGCNA_gsea_OLIG/universalBackground/result_textFiles/msigdbSetsStrictlyPruned.tsv.gz", sep="\t")
  df[(df$Set == "darkgrey") & (df$BH_AdjP<0.05),]
  
  
  df_plot = df %>%
    filter(Set == "darkgrey", BH_AdjP < 0.05) %>%
    mutate(minusLogP = -log10(pval)) %>%
    arrange(desc(minusLogP))
  
  df_plot$Pathway = gsub("^Gobp ", "", sapply(df_plot$name_full, function(x) strsplit(split="\\(", x)[[1]][1]))
  gseaOligPlot = ggplot(df_plot, aes(x = reorder(Pathway, minusLogP), y = minusLogP)) + geom_col(width = 0.75) +
    coord_flip() + scale_fill_npg() + labs(x = NULL, y = expression(-log[10]("P-value")), title = "Top enriched pathways — darkgrey module") +
    theme_classic(base_size = 12) + theme(axis.text.y = element_text(size = 10), plot.title = element_text(face = "bold"))
  
  pdf(file=file.path(outDir, "wgcna_gsea_olig.pdf"), width=10, height=4); print(gseaOligPlot); dev.off()
  
  ###
  # Table exports (not necessarily papers' tables)
  {
    dacMerged = do.call("rbind.data.frame", DAC_ENRICH)
    write.csv(dacMerged, file=file.path(outDir, "Table_WGCNA_DAC_all.csv"), row.names=F)
    
    magmaMerged = do.call("rbind.data.frame", MAGMA_ENRICH)
    table(p.adjust(magmaMerged$P < 0.05))
    table(magmaMerged$FDR < 0.05)
    magmaMerged$FDR_new = p.adjust(magmaMerged$P)
    write.csv(magmaMerged, file=file.path(outDir, "Table_WGCNA_MAGMA_all.csv"), row.names=F)
    
    eigenVectorCorrel = do.call("rbind.data.frame", EIGENCORREL)
    write.csv(eigenVectorCorrel, file=file.path(outDir, "Table_eigenVectorCorrel.csv"), row.names=F)
    
    eigenVectorCorrelSum = do.call("rbind.data.frame", EIGENCORRELSUM)
    write.csv(eigenVectorCorrelSum, file=file.path(outDir, "Table_eigenVectorCorrelSum.csv"), row.names=F)
    
    # Tables for export
    gene_modules_df = do.call("rbind.data.frame", lapply(names(NETs_LIST), function(ctype) {
      cbind.data.frame("Cell_type"=ctype, "Module_name"=labels2colors(NETs_LIST[[ctype]]$colors), "Ensembl_Gene_ID"=names(NETs_LIST[[ctype]]$colors))
    }))
    write.csv(gene_modules_df, file=file.path(outDir, "Table_eigenVectorCorrelSum.csv"), row.names=F)
    
    eigenVectorCorrelExport = eigenVectorCorrel[(eigenVectorCorrel$FDR<0.05),c("ctype", "module", "peak", "r", "FDR")]
    eigenVectorCorrelExport$r = round(eigenVectorCorrelExport$r, 3)
    eigenVectorCorrelExport$p = round(eigenVectorCorrelExport$r, 3)
    eigenVectorCorrelExport$FDR = round(eigenVectorCorrelExport$FDR, 3)
    colnames(eigenVectorCorrelExport) = c("Cell_type", "Module_name", "OCR ID", "R", "P-value", "FDR")
    write.csv(eigenVectorCorrelExport, file=file.path(outDir, "Table_coexpress_OCR_regulating.csv"), row.names=F)
    
    dacMerged$comb = paste0(dacMerged$ctype, "_", dacMerged$Module)
    colnames(dacMerged) = paste0("dacEnrich_", colnames(dacMerged))
    magmaMerged$comb = paste0(magmaMerged$ctype, "_", magmaMerged$VARIABLE)
    colnames(magmaMerged) = paste0("magmaEnrich_", colnames(magmaMerged))
    
    outputDf = cbind.data.frame(magmaMerged, dacMerged[match(magmaMerged$magmaEnrich_comb, dacMerged$dacEnrich_comb),])
    write.csv(file=file.path(outDir, "Table_coexpress_dacEnrich.csv"), row.names=F, outputDf)
    write.csv(gene_modules_df, file=file.path(outDir, "Table_eigenVectorCorrelSum.csv"), row.names=F)
  }
  
  #save(list = ls(all = T),file=file.path(outDir, "gene_module_analysis.RData"),envir=environment())
  
}
