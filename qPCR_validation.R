library(tidyverse)    # includes ggplot2, dplyr, tidyr, etc.
library(ggpubr)       # for easy p-value annotation on plots
library(ggsci)
library(ggplot2)
library(dplyr)
library(patchwork)

################################################################################
##### CONFIG ###################################################################

{
  ROOT = "~/Desktop/molecular_profiling/molecular_profiling/" # !!! FIXME: SET TO YOUR CUSTOM DIRECTORY !!
  dir.create(file.path(ROOT, "outputs"))
}

################################################################################
##### DATA PREPARATION #########################################################

{
  # Load qPCR raw input data
  df_raw = read.csv(file.path(ROOT, "inputs", "qPCR_raw_input.csv"))  
  df_raw$combo = paste0(df_raw$Donor_ID, "_", df_raw$Transcript_ID)   # Helper column (unique combination of Sample/Donor - Transcript/Probe)
  
  # Load differential analysis results from RNA-seq transcript-level analysis for those 8 qPCR-validates transcripts
  det_orig = read.csv(file.path(ROOT, "inputs", "qPCR_DET.csv")) 
  det_orig$SE =  det_orig$logFC / det_orig$t  # Calculate standard error
  det_orig$method = "RNA-seq"                 # Add method name for compatibility (and discrimination from qPCR in plots)
  
  # Load differential analysis results from RNA-seq transcript-level analysis for those 8 qPCR-validates transcripts
  exprMx_orig = read.csv(file.path(ROOT, "inputs", "qPCR_exprMx.csv"))
  rownames(exprMx_orig) = exprMx_orig[,1]
  exprMx_orig = exprMx_orig[,2:ncol(exprMx_orig)]
  
  # Shrink input qPCR measurement - keep only one value (median CT) for each triplicate
  df_summary <- df_raw %>%
    group_by(combo) %>%
    summarise(
      mean_CT = median(CT, na.rm = TRUE),
      sd_CT = sd(CT, na.rm = TRUE),
      Donor_ID = first(Donor_ID),
      Transcript_ID = first(Transcript_ID),
      Pair = first(Pair),
      Condition = first(Condition),
      Sex = first(Sex),
      RIN = first(RIN),
      pH = first(pH),
      Age.of.death = first(Age.of.death),
      .groups = "drop"
    )
  
  # Separate qPCR measurements for targets and reference 
  targets = data.frame(df_summary[(df_summary$Transcript_ID != "GAPDH"),])
  gapdh = data.frame(df_summary[(df_summary$Transcript_ID == "GAPDH"),])
  
  # ΔCT calculation
  targets$del_ct = sapply(1:nrow(targets), function(i) {
    val = (targets[i,"mean_CT"] - gapdh[which(gapdh$Donor_ID == targets[i,"Donor_ID"]),"mean_CT"])
    ifelse(length(val) == 0, NA, val)
  })
  targets = targets[!is.na(targets$del_ct),]
  
  # LogFC calculation
  fc_results <- targets %>%
    group_by(Transcript_ID, Condition) %>%
    summarise(
      mean_deltaCT = mean(del_ct, na.rm = TRUE),
      se_deltaCT = sd(del_ct, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Condition,
      values_from = c(mean_deltaCT, se_deltaCT),
      names_sep = "_"
    ) %>%
    mutate(
      logFC = mean_deltaCT_Control - mean_deltaCT_SCZ,
      fold_change = 2^(logFC),
      SE = sqrt(se_deltaCT_Control^2 + se_deltaCT_SCZ^2)
    )
  fc_results$method = "qPCR"
  
  # ---------------------------------------------------------------------------
  # NEW: Paired t-tests at the qPCR level (one per transcript)
  # ---------------------------------------------------------------------------
  # Pivot ΔCt to wide format so each row is a matched pair (SCZ vs. Control)
  targets_wide <- targets %>%
    select(Pair, Transcript_ID, Condition, del_ct) %>%
    pivot_wider(names_from = Condition, values_from = del_ct)
  
  # Run paired t-test for each transcript and extract p-value + t-statistic
  # (compute t.test once per group using group_map to avoid calling it twice in summarise)
  qpcr_stats <- targets_wide %>%
    group_by(Transcript_ID) %>%
    group_modify(~ {
      scz  <- .x$SCZ[!is.na(.x$SCZ) & !is.na(.x$Control)]
      ctrl <- .x$Control[!is.na(.x$SCZ) & !is.na(.x$Control)]
      n    <- length(scz)
      tt   <- tryCatch(t.test(scz, ctrl, paired = TRUE), error = function(e) NULL)
      data.frame(
        n_pairs = n,
        t_stat  = if (!is.null(tt)) as.numeric(tt$statistic) else NA_real_,
        p_value = if (!is.null(tt)) tt$p.value               else NA_real_
      )
    }) %>%
    ungroup() %>%
    # FDR correction across the 8 transcripts tested by qPCR
    mutate(adj_p_value = p.adjust(p_value, method = "BH"))
  
  # Add significance stars (used for plot annotations)
  sig_stars <- function(p) {
    case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  }
  qpcr_stats$sig_label     <- sig_stars(qpcr_stats$p_value)
  qpcr_stats$sig_adj_label <- sig_stars(qpcr_stats$adj_p_value)
  
  # Attach qPCR stats to fc_results
  fc_results <- fc_results %>% left_join(qpcr_stats, by = "Transcript_ID")
  
  # ---------------------------------------------------------------------------
  # Build combined dataframe for Panel D (logFC + SE from both methods)
  # ---------------------------------------------------------------------------
  df = rbind.data.frame(
    fc_results[, c("Transcript_ID", "logFC", "method", "SE")],
    det_orig[,   c("Transcript_ID", "logFC", "method", "SE")]
  )
  
  # ---------------------------------------------------------------------------
  # NEW: Build a tidy significance comparison table (for new Panel E)
  # Columns: Transcript_ID | method | p_value | adj_p_value | sig_label
  # ---------------------------------------------------------------------------
  rnaseq_sig <- det_orig %>%
    select(Transcript_ID, p_value = P.Value, adj_p_value = adj.P.Val) %>%
    mutate(
      method    = "RNA-seq",
      sig_label = sig_stars(p_value),
      sig_adj_label = sig_stars(adj_p_value)
    )
  
  qpcr_sig <- qpcr_stats %>%
    select(Transcript_ID, p_value, adj_p_value) %>%
    mutate(
      method    = "qPCR",
      sig_label = sig_stars(p_value),
      sig_adj_label = sig_stars(adj_p_value)
    )
  
  sig_combined <- bind_rows(rnaseq_sig, qpcr_sig)
}

################################################################################
##### PANEL "A" :: Demographic and technical characteristics of 15 matched     #
################## SCZ–control pairs ###########################################

{
  targets$Condition = ordered(targets$Condition, levels=c("SCZ", "Control"))
  ph_plot = ggplot(targets[!duplicated(targets$Donor_ID),], aes(x = Condition, y = pH, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = "pH Distribution by Condition", x = "", y = "pH") + theme_minimal() + scale_fill_nejm()  
  
  rin_plot = ggplot(targets[!duplicated(targets$Donor_ID),], aes(x = Condition, y = RIN, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = "RIN Distribution by Condition", x = "", y = "RIN") + theme_minimal() + scale_fill_nejm()
  
  age_plot = ggplot(targets[!duplicated(targets$Donor_ID),], aes(x = Age.of.death, color = Condition, fill = Condition)) + geom_density(alpha = 0.3) + 
    facet_wrap(~ Sex) + labs(title = "Age.of.death Distribution by Condition and Sex", x = "Age at Death", y = "Density") +
    theme_minimal() + scale_fill_nejm() + scale_color_nejm()
  
  combined_plot = ph_plot + rin_plot + age_plot + plot_layout(ncol = 3)
  pdf(file=file.path(ROOT, "outputs", "Fig_S10a.pdf"), width=10, height=4); print(combined_plot); dev.off()
}

################################################################################
##### PANEL "B" :: Summary of selected transcripts for qPCR experiment #########

# This panel was manually created in graphic editor by editing (adding arrows and comments) to Fig. 5B

################################################################################
##### PANEL "C" :: Concordance between RNA-seq - qPCR measurements of overall ##
################## transcript abundance ########################################

{
  donors_scz = paste0("Sample_", unique(targets[(targets$Condition == "SCZ"), "Donor_ID"]), ".Olig.RNA")
  donors_ctrl = paste0("Sample_", unique(targets[(targets$Condition == "Control"), "Donor_ID"]), ".Olig.RNA")
  cacna1c_transcripts = c("ENST00000483136", "ENST00000496818", "ENST00000492150", "ENST00000465278")
  trim2_transcripts = c("ENST00000437508", "ENST00000502281", "ENST00000460908", "ENST00000338700")
  
  # RNA-seq :: Generate dataframe, each row one combination of Donor_ID & Transcript_ID & Expression
  rnaseq_expr_long <- as.data.frame(exprMx_orig) %>%
    tibble::rownames_to_column("Transcript_ID") %>%
    pivot_longer(-Transcript_ID, names_to = "Donor", values_to = "Expression") %>%
    mutate(Group = case_when(
      Donor %in% donors_scz ~ "SCZ",
      Donor %in% donors_ctrl ~ "CTRL",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Group)) %>%
    mutate(Gene = case_when(
      Transcript_ID %in% cacna1c_transcripts ~ "CACNA1C",
      Transcript_ID %in% trim2_transcripts ~ "TRIM2"
    ))
  
  # RNA-seq :: Shrink donor+transcript expression to just diagnosis+transcript expression (averaged)
  rnaseq_expr_summary <- rnaseq_expr_long %>%
    group_by(Gene, Transcript_ID) %>%
    summarise(
      mean_expr = mean(Expression, na.rm = TRUE),
      se_expr = sd(Expression, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # qPCR :: Generate dataframe, each row diagnosis+transcript expression (averaged)
  qPCR_expr_summary = targets %>%
    group_by(Transcript_ID) %>%
    summarise(
      mean_del_ct = mean(del_ct, na.rm = TRUE),
      sd_del_ct = sd(del_ct, na.rm = TRUE),
      n = n()
    )
  
  # Top plot: RNA-seq
  rnaseq_expr_summary$Transcript_ID = ordered(rnaseq_expr_summary$Transcript_ID, levels=c("ENST00000496818", "ENST00000465278", "ENST00000483136", "ENST00000492150", "ENST00000437508", "ENST00000502281", "ENST00000460908", "ENST00000338700"))
  rnaseq_avgExprPlot = ggplot(rnaseq_expr_summary, aes(x = Transcript_ID, y = mean_expr)) +
    geom_col(fill = pal_nejm("default")(8)[4]) +
    geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr), width = 0.2) +
    labs(y = "logCPM (RNA-seq)", x = NULL) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_fill_manual(values = pal_nejm("default")(8)[3:8])
  
  # Bottom plot: qPCR
  qPCR_expr_summary$Transcript_ID = ordered(qPCR_expr_summary$Transcript_ID, levels=c("ENST00000496818", "ENST00000465278", "ENST00000483136", "ENST00000492150", "ENST00000437508", "ENST00000502281", "ENST00000460908", "ENST00000338700"))
  qPCR_avgExprPlot = ggplot(qPCR_expr_summary, aes(x = Transcript_ID, y = -mean_del_ct)) +
    geom_col(fill = pal_nejm("default")(8)[3]) +
    geom_errorbar(aes(ymin = -mean_del_ct - sd_del_ct, ymax = -mean_del_ct + sd_del_ct), width = 0.2) +
    labs(y = "delta(CT) [qPCR]", x = "Transcript Target") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_fill_manual(values = pal_nejm("default")(8)[4:8])
  qPCR_avgExprPlot
  # Combine plots
  finalPlot = rnaseq_avgExprPlot / qPCR_avgExprPlot + plot_layout(heights = c(1, 1))
  pdf(file=file.path(ROOT, "outputs", "Fig_S10c.pdf"), width=10, height=6); print(finalPlot); dev.off()
}

################################################################################
##### PANEL "D" :: Comparison of SCZ versus control transcript-level           #  
################## differences measured by RNA-seq and qPCR ####################
# FIXED: added position = position_dodge(width = 0.7) to geom_errorbar so that
#        error bars are correctly aligned with their respective dodged bars.
#        Without this, all error bars sat at the group centre x-position,
#        making the 3rd and 4th transcript bars appear mislabelled.
# NEW:   added significance stars above each qPCR bar (raw p-value from
#        paired t-test). Stars are positioned just above the top of each bar.

{
  # Transcript display order
  tx_order = c("ENST00000496818", "ENST00000465278", "ENST00000483136", "ENST00000492150",
               "ENST00000437508", "ENST00000502281", "ENST00000460908", "ENST00000338700")
  
  df$Transcript_ID = ordered(df$Transcript_ID, levels = tx_order)
  df$method        = ordered(df$method, levels = c("RNA-seq", "qPCR"))
  
  # Prepare star annotation data for qPCR bars only
  # Position the star just above the top of the error bar
  qpcr_annotations <- fc_results %>%
    select(Transcript_ID, logFC, SE, sig_label) %>%
    mutate(
      Transcript_ID = ordered(Transcript_ID, levels = tx_order),
      method        = ordered("qPCR", levels = c("RNA-seq", "qPCR")),
      # Place label above the bar + 1 SE, with a small extra gap
      label_y       = ifelse(logFC >= 0, logFC + SE + 0.07, logFC - SE - 0.07)
    )
  
  dodge_width = 0.7   # keep consistent with bar width
  
  logfc_plot = ggplot(df, aes(x = Transcript_ID, y = logFC, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = dodge_width), width = 0.6) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
    
    # FIX: position_dodge added so error bars align with their bars
    geom_errorbar(
      aes(ymin = logFC - SE, ymax = logFC + SE),
      position = position_dodge(width = dodge_width),
      width = 0.25
    ) +
    
    # NEW: significance stars above qPCR bars
    geom_text(
      data    = qpcr_annotations,
      aes(x   = Transcript_ID, y = label_y, label = sig_label, group = method),
      position = position_dodge(width = dodge_width),
      size    = 4,
      vjust   = 0,
      inherit.aes = FALSE
    ) +
    
    labs(
      title    = "Comparison of log2 Fold Change (logFC) with qPCR significance",
      subtitle = "Stars above qPCR bars: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant (paired t-test)",
      x        = "Target",
      y        = "log2 Fold Change",
      fill     = "Method"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = pal_nejm("default")(8)[c(4, 3)])
  
  logfc_plot
  pdf(file = file.path(ROOT, "outputs", "Fig_S10d.pdf"), width = 10, height = 6)
  print(logfc_plot)
  dev.off()
}

################################################################################
##### PANEL "E" (NEW) :: Side-by-side significance comparison                  #
#####                    RNA-seq vs qPCR p-values per transcript                #
################################################################################
# This panel directly addresses the reviewer's request to compare significance
# levels from both methods. It shows -log10(p) for both methods as a dot plot
# with the significance threshold lines marked, allowing easy visual assessment
# of concordance in statistical support.

{
  tx_order = c("ENST00000496818", "ENST00000465278", "ENST00000483136", "ENST00000492150",
               "ENST00000437508", "ENST00000502281", "ENST00000460908", "ENST00000338700")
  
  sig_combined$Transcript_ID = ordered(sig_combined$Transcript_ID, levels = tx_order)
  sig_combined$method        = ordered(sig_combined$method, levels = c("RNA-seq", "qPCR"))
  sig_combined$neg_log10_p   = -log10(sig_combined$p_value)
  
  # ----- Panel E1: dot plot of -log10(p) for both methods --------------------
  sig_dot_plot = ggplot(sig_combined, aes(x = Transcript_ID, y = neg_log10_p,
                                          color = method, shape = method)) +
    geom_point(size = 4, alpha = 0.9) +
    # p < 0.05 threshold line
    geom_hline(yintercept = -log10(0.05),  linetype = "dashed", color = "gray40", linewidth = 0.5) +
    # p < 0.01 threshold line
    geom_hline(yintercept = -log10(0.01),  linetype = "dotted", color = "gray60", linewidth = 0.5) +
    annotate("text", x = 0.6, y = -log10(0.05) + 0.08, label = "p = 0.05", size = 3, hjust = 0, color = "gray40") +
    annotate("text", x = 0.6, y = -log10(0.01) + 0.08, label = "p = 0.01", size = 3, hjust = 0, color = "gray60") +
    labs(
      title  = "Significance comparison: RNA-seq vs. qPCR",
      x      = "Transcript",
      y      = expression(-log[10](p-value)),
      color  = "Method",
      shape  = "Method"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(values = pal_nejm("default")(8)[c(4, 3)])
  
  # ----- Panel E2: table-style heatmap of significance labels ----------------
  # This gives a compact at-a-glance concordance view (as a tile matrix)
  sig_tile_plot = ggplot(sig_combined, aes(x = Transcript_ID, y = method, fill = sig_label)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sig_label), size = 4, fontface = "bold") +
    scale_fill_manual(
      values = c("***" = "#2166AC", "**" = "#74ADD1", "*" = "#ABD9E9", "ns" = "#F0F0F0"),
      name   = "Significance"
    ) +
    labs(
      title = "Significance concordance (nominal p-value)",
      x     = "Transcript",
      y     = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid  = element_blank()
    )
  
  # Combine into one figure
  sig_fig = sig_dot_plot / sig_tile_plot + plot_layout(heights = c(2, 1))
  
  pdf(file = file.path(ROOT, "outputs", "Fig_S10e.pdf"), width = 10, height = 8)
  print(sig_fig)
  dev.off()
  
  # Also print the qPCR stats table for inspection / Table S16 supplement
  qpcr_stats_export <- qpcr_stats %>%
    left_join(det_orig %>% select(Transcript_ID, rnaseq_P = P.Value, rnaseq_adjP = adj.P.Val),
              by = "Transcript_ID") %>%
    select(Transcript_ID, n_pairs, t_stat, qpcr_P = p_value, qpcr_adjP = adj_p_value,
           qpcr_sig = sig_label, rnaseq_P, rnaseq_adjP)
  
  write.csv(qpcr_stats_export,
            file = file.path(ROOT, "outputs", "qPCR_significance_summary.csv"),
            row.names = FALSE)
  
  print(qpcr_stats_export)
}
