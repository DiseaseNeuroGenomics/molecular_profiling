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

# Load qPCR raw input data
df = read.csv(file.path(ROOT, "inputs", "qPCR_raw_input.csv"))  
df$combo = paste0(df$Donor_ID, "_", df$Transcript_ID)           # Helper column (unique combination of Sample/Donor - Transcript/Probe)

# Load differential analysis results from RNA-seq transcript-level analysis for those 8 qPCR-validates transcripts
det_orig = read.csv(file.path(ROOT, "inputs", "qPCR_DET.csv")) 
det_orig$SE =  det_orig$logFC / det_orig$t  # Calculate standard error
det_orig$method = "RNA-seq"                 # Add method name for compatibility (and discrimination from qPCR in plots)

# Load differential analysis results from RNA-seq transcript-level analysis for those 8 qPCR-validates transcripts
exprMx_orig = read.csv(file.path(ROOT, "inputs", "qPCR_exprMx.csv"))
rownames(exprMx_orig) = exprMx_orig[,1]
exprMx_orig = exprMx_orig[,2:ncol(exprMx_orig)]

# Shrink input qPCR measurement - keep only one value (median CT) for each triplicate
df_summary <- df %>%
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

# Merge qPCR and DET data
df = (rbind.data.frame(fc_results[,c("Transcript_ID", "logFC", "method", "SE")], det_orig[,c("Transcript_ID", "logFC", "method", "SE")]))

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

df$Transcript_ID = ordered(df$Transcript_ID, levels=c("ENST00000496818", "ENST00000465278", "ENST00000483136", "ENST00000492150", "ENST00000437508", "ENST00000502281", "ENST00000460908", "ENST00000338700"))
df$method = ordered(df$method, levels=c("RNA-seq", "qPCR"))
logfc_plot = ggplot(df, aes(x = Transcript_ID, y = logFC, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
  labs(
    title = "Comparison of log2 Fold Change (logFC)",
    x = "Target",
    y = "log2 Fold Change",
    fill = "Method"
  ) +
  geom_errorbar(aes(ymin = logFC - SE, ymax = logFC + SE), width = 0.3) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +   scale_fill_manual(values = pal_nejm("default")(8)[c(4,3)])
logfc_plot

pdf(file=file.path(ROOT, "outputs", "Fig_S10d.pdf"), width=10, height=6); print(logfc_plot); dev.off()
