# MolecularProfiling

This relates to the publication: **The cell type specific molecular landscape of schizophrenia**

doi: TBA

### main.R

This script when run as:

```bash
Rscript main.R
```

produces the following output files:

- outputs/Fig_2_a.pdf     # Numbers of differentially accessible OCRs (BH-adjusted P-value < 0.05) stratified by cell type and direction of change.
- outputs/Fig_4_a.pdf     # Numbers of differentially expressed genes (DEGs) per cell type.
- outputs/Fig_5_a.pdf.pdf # Numbers of dysregulated transcripts as well as genes with at least one differentially expressed transcript (FDR < 0.05) stratified by cell type. 
- outputs/Fig_5_d.pdf.pdf # Comparison of ENST00000465278 and ENST00000483136 expression in the study profiling OPC and mature oligodendrocytes (MO) in infants and adults.
- outputs/Fig_5_f.pdf.pdf # Comparison of ENST00000338700 and ENST00000460908 expression in the study profiling OPC and MO in infants and adults.
- outputs/Fig_S1_*.pdf    # Quality control for RNA-seq and ATAC-seq data. 
- outputs/Fig_S3_coleman__[celltype1]__[celltype2].pdf # Correlation of log2(cpm+1) counts between our RNA-seq data and external RNA-seq dataset of 3 cell types from parahippocampal gyrus.
- outputs/Fig_S3_hauberg__[celltype1]__[celltype2].pdf # Correlation of log2(cpm+1) counts between our ATAC-seq data and external ATAC-seq data from 4 cell types from the prefrontal cortex.

### qtl_analysis

Folder contains scripts for prepping files for qtl detection with MMQTL.

- 1_peer_factors.R: generates PEER factors using expression matrix
- 2_prep_mmqtl_files.R: generates residualized matrix using PEER factors and bed file containing gene information
- 3_generate_mmqtl_run_script.sh: generates scritps for running MMQTL per chromosome

