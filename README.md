# hackathon_TSG_MEN1
This repository investigates the correlation of expression of the MEN1 gene with a curated list of DNA repair and remodelling genes. The curated list is saved in "DDR_epi_genes_names.rds".

Raw counts of gene expression for various different cancers were extracted from the TCGA website using the scripts provided in the folder 'read_extraction'. The resulting raw counts for the curated list of genes from cancer tissues (TCGA) can be found in the folder 'cancer-data-all-files'. The raw counts from healthy tissues (GTEx) were processed locally using the script provided ("x") and not uploaded to the GitHub.

For normalisation of the raw counts, the MOR/MRN method was used as described by Perez et al (2021). Spearman correlation analysis was performed for every gene w.r.t. Men1 and the top scoring 10 % were compared across different tissues. Finally, genes with high correlation across different tissues were compared in health (GTEx) and disease (TCGA). The according scripts and resulting candidate genes can be found in folder 'correlation_analysis'. 

As an extension, the absolute difference of correlation with Men1 expression in healthy vs cancer tissue was calculated and visualized as a heatmap ("correlation_analysis/visualiser.py") for five different tissues.
