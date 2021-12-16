# hackathon_TSG_MEN1
This repository investigates the correlation of expression of the MEN1 gene with a curated list of DNA repair and remodelling genes. The curated list is saved in "DDR_epi_genes_names.rds".

Raw counts of gene expression for different tissues were extracted from the GTEx and TCGA websites. For normalisation, the MOR/MRN method was used as described by Perez et al (2021). Correlation analysis was performed for every gene w.r.t. MEN1 and the top scoring 10 % were compared across different tissues. Finally, genes with high correlation across different tissues were compared in health (GTEx) and disease (TCGA).
