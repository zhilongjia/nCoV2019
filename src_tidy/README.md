
# merge data and gene symbols
./0.1_preprocessing.R  

# Differential Expression Analysis for three comparsion pairs
./1.2_DE_pneu.R

# norm all samples for PCA/Heatmap only. 
./1.1_norm3groups.R  

# Pathway analysis
./2.0_KEGG.R

# cogena analysis
./2.1_cogena.R


# heatmap and venn
./3.3_DEG_heatmap.R 
./3.1_drug_pathway.R
