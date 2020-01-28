
# merge data and gene symbols
./0.1_preprocessing.R  

# Differential Expression Analysis for three comparsion pairs
./1.0_DE.R 

# norm all samples for PCA/Heatmap only. 
./1.1_norm3groups.R  

# Pathway analysis
./2.0_KEGG.R

# cogena analysis
./2.1_cogena.R
