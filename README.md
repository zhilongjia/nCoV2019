# nCoV2019

## Transcriptonal Analysis of 2019-nCoV

## todo
    filter other genes except protein-coding genes
    CMap/LINCS-based Drug repositioning analysis
    Key Pathways
    compare HIV drug genes and pathways with nCoV's

## Files: 

data: 

    final: raw read counts

documents: 
    
    slides

src/

    # mapping
    ./0.0.1_mapping_cmd.sh
    
    # merge data and gene symbols
    ./0.1_preprocessing.R  

    # Differential Expression Analysis for three comparsion pairs
    ./1.0_DE.R 
    ./1.0_DE_q0.001.R

    # norm all samples for PCA/Heatmap only. 
    ./1.1_norm3groups.R  

    # Pathway analysis
    ./2.0_KEGG.R

    # cogena analysis
    ./2.1_cogena.R
    ./2.3_cmap_lincs.R

# Citation

Zhou, Z., Ren, L., Zhang, L., Zhong, J., Xiao, Y., Jia, Z., Guo, L., Yang, J., Wang, C., Jiang, S. and Yang, D., 2020. Heightened innate immune responses in the respiratory tract of COVID-19 patients. Cell Host & Microbe. https://doi.org/10.1016/j.chom.2020.04.017

# log
Created: Jan 28, 2020
