

load("../results/DEA_list.RData")
load("../results/DEA_pneu_list.RData")

################################################################################
limma_res_list <- c(DEA_list[["limma_DEA"]], DEA_pneu_list[["limma_DEA"]])



limma_res_list_df <- list()

for (comp_i in names(limma_res_list) ) {
    # comp_i <- "nCoV_Heal"
    res_df_i <- limma_res_list[[comp_i]][,1:2]
    colnames(res_df_i) <- c("symbol", comp_i)
    
    limma_res_list_df[[comp_i]] <- res_df_i
    
}


comp_logFC <- Reduce(dplyr::full_join, limma_res_list_df) %>% as.data.frame()



nCoV_pneu_comp_logFC <- comp_logFC[,c("symbol", "nCoV_Heal", "pneu_Heal", "pneuVir_Heal", "pneuBac_Heal")]


readr::write_tsv(nCoV_pneu_comp_logFC, path="../results/nCoV_pneu_comp_logFC.txt")


save(nCoV_pneu_comp_logFC, file="../results/nCoV_pneu_comp_logFC.RData")

# add annotation
source("../../../SynologyDrive/zzz.R")
DEGs_anno <- gene2KEGG( nCoV_pneu_comp_logFC$symbol )
nCoV_pneu_comp_logFC_anno <- dplyr::left_join(nCoV_pneu_comp_logFC, DEGs_anno, by=c("symbol"="SYMBOL") )
readr::write_tsv(nCoV_pneu_comp_logFC_anno, path="../results/nCoV_pneu_comp_logFC_anno.txt")




