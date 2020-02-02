

load("../results/DEA_list.RData")
load("../results/DEA_pneu_list.RData")

################################################################################
limma_res_list <- c(DEA_list[["limma_res"]], DEA_pneu_list[["limma_res"]])


limma_res_list_df <- list()

for (comp_i in names(limma_res_list) ) {
    # comp_i <- "nCoV_Heal"
    res_df_i <- dplyr::filter(limma_res_list[[comp_i]], type=="protein_coding")[,1:2]
    
    colnames(res_df_i) <- c("symbol", comp_i)
    limma_res_list_df[[comp_i]] <- res_df_i
    
}


Allgene_logFC <- Reduce(dplyr::full_join, limma_res_list_df) %>% as.data.frame()


################################################################################
# abs(nCoV_Heal)>=6 | abs(pneu_Heal)>=6
DEG_logFC <- dplyr::select(Allgene_logFC, symbol, nCoV_Heal, pneu_Heal, pneuVir_Heal, pneuBac_Heal) %>% 
    dplyr::filter(abs(nCoV_Heal)>=6 | abs(pneu_Heal)>=6 )
nCoV_pneu_logFC_df <- as.data.frame(DEG_logFC[,-1])
rownames(nCoV_pneu_logFC_df) <- DEG_logFC$symbol

gplots::heatmap.2(as.matrix(nCoV_pneu_logFC_df), trace = "none", 
                  col="bluered", scale="none", 
                  srtCol=50, keysize=1, Colv=T,
                  ylab=paste(nrow(nCoV_pneu_logFC_df), " Genes"), 
                  margins = c(10, 7))
################################################################################
# abs(nCoV_Heal)>=4 & abs(nCoV_pneu)>=4
DEG_logFC <- dplyr::select(Allgene_logFC, symbol, nCoV_Heal, pneu_Heal, pneuVir_Heal, nCoV_pneu) %>% 
    dplyr::filter(abs(nCoV_Heal)>=4 & abs(nCoV_pneu)>=4 ) %>% 
    dplyr::arrange(desc(nCoV_Heal))
nCoV_pneu_logFC_df <- as.data.frame(DEG_logFC[,-1])
rownames(nCoV_pneu_logFC_df) <- DEG_logFC$symbol


gplots::heatmap.2(as.matrix(nCoV_pneu_logFC_df), trace = "none", 
                  col="bluered", scale="none", Rowv=F,
                  srtCol=50, keysize=1, Colv=T,
                  ylab=paste(nrow(nCoV_pneu_logFC_df), " Genes"), 
                  margins = c(10, 7))

## abs(nCoV_Heal)>=5 & abs(nCoV_pneuVir)>=5
DEG_logFC <- dplyr::select(Allgene_logFC, symbol, nCoV_Heal, pneu_Heal, pneuVir_Heal, nCoV_pneuVir) %>% 
    dplyr::filter(abs(nCoV_Heal)>=4 & abs(nCoV_pneuVir)>=4 ) %>% 
    dplyr::arrange(desc(nCoV_Heal))

nCoV_pneu_logFC_df <- as.data.frame(DEG_logFC[,-1])
rownames(nCoV_pneu_logFC_df) <- DEG_logFC$symbol

gplots::heatmap.2(as.matrix(nCoV_pneu_logFC_df), trace = "none", 
                  col="bluered", scale="none", Rowv=F,
                  srtCol=50, keysize=1, Colv=T,
                  ylab=paste(nrow(nCoV_pneu_logFC_df), " Genes"), 
                  margins = c(10, 8))

################################################################################
# logFC6 EXP heatmap
load("../results/nCoV_pneu_Heal_norm_symbol_GE.RData")

iid <- c(nCoV_samp, Heal_samp, pneuVir_samp)
DEG_expr <- nCoV_pneu_Heal_norm_symbol_GE[DEG_logFC$symbol, iid]
gplots::heatmap.2(as.matrix(DEG_expr), trace = "none", 
                  col="bluered", scale="row", 
                  srtCol=50, keysize=1, Colv=T,
                  ylab=paste(nrow(DEG_expr), " Genes"), 
                  # lmat=rbind(c(1,3),c(2,4)),
                  margins = c(10, 8))


