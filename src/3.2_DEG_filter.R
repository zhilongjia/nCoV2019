

################################################################################
load("../results/DEA_list.RData")
load("../results/DEA_pneu_list.RData")
DEA_list[["limma_DEA"]] <- c(DEA_list[["limma_DEA"]], DEA_pneu_list[["limma_DEA"]])
DEA_list[["limma_DE"]] <- c(DEA_list[["limma_DE"]], DEA_pneu_list[["limma_DE"]])

################################################################################
# filter lower logfc genes (nCoV_Heal)
nCoV_logFC <- dplyr::filter(DEA_list[["limma_DEA"]][["nCoV_Heal"]], abs(logFC)>=4) 

readr::write_csv(nCoV_logFC, path="../results/nCoV_logFC4.csv")
readr::write_csv(nCoV_logFC, path="../results/nCoV_logFC6.csv")

#
pneu_logFC <- dplyr::filter(DEA_list[["limma_DEA"]][["pneu_Heal"]], abs(logFC)>=3)
readr::write_csv(pneu_logFC, path="../results/pneu_logFC3.csv")

#
pneuVir_logFC <- dplyr::filter(DEA_list[["limma_DEA"]][["pneuVir_Heal"]], abs(logFC)>=3)
readr::write_csv(pneuVir_logFC, path="../results/pneuVir_logFC3.csv")


nCoV_Heal_DE <- DEA_list[["limma_DE"]][["nCoV_Heal"]][nCoV_logFC$nCoV_Heal,]
nCoV_Heal_pheno <- DEA_list[["subsample_pheno"]][["nCoV_Heal"]] 


################################################################################

nCoV_logFC_top100 <- dplyr::top_n(DEA_list[["limma_DEA"]][["nCoV_Heal"]], 100, abs(logFC)) 
pneu_logFC_top100 <- dplyr::top_n(DEA_list[["limma_DEA"]][["pneu_Heal"]], 100, abs(logFC)) 

nCoV_pneu_logFC_list <- list(nCoV_logFC_top100=nCoV_logFC_top100$nCoV_Heal, 
     pneu_logFC_top100=pneu_logFC_top100$pneu_Heal)

venn::venn(nCoV_pneu_logFC_list, ilab=TRUE, zcolor = "style", 
           box=FALSE, ilcs=1.5, sncs=1.2, cexil=10, cexsn=10)

load("../results/nCoV_pneu_comp_logFC.RData")

################################################################################
nCoV_logFC <- dplyr::filter(DEA_list[["limma_DEA"]][["nCoV_Heal"]], abs(logFC)>=6) 

# filter gene and heatmap
nCoV_pneu_logFC <- dplyr::filter(nCoV_pneu_comp_logFC, symbol %in% c(nCoV_pneu_logFC_list$nCoV_logFC_top100, nCoV_pneu_logFC_list$pneu_logFC_top100) ) 
nCoV_pneu_logFC_df <- as.data.frame(nCoV_pneu_logFC[,-1])
rownames(nCoV_pneu_logFC_df) <- nCoV_pneu_logFC$symbol
nCoV_pneu_logFC_df[is.na(nCoV_pneu_logFC_df)] <- 0


gplots::heatmap.2(as.matrix(nCoV_pneu_logFC_df), trace = "none", 
                  col="bluered", scale="none", 
                  srtCol=50, keysize=1, Colv=T,
                  ylab=paste(nrow(nCoV_pneu_logFC_df), " Genes"), 
                  margins = c(10, 6))

################################################################################
# filter gene and heatmap 
nCoV_pneu_logFC <- dplyr::filter(nCoV_pneu_comp_logFC, symbol %in% nCoV_logFC$nCoV_Heal ) 
nCoV_pneu_logFC_df <- as.data.frame(nCoV_pneu_logFC[,-1])
rownames(nCoV_pneu_logFC_df) <- nCoV_pneu_logFC$symbol
nCoV_pneu_logFC_df[is.na(nCoV_pneu_logFC_df)] <- 0

gplots::heatmap.2(as.matrix(nCoV_pneu_logFC_df), trace = "none", 
                  col="bluered", scale="none", 
                  srtCol=50, keysize=1, Colv=T,
                  ylab=paste(nrow(nCoV_pneu_logFC_df), " Genes"), 
                  margins = c(10, 6))
