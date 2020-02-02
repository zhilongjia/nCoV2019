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


g3 <- c("nCoV_Heal", "pneu_Heal", "pneuVir_Heal")

################################################################################
# comapred with innatedb_curated_genes
innatedb_curated_genes <- readxl::read_xls("../data/innatedb_curated_genes.xls", sheet="Sheet2")

limma_DEG_list <- c(DEA_list[["limma_DEG"]], DEA_pneu_list[["limma_DEG"]])

nCoV_pneu_innate_list <- limma_DEG_list[g3]

nCoV_pneu_innate_list[["innatedb"]] <- innatedb_curated_genes$symbol


venn::venn(nCoV_pneu_innate_list, ilab=TRUE, zcolor = "style", 
           box=FALSE, ilcs=1.5, sncs=1.2, cexil=10, cexsn=10)


################
# innate heatmap
overlapped_Innate <- intersect(c(nCoV_pneu_innate_list$nCoV_Heal, nCoV_pneu_innate_list$pneu_Heal, nCoV_pneu_innate_list$pneuVir_Heal), 
                               innatedb_curated_genes$symbol)


nCoV_pneu_innate_logFC <- dplyr::filter(Allgene_logFC, symbol %in%  overlapped_Innate )

nCoV_pneu_innate_logFC_df <- nCoV_pneu_innate_logFC[,-1]
rownames(nCoV_pneu_innate_logFC_df) <- nCoV_pneu_innate_logFC$symbol

# iid <- which(rowSums(is.na(nCoV_pneu_innate_logFC_df)) != 4)
# nCoV_pneu_innate_logFC_df <- nCoV_pneu_innate_logFC_df[iid,]
nCoV_pneu_innate_logFC_df[is.na(nCoV_pneu_innate_logFC_df)] <-0

nCoV_pneu_innate_logFC_df <- nCoV_pneu_innate_logFC_df[order(-nCoV_pneu_innate_logFC_df$nCoV_Heal),g3]


gplots::heatmap.2(as.matrix(nCoV_pneu_innate_logFC_df), trace = "none", 
                  col="bluered", scale="none", 
                  srtCol=50, keysize=1, Colv=T, Rowv=F,
                  ylab=paste(nrow(nCoV_pneu_innate_logFC_df), " Genes"), 
                  margins = c(10, 6))

################################################################################
# SARS heatmap

SARS_gene <- readxl::read_xlsx("../data/SARS_up_comparison_result.xlsx")

limma_DEG_list <- c(DEA_list[["limma_DEG"]], DEA_pneu_list[["limma_DEG"]])

nCoV_pneu_SARS_list <- limma_DEG_list[g3]

nCoV_pneu_SARS_list[["SARS"]] <- SARS_gene$t_name


venn::venn(nCoV_pneu_SARS_list, ilab=TRUE, zcolor = "style", 
           box=FALSE, ilcs=1.5, sncs=1.2, cexil=10, cexsn=10)


################
# innate heatmap
overlapped_SARS <- intersect(c(nCoV_pneu_SARS_list$nCoV_Heal, nCoV_pneu_SARS_list$pneu_Heal, nCoV_pneu_SARS_list$pneuVir_Heal), 
                               SARS_gene$t_name)


nCoV_pneu_SARS_logFC <- dplyr::filter(Allgene_logFC, symbol %in%  overlapped_SARS )

nCoV_pneu_SARS_logFC_df <- nCoV_pneu_SARS_logFC[,-1]
rownames(nCoV_pneu_SARS_logFC_df) <- nCoV_pneu_SARS_logFC$symbol

# iid <- which(rowSums(is.na(nCoV_pneu_innate_logFC_df)) != 4)
# nCoV_pneu_innate_logFC_df <- nCoV_pneu_innate_logFC_df[iid,]
nCoV_pneu_SARS_logFC_df[is.na(nCoV_pneu_SARS_logFC_df)] <- 0

nCoV_pneu_SARS_logFC_df <- nCoV_pneu_SARS_logFC_df[order(-nCoV_pneu_SARS_logFC_df$nCoV_Heal),g3]


gplots::heatmap.2(as.matrix(nCoV_pneu_SARS_logFC_df), trace = "none", 
                  col="bluered", scale="none", 
                  srtCol=50, keysize=1, Colv=T, Rowv=F,
                  ylab=paste(nrow(nCoV_pneu_SARS_logFC_df), " Genes"), 
                  margins = c(12, 6))


