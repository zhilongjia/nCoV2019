

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


# add annotation
source("../../../SynologyDrive/zzz.R")
DEGs_anno <- gene2KEGG( nCoV_pneu_comp_logFC$symbol )
nCoV_pneu_comp_logFC_anno <- dplyr::left_join(nCoV_pneu_comp_logFC, DEGs_anno, by=c("symbol"="SYMBOL") )
readr::write_tsv(nCoV_pneu_comp_logFC_anno, path="../results/nCoV_pneu_comp_logFC_anno.txt")


################################################################################
# filter gene and heatmap
nCoV_pneu_logFC <- dplyr::filter(nCoV_pneu_comp_logFC, abs(nCoV_Heal)>=4) 
nCoV_pneu_logFC_df <- as.data.frame(nCoV_pneu_logFC[,-1])
rownames(nCoV_pneu_logFC_df) <- nCoV_pneu_logFC$symbol

gplots::heatmap.2(as.matrix(nCoV_pneu_logFC_df), trace = "none", 
                  col="bluered", scale="none", 
                  srtCol=50, keysize=1, Colv=T,
                  ylab=paste(nrow(nCoV_pneu_logFC_df), " Genes"), 
                  margins = c(10, 6))

################################################################################
# comapred with innatedb_curated_genes
innatedb_curated_genes <- readxl::read_xls("../data/innatedb_curated_genes.xls", sheet="Sheet2")

limma_DEG_list <- c(DEA_list[["limma_DEG"]], DEA_pneu_list[["limma_DEG"]])

nCoV_pneu_innate_list <- limma_DEG_list[c("nCoV_Heal", "pneu_Heal", "pneuVir_Heal")]
nCoV_pneu_innate_list[["innatedb"]] <- innatedb_curated_genes$symbol


venn::venn(nCoV_pneu_innate_list, ilab=TRUE, zcolor = "style", 
           box=FALSE, ilcs=1.5, sncs=1.2, cexil=10, cexsn=10)


################################################################################
# innate heatmap
overlapped_Innate <- intersect(nCoV_pneu_comp_logFC$symbol, innatedb_curated_genes$symbol)


nCoV_pneu_innate_logFC <- dplyr::filter(nCoV_pneu_comp_logFC, symbol %in%  overlapped_Innate )

nCoV_pneu_innate_logFC_df <- nCoV_pneu_innate_logFC[,-1]
rownames(nCoV_pneu_innate_logFC_df) <- nCoV_pneu_innate_logFC$symbol

iid <- which(rowSums(is.na(nCoV_pneu_innate_logFC_df)) != 4)
nCoV_pneu_innate_logFC_df <- nCoV_pneu_innate_logFC_df[iid,]
nCoV_pneu_innate_logFC_df[is.na(nCoV_pneu_innate_logFC_df)] <-0

gplots::heatmap.2(as.matrix(nCoV_pneu_innate_logFC_df), trace = "none", 
                  col="bluered", scale="none", 
                  srtCol=50, keysize=1, Colv=T,
                  ylab=paste(nrow(nCoV_pneu_innate_logFC_df), " Genes"), 
                  margins = c(10, 6))







