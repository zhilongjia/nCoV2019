

# load("../results/DEA_list.RData")
load("../results/DEA_pneu_list.RData")

################################################################################
# limma_res_list <- c(DEA_list[["limma_res"]], DEA_pneu_list[["limma_res"]])
limma_res_list <- DEA_list[["limma_res"]]

limma_res_list_df <- list()

for (comp_i in names(limma_res_list) ) {
    # comp_i <- "nCoV_Heal"
    res_df_i <- dplyr::filter(limma_res_list[[comp_i]], type=="protein_coding")[,1:2]
    
    colnames(res_df_i) <- c("symbol", comp_i)
    limma_res_list_df[[comp_i]] <- res_df_i
    
}


Allgene_logFC <- Reduce(dplyr::full_join, limma_res_list_df) %>% as.data.frame()

################################################################################
## abs(nCoV_Heal)>=5 & DEG_nCoV_Vir

DEG_nCoV_Vir <- DEA_list[["limma_DEG"]][["nCoV_Vir"]]

DEG_logFC <- dplyr::select(Allgene_logFC, symbol, nCoV_Heal, Vir_Heal, Others_Heal, nCoV_Vir) %>% 
    dplyr::filter(abs(nCoV_Heal)>=5 & (symbol %in% DEG_nCoV_Vir) ) %>% 
    dplyr::arrange(desc(nCoV_Heal))

DEG_logFC[is.na(DEG_logFC)] <- 0

nCoV_pneu_logFC_df <- as.data.frame(DEG_logFC[,-1])
rownames(nCoV_pneu_logFC_df) <- DEG_logFC$symbol

gplots::heatmap.2(as.matrix(nCoV_pneu_logFC_df), trace = "none", density.info="none",
                  col=pal$s, scale="none", Rowv=F,
                  srtCol=50, keysize=1, Colv=T,
                  ylab=paste(nrow(nCoV_pneu_logFC_df), " Genes"), 
                  margins = c(10, 6))

################################################################################
# comapred with innatedb_curated_genes
g3 <- c("nCoV_Heal", "Vir_Heal", "Others_Heal")

innatedb_curated_genes <- readxl::read_xls("../data/innatedb_curated_genes.xls", sheet="Sheet2")

limma_DEG_list <- DEA_list[["limma_DEG"]]

nCoV_pneu_innate_list <- limma_DEG_list[g3]

nCoV_pneu_innate_list[["innatedb"]] <- innatedb_curated_genes$symbol


venn::venn(nCoV_pneu_innate_list, ilab=TRUE, zcolor = "style", 
           box=FALSE, ilcs=1.5, sncs=1.2, cexil=10, cexsn=10)

# output veen res
venn_res <- gplots::venn(nCoV_pneu_innate_list)
aa <- attr(venn_res,"intersections")
aa_df <- as.data.frame(sapply(aa, "length<-", max(lengths(aa))))
readr::write_csv(aa_df, path="../results/venn/nCoV_pneu_innate_venn_res.csv", na="")


################
# innate heatmap
overlapped_Innate <- intersect(c(nCoV_pneu_innate_list$nCoV_Heal, nCoV_pneu_innate_list$Vir_Heal, nCoV_pneu_innate_list$Others_Heal), 
                               innatedb_curated_genes$symbol)


nCoV_pneu_innate_logFC <- dplyr::filter(Allgene_logFC, symbol %in%  overlapped_Innate )

nCoV_pneu_innate_logFC_df <- nCoV_pneu_innate_logFC[,-1]
rownames(nCoV_pneu_innate_logFC_df) <- nCoV_pneu_innate_logFC$symbol

# iid <- which(rowSums(is.na(nCoV_pneu_innate_logFC_df)) != 4)
# nCoV_pneu_innate_logFC_df <- nCoV_pneu_innate_logFC_df[iid,]
nCoV_pneu_innate_logFC_df[is.na(nCoV_pneu_innate_logFC_df)] <-0

nCoV_pneu_innate_logFC_df <- nCoV_pneu_innate_logFC_df[order(-nCoV_pneu_innate_logFC_df$nCoV_Heal),g3]


gplots::heatmap.2(as.matrix(nCoV_pneu_innate_logFC_df), trace = "none", density.info="none",
                  col=pal$s, scale="none", labRow=NA,
                  srtCol=50, keysize=1, Colv=T, Rowv=F,
                  ylab=paste(nrow(nCoV_pneu_innate_logFC_df), " Genes"), 
                  margins = c(10, 3))

######
# exp
load("../results/nCoV_pneu_Heal_norm_symbol_GE.RData")

nCoV_samp <- dplyr::filter(sample_meta, Diagnostics=="nCoV pneumonia")  %>% 
    dplyr::select(NewID) %>% unlist(use.names = FALSE)

innated_pneu_gene <- intersect(overlapped_Innate, rownames(nCoV_pneu_Heal_norm_symbol_GE) )

nCoV_exp <- nCoV_pneu_Heal_norm_symbol_GE[innated_pneu_gene, c(nCoV_samp,Heal_samp)]
gplots::heatmap.2(as.matrix(nCoV_exp), trace = "none", density.info="none",
                  col=pal$s, scale="row", labRow=NA,
                  srtCol=50, keysize=1, Colv=T, Rowv=F,
                  ylab=paste(nrow(nCoV_exp), " Genes"), 
                  margins = c(10, 3))


################################################################################
################################################################################
################################################################################
# below to be delete


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


