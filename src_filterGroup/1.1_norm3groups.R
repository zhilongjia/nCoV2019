
load("../data/genelist_filter_lz.RData")


load("../results/sample_count_df_symbol.RData")
library(edgeR)
library(limma)
################################################################################

g3_names <- c("Healthy", "Others", "Viral-like", "nCoV" )
subsample_pheno <- dplyr::filter(sample_meta, Types %in% g3_names )
subsample_pheno$Types <- factor(subsample_pheno$Types, levels = g3_names )
Expdesign <- model.matrix(~subsample_pheno$Types)

S_raw <- as.matrix(dplyr::select(sample_count_df_symbol, -SYMBOL) )
rownames(S_raw) <- sample_count_df_symbol$SYMBOL
dge <- DGEList(counts=S_raw)

subsample_dge <- dge[,subsample_pheno$NewID]

# filter and norm
keep <- filterByExpr(subsample_dge, Expdesign)
# use filter genes from Li Zhang (in gene.list)
keep[names(keep)] <- FALSE

keep.gene <- union(union(rownames(gene.list$nCoV_Heal), rownames(gene.list$vir_Heal) ), rownames(gene.list$Others_Heal) )

keep[keep.gene] <- TRUE


subsample_dge <- subsample_dge[keep,,keep.lib.sizes=FALSE]
subsample_dge <- calcNormFactors(subsample_dge)

v <- voom(subsample_dge, Expdesign, plot=FALSE, normalize="quantile")

nCoV_pneu_Heal_norm_symbol_GE <- v$E


################################################################################
#PCA 3 groups
load("../results/DEA_pneu_list.RData")

DEG_nCoV_Heal <- rownames(DEA_list[["limma_DE"]][["nCoV_Heal"]])
DEG_pneuVir_Heal <- rownames(DEA_list[["limma_DE"]][["Vir_Heal"]])
DEG_pneuBac_Heal <- rownames(DEA_list[["limma_DE"]][["Others_Heal"]])
DEG_nCoV_pneuVir <- rownames(DEA_list[["limma_DE"]][["nCoV_Vir"]])

DEG_united <- unique(c(DEG_nCoV_Heal, DEG_pneuVir_Heal, DEG_pneuBac_Heal))

library(ggfortify)
# all union genes
autoplot(prcomp(t(v$E[intersect(DEG_united, rownames(v$E) ),]), scale=F), 
         data=subsample_pheno, colour = "Types", 
         size = 5, label = F, label.colour="black", ts.colour="black" )

autoplot(prcomp(t(v$E[intersect(DEG_united, rownames(v$E) ),]), scale=F), x=1,y=3,
         data=subsample_pheno, colour = "Types", 
         size = 5, label = F, label.colour="black", ts.colour="black" )

autoplot(prcomp(t(v$E[intersect(DEG_united, rownames(v$E) ),]), scale=F), x=1,y=2,
         data=subsample_pheno, colour = "virus", 
         size = 5, label = F, label.colour="black", ts.colour="black" )


readr::write_tsv(tibble::rownames_to_column(as.data.frame(v$E)), path="../results/nCoV_pneu_Heal_norm_symbol_GE.tsv")

save.image("../results/1.1_norm3groups.RData")
save(nCoV_pneu_Heal_norm_symbol_GE, sample_meta, file="../results/nCoV_pneu_Heal_norm_symbol_GE.RData")


