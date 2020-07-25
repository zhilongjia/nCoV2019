
# load("../results/0.1_preprocessing.RData")
load("../results/sample_count_df_symbol.RData")
library(edgeR)
library(limma)
################################################################################
# consider 3 groups as 2 groups
# subsample_pheno <- dplyr::filter(sample_meta, Group2 %in% c("Healthy", "Disease") )
# subsample_pheno$Group2 <- factor(subsample_pheno$Group2, levels=c("Healthy", "Disease") )
# Expdesign <- model.matrix(~subsample_pheno$Group2)

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
subsample_dge <- subsample_dge[keep,,keep.lib.sizes=FALSE]
subsample_dge <- calcNormFactors(subsample_dge)

v <- voom(subsample_dge, Expdesign, plot=FALSE, normalize="quantile")

nCoV_pneu_Heal_norm_symbol_GE <- v$E

# par(mar=c(10,2,1,1))
# boxplot(nCoV_pneu_Heal_norm_symbol_GE[,c(nCoV_samp, Heal_samp)], las=2)
# 
# boxplot(nCoV_pneu_Heal_norm_symbol_GE[,c(nCoV_samp, Heal_samp, pneuVir_samp)], las=2)

################################################################################
#PCA 3 groups
# load("../results/DEA_list.RData")
load("../results/DEA_pneu_list.RData")
DEA_list <- DEA_pneu_list

DEG_nCoV_Heal <- rownames(DEA_list[["limma_DE"]][["nCoV_Heal"]])
DEG_pneuVir_Heal <- rownames(DEA_list[["limma_DE"]][["Vir_Heal"]])
DEG_pneuBac_Heal <- rownames(DEA_list[["limma_DE"]][["Others_Heal"]])
# DEG_nCoV_pneuBac <- rownames(DEA_list[["limma_DE"]][["nCoV_pneuBac"]])
DEG_nCoV_pneuVir <- rownames(DEA_list[["limma_DE"]][["nCoV_Vir"]])


# DEG_pneu_Heal <- rownames(DEA_list[["limma_DE"]][["pneu_Heal"]])
# DEG_nCoV_pneu <- rownames(DEA_list[["limma_DE"]][["nCoV_pneu"]])

DEG_united <- unique(c(DEG_nCoV_Heal, DEG_pneuVir_Heal, DEG_pneuBac_Heal))

library(ggfortify)
# only DEG_nCoV_Heal
# autoplot(prcomp(t(v$E[DEG_nCoV_Heal,]), scale=T), x=1,y=2,
#          data=subsample_pheno, colour = "Types", 
#          size = 5, label = F, label.colour="black", ts.colour="black" )
# 
# autoplot(prcomp(t(v$E[DEG_nCoV_Heal,]), scale=F), x=1,y=3,
#          data=subsample_pheno, colour = "Group1", 
#          size = 5, label = F, label.colour="black", ts.colour="black" )
# 
# autoplot(prcomp(t(v$E[DEG_nCoV_Heal,]), scale=F), x=1,y=2,
#          data=subsample_pheno, colour = "Types", 
#          size = 5, label = F, label.colour="black", ts.colour="black" )
# 
# autoplot(prcomp(t(v$E[DEG_nCoV_Heal,]), scale=F), x=1,y=2,
#          data=subsample_pheno, colour = "FTD", 
#          size = 5, label = F, label.colour="black", ts.colour="black" )
# 
# autoplot(prcomp(t(v$E[DEG_nCoV_Heal,]), scale=F), 
#          data=subsample_pheno, colour = 'Abundance_among_ABFV', 
#          size = 8, label = TRUE, label.colour="black", ts.colour="black" )

# all union genes
autoplot(prcomp(t(v$E[intersect(DEG_united, rownames(v$E) ),]), scale=T), 
         data=subsample_pheno, colour = "Types", 
         size = 5, label = F, label.colour="black", ts.colour="black" )

autoplot(prcomp(t(v$E[intersect(DEG_united, rownames(v$E) ),]), scale=F), x=1,y=3,
         data=subsample_pheno, colour = "Types", 
         size = 5, label = F, label.colour="black", ts.colour="black" )

autoplot(prcomp(t(v$E[intersect(DEG_united, rownames(v$E) ),]), scale=F), x=1,y=2,
         data=subsample_pheno, colour = "virus", 
         size = 5, label = F, label.colour="black", ts.colour="black" )

# gplots::heatmap.2(v$E[DEG_nCoV_Heal,], trace = "none", 
#                   col="bluered", scale="row", labRow=FALSE,
#                   srtCol=60, keysize=1, Colv=TRUE,
#                   ylab=paste(nrow(v$E[DEG_nCoV_Heal,]), "DE Genes"))
# 
# 
# gplots::heatmap.2(v$E[intersect(DEG_united, rownames(v$E) ),], trace = "none", 
#                   col="bluered", scale="row", labRow=FALSE,
#                   srtCol=60, keysize=1, Colv=TRUE,
#                   ylab=paste(nrow(v$E[intersect(DEG_united, rownames(v$E) ),]), "DE Genes"), 
#                   margins = c(10, 8))

readr::write_tsv(tibble::rownames_to_column(as.data.frame(v$E)), path="../results/nCoV_pneu_Heal_norm_symbol_GE.tsv")

save.image("../results/1.1_norm3groups.RData")
save(nCoV_pneu_Heal_norm_symbol_GE, sample_meta, file="../results/nCoV_pneu_Heal_norm_symbol_GE.RData")


