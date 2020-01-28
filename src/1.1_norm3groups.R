
load("../results/0.1_preprocessing.RData")
library(edgeR)
library(limma)
################################################################################
# 3 groups
sample_meta$Group2 <- sample_meta$Group
sample_meta$Group2 <- gsub("nCoV","pneumonia", sample_meta$Group2)
subsample_pheno <- dplyr::filter(sample_meta, Group2 %in% c("Healthy", "pneumonia") )
subsample_pheno$Group2 <- factor(subsample_pheno$Group2, levels=c("Healthy", "pneumonia") )
Expdesign <- model.matrix(~subsample_pheno$Group2)

S_raw <- as.matrix(dplyr::select(sample_count_df_symbol, -SYMBOL) )
rownames(S_raw) <- sample_count_df_symbol$SYMBOL
dge <- DGEList(counts=S_raw)

subsample_dge <- dge[,subsample_pheno$sampleID]

# filter and norm
keep <- filterByExpr(subsample_dge, Expdesign)
subsample_dge <- subsample_dge[keep,,keep.lib.sizes=FALSE]
subsample_dge <- calcNormFactors(subsample_dge)

v <- voom(subsample_dge, Expdesign, plot=FALSE, normalize="quantile")

nCoV_pneu_Heal_norm_symbol_GE <- v$E
boxplot(nCoV_pneu_Heal_norm_symbol_GE)

################################################################################
#PCA 3 groups
load("../results/DEA_list.RData")
DEG_nCoV_Heal <- rownames(DEA_list[["limma_DE"]][["nCoV_Heal"]])
DEG_pneu_Heal <- rownames(DEA_list[["limma_DE"]][["pneu_Heal"]])
DEG_nCoV_pneu <- rownames(DEA_list[["limma_DE"]][["nCoV_pneu"]])

DEG_united <- unique(c(DEG_nCoV_Heal, DEG_pneu_Heal, DEG_nCoV_pneu))

library(ggfortify)
# only DEG_nCoV_Heal
autoplot(prcomp(t(v$E[DEG_nCoV_Heal,]), scale=T), 
         data=subsample_pheno, colour = "Group", 
         size = 8, label = TRUE, label.colour="black", ts.colour="black" )

# all union genes
autoplot(prcomp(t(v$E[intersect(DEG_united, rownames(v$E) ),]), scale=T), 
         data=subsample_pheno, colour = "Group", 
         size = 8, label = TRUE, label.colour="black", ts.colour="black" )

readr::write_tsv(as.data.frame(v$E), path="../results/nCoV_pneu_Heal_norm_symbol_GE.tsv")

save.image("../results/1.1_norm3groups.RData")
save(nCoV_pneu_Heal_norm_symbol_GE, sample_meta, file="../results/nCoV_pneu_Heal_norm_symbol_GE.RData")

