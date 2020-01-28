
load("../results/0.1_preprocessing.RData")

################################################################################
# Normlization
library(edgeR)
library(limma)

# All gene count noramlization
S_raw <- as.matrix(dplyr::select(sample_count_df_symbol, -SYMBOL) )
rownames(S_raw) <- sample_count_df_symbol$SYMBOL
dge <- DGEList(counts=S_raw)

################################################################################
comp2groups <- list(nCoV_Heal=c("Healthy", "nCoV"), 
                    pneu_Heal=c("Healthy", "pneumonia"), 
                    nCoV_pneu=c("pneumonia", "nCoV") )

DEA_list <- list()

for (comp_var in names(comp2groups) ) {
    # comp_var <- names(comp2groups)[[1]]
    comp_i <- comp2groups[[comp_var]]
    print (comp_i)
    subsample_pheno <- dplyr::filter(sample_meta, Group %in% comp_i )
    subsample_pheno$Group <- factor(subsample_pheno$Group, levels=comp_i )
    Expdesign <- model.matrix(~subsample_pheno$Group)
    
    subsample_dge <- dge[,subsample_pheno$sampleID]
    
    # DEA
    keep <- filterByExpr(subsample_dge, Expdesign)
    subsample_dge <- subsample_dge[keep,,keep.lib.sizes=FALSE]
    subsample_dge <- calcNormFactors(subsample_dge)
    
    v <- voom(subsample_dge, Expdesign, plot=FALSE, normalize="quantile")
    
    # boxplot(v$E[,nCoV_samp])
    # boxplot(v$E)
    
    Expfit1 <- lmFit(v, Expdesign)
    Expfit2 <- eBayes(Expfit1)
    limma_res <- topTable(Expfit2, coef=tail(colnames(Expdesign), 1), number=Inf) %>% 
        tibble::rownames_to_column(var=comp_var )
    
    limma_DEA <- dplyr::filter(limma_res, adj.P.Val<=0.05, abs(logFC)>=2 )
    readr::write_tsv(limma_DEA, paste0("../results/", comp_var, "_limma_res.tsv") )
    
    limma_UP <- dplyr::filter(limma_DEA, logFC>0)[[comp_var]]
    limma_DN <- dplyr::filter(limma_DEA, logFC<0)[[comp_var]]
    
    limma_DEG <- limma_DEA[[comp_var]]
    print(length(limma_DEG))
    limma_DE <- v$E[limma_DEG,]
    
    
    DEA_list[["limma_res"]][[comp_var]] <- limma_res
    DEA_list[["limma_DEA"]][[comp_var]] <- limma_DEA
    DEA_list[["limma_DEG"]][[comp_var]] <- limma_DEG
    DEA_list[["limma_UPDN"]][[paste0(comp_var, "_UP")]] <- limma_UP
    DEA_list[["limma_UPDN"]][[paste0(comp_var, "_DN")]] <- limma_DN
    DEA_list[["limma_DE"]][[comp_var]] <- limma_DE
    DEA_list[["subsample_pheno"]][[comp_var]] <- subsample_pheno
    
}

save.image("../results/1.0_DE.RData")

################################################################################
library(ggfortify)
autoplot(prcomp(t(DEA_list[["limma_DE"]][["nCoV_pneu"]]), scale=F), 
         data=DEA_list[["subsample_pheno"]][["nCoV_pneu"]], colour = "Group", 
         size = 8, label = TRUE, label.colour="black", ts.colour="black" )

library(venn)
venn(DEA_list[["limma_UPDN"]][c("nCoV_Heal_UP", "nCoV_Heal_DN", "pneu_Heal_UP", "pneu_Heal_DN")], ilab=TRUE, zcolor = "style", cexil=10, cexsn=10)

save(S_raw, DEA_list, sample_meta, file="../results/DEA_list.RData")


