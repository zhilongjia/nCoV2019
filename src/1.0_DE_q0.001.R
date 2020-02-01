
# load("../results/0.1_preprocessing.RData")
load("../results/sample_count_df_symbol.RData")
################################################################################
# Normlization
library(edgeR)
library(limma)

# All gene count noramlization
S_raw <- as.matrix(dplyr::select(sample_count_df_symbol, -SYMBOL) )
rownames(S_raw) <- sample_count_df_symbol$SYMBOL
dge <- DGEList(counts=S_raw)

################################################################################
mRNA_lncRNA <- readr::read_tsv("../data/ensembl_symbol_type.tsv",col_names = FALSE)
colnames(mRNA_lncRNA) <- c("ensembl", "symbol", "type")
mRNA_lncRNA <- unique(mRNA_lncRNA[,-1])
################################################################################
comp2groups <- list(nCoV_Heal=c("Healthy", "nCoV"), 
                    pneu_Heal=c("Healthy", "pneumonia"), 
                    nCoV_pneu=c("pneumonia", "nCoV") )

DEA_list <- list()

for (comp_var in names(comp2groups) ) {
    # comp_var <- names(comp2groups)[[1]]
    comp_i <- comp2groups[[comp_var]]
    print (comp_i)
    subsample_pheno <- dplyr::filter(sample_meta, Group1 %in% comp_i )
    subsample_pheno$Group1 <- factor(subsample_pheno$Group1, levels=comp_i )
    Expdesign <- model.matrix(~subsample_pheno$Group1)
    
    # DEA
    subsample_dge <- dge[,subsample_pheno$sampleID]
    keep <- filterByExpr(subsample_dge, Expdesign)
    subsample_dge <- subsample_dge[keep,,keep.lib.sizes=FALSE]
    subsample_dge <- calcNormFactors(subsample_dge)
    
    v <- voom(subsample_dge, Expdesign, plot=FALSE, normalize="quantile")
    
    # boxplot test
    # par(mar=c(9,2,1,1))
    # boxplot(v$E[,nCoV_samp], las=2)
    # boxplot(v$E, las=2)
    
    Expfit1 <- lmFit(v, Expdesign)
    Expfit2 <- eBayes(Expfit1)
    limma_res <- topTable(Expfit2, coef=tail(colnames(Expdesign), 1), number=Inf) %>% 
        tibble::rownames_to_column() %>% 
        dplyr::left_join(mRNA_lncRNA, by=c("rowname"="symbol")) 
    colnames(limma_res)[1] <- comp_var
    
    # only protein_coding genes
    limma_DEA <- dplyr::filter(limma_res, adj.P.Val<=0.001, abs(logFC)>=2, type=="protein_coding" )
    readr::write_tsv(limma_DEA, paste0("../results/", comp_var, "_limma_q0.001_res.tsv") )
    
    limma_UP <- dplyr::filter(limma_DEA, logFC>0)[[comp_var]]
    limma_DN <- dplyr::filter(limma_DEA, logFC<0)[[comp_var]]
    
    limma_DEG <- limma_DEA[[comp_var]]
    print(length(limma_DEG))
    limma_DE <- v$E[limma_DEG,]
    
    # limma_res include lncRNA, other objs exclude lncRNA
    DEA_list[["limma_res"]][[comp_var]] <- limma_res
    DEA_list[["limma_DEA"]][[comp_var]] <- limma_DEA
    DEA_list[["limma_DEG"]][[comp_var]] <- limma_DEG
    DEA_list[["limma_UPDN"]][[paste0(comp_var, "_UP")]] <- limma_UP
    DEA_list[["limma_UPDN"]][[paste0(comp_var, "_DN")]] <- limma_DN
    DEA_list[["limma_DE"]][[comp_var]] <- limma_DE
    DEA_list[["subsample_pheno"]][[comp_var]] <- subsample_pheno
    
}

save.image("../results/1.0_DE_q0.001.RData")

################################################################################
library(ggfortify)
autoplot(prcomp(t(DEA_list[["limma_DE"]][["nCoV_pneu"]]), scale=F), 
         data=DEA_list[["subsample_pheno"]][["nCoV_pneu"]], colour = "Group1", 
         size = 5, label = TRUE, label.colour="black", ts.colour="black" )

library(venn)

venn(DEA_list[["limma_DEG"]], ilab=TRUE, zcolor = "style", 
     box=FALSE, ilcs=1.5, sncs=1.2, cexil=10, cexsn=10)

venn(DEA_list[["limma_UPDN"]][c("nCoV_Heal_UP", "nCoV_Heal_DN", "pneu_Heal_UP", "pneu_Heal_DN")], 
     ilab=TRUE, zcolor = "style", box=FALSE, ilcs=1.5, sncs=1.2, cexil=10, cexsn=10)

save(S_raw, DEA_list, sample_meta, file="../results/DEA_list_q0.001.RData")


