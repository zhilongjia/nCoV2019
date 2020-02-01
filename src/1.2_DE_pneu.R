
# load("../results/0.1_preprocessing.RData")
load("../results/sample_count_df_symbol.RData")
################################################################################
# Normlization
library(limma)
library(edgeR)

# All gene count noramlization
S_raw <- as.matrix(dplyr::select(sample_count_df_symbol, -SYMBOL) )
rownames(S_raw) <- sample_count_df_symbol$SYMBOL
dge <- DGEList(counts=S_raw)

sample_meta <- readr::read_csv("../data/sample_meta-v2.csv")

################################################################################
mRNA_lncRNA <- readr::read_tsv("../data/ensembl_symbol_type.tsv",col_names = FALSE)
colnames(mRNA_lncRNA) <- c("ensembl", "symbol", "type")
mRNA_lncRNA <- unique(mRNA_lncRNA[,-1])
################################################################################
comp2groups <- list(pneuVir_Heal=c("Healthy", "pneumonia-vir"), 
                    pneuBac_Heal=c("Healthy", "pneumonia-bca"), 
                    pneuBac_pneuVir=c("pneumonia-bca", "pneumonia-vir"),
                    nCoV_pneuBac=c("pneumonia-bca", "nCoV"),
                    nCoV_pneuVir=c("pneumonia-vir", "nCoV") )

DEA_pneu_list <- list()

for (comp_var in names(comp2groups) ) {
    # comp_var <- names(comp2groups)[[1]]
    comp_i <- comp2groups[[comp_var]]
    print (comp_i)
    subsample_pheno <- dplyr::filter(sample_meta, Group3 %in% comp_i )
    subsample_pheno$Group3 <- factor(subsample_pheno$Group3, levels=comp_i )
    Expdesign <- model.matrix(~subsample_pheno$Group3)
    
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
        tibble::rownames_to_column() %>% 
        dplyr::left_join(mRNA_lncRNA, by=c("rowname"="symbol")) 
    colnames(limma_res)[1] <- comp_var
    
    
    # only protein_coding genes
    limma_DEA <- dplyr::filter(limma_res, adj.P.Val<=0.05, abs(logFC)>=2, type=="protein_coding" )
    readr::write_tsv(limma_DEA, paste0("../results/", comp_var, "_limma_res.tsv") )
    
    limma_UP <- dplyr::filter(limma_DEA, logFC>0)[[comp_var]]
    limma_DN <- dplyr::filter(limma_DEA, logFC<0)[[comp_var]]
    
    limma_DEG <- limma_DEA[[comp_var]]
    print(length(limma_DEG))
    limma_DE <- v$E[limma_DEG,]
    
    # limma_res include lncRNA, other objs exclude lncRNA
    DEA_pneu_list[["limma_res"]][[comp_var]] <- limma_res
    DEA_pneu_list[["limma_DEA"]][[comp_var]] <- limma_DEA
    DEA_pneu_list[["limma_DEG"]][[comp_var]] <- limma_DEG
    DEA_pneu_list[["limma_UPDN"]][[paste0(comp_var, "_UP")]] <- limma_UP
    DEA_pneu_list[["limma_UPDN"]][[paste0(comp_var, "_DN")]] <- limma_DN
    DEA_pneu_list[["limma_DE"]][[comp_var]] <- limma_DE
    DEA_pneu_list[["subsample_pheno"]][[comp_var]] <- subsample_pheno
    
}

save.image("../results/1.0_DE_pneu.RData")

################################################################################
library(ggfortify)
autoplot(prcomp(t(DEA_pneu_list[["limma_DE"]][["pneuBac_pneuVir"]]), scale=F), 
         data=DEA_pneu_list[["subsample_pheno"]][["pneuBac_pneuVir"]], colour = "Group3", 
         size = 8, label = F, label.colour="black", ts.colour="black" )

autoplot(prcomp(t(DEA_pneu_list[["limma_DE"]][["nCoV_pneuBac"]]), scale=F), 
         data=DEA_pneu_list[["subsample_pheno"]][["nCoV_pneuBac"]], colour = "Group3", 
         size = 8, label = TRUE, label.colour="black", ts.colour="black" )

autoplot(prcomp(t(DEA_pneu_list[["limma_DE"]][["nCoV_pneuVir"]]), scale=F), 
         data=DEA_pneu_list[["subsample_pheno"]][["nCoV_pneuVir"]], colour = "Group3", 
         size = 8, label = TRUE, label.colour="black", ts.colour="black" )



library(venn)

venn::venn(DEA_pneu_list[["limma_DEG"]], ilab=TRUE, zcolor = "style", 
           box=FALSE, ilcs=1.5, sncs=1.2, cexil=10, cexsn=10)

tmp <- c("pneuVir_Heal","pneuBac_Heal")

venn::venn(DEA_pneu_list[["limma_DEG"]][c("pneuVir_Heal","pneuBac_Heal")], ilab=TRUE, zcolor = "style", 
           box=FALSE, ilcs=1.5, sncs=1.2, cexil=10, cexsn=10)


venn::venn(DEA_pneu_list[["limma_UPDN"]][c("pneuVir_Heal_UP", "pneuVir_Heal_DN", 
                                      "pneuBac_Heal_UP", "pneuBac_Heal_DN")], 
           box=FALSE, ilcs=1.5, sncs=1.5, ilab=TRUE, zcolor = "style", cexil=5, cexsn=10)

venn::venn(DEA_pneu_list[["limma_UPDN"]][c("nCoV_pneuVir_UP", "nCoV_pneuVir_DN", 
                                      "nCoV_pneuBac_UP", "nCoV_pneuBac_DN")], 
           box=FALSE, ilcs=1.5, sncs=1.5, ilab=TRUE, zcolor = "style", cexil=5, cexsn=10)

save(S_raw, DEA_pneu_list, sample_meta, file="../results/DEA_pneu_list.RData")


