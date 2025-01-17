
load("../data/genelist_filter_lz.RData")
names(gene.list)[2] <- "Vir_Heal"

load("../results/sample_count_df_symbol.RData")
################################################################################
# Normlization
library(limma)
library(edgeR)

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
                    Vir_Heal=c("Healthy", "Viral-like"), 
                    Others_Heal=c("Healthy", "Others"), 
                    # pneuBac_pneuVir=c("Others", "Viral-like"),
                    # nCoV_pneuBac=c("Others", "nCoV"),
                    nCoV_Vir=c("Viral-like", "nCoV") )



DEA_pneu_list <- list()

for (comp_var in names(comp2groups) ) {
    # comp_var <- names(comp2groups)[[1]]
    comp_i <- comp2groups[[comp_var]]
    print (comp_i)
    subsample_pheno <- dplyr::filter(sample_meta, Types %in% comp_i )
    subsample_pheno$Types <- factor(subsample_pheno$Types, levels=comp_i )
    Expdesign <- model.matrix(~subsample_pheno$Types)
    
    subsample_dge <- dge[,subsample_pheno$NewID]
    
    # DEA
    keep <- filterByExpr(subsample_dge, Expdesign)
    # use filter genes from Li Zhang (in gene.list)
    keep[names(keep)] <- FALSE
    keep[rownames(gene.list[[comp_var]])] <- TRUE
    
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
    readr::write_csv(limma_DEA, paste0("../results/limma/", comp_var, "_limma_res.csv") )
    
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
library(venn)

venn::venn(DEA_pneu_list[["limma_DEG"]][c("nCoV_Heal","Vir_Heal","Others_Heal")], 
           ilab=TRUE, zcolor = "style", 
           box=FALSE, ilcs=1.5, sncs=1.2, cexil=10, cexsn=10)

venn::venn(DEA_pneu_list[["limma_UPDN"]][c("nCoV_Heal_UP", "nCoV_Heal_DN", 
                                      "Vir_Heal_UP", "Vir_Heal_DN")], 
           box=FALSE, ilcs=1.5, sncs=1.2, ilab=TRUE, zcolor = "style", cexil=5, cexsn=10)

# output veen res
venn_res <- gplots::venn(DEA_pneu_list[["limma_DEG"]][c("nCoV_Heal","Vir_Heal","Others_Heal")])
aa <- attr(venn_res,"intersections")
aa_df <- as.data.frame(sapply(aa, "length<-", max(lengths(aa))))
readr::write_csv(aa_df, path="../results/venn/DEG_venn_res.csv", na="")

DEA_list <- DEA_pneu_list

pal <- list()
pal$d <- scales::hue_pal()(4)
pal$s <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)

save(pal, S_raw, DEA_list, sample_meta, file="../results/DEA_pneu_list.RData")
