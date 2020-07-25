
library(dplyr)
gene_count_rpkm_fn <- dir("../data/final/",
                          pattern="Counts", full.names=TRUE)

gene_count_list <- list()
for (fn_i in gene_count_rpkm_fn ) {
    gene_count_rpkm <- readr::read_tsv(fn_i) %>% 
        # dplyr::select(Symbol, expected_count) %>% 
        dplyr::group_by(Geneid) %>%
        summarise_all(funs(sum))
    
    sample_id <- sub("_Counts.txt","",basename(fn_i) )
    colnames(gene_count_rpkm)[2] <- sample_id
    
    gene_count_list[[sample_id]] <- gene_count_rpkm
}

sample_count_df <- Reduce(dplyr::inner_join, gene_count_list )


readr::write_tsv(sample_count_df, "../results/all_genecounts.tsv")

sample_meta <- readr::read_csv("../data/sample_meta.csv")

################################################################################
#SE and PE evaluation

samp_PE_SE <- c("Paried_x19001-R", "Paried_x19002-R", "Paried_x19003-R", "Paried_x19004-R", "x19001-R", "x19002-R", "x19003-R", "x19004-R")

rep_G <- as.data.frame(sample_count_df[, samp_PE_SE])
rownames(rep_G) <- sample_count_df$Geneid
corrplot::corrplot.mixed(cor(rep_G, method="spearman"), order="hclust")
# corrplot::corrplot.mixed(cor(rep_G, method="pearson"), order="hclust")

################################################################################
# sample gorouped
nCoV_samp <- dplyr::filter(sample_meta, Group=="nCoV")  %>% 
    select(sampleID) %>% unlist(use.names = FALSE)  
Heal_samp <- dplyr::filter(sample_meta, Group=="Healthy")  %>% 
    select(sampleID) %>% unlist(use.names = FALSE)  
pneu_samp <- dplyr::filter(sample_meta, Group=="pneumonia")  %>% 
    select(sampleID) %>% unlist(use.names = FALSE)  

################################################################################
boxplot(sample_count_df[,c(nCoV_samp)])



# sample_count_mat <- sample_count_mat[names(which(rowSums(sample_count_mat)!=0)),]

# 54.92% of input gene IDs are fail to map
symbol_ENTREZID <- clusterProfiler::bitr(gsub("\\..+", "", rownames(sample_count_mat) ), 
                                         fromType="ENSEMBL", toType=c("ENSEMBL", "SYMBOL"), OrgDb="org.Hs.eg.db", drop=TRUE)


sample_count_df_symbol <- dplyr::select(sample_count_df, c("Geneid", nCoV_samp, Heal_samp, pneu_samp)) %>% 
    dplyr::filter(Geneid %in% nlq_ensembl) %>% 
    dplyr::mutate(ENSEMBL=gsub("\\..+", "", Geneid )) %>% 
    dplyr::left_join(symbol_ENTREZID ) %>% 
    dplyr::filter(!is.na(SYMBOL) ) %>% 
    dplyr::select("SYMBOL",nCoV_samp, Heal_samp, pneu_samp) %>% 
    dplyr::group_by(SYMBOL) %>% 
    dplyr::summarise_all(sum)

save.image("../results/1.1_DE.RData")

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
    
    limma_DEG <- dplyr::filter(limma_res, adj.P.Val<=0.05, abs(logFC)>=2 )
    readr::write_tsv(limma_DEG, paste0("../results/", comp_var, "_limma_res.tsv") )
    
    limma_DEG <- limma_DEG[[comp_var]]
    
    print(length(limma_DEG))
    limma_DE <- v$E[limma_DEG,]
    
    
    DEA_list[["limma_res"]][[comp_var]] <- limma_res
    DEA_list[["limma_DEG"]][[comp_var]] <- limma_DEG
    DEA_list[["limma_DE"]][[comp_var]] <- limma_DE
    DEA_list[["subsample_pheno"]][[comp_var]] <- subsample_pheno
    
}


save.image("../results/1.2_DE_stage2.RData")


################################################################################
gplots::heatmap.2(nCoV_Heal_DE, trace = "none", col="greenred", scale="row",
                  srtCol=60, keysize=1, cexCol=1.3, cexRow=1.1, Colv=TRUE, labRow=FALSE,
                  ylab=paste(nrow(nCoV_Heal_DE), "DE Genes"), margins = c(10, 2))


library(ggfortify)
autoplot(prcomp(t(DEA_list[["limma_DE"]][["nCoV_pneu"]]), scale=F), 
         data=DEA_list[["subsample_pheno"]][["nCoV_pneu"]], colour = "Group", 
         size = 8, label = TRUE, label.colour="black", ts.colour="black" )

library(venn)
venn(DEA_list[["limma_DEG"]], ilab=TRUE, zcolor = "style", cexil=10, cexsn=10)


DEG_list_raw <- readr::read_csv("../results/DEG_list.csv")

DEG_list <- lapply(DEG_list_raw, function(x) {tmp = na.omit(x); as.vector(tmp)} )


DEG_list_venn <- DEG_list[c("pneu_Heal_UP", "pneu_Heal_DN", "nCoV_Heal_UP", "nCoV_Heal_DN")]
venn(DEG_list_venn, ilab=TRUE, zcolor = "style", cexil=10, cexsn=10)

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

# DEA
keep <- filterByExpr(subsample_dge, Expdesign)
subsample_dge <- subsample_dge[keep,,keep.lib.sizes=FALSE]
subsample_dge <- calcNormFactors(subsample_dge)

v <- voom(subsample_dge, Expdesign, plot=FALSE, normalize="quantile")


DEG <- rownames(DEA_list[["limma_DE"]][["nCoV_Heal"]])

library(ggfortify)
    autoplot(prcomp(t(v$E[DEG,]), scale=T), 
         data=subsample_pheno, colour = "Group", 
         size = 8, label = TRUE, label.colour="black", ts.colour="black" )

