
library(dplyr)
gene_count_rpkm_fn <- dir("../data/final/", pattern="Counts", full.names=TRUE)

gene_count_list <- list()
for (fn_i in gene_count_rpkm_fn ) {
    gene_count_rpkm <- readr::read_tsv(fn_i) %>% 
        dplyr::group_by(Geneid) %>%
        summarise_all(funs(sum))
    
    sample_id <- sub("_Counts.txt","",basename(fn_i) )
    colnames(gene_count_rpkm)[2] <- sample_id
    
    gene_count_list[[sample_id]] <- gene_count_rpkm
}

sample_count_df <- Reduce(dplyr::inner_join, gene_count_list )

readr::write_tsv(sample_count_df, "../results/all_genecounts.tsv")

sample_meta <- readr::read_csv("../data/sample_meta.csv")


save.image("../results/0.1_preprocessing_stage1.RData")

################################################################################
# sample grouped
nCoV_samp <- dplyr::filter(sample_meta, Group=="nCoV")  %>% 
    select(sampleID) %>% unlist(use.names = FALSE)  
Heal_samp <- dplyr::filter(sample_meta, Group=="Healthy")  %>% 
    select(sampleID) %>% unlist(use.names = FALSE)  
pneu_samp <- dplyr::filter(sample_meta, Group=="pneumonia")  %>% 
    select(sampleID) %>% unlist(use.names = FALSE)  


################################################################################
# all non-zeros genes (to filter)
sample_count_mat <- as.data.frame(sample_count_df[,c(nCoV_samp, Heal_samp, pneu_samp)])
rownames(sample_count_mat) <- sample_count_df$Geneid
nlq_ensembl <- names(which(rowSums(sample_count_mat)!=0))

# 56.39% of input gene IDs are fail to map
symbol_ENTREZID <- clusterProfiler::bitr(gsub("\\..+", "", rownames(sample_count_mat) ), 
                                         fromType="ENSEMBL", toType=c("ENSEMBL", "SYMBOL"), OrgDb="org.Hs.eg.db", drop=TRUE)

# get symbol gene expression
sample_count_df_symbol <- dplyr::select(sample_count_df, c("Geneid", nCoV_samp, Heal_samp, pneu_samp)) %>% 
    dplyr::filter(Geneid %in% nlq_ensembl) %>% 
    dplyr::mutate(ENSEMBL=gsub("\\..+", "", Geneid )) %>% 
    dplyr::left_join(symbol_ENTREZID ) %>% 
    dplyr::filter(!is.na(SYMBOL) ) %>% 
    dplyr::select("SYMBOL",nCoV_samp, Heal_samp, pneu_samp) %>% 
    dplyr::group_by(SYMBOL) %>% 
    dplyr::summarise_all(sum)

save.image("../results/0.1_preprocessing.RData")





