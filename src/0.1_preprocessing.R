
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


# combine PE and SE samples
# ref: https://www.biostars.org/p/252552/
sample_count_df$`x19001-R` <- sample_count_df$`x19001-R` + sample_count_df$`Paried_x19001-R`
sample_count_df$`x19002-R` <- sample_count_df$`x19002-R` + sample_count_df$`Paried_x19002-R`
sample_count_df$`x19003-R` <- sample_count_df$`x19003-R` + sample_count_df$`Paried_x19003-R`
sample_count_df$`x19004-R` <- sample_count_df$`x19004-R` + sample_count_df$`Paried_x19004-R`


readr::write_tsv(sample_count_df, "../results/all_genecounts.tsv")

# sample_meta <- readr::read_csv("../data/sample_meta-v1.1.csv") # ALL
sample_meta <- readr::read_csv("../data/sample_meta-v2.csv") # SE
# sample_meta <- readr::read_csv("../data/sample_meta-v3.csv") # PE


save.image("../results/0.1_preprocessing_stage1.RData")

################################################################################
# sample grouped
nCoV_samp <- dplyr::filter(sample_meta, Group1=="nCoV")  %>% 
    dplyr::select(sampleID) %>% unlist(use.names = FALSE)  
Heal_samp <- dplyr::filter(sample_meta, Group1=="Healthy")  %>% 
    dplyr::select(sampleID) %>% unlist(use.names = FALSE)  
pneu_samp <- dplyr::filter(sample_meta, Group1=="pneumonia")  %>% 
    dplyr::select(sampleID) %>% unlist(use.names = FALSE)  

pneuVir_samp <- dplyr::filter(sample_meta, Group3=="pneumonia-vir")  %>% 
    dplyr::select(sampleID) %>% unlist(use.names = FALSE)  

pneuBac_samp <- dplyr::filter(sample_meta, Group3=="pneumonia-bac")  %>% 
    dplyr::select(sampleID) %>% unlist(use.names = FALSE)  


################################################################################
# all non-zeros genes (to filter)
# sample_count_mat <- as.data.frame(sample_count_df[,c(nCoV_samp, Heal_samp, pneu_samp)])
# rownames(sample_count_mat) <- sample_count_df$Geneid
# nlq_ensembl <- names(which(rowSums(sample_count_mat)!=0))

# 56.39% of input gene IDs are fail to map
symbol_ENTREZID <- clusterProfiler::bitr(gsub("\\..+", "", sample_count_df$Geneid ), 
                                         fromType="ENSEMBL", toType=c("ENSEMBL", "SYMBOL"), OrgDb="org.Hs.eg.db", drop=TRUE)

# get symbol gene expression
sample_count_df_symbol <- dplyr::select(sample_count_df, c("Geneid", all_of(c(nCoV_samp, Heal_samp, pneu_samp)))) %>% 
    # dplyr::filter(Geneid %in% nlq_ensembl) %>% 
    dplyr::mutate(ENSEMBL=gsub("\\..+", "", Geneid )) %>% 
    dplyr::left_join(symbol_ENTREZID ) %>% 
    dplyr::filter(!is.na(SYMBOL) ) %>% 
    dplyr::select("SYMBOL",nCoV_samp, Heal_samp, pneu_samp) %>%
    dplyr::group_by(SYMBOL) %>% 
    dplyr::summarise_all(sum)

save.image("../results/0.1_preprocessing.RData")


save(sample_count_df_symbol, sample_meta, file="../results/sample_count_df_symbol.RData")


