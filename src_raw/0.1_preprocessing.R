

# sample_meta <- readr::read_csv("../data/sample_meta-v1.1.csv") # ALL
# sample_meta <- readr::read_csv("../data/sample_meta-v2.csv") # SE
# sample_meta <- readr::read_csv("../data/sample_meta-v3.csv") # PE
sample_meta <- readxl::read_xlsx("../data/sample_meta_v8.xlsx") # SE

sampID_newID <- sample_meta$NewID
names(sampID_newID) <- sample_meta$sampleID


# RNA-Seq data
library(dplyr)
gene_count_rpkm_fn <- dir("../data/final/", pattern="Counts", full.names=TRUE)

gene_count_list <- list()
for (fn_i in gene_count_rpkm_fn ) {
    # fn_i <- gene_count_rpkm_fn[1]
    gene_count_rpkm <- readr::read_tsv(fn_i) %>% 
        dplyr::group_by(Geneid) %>%
        summarise_all(funs(sum))
    
    sample_id <- sub("_Counts.txt","",basename(fn_i) )
    
    #  use NewID as sample names
    if (!is.na(sampID_newID[sample_id]) ) {
        colnames(gene_count_rpkm)[2] <- sampID_newID[sample_id]
        gene_count_list[[sample_id]] <- gene_count_rpkm
    }
    
}

sample_count_df <- Reduce(dplyr::inner_join, gene_count_list )


# combine PE and SE samples
# ref: https://www.biostars.org/p/252552/
sample_count_df$C1 <- sample_count_df$C1_1 + sample_count_df$C1
sample_count_df$C2 <- sample_count_df$C2_1 + sample_count_df$C2
sample_count_df$C3 <- sample_count_df$C3_1 + sample_count_df$C3
sample_count_df$C4 <- sample_count_df$C4_1 + sample_count_df$C4

sample_count_df <- dplyr::select(sample_count_df, -one_of(c("C1_1", "C2_1", "C3_1", "C4_1")) )


readr::write_tsv(sample_count_df, "../results/all_genecounts.tsv")

################################################################################


save.image("../results/0.1_preprocessing_stage1.RData")

################################################################################
# sample grouped
nCoV_samp <- dplyr::filter(sample_meta, Diagnostics=="nCoV pneumonia")  %>% 
    dplyr::select(NewID) %>% unlist(use.names = FALSE)  
Heal_samp <- dplyr::filter(sample_meta, Diagnostics=="Healthy")  %>% 
    dplyr::select(NewID) %>% unlist(use.names = FALSE)  
pneu_samp <- dplyr::filter(sample_meta, Diagnostics=="CAP")  %>% 
    dplyr::select(NewID) %>% unlist(use.names = FALSE)  

pneuVir_samp <- dplyr::filter(sample_meta, Types=="Viral-like")  %>% 
    dplyr::select(NewID) %>% unlist(use.names = FALSE)  

pneuBac_samp <- dplyr::filter(sample_meta, Types=="Others")  %>% 
    dplyr::select(NewID) %>% unlist(use.names = FALSE)  


################################################################################
# all non-zeros genes (to filter)
# sample_count_mat <- as.data.frame(sample_count_df[,c(nCoV_samp, Heal_samp, pneu_samp)])
# rownames(sample_count_mat) <- sample_count_df$Geneid
# nlq_ensembl <- names(which(rowSums(sample_count_mat)!=0))

# 56.39% of input gene IDs are fail to map
symbol_ENTREZID <- clusterProfiler::bitr(gsub("\\..+", "", sample_count_df$Geneid ), 
                                         fromType="ENSEMBL", toType=c("ENSEMBL", "SYMBOL"), OrgDb="org.Hs.eg.db", drop=TRUE)

# get symbol gene expression
sample_count_df_symbol <- dplyr::mutate(sample_count_df, ENSEMBL=gsub("\\..+", "", Geneid )) %>% 
    dplyr::left_join(symbol_ENTREZID ) %>% 
    dplyr::filter(!is.na(SYMBOL) ) %>% 
    dplyr::select("SYMBOL",nCoV_samp, Heal_samp, pneu_samp) %>%
    dplyr::group_by(SYMBOL) %>% 
    dplyr::summarise_all(sum)



sample_meta <- dplyr::filter(sample_meta, Diagnostics %in% c("nCoV pneumonia",
                                                             "Healthy",
                                                             "CAP") )


save.image("../results/0.1_preprocessing.RData")


save(sample_count_df_symbol, sample_meta, 
     nCoV_samp, Heal_samp, pneu_samp, pneuVir_samp, pneuBac_samp,
     file="../results/sample_count_df_symbol.RData")


