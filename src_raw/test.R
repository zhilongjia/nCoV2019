
# gene symbols mapping test
library(org.Hs.eg.db)
ensembl_symbol1 <- dplyr::full_join( toTable(org.Hs.egENSEMBL), toTable(org.Hs.egSYMBOL) )[,c("ensembl_id", "symbol")]

readr::write_tsv(ensembl_symbol1, path="../results/ensembl_symbol1.tsv")

ensembl_symbol2 <- readr::read_tsv("../documents/降帅/ensembl_symbol_type.tsv", col_names=FALSE)
colnames(ensembl_symbol2) <- c("ENSEMBL", "SYMBOL1", "Type")



###
symbol_ENTREZID <- clusterProfiler::bitr(gsub("\\..+", "", rownames(sample_count_mat) ), 
                                         fromType="ENSEMBL", toType=c("ENSEMBL", "SYMBOL"), OrgDb="org.Hs.eg.db", drop=F)


ensembl_ids <- gsub("\\..+", "", rownames(sample_count_mat) ) %>% as.data.frame()
colnames(ensembl_ids) <- "Ensembl gene ID"

###
HGCN_raw <- readr::read_tsv("../data/HGCN_genesymbols.txt")
tmp0 <- dplyr::left_join(symbol_ENTREZID, ensembl_symbol2) %>% 
    dplyr::left_join(HGCN_raw, by=c("ENSEMBL"="Ensembl gene ID"))


tmp1 <- dplyr::filter(tmp0, is.na(SYMBOL) )
tmp3 <- dplyr::filter(tmp0, !is.na(SYMBOL) )

tmp4 <- dplyr::filter(tmp3, Type=="lncRNA")


tmp2 <- dplyr::left_join(ensembl_ids, HGCN_raw)


tmp5 <- dplyr::left_join(DEA_list[["limma_DEA"]][["nCoV_Heal"]], ensembl_symbol2, by=c("nCoV_Heal"="SYMBOL1")) %>% 
    dplyr::left_join(HGCN_raw, by=c("ENSEMBL"="Ensembl gene ID"))


tmp6 <- dplyr::left_join(DEA_list[["limma_DEA"]][["pneu_Heal"]], ensembl_symbol2, by=c("pneu_Heal"="SYMBOL1")) %>% 
    dplyr::left_join(HGCN_raw, by=c("ENSEMBL"="Ensembl gene ID"))

###
library("biomaRt")
mart = useMart('ensembl')
listDatasets(mart)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# attributes = listAttributes(ensembl)
ensembl_ids <- ensembl_ids$`Ensembl gene ID`
ensembl_symbol <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
                        filters = 'ensembl_gene_id', 
                        values = ensembl_ids, 
                        mart = ensembl)

c2 <- readr::read_tsv("~/Downloads/c2.cp.kegg.v7.0.symbols.tsv", col_names = FALSE)

c2_type <- dplyr::left_join(c2, ensembl_symbol2, by=c("X1"="SYMBOL1"))


################################################################################
# tmp <- colSums(sample_count_df[,-1])
# names(tmp) <- names(sample_count_df[,-1])
# View(tibble::rownames_to_column(as.data.frame(tmp)))
# 
# 

