
load("../results/DEA_list.RData")
load("../results/DEA_pneu_list.RData")

symbol2entrezID <- function(gene_symbols) {
    symbol_ENTREZID <- clusterProfiler::bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    return(symbol_ENTREZID$ENTREZID)
}


DEA_list <- c(DEA_list[["limma_DEG"]], DEA_pneu_list[["limma_DEG"]])

c1 <- c("nCoV_Heal", "pneu_Heal", "nCoV_pneu", 
        "pneuVir_Heal", "pneuBac_Heal",  
        "pneuBac_pneuVir", "nCoV_pneuBac", "nCoV_pneuVir")

c1 <- c("nCoV_Heal", "pneu_Heal", "pneuVir_Heal")


################################################################################
# Pathway & GO analysis of DEG
library(clusterProfiler)
genesEntrezID_3g <- sapply(DEA_list[c1], symbol2entrezID )
sapply(genesEntrezID_3g, length)

genesEntrezID_3g_KEGG <- compareCluster(genesEntrezID_3g, fun='enrichKEGG')
dotplot(genesEntrezID_3g_KEGG, showCategory=10)

genesEntrezID_3g_BP <- compareCluster(genesEntrezID_3g, fun='enrichGO', OrgDb='org.Hs.eg.db', ont="BP")
dotplot(genesEntrezID_3g_BP, showCategory=20)


save.image("../results/2.0_KEGG_DEG.RData")
