
load("../results/DEA_list.RData")
source("~/Dropbox/zzz.R")

################################################################################
# Pathway & GO analysis of DEG at each timepoint
library(clusterProfiler)
genesEntrezID_3g <- sapply(DEA_list[["limma_UPDN"]], symbol2entrezID )

genesEntrezID_3g_KEGG <- compareCluster(genesEntrezID_3g, fun='enrichKEGG')
dotplot(genesEntrezID_3g_KEGG, showCategory=30)

genesEntrezID_3g_GO <- compareCluster(genesEntrezID_3g, fun='enrichGO', OrgDb='org.Hs.eg.db')
dotplot(genesEntrezID_3g_GO, showCategory=30)

save.image("../results/2.0_KEGG.RData")
