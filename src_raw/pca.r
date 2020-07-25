




library(factoextra)
library(FactoMineR)
load("../results/nCoV_pneu_Heal_norm_symbol_GE.RData")
load("../results/DEA_pneu_list.RData")

dd=nCoV_pneu_Heal_norm_symbol_GE[DEA_list$limma_DEG$nCoV_Hea,rownames(map)]
#dd=vsd_q[filter(nCoV_hea,abs(logFC)>3)[,1] %>% as.vector(),rownames(ann_col)[-107]]
pp <- PCA(t(dd %>% na.omit), scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(pp,col.ind =map$Group,
             mean.point=F,geom.ind = "point",
             palette = "jco",pointsize=3,pointshape=16,
             axes = c(1,3),#PC1 PC2 PC3
             legend.title = "Groups",title="")

load("../results/nCoV_pneu_Heal_norm_symbol_GE.RData")
load("../results/DEA_pneu_list.RData")

library(ggfortify)

DEG_nCoV_Heal <- rownames(DEA_list[["limma_DE"]][["nCoV_Heal"]])
DEG_pneuVir_Heal <- rownames(DEA_list[["limma_DE"]][["Vir_Heal"]])
DEG_pneuBac_Heal <- rownames(DEA_list[["limma_DE"]][["Others_Heal"]])

DEG_united <- unique(c(DEG_nCoV_Heal, DEG_pneuVir_Heal, DEG_pneuBac_Heal))
autoplot(prcomp(t(nCoV_pneu_Heal_norm_symbol_GE[intersect(DEG_united, rownames(nCoV_pneu_Heal_norm_symbol_GE) ),]), scale=F), 
         data=sample_meta, colour = "Types", 
         size = 5, label = F, label.colour="black", ts.colour="black" )
