

load("../results/DEA_list.RData")


################################################################################
# filter lower logfc genes
nCoV_logFC <- dplyr::filter(DEA_list[["limma_DEA"]][["nCoV_Heal"]], abs(logFC)>=4) 


nCoV_Heal_DE <- DEA_list[["limma_DE"]][["nCoV_Heal"]][nCoV_logFC$nCoV_Heal,]
nCoV_Heal_pheno <- DEA_list[["subsample_pheno"]][["nCoV_Heal"]] 


sampleLabel <- nCoV_Heal_pheno$Group1
names(sampleLabel) <- nCoV_Heal_pheno$sampleID
################################################################################
devtools::load_all("../../cogena")
# library(cogena)

nClust <- 2:10
ncore <- 7
clMethods <- c("hierarchical","kmeans","diana","fanny","som","sota","pam","clara","agnes")
clMethods <- c("hierarchical","kmeans","pam")

genecl_result <- coExp(nCoV_Heal_DE, nClust=nClust, 
                       clMethods=clMethods, 
                       metric="correlation", 
                       method="complete", 
                       ncore=ncore, 
                       verbose=TRUE)


###########################################
##########################################
annoGMT <- "c2.cp.kegg.v7.0.symbols.gmt"


annofile <- system.file("extdata", annoGMT, package="cogena")


cogena_result <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel, ncore=ncore)

summary(cogena_result)

heatmapCluster(cogena_result, "h", "3", maintitle="nCoV_Heal", add2=FALSE,
               cexCol=1.2 )

heatmapPEI(cogena_result, "h", "3", maintitle="nCoV_Heal", add2=FALSE,
           CutoffNumGeneset=20)


save.image("../results/2.1_cogena_logFC3_KEGG.RData")

#########################################################################
annoGMT <- "CmapUp100.gmt.xz";
annofile <- system.file("extdata", annoGMT, package="cogena")
drugUP_result <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel, ncore=ncore)

heatmapPEI(drugUP_result, "h", "3", maintitle="nCoV_Heal", add2=TRUE,
           CutoffNumGeneset=20, orderMethod = "2", printGS = TRUE)

heatmapPEI(drugUP_result, "h", "3", maintitle="nCoV_Heal", add2=TRUE,
           CutoffNumGeneset=20, orderMethod = "3", printGS = TRUE)

heatmapPEI(drugUP_result, "h", "3", maintitle="nCoV_Heal", add2=TRUE,
           CutoffNumGeneset=20, orderMethod = "Down", printGS = TRUE)

save.image("../results/2.1_cogena_CmapUP.RData")


#########################################################################
annoGMT <- "CmapDn100.gmt.xz"; 
annofile <- system.file("extdata", annoGMT, package="cogena")
drugDN_result <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel, ncore=ncore)

heatmapPEI(drugDN_result, "h", "3", maintitle="nCoV_Heal", add2=TRUE,
           CutoffNumGeneset=20, orderMethod = "1", printGS = TRUE)

heatmapPEI(drugDN_result, "h", "3", maintitle="nCoV_Heal", add2=TRUE,
           CutoffNumGeneset=20, orderMethod = "Up", printGS = TRUE)

save.image("../results/2.1_cogena_CmapDN.RData")

#########################################################################

