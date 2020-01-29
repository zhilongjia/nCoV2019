

load("../results/DEA_list.RData")

nCoV_Heal_DE <- DEA_list[["limma_DE"]][["nCoV_Heal"]]
nCoV_Heal_pheno <- DEA_list[["subsample_pheno"]][["nCoV_Heal"]] 


sampleLabel <- nCoV_Heal_pheno$Group
names(sampleLabel) <- nCoV_Heal_pheno$sampleID
################################################################################
# devtools::load_all("~/data/cogena")
library(cogena)

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
# annoGMT <- "c5.all.v6.1.symbols.gmt"; nClust =3; clMethods=c("kmeans")
annoGMT <- "CmapDn100.gmt.xz"; nClust =5; clMethods=c("kmeans")
annoGMT <- "CmapUp100.gmt.xz"; nClust =3; clMethods=c("kmeans")


annofile <- system.file("extdata", annoGMT, package="cogena")


cogena_result <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel, ncore=ncore)

summary(cogena_result)

heatmapCluster(cogena_result, "kmeans", "7", maintitle="nCoV_Heal", add2=FALSE,
               cexCol=1.2 )

heatmapPEI(cogena_result, "k", "3", maintitle="nCoV_Heal", add2=FALSE,
           CutoffNumGeneset=20)

heatmapPEI(cogena_result, "k", "3", maintitle="nCoV_Heal", add2=FALSE,
           CutoffNumGeneset=20, orderMethod = "1", printGS = TRUE)
heatmapPEI(cogena_result, "k", "5", maintitle="nCoV_Heal", add2=FALSE,
           CutoffNumGeneset=20, orderMethod = "5", printGS = TRUE)


heatmapPEI(cogena_result, "k", "5", maintitle="nCoV_Heal", add2=FALSE,
           CutoffNumGeneset=20, orderMethod = "2", printGS = TRUE)

save.image("../results/2.1_cogena_KEGG.RData")
save.image("../results/2.1_cogena_CmapUP.RData")
save.image("../results/2.1_cogena_CmapDN.RData")
