

CmapUp_raw <- cogena::gmt2list("../data/CmapUp100.gmt")


drug_DEA_list <- list()
# saquinavir@MCF7#5.2E-06M_6246
drug_DEA_list[["Upgene_sauinavir"]] <- CmapUp_raw[[4910]]

# ribavirin@PC3#1.64E-05M_7316
drug_DEA_list[["Upgene_ribavirin"]] <- CmapUp_raw[[5858]]


load("../results/DEA_list.RData")
load("../results/DEA_pneu_list.RData")


drug_DEA_list[["nCoV_Heal_DN"]] <- DEA_list[["limma_UPDN"]][["nCoV_Heal_DN"]]


library(venn)

venn::venn(drug_DEA_list, ilab=TRUE, zcolor = "style", 
           box=FALSE, ilcs=1.5, sncs=1.2, cexil=10, cexsn=10)


drug_gene_venn_res <- gplots::venn(drug_DEA_list, show.plot=FALSE)

aa <- attr(drug_gene_venn_res,"intersections")
aa_df <- as.data.frame(sapply(aa, "length<-", max(lengths(aa))))

readr::write_tsv(aa_df, path="../results/drug_gene_venn_res.tsv", na="")

save.image("../results/3.1_drug_pathway.RData")
