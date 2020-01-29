
load("../results/DEA_list_q0.001.RData")
source("~/Dropbox/zzz.R")

################################################################################
# Convert gene symbols to probes in HGU133a.
symbol2Probe <- function(gs){
    library(hgu133a.db)
    p <- AnnotationDbi::select(hgu133a.db, gs, "PROBEID", "SYMBOL")$PROBEID
    p <- unique(p[which(!is.na(p))])
    print (length(p))
    p
}


limma_UPDN_df <- sapply(DEA_list$limma_UPDN, "length<-", max(lengths(DEA_list$limma_UPDN)))
write.table(limma_UPDN_df, "../results/drp/limma_UPDN_df.tsv", sep="\t", na="",row.names = FALSE,  quote=FALSE )

# for CMap
limma_UPDN_probe <- sapply(DEA_list$limma_UPDN, symbol2Probe)

for (i in names(limma_UPDN_probe)) {
    write.table(limma_UPDN_probe[[i]], file=paste0("../results/drp/", i, ".grp"), quote=F, col.names = F, row.names = F)
}

