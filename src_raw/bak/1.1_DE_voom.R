

dge <- DGEList(counts=S_raw)

keep <- filterByExpr(subsample_dge, Expdesign)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
