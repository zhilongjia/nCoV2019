


library(GEOquery)


SARS_geo_raw <- getGEO("GSE1739", destdir="../data/")



SARS_gpl_raw <- getGEO("GSE1739", destdir="../data/")


gset <- SARS_geo_raw$GSE1739_series_matrix.txt.gz

SARS_lablel <- pData(gset)

# log2 transform
ex <- exprs(gset)
if (any(is.na(ex)) | nrow(ex)==0 | ncol(ex)==0 ) {print ("NA in gene expresion."); next }
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
print(LogC)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex<- exprs(gset) <- log2(ex)
}


