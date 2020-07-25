
filterByExpr <- function(y, ...)
    UseMethod("filterByExpr")

filterByExpr.DGEList <- function(y, design=NULL, group=NULL, lib.size=NULL, ...)
{
    
    filterByExpr.default(y$counts, design=design, group=group, lib.size=lib.size, ...)
}

filterByExpr.default <- function(y, design=NULL, group=NULL, lib.size=NULL, min.count=10, min.total.count=15, large.n=10, min.prop=0.7, ...)
    #	Filter low expressed genes given count matrix
    #	Computes TRUE/FALSE index vector indicating which rows to keep
    #	Gordon Smyth
    #	Created 13 Nov 2017. Last revised 4 August 2019.
{
    if(is.null(design) && is.null(group)) {
        design <- y$design
        if(is.null(design)) group <- y$samples$group
    }
    if(is.null(lib.size)) lib.size <- y$samples$lib.size * y$samples$norm.factors
    y <- y$counts
    
    y <- as.matrix(y)
    if(mode(y) != "numeric") stop("y is not a numeric matrix")
    if(is.null(lib.size)) lib.size <- colSums(y)
    if(is.null(design) && is.null(group)) group <- rep_len(1L,ncol(y))
    
    #	Minimum effect sample sample size for any of the coefficients
    if(is.null(group)) {
        h <- hat(design)
        MinSampleSize <- 1/max(h)
    } else {
        group <- as.factor(group)
        n <- tabulate(group)
        MinSampleSize <- min(n[n > 0L])
    }
    if(MinSampleSize > large.n) MinSampleSize <- large.n + (MinSampleSize-large.n)*min.prop
    
    #	CPM cutoff
    MedianLibSize <- median(lib.size)
    CPM.Cutoff <- min.count/MedianLibSize*1e6
    CPM <- cpm(y,lib.size=lib.size)
    tol <- 1e-14
    keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - tol)
    
    #	Total count cutoff
    keep.TotalCount <- (rowSums(y) >= min.total.count - tol)
    
    keep.CPM & keep.TotalCount
}


################################################################################
y <- subsample_dge
design <- Expdesign
group <- subsample_pheno[,c("NewID", "Types")]

filterByGroup(y, design = NULL, group = NULL, min.count = 10, 
              min.total.count = 15, large.n = 10, min.prop = 0.7 ) {
    
    tmp <- y$counts[,which(group == levels(group)[1])]
    
    min_group <- 
        
        
        gene_sample_group_count <- tibble::rownames_to_column(as.data.frame(y$counts)) %>% 
        tidyr::gather("NewID", "counts", 2:(ncol(y$counts)+1) ) %>% 
        dplyr::left_join(group)
    
    nonZero=function(x){length(which(x!=0)) }
    Num_pergroup <- function(x) {
        table(subsample_pheno$Types)[x]
    }
    
    
    gene_group_count_stat <- dplyr::group_by(gene_sample_group_count, rowname, Types)  %>% 
        dplyr::summarise(countALL=sum(counts), CountNonZero=nonZero(counts)  ) %>% 
        dplyr::mutate(NumperGroup=Num_pergroup(Types), pent=CountNonZero/NumperGroup )
    
    filter_gene <- dplyr::filter(gene_group_count_stat, pent > 0.5,  countALL>15 )
    
    
    if(is.null(lib.size)) lib.size <- y$samples$lib.size * y$samples$norm.factors
    
}

