

hisat2  -p 45  -t -x /pnas/limk_group/zhongjx/megablast/rnaseq/hg38_hisat2/hg38 -1 /pnas/limk_group/zhongjx/rna/clean_fq/$file  -2 /pnas/limk_group/zhongjx/rna/clean_fq/$file2 -S /pnas/limk_group/zhongjx/rna/tmp/${file%%.*}".sam" --summary-file /pnas/limk_group/zhongjx/rna/tmp/${file%%.*}"_hisat2.summary" "  >>job.pbs

featureCounts -T 45 -p -t exon -g gene_id  -a /pnas/limk_group/zhongjx/megablast/rnaseq/gencode.v32.annotation.gtf -o /pnas/limk_group/zhongjx/rna/counts/${file%%.*}"_featureCounts.txt" /pnas/limk_group/zhongjx/rna/bam/${file%%.*}".bam"" >>job.pbs

echo "cut -f1,7 /pnas/limk_group/zhongjx/rna/counts/${file%%.*}"_featureCounts.txt"|grep -v '^#'> /pnas/limk_group/zhongjx/rna/final/${file%%.*}"_Counts.txt" " >>job.pbs

