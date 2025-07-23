suppressMessages(library(QDNAseq))
suppressMessages(library(DNAcopy))

plot_cnv_from_bam_DNAcopy <- function(bam, output_file = NULL, utils_path= NULL, makeplot = TRUE,
                                      lines_only = FALSE,
                                      binsize = 1e6){
  #' Plot CNV
  #'
  #' PLOT CNVs from the given bam file
  #' Change this to a config getter
  SOURCE_DIR = utils_path
  cnvplot_tmpfile <- paste0(SOURCE_DIR, "/bins_sample_counts_tmp.bed")

  # Collect binned read counts and store it in a bed file
  readCounts <- QDNAseq::binReadCounts(readRDS(paste0(SOURCE_DIR, "/chm13_v2.0_bins.rds")), bamfiles = bam,
                                       isNotPassingQualityControls = NA,
                                       minMapq = 2, isDuplicate = NA, isSecondaryAlignment = NA)
  QDNAseq::exportBins(readCounts, file = cnvplot_tmpfile, format = "bed", filter = FALSE, logTransform = FALSE)

  # Prep reference file
  reference <- read.table(paste0(SOURCE_DIR, "/bins_export.bed"), skip = 1, stringsAsFactors = FALSE)
  colnames(reference) <- c("chrom", "start", "end", "name", "coverage", "orientation")
  # Prep sample file

  sample <- read.table(cnvplot_tmpfile, skip = 1, stringsAsFactors = FALSE)
  system(paste("rm", cnvplot_tmpfile))
  colnames(sample) <- c("chrom", "start", "end", "name", "coverage", "orientation")

  sum_reads_title <- sum(sample$coverage)                     # To be used in a plot title downstream

  sample$coverage <- as.numeric(sample$coverage) + 0.001      # Add to coverage to ensure none are absolute 0
  sample$refcov <- reference$coverage
  sample <- sample[sample$chrom != "Y" & sample$chrom != "X", ]
  sample$reltoref <- sample$coverage / sample$refcov
  exclusionListBins <- sample[sample$refcov < 500, ]
  exclusionListBins <- row.names(exclusionListBins)
  sample <- sample[!row.names(sample) %in% exclusionListBins, ]
  mn <- mean(sample$reltoref)

  sample$logtr <- log(sample$reltoref / mn, 2)
  sample$logtr <- ifelse(sample$logtr > 3, 3, ifelse(sample$logtr<(-3), -3, sample$logtr))
  sample$chromname <- ifelse(nchar(sample$chrom) == 1, paste0("chr0", sample$chrom), paste0("chr", sample$chrom))

  bam_cna <- CNA(sample$logtr, chrom = sample$chromname, maploc = sample$start + binsize,
                 data.type = "logratio", sampleid = "cnvplot")
  smoothed.CNA.object <- smooth.CNA(bam_cna)
  segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
  pointcolorframe <- data.frame(pointx = 1:length(bam_cna$cnvplot), pointy = bam_cna$cnvplot, pointcolor = "black")

  for(x in 1:nrow(segment.smoothed.CNA.object$segRows)){
    if(segment.smoothed.CNA.object$output[x, "seg.mean"]>0.5){
      pointcolorframe[segment.smoothed.CNA.object$segRows[x,1]:segment.smoothed.CNA.object$segRows[x,2],"pointcolor"]="green"
    }
    if(segment.smoothed.CNA.object$output[x, "seg.mean"]<(-0.5)){
      pointcolorframe[segment.smoothed.CNA.object$segRows[x,1]:segment.smoothed.CNA.object$segRows[x,2],"pointcolor"]="blue"
    }    }
  if(makeplot){
    pdf(file=paste0(output_file, ".pdf"))
    a <- dev.cur()
    png(file=paste0(output_file, ".png"), width = 1200, height = 1200, res = 120 )
    dev.control('enable')

    plot("", xlim=c(0,length(bam_cna$cnvplot)), ylim=c(-3,3), main=paste0("CNVs nreads=",sum_reads_title), xaxt="n", xlab="genomic pos",
         ylab="log ratio")
    if(lines_only==F){
      points(x=pointcolorframe$pointx, y=pointcolorframe$pointy, pch=16, cex=0.5, col=pointcolorframe$pointcolor)
    }
    for(x in 1:nrow(segment.smoothed.CNA.object$segRows)){
      segments(x0=  segment.smoothed.CNA.object$segRows[x,1], x1=segment.smoothed.CNA.object$segRows[x,2],
               y0=segment.smoothed.CNA.object$output[x,"seg.mean"], y1=segment.smoothed.CNA.object$output[x,"seg.mean"],
               col="red", lwd=2)
    }
    fr = data.frame(chrom=smoothed.CNA.object$chrom, loc=smoothed.CNA.object$maploc, val=smoothed.CNA.object$cnvplot)
    nr = 1
    for(i in 1:(nrow(fr)-1)){
      if(fr$chrom[i]!=fr$chrom[i+1]){
        abline(v=as.numeric(row.names(fr)[i])+0.5, col="grey")
        chrname= gsub(pattern = "chr", replacement = "", x = fr[i,"chrom"])
        ypos = ifelse(nr %% 2 == 1, 2, 1.8)
        nr=nr+1
        text(x = i, adj=c(1,0), labels = chrname, y = ypos)
      }
    }
    if(lines_only==F){
      relevant_genes = read.table(paste0(SOURCE_DIR, "/relevant_genes_with_chm13v2_1Mb_bin_nrs.bed"), header=T,
                                  stringsAsFactors = FALSE)
      fr$start=fr$loc-binsize
      for(i in 1:nrow(relevant_genes)){
        genestart = relevant_genes[i,"start"]
        genechrom = relevant_genes[i,"chrom"]
        genechrom = strsplit(genechrom, split="hr")[[1]][2]
        genechrom=ifelse(nchar(genechrom)==1, paste0("chr0", genechrom), paste0("chr", genechrom))
        binline = fr[fr$chrom==genechrom&fr$start<genestart&(fr$start+binsize*2)>genestart,]
        ycoord = binline$val
        text(x=as.numeric(row.names(binline)), y=ycoord-0.1+relevant_genes[i,"yoffset"], cex=0.5, srt=-90,adj=c(0.5,1), labels = relevant_genes[i, "name"], pos=1)
        points(x=as.numeric(row.names(binline)), y=ycoord, col="red")
      }}

    dev.copy(which=a)
    dev.off()
    dev.off()
  }
  output = list(pointdata=pointcolorframe, segdata=segment.smoothed.CNA.object)
  return(output)
}