#!/usr/bin/env Rscript
# peroperative sequencing functions

suppressMessages(library(QDNAseq))
suppressMessages(library(DNAcopy))
suppressMessages(library(yaml))
suppressMessages(library(argparser))

SOURCE_DIR <- "/opt/sturgeon/R_scripts/"
# SOURCE_DIR <- "/opt/docker/R_scripts/"
# SOURCE_DIR <- "/Users/a.janse-3/Documents/sturgeon/R_scripts/"
CONFIG <- yaml::read_yaml(paste0(SOURCE_DIR, "/config.yml"))$peroperative_seq_functions_shorter
COLOR_TRANSLATION <- readRDS(paste0(SOURCE_DIR, "/color_translation.rds"))   #rds file containing a color list

# wrappert <- function(main_folder, fast5, iteration, bcoverride=F){
#   #this wrapper runs megalodon on a fast5 file
#   #then uses qcat to split out the most frequent barcode
#   #then runs sturgeon
#   #fast5 = "~/nanocns/data/example_dataII/FAR96725_b5b5d4b3_0.fast5"
#   out_folder = paste0(main_folder,"/iteration_",iteration)
#
#   system(paste0("mkdir ",out_folder))
#   system(paste0("cp ", fast5," ", out_folder))
#
#   mega_command = paste0("/home/sturgeon/miniconda3/envs/megalodon/bin/megalodon ", out_folder, "  --outputs mods basecalls mappings --mappings-format bam ",
#                         "--reference ", CONFIG$refgenome, " --write-mods-text --mod-motif m CG 0 --processes 10 ",
#                         "--guppy-server-path /home/sturgeon/nanocns/software/guppy_v5/bin/guppy_basecall_server ",
#                         "--guppy-params \"-d /home/sturgeon/nanocns/software/rerio/basecall_models/ --num_callers 2 --ipc_threads 3\"",
#                         " --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg --devices cuda:0 --overwrite --output-directory ", out_folder,"/meg_out/ ",
#                         "--suppress-progress-bars --suppress-queues-status")
#
#   print(mega_command)
#   system(mega_command, ignore.stdout = T, ignore.stderr=T)
#
#   system(paste0("mv ", out_folder,"/meg_out/* ", out_folder))
#   #run qcat barcoding
#
#   system(paste0("/home/sturgeon/miniconda3/bin/qcat --tsv -k RBK004 -f ",
#                 out_folder,"/basecalls.fastq > ",out_folder,"/barcoding.tsv"))
#   #select reads from the most common barcode
#   bcds = read.table(paste0(out_folder,"/barcoding.tsv"),header=T)
#   if(bcoverride==F){
#
#     mostfreq = names(sort(table(bcds$barcode), decreasing = T, na.last = T)[1])
#     print(paste("using most frequent barcode, barcode nr",mostfreq))
#     bcds = bcds[bcds$barcode==mostfreq,"name"]
#     print(paste("total reads:",length(bcds)))}
#
#   if(bcoverride!=F){
#     print(paste("using override, barcode nr",bcoverride))
#     bcds = bcds[bcds$barcode==bcoverride,"name"]
#     print(paste("total reads:",length(bcds)))
#   }
#
#
#   #filter read calls
#   methfile= read.table(paste0(out_folder,"/per_read_modified_base_calls.txt"),header=T)
#   system(paste0("mv ",out_folder,"/per_read_modified_base_calls.txt ",out_folder,"/meg_out/"))
#   methfile=methfile[methfile$read_id%in%bcds,]
#   write.table(methfile, file=paste0(main_folder,"/per_read_modified_base_calls_it",iteration,".txt"),quote=F, row.names=F, sep="\t")
# }

fix_name=function(x){
  return(gsub(gsub(x = x, pattern = "\\.\\.\\.", replacement = " - "), pattern = "\\.", replacement = " "))}

add_and_plot=function(merged_data, output_folder, iteration){
  result_file <- paste0(output_folder, "/iteration_", iteration,"/merged_probes_methyl_calls_", CONFIG$modelname, "_iteration_",iteration,".csv")
  if(exists("COLOR_TRANSLATION")==F){COLOR_TRANSLATION=readRDS(paste0(SOURCE_DIR, "/color_translation.rds"))}
  reslist_2 = t(data.frame(read.table(result_file, sep=",", stringsAsFactors = FALSE)))
  reslist_2=reslist_2[2:nrow(reslist_2),]
  colnames(reslist_2)=c("class", "score")
  tm = Sys.time()
  tm= paste(format(as.POSIXct(tm), format = "%H:%M"))
  reslist_2=rbind(reslist_2,data.frame(class="TIME", score=tm))

  if(nrow(merged_data)==0){
    return(reslist_2)
  }else{

    mgd = merge(merged_data, reslist_2, by="class")
    colnames(mgd)[ncol(mgd)]=paste0("iteration_", ncol(mgd)-1)

    # plot('', xlim=c(0,ifelse(ncol(mgd)>5, ncol(mgd)+5, 10)), ylim=c(0,1), main="confidence over time", xlab="iteration", ylab="confidence",
    #      xaxt="n")
    # for(i in 1:(nrow(mgd)-1)){
    #   clr = unlist(unname(COLOR_TRANSLATION[mgd[i,"class"]]))
    #   lines(x = 0:(ncol(mgd)-2), y= as.numeric(mgd[i,2:ncol(mgd)]), col=clr)
    #   if(as.numeric(mgd[i,ncol(mgd)])>0.5){text(x=ncol(mgd), y=as.numeric(mgd[i,ncol(mgd)])-0.02, labels = mgd[i,"class"])}
    # }
    # maxtics= ifelse(ncol(mgd)>5,ncol(mgd), ncol(mgd)+5)
    # axis(side = 1, at = 0:maxtics, labels = 1:(maxtics+1))
    # axis(side = 1, at = 0:(ncol(mgd)-2), labels = mgd[mgd$class=="TIME",2:ncol(mgd)], line=1, tick=F, cex.axis=0.7)
    # abline(h=0.95, col="red")
    # abline(h=0.80, col="orange")
    return(mgd)
  }
}

confidence_over_time_plot=function(merged_data, output_file, color_translation= COLOR_TRANSLATION){
  pdf(file=paste0(output_file, ".pdf"))
  a<-dev.cur()
  png(file=paste0(output_file, ".png"))
  dev.control("enable")

  mgd = merged_data
  colnames(mgd)[ncol(mgd)]=paste0("iteration_", ncol(mgd)-1)
  plot('', xlim=c(0,ifelse(ncol(mgd)>5, ncol(mgd)+5, 10)), ylim=c(0,1), main="confidence over time", xlab="iteration", ylab="confidence",
       xaxt="n")
  for(i in 1:(nrow(mgd)-1)){
    tumor_class <- gsub("/", " ", mgd[i,"class"])
    clr = unlist(unname(color_translation[tumor_class]))
    lines(x = 0:(ncol(mgd)-2), y= as.numeric(mgd[i,2:ncol(mgd)]), col=clr)
    if(as.numeric(mgd[i,ncol(mgd)])>0.5){text(x=ncol(mgd), y=as.numeric(mgd[i,ncol(mgd)])-0.02, labels = mgd[i,"class"])}
  }
  maxtics= ifelse(ncol(mgd)>5,ncol(mgd), ncol(mgd)+5)
  axis(side = 1, at = 0:maxtics, labels = 1:(maxtics+1))
  axis(side = 1, at = 0:(ncol(mgd)-2), labels = mgd[mgd$class=="TIME",2:ncol(mgd)], line=1, tick=F, cex.axis=0.7)
  abline(h=0.95, col="red")
  abline(h=0.80, col="orange")

  dev.copy(which=a)
  dev.off()
  dev.off()
}

#still contains hardcoded paths
wrappert_guppy_R10_guppy6.5 <- function(main_folder, fast5, iteration, bcoverride=F, include_unclassified=F){
  #this wrapper runs guppy on a fast5 file

  out_folder = paste0(main_folder, "/iteration_", iteration)

  system(paste("mkdir", out_folder))
  system(paste("cp", fast5, out_folder))

  guppy_command = paste0("guppy_basecaller -i ", out_folder, " -s ", out_folder, "/guppy_out ",
                         "-c ", CONFIG$guppy, "/data/dna_r10.4.1_e8.2_400bps_5khz_modbases_5mc_cg_hac_mk1c.cfg --device auto --bam_out ",
                         "--disable_qscore_filtering --min_score_barcode_front 6 ",
                         "--align_ref ", CONFIG$refgenome, " --barcode_kits SQK-RBK114-24 ",
                         "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4 --allow_inferior_barcodes")
  print(paste0("FLAG: starting guppy for iteration_", iteration))
  print(paste0("DEBUG: guppy command: ", guppy_command))
  system(guppy_command)

  barcode = ifelse(nchar(bcoverride)==1, paste0("0", bcoverride), bcoverride)
  system(paste0("mv ", out_folder,"/guppy_out/barcode",barcode,"/*.bam ",main_folder,"/guppy_output_iteration_",iteration,".bam"))

  #include unclassified reads yes or no?
  if(include_unclassified==T){
    system(paste0("mv ", out_folder,"/guppy_out/unclassified/*.bam ",main_folder,"/guppy_output_iteration_",iteration,"_unclassified.bam"))}

  print(paste0("FLAG: starting sturgeon for iteration_", iteration))
  system(paste(paste0(SOURCE_DIR, "/sturgeon_guppy.sh"), main_folder, paste0(SOURCE_DIR, CONFIG$probes), CONFIG$model))

  system(paste0("cp ", main_folder, "/merged_probes_methyl_calls_", CONFIG$modelname,".pdf ", out_folder,
                "/merged_probes_methyl_calls_", CONFIG$modelname, "_iteration_", iteration, ".pdf"))

  system(paste0("cp ", main_folder, "/merged_probes_methyl_calls_", CONFIG$modelname,".png ", out_folder,
                "/merged_probes_methyl_calls_", CONFIG$modelname, "_iteration_", iteration, ".png"))

  system(paste0("cp ", main_folder, "/merged_probes_methyl_calls_", CONFIG$modelname,".csv ", out_folder,
                "/merged_probes_methyl_calls_", CONFIG$modelname, "_iteration_", iteration, ".csv"))

}

plot_cnv_from_bam_DNAcopy <- function(bam, output_file = NULL, makeplot = TRUE,
                                      lines_only = FALSE,
                                      binsize = 1e6){
  #' Plot CNV
  #'
  #' PLOT CNVs from the given bam file

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
  blacklistbins <- sample[sample$refcov < 500, ]
  blacklistbins <- row.names(blacklistbins)
  sample <- sample[!row.names(sample) %in% blacklistbins, ]
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
    png(file=paste0(output_file, ".png"))
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

plot_cnv_from_live_dir_DNAcopy=function(directory, output_file){
  #' Merge BAMs and plot CNV
  #'
  #' Merges all BAM files and call function to plot the CNVs
  #' @param directory Directory to write to
  # Create output directory if it doesn't exist else cleanup previous merged bam file
  if(!dir.exists(paste0(directory,"/merged_bams/"))) {
    system(paste0("mkdir ", directory, "/merged_bams/"))
  } else {
    system(paste0("rm ", directory, "/merged_bams/merged_bam.bam"))
  }
  # Merge all bam files
  system(paste0("samtools merge -@ 10 -O BAM -o ",directory,"/merged_bams/merged_bam.bam ", directory,"/*.bam" ))
  # Make CNV plot based on the merged bam file
  output <- plot_cnv_from_bam_DNAcopy(paste0(directory,"/merged_bams/merged_bam.bam"), output_file)
}

