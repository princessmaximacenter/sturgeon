#!/usr/bin/env Rscript
#live processing with guppy R10

suppressMessages(library(argparser))
suppressMessages(library(docstring))

# Source R script to be able to call functions from it
source_dir <- "/opt/sturgeon/R_scripts/"
source(paste0(source_dir, "/peroperative_seq_functions_shorter.R"))

main <- function(args = NULL){
  args <- getArgs(args)                                                       # Collect User Args
  print("Starting live processing")
  if(file.exists(paste0(args$wrapper, "/wrapper_stop.txt"))){
    system(paste0("rm -f ", args$wrapper, "/wrapper_running.txt ", args$wrapper, "/wrapper_stop.txt"))
  } else {
    if (file.exists(paste0(args$wrapper, "/wrapper_running.txt"))) {
      stop(paste0("ERROR: The ", args$wrapper, "/wrapper_running.txt file exists. This could mean sturgeon is ",
      "already running. Only one sturgeon can be running at the time! ",
      "If this is not the case delete this file and try again."))
    }
  }

  iteration <- 1                                                              # Counter to decide what to run when
  system(paste("mkdir -p", args$output))                                      # Make output dir
  running <- T                                        # TODO: find more elegant way to keep running and kill the script
  progress <- data.frame()
  setwd(args$output)
  
  # Create progress file
  file.create("processed_files.txt")
  write("start", file = "processed_files.txt", append = T)
  
  while(running == T) {

    if(file.exists(paste0(args$wrapper, "/wrapper_stop.txt"))){
      running = FALSE
      break
    }

    # Keep track of processed files
    processed_fast5s <- read.table("processed_files.txt")
    made_fast5s <- list.files(args$input)

    # Check if there are new files, if not, wait 30 seconds and check again
    new_fast5s <- made_fast5s[!made_fast5s%in%processed_fast5s$V1]
    if(length(new_fast5s) == 0) {
      print("FLAG: No new files, waiting 30 seconds...")
      Sys.sleep(30)
      next    # Ignores the rest of the code and enters the next loop
    }

    #find out if they are ready for processing ie. if it contains 4000 reads
    fast5_cpy <- get_fast5_copy(new_fast5s, args$input, args$sizePod5)

    #if it is not complete, wait 30 seconds and start over.
    if(length(fast5_cpy) == 0) {
      print("No complete files found to copy, waiting 30 seconds...")
      Sys.sleep(30)
      next
    }
    for(f5 in fast5_cpy){
      if(file.exists(paste0(args$wrapper, "/wrapper_stop.txt"))){
      running = FALSE
      break
    }
      #process files
      print(paste0("FLAG: starting processing of iteration_", iteration))
      progress <- process_file(f5, iteration, progress, args$output, 
                               args$input, !args$useClassifiedBarcode, args$barcode, args$wrapper)
      if(iteration %% args$cnvFreq == 0){
        plot_cnv(args$output, iteration)
      }
      print(paste0("FLAG: iteration_", iteration, " completed!"))
      iteration = iteration + 1
    }
  }
  print("Live processing has stopped.")
}

getArgs <- function(args = NULL){
  #' Collect Arguments
  #' 
  #' Function to capture and store user arguments in an object
  #' @param args vector of strings with each string being a key and value space separated. 
  #' Keep null to collect args from command line
  #' @returns Object with all keys and their values
  p <- argparser::arg_parser(name = "Sturgeon: Live Processing (mode: Guppy)",
                             description = "Live processing pod5 files into bam files and figures. 
                             Expects guppy for variant calling and pod5 python package to be installed on the system")
  p <- argparser::add_argument(parser = p, arg = "--input", 
                               help = "Folder where the pod5 files are written to by minKnow.",
                               default = "/home/docker/input/")
  p <- argparser::add_argument(parser = p, arg = "--output", 
                               help = "Folder where the result files are written to. Folder does not have to exist yet",
                               default = "/home/docker/output/")
  p <- argparser::add_argument(parser = p, arg = "--wrapper",
                               help = "Path to the location where the wrapper_running.txt is written to. 
                               When a megalodon instance is running, this file will be created and removed once megalodon is done. 
                               If this file already exists, this script will wait until it is gone before starting a new megalodon, to prevent a system crash. 
                               Remove it manually, if no other version of this script is running or else it will wait forever.",
                               default = "/home/docker/wrapper/")
  p <- argparser::add_argument(parser = p, arg = "--barcode",
                               help = "Barcode used in the library preparation",
                               type = "numeric", default = 18)
  p <- argparser::add_argument(parser = p, arg = "--useClassifiedBarcode",
                               help = "Use classified barcodes", flag = TRUE)
  p <- argparser::add_argument(parser = p, arg = "--cnvFreq",
                               help = "The number of iterations before merging bams and creating a new CNV plot",
                               default = 5, type = "numeric")
  p <- argparser::add_argument(parser = p, arg = "--sizePod5",
                               help = "Number of reads in the pod5 files. (Recommended to use default)",
                               default = 4000, type = "numeric")
  return(if (is.null(args)) argparser::parse_args(p) else argparser::parse_args(p, args))
}

get_fast5_copy <- function(new_fast5s, input_folder, nr_reads){
  #' Collect copy ready fast5/pod5 files
  #' 
  #' Function to check if the new fast5 files are ready to be processed based on the number of reads present
  #' @param new_fast5s list of new fast5/pod5 file names
  #' @param input_folder The location of the fast5/pod5 files
  #' @param nr_reads Number of expected reads when completed
  #' @returns A list of file names that are ready to be processed
  fast5_cpy <- c()
  for(fast5 in new_fast5s){
    tryCatch(len <- system(paste0("pod5 view ", input_folder, "/", fast5, " | wc -l"), intern=T), 
             error=function(x){ len <- 0 })
    len <- as.numeric(len)-1
    
    if(len == nr_reads){ fast5_cpy <- c(fast5_cpy, fast5) }
    if(len != 0 & len < nr_reads){ print(paste0("incomplete file found: ", fast5)) }
  }
  return(fast5_cpy)
}

process_file <- function(f5, iteration, progress, output_folder, 
                         input_folder, include_unclassified, barcode, wrapper){
  #' Process a Fast5/Pod5 file
  #' 
  #' Function to call guppy and perform analyses on the output on the given file
  #' @param f5 File name of file to process
  #' @param iteration Number of the current iteration
  #' @param progress table with processed files
  #' @param output_folder folder to write to
  #' @param input_folder folder with the input files
  #' @param include_unclassified Boolean whether to use unclassified barcodes or not
  #' @param barcode Use barcode
  #' @param wrapper file location of the wrapper flag file to make sure only one script is running at the time
  #' @returns Updated progress table
  #make sure no wrapper is running already, this will overload the GPU RAM
  holdup <- T
  while(holdup){
    if(file.exists(paste0(wrapper, "/wrapper_running.txt"))){ Sys.sleep(1) } else{ holdup <- F }
  }
  #make a file preventing someone else from starting a megalodon run
  system(paste0("touch ", wrapper, "/wrapper_running.txt"))
  
  print(paste("starting wrapper for file:", f5))
  #wrapper starts guppy and extracts reads from the correct barcode
  wrappert_guppy_R10_guppy6.5(output_folder, paste0(input_folder,"/", f5), iteration = iteration, 
                              bcoverride = barcode, include_unclassified = include_unclassified)
  system(paste0("rm ", wrapper, "/wrapper_running.txt"))
  
  write(f5, file="processed_files.txt", append=T)
  
  #add iteration results to the progress dataframe, and make a plot
  progress <- plot_process(progress, output_folder, iteration)
  Sys.sleep(5)
  
  return(progress)
}

plot_process <- function(progress, output_folder, iteration){
  progress <- add_and_plot(progress, output_folder)
  write.table(progress, file=paste0(output_folder, "/iteration_", iteration, "/classifier_progress_iteration_", iteration, ".tsv"), quote=F, row.names = F, sep = "\t")

  dev.new()
  print(paste0("FLAG: making confidence over time plot for iteration_", iteration))
  confidence_over_time_plot(progress)

  dev.copy2pdf(file=paste0(output_folder, "/iteration_", iteration, "/confidence_over_time_plot_iteration_", iteration, ".pdf"))
  dev.copy(png, file=paste0(output_folder, "/iteration_", iteration, "/confidence_over_time_plot_iteration_", iteration, ".png"))
  dev.off()
  return(progress)
}

plot_cnv <- function(output_folder, iteration){
  dev.new()
  print(paste0("FLAG: creating CNV plot for iteration_", iteration))
  plot_cnv_from_live_dir_DNAcopy(output_folder)
  dev.copy2pdf(file=paste0(output_folder, "/iteration_", iteration, "/CNV_plot_iteration_", iteration, ".pdf"))
  dev.copy(png, file=paste0(output_folder, "/iteration_", iteration, "/CNV_plot_iteration_", iteration, ".png"))
  dev.off()
}

main()
