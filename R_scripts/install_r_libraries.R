################
# Installation #
# settings     #
################

# Cran mirror to get R packages
cran_mirror <- "https://cran.r-project.org/"

# The proportion of available cpu to be used for installing R packages
# 1 means use all available CPUs
proportion_cpu_used <- 1

## list of libraries/packages and version requirements
packages <-  list(yaml = "2.3.9",
                  argparser = "0.7.2",
                  docstring = "1.0.0",
                  rstudioapi = "0.16.0")

# Biocoductor libraries are version bound to the BiocManager version
# so we only need BiocManager version and names of libraries to install
bioconductor_version <- "3.19"
bio_packages <- c("QDNAseq", "DNAcopy")


################
# Installation #
################

# Number of threads/cpu (round up)
ncpus_for_installation <- ceiling(parallel::detectCores() * proportion_cpu_used)

# Install devtools to allow version specific installation
install.packages("devtools", dependencies = TRUE, repos = cran_mirror,
                 Ncpus = ncpus_for_installation)

# Install CRAN libraries/packages
sapply(names(packages), FUN = function(pkg_name) {
  if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
    devtools::install_version(
      pkg_name, version = packages[pkg_name],
      dependencies = TRUE, repos = cran_mirror,
      threads = ncpus_for_installation
    )
  }
})

# Install Biocoductor libraries/packages
lapply(bio_packages, FUN = function(pkg) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = cran_mirror,
                     Ncpus = ncpus_for_installation)
  BiocManager::install(pkg, version = bioconductor_version,
                       Ncpus = ncpus_for_installation)
})

# remove devtools/biocmanager from the final image
remove.packages(c("devtools", "BiocManager"))


################
# Validate     #
# installation #
################

for (pkg in c(names(packages), bio_packages)) {
  print(pkg)
  stopifnot(pkg %in% installed.packages()[,'Package'])
}