# dhsrmethtools/R/functions.R

# Import necessary functions
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table fread
#' @importFrom HDF5Array writeHDF5Array
#' @importFrom DSS DMLtest

#' dragen2bsseq
#' 
#' @param file A DRAGEN bed file
#' @param samplename basename(file)
#' @return a bsseq object
#' @examples
#' dragen2bsseq(file=system.file("extdata", "test.bed.gz", package = "dhsrmethtools"),samplename="test")
#' @export
dragen2bsseq <- function(file,samplename=NULL){
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  # check for bsseq and HDF5Array
  if (!requireNamespace("HDF5Array", quietly = TRUE)) {
    stop("Package 'HDF5Array' is required but not installed.")
  }
  if (!requireNamespace("bsseq", quietly = TRUE)) {
    stop("Package 'bsseq' is required but not installed.")
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed.")
  } 

  if (is.null(samplename))
    samplename <- gsub(".CX_report.txt(.gz)?", "", basename(file))

  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    data <- data.table::fread(cmd = paste("gunzip -c ", file, " | awk -F'\t' '$1~/^chr[0-9XYM]+$/ && $6==\"CG\"'"), header = FALSE, col.names=c("chr","pos","strand","meth","unmeth","context","seq"))
  } else {
    # Use awk to filter regular file
    data <- data.table::fread(cmd = paste("awk -F'\t' '$1~/^chr[0-9XYM]+$/ && $6==\"CG\"'"), header = FALSE, col.names=c("chr","pos","strand","meth","unmeth","context","seq"))
  }
  chromosomes <- paste0("chr",c(1:22,"X","Y"))
  data <- data[chr %in% chromosomes]
  data$cov <- data$meth + data$unmeth
  hdf5_M <- HDF5Array::writeHDF5Array(as.matrix(data$meth))
  hdf5_C <- HDF5Array::writeHDF5Array(as.matrix(data$cov))
  bsOut <- strandCollapse(bsseq::BSseq(M = hdf5_M,Cov = hdf5_C,gr=GRanges(data$chr,IRanges(data$pos,data$pos),strand=data$strand),sampleNames=samplename),shift=TRUE)
  return(bsOut)

}

#' bed2bsseq
#'
#' @param file A simplfied bed file with methylation: chr\\tstart\\tend\\tmethylation\\tcoverage
#' @param samplename basename(file)
#' @return a bsseq object
#' @examples
#' bed2bsseq(file=system.file("extdata", "test.bed.gz", package = "dhsrmethtools"),samplename="test")
#' @export
bed2bsseq <- function(file,samplename=basename(file),hdf5=TRUE){
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  # check for bsseq and HDF5Array
  if (!requireNamespace("HDF5Array", quietly = TRUE)) {
    stop("Package 'HDF5Array' is required but not installed.")
  }
  if (!requireNamespace("bsseq", quietly = TRUE)) {
    stop("Package 'bsseq' is required but not installed.")
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed.")
  }

  samplename <- gsub(".bed(.gz)?", "", samplename)
  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    data <- data.table::fread(cmd = paste("gunzip -c ", file, " | awk -F'\t' '$1~/^chr[0-9XYM]+$/'"), header = FALSE, col.names=c("chr","start","end","meth","cov"))
  } else {
    # Use awk to filter regular file
    data <- data.table::fread(cmd = paste("awk -F'\t' '$1~/^chr[0-9XYM]+$/' ", file), header = FALSE, col.names=c("chr","start","end","meth","cov"))
  }
  chromosomes <- paste0("chr",c(1:22,"X","Y"))
  data <- data[chr %in% chromosomes]
  data$meth <- data$meth * data$cov
  if (hdf5){
    hdf5_M <- HDF5Array::writeHDF5Array(as.matrix(data$meth))
    hdf5_C <- HDF5Array::writeHDF5Array(as.matrix(data$cov))
    bsOut <- bsseq::BSseq(M = hdf5_M,Cov = hdf5_C,gr=GRanges(data$chr,IRanges(data$end,data$end)),sampleNames=samplename)
  } else {
    bsOut <- bsseq::BSseq(M = as.matrix(data$meth),Cov = as.matrix(data$cov),gr=GRanges(data$chr,IRanges(data$end,data$end)),sampleNames=samplename)
  }
  return(bsOut)
}

#' modkit2bsseq
#'
#' @param file An ONT bedmethyl file
#' @param samplename basename(file)
#' @return a bsseq object
#' @examples
#' modkit2bsseq(file=system.file("extdata", "test.basemods.bedmethyl.gz", package = "dhsrmethtools"),samplename="test")
#' @export
modkit2bsseq <- function(file,samplename=basename(file)){
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  # check for bsseq and HDF5Array
  if (!requireNamespace("HDF5Array", quietly = TRUE)) {
    stop("Package 'HDF5Array' is required but not installed.")
  }
  if (!requireNamespace("bsseq", quietly = TRUE)) {
    stop("Package 'bsseq' is required but not installed.")
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed.")
  }

  samplename <- gsub(".basemods.bedmethyl(.gz)?", "", samplename)
  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    # include column names to match the ONT bedmethyl file from modkit
    data <- data.table::fread(cmd = paste("gunzip -c ", file, "| awk -F'\t' '$4 == \"m\"'"), header = FALSE, col.names=c("chr","start","end","type","score","strand","start2","end2","color","valid_cov","fraction_mod","n_mod","n_canonical","n_othermod","n_deleted","n_fail","n_diff","n_nocall"))
  } else {
    # Use awk to filter regular file
    data <- data.table::fread(cmd = paste("awk -F'\t' '$4 == \"m\"' ", file), header = FALSE, col.names=c("chr","start","end","type","score","strand","start2","end2","color","valid_cov","fraction_mod","n_mod","n_canonical","n_othermod","n_deleted","n_fail","n_diff","n_nocall"))
  }
  chromosomes <- paste0("chr",c(1:22,"X","Y"))
  data <- data[chr %in% chromosomes]
  hdf5_M <- HDF5Array::writeHDF5Array(as.matrix(data$n_mod))
  hdf5_C <- HDF5Array::writeHDF5Array(as.matrix(data$valid_cov))
  bsOut <- bsseq::BSseq(M = hdf5_M,Cov = hdf5_C,gr=GRanges(data$chr,IRanges(data$end,data$end)),sampleNames=samplename)
  return(bsOut)
}

#' pbcpg2bsseq
#'
#' @param file A pb-CpG tools bedmethyl file
#' @param samplename basename(file)
#' @return a bsseq object
#' @examples
#' pbcpg2bsseq(file=system.file("extdata", "test.pbcpg.bedmethyl.gz", package = "dhsrmethtools"),samplename="test")
#' @export
pbcpg2bsseq <- function(file,samplename=basename(file)){
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  # check for bsseq and HDF5Array
  if (!requireNamespace("HDF5Array", quietly = TRUE)) {
    stop("Package 'HDF5Array' is required but not installed.")
  }
  if (!requireNamespace("bsseq", quietly = TRUE)) {
    stop("Package 'bsseq' is required but not installed.")
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed.")
  }

  samplename <- gsub(".basemods.bedmethyl(.gz)?|.bed|.bed.gz|.5mC.bed.gz", "", samplename)
  data <- NULL
  # Check if the file is gzipped
  if (grepl("\\.gz$", file)) {
    # Use zcat and awk to filter gzipped file
    data <- data.table::fread(cmd = paste("gunzip -c ", file), header=FALSE,col.names=c("chr","start","end","total_meth","type","coverage","mod","unmod","eff_meth"))
  } else {
    # Use awk to filter regular file
    data <- data.table::fread(file, header=FALSE,col.names=c("chr","start","end","total_meth","type","coverage","mod","unmod","eff_meth"))
  }
  chromosomes <- paste0("chr",c(1:22,"X","Y"))
  data <- data[chr %in% chromosomes]
  hdf5_M <- HDF5Array::writeHDF5Array(data.matrix(data$mod))
  hdf5_C <- HDF5Array::writeHDF5Array(data.matrix(data$coverage))
  bsOut <- bsseq::BSseq(M = hdf5_M,Cov = hdf5_C,gr=GRanges(data$chr,IRanges(data$end,data$end)),sampleNames=samplename)
  return(bsOut)
}

#' hapBSseq
#'
#' Create a BSseq object from three input files (hap1, hap2, combined) based on the platform type.
#'
#' @param hap1_file Path to the haplotype 1 bedmethyl file.
#' @param hap2_file Path to the haplotype 2 bedmethyl file.
#' @param combined_file Path to the combined bedmethyl file.
#' @param sample_name Sample name. Default is derived from the combined file basename.
#' @param platform Platform type ('ont' or 'hifi').
#' @param results_path Base path to save the BSseq results. If NULL, the function returns the BSseq object.
#' @param replace Logical, whether to replace the existing directory. Default is FALSE.
#' @return A message indicating the success or failure of the operation.
#' @examples
#' hapBSseq(
#'   hap1_file = "path/to/hap1_file",
#'   hap2_file = "path/to/hap2_file",
#'   combined_file = "path/to/combined_file",
#'   sample_name = "sample1",
#'   platform = "ont",
#'   results_path = "path/to/results",
#'   replace = FALSE
#' )
#' @export
hapBSseq <- function(hap1_file, hap2_file, combined_file, sample_name = NULL, platform, results_path = NULL, replace = FALSE) {
  if (!platform %in% c("ont", "hifi")) {
    stop("Unknown platform type. Supported types are 'ont' and 'hifi'.")
  }
  
  # Extract sample name from combined file if not provided
  if (is.null(sample_name)) {
    sample_name <- sub("\\..*", "", basename(combined_file))
  }
  
  # Construct the full path to save the BSseq object
  if (!is.null(results_path)) {
    final_results_path <- file.path(results_path, sample_name)
    
    # Check if the directory already exists and handle according to the 'replace' parameter
    if (dir.exists(final_results_path)) {
      if (replace) {
        unlink(final_results_path, recursive = TRUE)
      } else {
        stop(paste("Directory", final_results_path, "already exists. Use 'replace = TRUE' to replace it. Its content will be lost!"))
      }
    }
  }
  
  bsList <- list()
  
  if (platform == "ont") {
    bsList[[paste0(sample_name, ".hap1")]] <- modkit2bsseq(hap1_file, paste0(sample_name, ".hap1"))
    bsList[[paste0(sample_name, ".hap2")]] <- modkit2bsseq(hap2_file, paste0(sample_name, ".hap2"))
    bsList[[paste0(sample_name, ".combined")]] <- modkit2bsseq(combined_file, paste0(sample_name, ".combined"))
  } else if (platform == "hifi") {
    bsList[[paste0(sample_name, ".hap1")]] <- pbcpg2bsseq(hap1_file, paste0(sample_name, ".hap1"))
    bsList[[paste0(sample_name, ".hap2")]] <- pbcpg2bsseq(hap2_file, paste0(sample_name, ".hap2"))
    bsList[[paste0(sample_name, ".combined")]] <- pbcpg2bsseq(combined_file, paste0(sample_name, ".combined"))
  }
  
  setAutoRealizationBackend(BACKEND = "HDF5Array")
  bs <- bsseq::combineList(bsList, BACKEND = "HDF5Array")
  
  if (is.null(results_path)) {
    return(bs)
  } else {
    saveHDF5SummarizedExperiment(bs, dir = final_results_path)
    return(paste("BSseq object saved successfully to:", final_results_path))
  }
}

#' callHapDmrs
#'
#' Perform differential methylation analysis on a BSseq object for two haplotypes in the same sample and save or return the DMRs.
#'
#' @param bs BSseq object or path to a BSseq HDF5Array directory.
#'           If a path is provided, the BSseq object will be loaded from the specified directory.
#' @param group1 Character vector, sample names for group 1. Default is samples with suffix ".hap1".
#' @param group2 Character vector, sample names for group 2. Default is samples with suffix ".hap2".
#' @param cpg_smoothing Logical, whether to perform smoothing. Default is TRUE.
#' @param ncores Integer, number of cores to use. Default is 1.
#' @param dmr_delta Numeric, delta value for callDMR function. Default is 0.35.
#' @param dmr_p_threshold Numeric, p-value threshold for callDMR function. Default is 0.05.
#' @param min.cpg.coverage Integer, minimum coverage for a CpG site to be considered. Default is 2.
#' @param min.cpgs Integer, minimum number of CpG sites in a region to be considered. Default is 5.
#' @param min.length Integer, minimum length of a DMR region to be considered. Default is 100.
#' @param output Character, path to save the output. If NULL, the function returns the DMR object.
#'               The file extension determines the output format:
#'               - "rda": Saves the DMRs as an R data file (.rda).
#'               - "bed": Saves the DMRs as a BED file (.bed) with no header, including columns "chr", "start", "end", "diff.Methy".
#'               - "tsv": Saves the DMRs as a tab-separated values file (.tsv) with a header.
#'               If the file extension is not one of these, the function stops and returns an error.
#' @return A message indicating the success or failure of the operation, or the DMR object if output is NULL.
#' @examples
#' callHapDmrs(
#'   bs = "path/to/bsseq_object",
#'   output = "path/to/output/DMRs.sample_name.hap.rda"
#' )
#' @export
callHapDmrs <- function(bs, 
                        group1 = sampleNames(bs)[grep("\\.hap1$", sampleNames(bs))], 
                        group2 = sampleNames(bs)[grep("\\.hap2$", sampleNames(bs))], 
                        cpg_smoothing = TRUE, 
                        ncores = 1, 
                        dmr_delta = 0.35, 
                        dmr_p_threshold = 0.05,
                        min.cpg.coverage = 2,
                        min.cpgs = 5,
                        min.length = 100,
                        output = NULL) {
  
  # Load BSseq object if a path is provided
  if (is.character(bs)) {
    bs <- HDF5Array::loadHDF5SummarizedExperiment(dir = bs)
  }
  
  # added a filter to get all positions with >=min.cpg.coverage reads. This is to
  # remove CpG positions with SNPs that eliminate a CpG. Parameterized in the function.
  bs <- bs[which(rowMins(getCoverage(bs,type="Cov",what="perBase"))>=min.cpg.coverage),]
  
  # Print log for group1 and group2
  message("Calling DML test for the following groups:")
  message("Group 1: ", paste(group1, collapse = ", "))
  message("Group 2: ", paste(group2, collapse = ", "))
  
  # Perform DML test
  dml_sm <- DMLtest(bs, group1 = group1, group2 = group2, smoothing = cpg_smoothing, ncores = ncores)
  
  # Call DMRs
  dmrs <- callDMR(dml_sm, delta = dmr_delta, p.threshold = dmr_p_threshold)

  # filter to retain hapDMRs with min.cpgs CpGs and min.length length (in bp)
  dmrs <- dmrs[which(dmrs$nCG >= min.cpgs & dmrs$length >= min.length),]
    
  # Save or return DMRs based on output
  if (!is.null(output)) {
    file_ext <- tools::file_ext(output)
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    switch(file_ext,
           rda = {
             save(dmrs, file = output)
             message(paste("DMRs saved as .rda to:", output))
           },
           bed = {
             write.table(dmrs[c("chr", "start", "end", "diff.Methy")], 
                         file = output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
             message(paste("DMRs saved as .bed to:", output))
           },
           tsv = {
             write.table(dmrs, file = output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
             message(paste("DMRs saved as .tsv to:", output))
           },
           stop("Unsupported file extension in output. Accepted formats are: rda, bed, tsv"))
  } else {
    return(dmrs)
  }
}

#' callDmrs
#'
#' Perform differential methylation analysis on a BSseq object for a sample compared to a normal and save or return the DMRs.
#'
#' @param sample_dir Directory containing the sample BSseq HDF5Array.
#' @param normal_dir Directory containing the normal BSseq HDF5Array.
#' @param sample_name Name of the sample column in the sample BSseq object. If NULL, it will be inferred assuming a .combined pattern.
#' @param normal_name Name of the normal column in the normal BSseq object. If NULL, it will be inferred assuming a .combined pattern.
#' @param cpg_smoothing Logical, whether to perform smoothing. Default is TRUE.
#' @param ncores Integer, number of cores to use. Default is 1.
#' @param dmr_delta Numeric, delta value for callDMR function. Default is 0.35.
#' @param dmr_p_threshold Numeric, p-value threshold for callDMR function. Default is 0.05.
#' @param min.cpg.coverage Integer, minimum coverage for a CpG site to be considered. Default is 2.
#' @param min.cpgs Integer, minimum number of CpG sites in a region to be considered. Default is 5.
#' @param min.length Integer, minimum length of a DMR region to be considered. Default is 100.
#' @param output Character, path to save the output. If NULL, the function returns the DMR object.
#'               The file extension determines the output format:
#'               - "rda": Saves the DMRs as an R data file (.rda).
#'               - "bed": Saves the DMRs as a BED file (.bed) with no header, including columns "chr", "start", "end", "diff.Methy".
#'               - "tsv": Saves the DMRs as a tab-separated values file (.tsv) with a header.
#'               If the file extension is not one of these, the function stops and returns an error.
#' @return A message indicating the success or failure of the operation, or the DMR object if output is NULL.
#' @examples
#' callDmrs(
#'   sample_dir = "path/to/sample_dir",
#'   normal_dir = "path/to/normal_dir",
#'   output = "path/to/output/DMRs.sample_name.vs.normal_name.rda"
#' )
#' @export
callDmrs <- function(sample_dir, 
                        normal_dir, 
                        sample_name = NULL, 
                        normal_name = NULL, 
                        cpg_smoothing = TRUE, 
                        ncores = 1, 
                        dmr_delta = 0.35, 
                        dmr_p_threshold = 0.05, 
                        min.cpg.coverage = 2,
                        min.cpgs = 5,
                        min.length = 100,
                        output = NULL) {
  # Check if the output directory exists
  if (!is.null(output)) {
    output_dir <- dirname(output)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      if (dir.exists(output_dir)) {
        message(paste("Created output directory:", output_dir))
      } else {
        stop(paste("Failed to create output directory:", output_dir))
      }
    }
  }
  
  if (!requireNamespace("HDF5Array", quietly = TRUE)) {
    stop("Package 'HDF5Array' is required but not installed.")
  }
  if (!requireNamespace("bsseq", quietly = TRUE)) {
    stop("Package 'bsseq' is required but not installed.")
  }
  if (!requireNamespace("DSS", quietly = TRUE)) {
    stop("Package 'DSS' is required but not installed.")
  }

  # Load the sample data
  sample <- HDF5Array::loadHDF5SummarizedExperiment(dir = sample_dir)
  
  if (is.null(sample_name)) {
    combined_names <- grep("\\.combined$", colnames(sample), value = TRUE)
    if (length(combined_names) != 1) {
      stop("No .combined pattern found in sample column names or multiple .combined patterns found.")
    }
    sample_name <- combined_names[1]
  } else {
    if (!sample_name %in% colnames(sample)) {
      stop("Provided sample_name does not exist in the sample data.")
    }
  }
  
  # Load the normal sample
  normal <- HDF5Array::loadHDF5SummarizedExperiment(dir = normal_dir)
  
  # Check for available normal names
  normal_names <- colnames(normal)
  if (is.null(normal_name)) {
    if (length(normal_names) != 1) {
      stop("Multiple or no sample names found in the normal sample. Please specify one.")
    }
    normal_name <- normal_names[1]
  } else {
    if (!normal_name %in% normal_names) {
      stop("Provided normal_name does not exist in the normal data.")
    }
  }
  
  # Create a list for bsseq objects
  bsList <- list()
  bsList[[sample_name]] <- sample[, sample_name]
  bsList[[normal_name]] <- normal[, normal_name]
  
  # Combine the bsseq objects
  bs_combined <- bsseq::combineList(bsList, BACKEND = "HDF5Array")
  
  # Filter for coverage
  bs_filtered <- bs_combined[which(rowMax(as.matrix(getCoverage(bs_combined, type = "Cov"))) > min.cpg.coverage), ]

  
  # Perform DML testing
  dml_results <- DSS::DMLtest(bs_filtered, group1 = sample_name, group2 = normal_name, smoothing = cpg_smoothing, ncores = ncores)
  
  # Call DMRs
  dmrs <- DSS::callDMR(dml_results, delta = dmr_delta, p.threshold = dmr_p_threshold)
  
  # filter to retain hapDMRs with min.cpgs CpGs and min.length length (in bp)
  dmrs <- dmrs[which(dmrs$nCG >= min.cpgs & dmrs$length >= min.length),]

  # Save or return DMRs based on output
  if (!is.null(output)) {
    file_ext <- tools::file_ext(output)
    switch(file_ext,
           rda = {
             save(dmrs, file = output)
             message(paste("DMRs saved as .rda to:", output))
           },
           bed = {
             write.table(dmrs[c("chr", "start", "end", "diff.Methy")], 
                         file = output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
             message(paste("DMRs saved as .bed to:", output))
           },
           tsv = {
             write.table(dmrs, file = output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
             message(paste("DMRs saved as .tsv to:", output))
           },
           stop("Unsupported file extension in output. Accepted formats are: rda, bed, tsv"))
  } else {
    return(dmrs)
  }
}


#' filterHapDmrs
#' Filter Haplotype-Specific Differentially Methylated Regions
#'
#' This function filters haplotype-specific differentially methylated regions (DMRs) based on various criteria.
#' It can compare methylation levels between two haplotypes and optionally against a normal sample.
#'
#' @param gr A GRanges object containing the regions to be filtered.
#' @param aml.bs A BSseq object containing methylation data for the AML sample.
#' @param normal.bs A BSseq object containing methylation data for the normal sample (optional).
#' @param normal.sample.name Character string specifying the name of the normal sample (optional).
#' @param vcffile Path to a VCF file for identifying divergent regions (optional).
#' @param svvcffile Path to SV VCF file for filtering out DMRs near SVs (optional).
#' @param sv.breakend.window bp to flank SV breakends for filtering DMRs (default: 1000).
#' @param min.diff Numeric value specifying the minimum methylation difference to consider (default: 0.4).
#' @param purity.tolerance Numeric value specifying the maximum allowable purity deviation (default: 0.4).
#' @param purity.tolerance.p.value One-sided P-value cutoff for binomial test of methylation counts against purity.tolerance or 1-purity.tolerance (default: 0.01).
#' @param excludeRegions A GRanges object specifying regions to exclude from the analysis (default: empty GRanges).
#' @param namehap1 Character string specifying the name of haplotype 1 in the aml.bs object (default: first sample name containing "hap1").
#' @param namehap2 Character string specifying the name of haplotype 2 in the aml.bs object (default: first sample name containing "hap2").
#' @param filter.sex.chromosomes Logical indicating whether to filter out sex chromosomes (default: TRUE).
#' @param run.fisher.test Logical indicating whether to run Fisher's exact test for differential methylation (default: FALSE).
#'
#' @return A filtered GRanges object containing the haplotype-specific DMRs that meet the specified criteria.
#'
#' @details
#' The function performs the following main steps:
#' 1. Checks input validity and prepares data.
#' 2. Optionally removes sex chromosomes.
#' 3. Optionally runs Fisher's exact test for differential methylation.
#' 4. Identifies and excludes divergent regions if a VCF file is provided.
#' 5. Filters regions based on various criteria (number of CpGs, length, coverage, methylation difference).
#' 6. If a normal sample is provided, further filters and classifies regions based on comparison to the normal sample.
#'
#' @examples
#' # Basic usage:
#' filtered_dmrs <- filterHapDmrs(gr, aml.bs)
#'
#' # With normal sample and VCF file:
#' filtered_dmrs <- filterHapDmrs(gr, aml.bs, normal.bs, normal.sample.name = "normal",
#'                                vcffile = "path/to/vcf.file", min.diff = 0.3)
#'
#' @importFrom methylKit calculateDiffMeth
#' @importFrom bsseq sampleNames getCoverage
#' @importFrom vcfR read.vcfR
#' @importFrom GenomicRanges GRanges
#'
#' @export
filterHapDmrs <- function(gr,aml.bs,normal.bs=NULL,normal.sample.name=NULL,vcffile=NULL,svvcffile=NULL,sv.breakend.window=1000,min.diff=.4,purity.tolerance=0.4,purity.tolerance.p.value=0.01,excludeRegions=GRanges(),namehap1=grep("hap1",sampleNames(aml.bs),value=T)[1],namehap2=grep("hap2",sampleNames(aml.bs),value=T)[1],filter.sex.chromosomes=TRUE,run.fisher.test=FALSE){
    require(methylKit)
    require(bsseq)
    require(vcfR)
    
    # if nlbmbs is not null check number of samples. If >1 then collapse
    if (!is.null(normal.bs) && ncol(normal.bs)>1){
        # exit with error that says 'nlbmbs must be a single sample'
        stop("normal.bs must be a single sample")
    } else if(!is.null(normal.bs) && ncol(normal.bs)==1){
        normal.sample.name <- sampleNames(normal.bs)[1]
    }

    # make sure namehap1 and namehap2 are in amlbs
    if (!(namehap1 %in% sampleNames(aml.bs)) || !(namehap2 %in% sampleNames(aml.bs))){
        # exit with error that says 'namehap1 and namehap2 must be in amlbs'
        stop("namehap1 and namehap2 must be in amlbs")
    }

    # function in methylKit package to do two sample fisher test on methylation counts
    bsseqDiffRegion <- function(bs,regions,group1,group2){   
        ids <- c(group1,group2)    
        cov <- getCoverage(bs[,c(group1,group2)],regions=regions,type="Cov",what = "perRegionTotal")
        colnames(cov) <- paste0("coverage",1:length(ids),sep="")    
        numCs <- getCoverage(bs[,c(group1,group2)],regions=regions,type="M",what = "perRegionTotal")
        colnames(numCs) <- paste0("numCs",1:length(ids),sep="")    
        numTs <- cov - numCs
        colnames(numTs) <- paste0("numTs",1:length(ids),sep="")    
        df <- as.data.frame(gr)[,c(1:3,5)]
        colnames(df) <- c("chr","start","end","strand")
        df <- cbind(df,cov[,1],numCs[,1],numTs[,1],cov[,2],numCs[,2],numTs[,2])  
        mb <- new("methylBase",as.data.frame(df),
                sample.ids=ids,
                coverage.index=grep("coverage",colnames(df)),
                numCs.index=grep("numCs",colnames(df)),
                numTs.index=grep("numTs",colnames(df)),
                assembly="hg38",context="CpG",
                destranded=TRUE,resolution="region",
                treatment=c(rep(0,length(group1)),rep(1,length(group2))))    
        mbdiff <- methylKit::calculateDiffMeth(mb)    
        return(mbdiff)    
    }
  
    # Function to import VCF and calculate heterozygous passing variants in 10kbp windows
    getDivergentRegions <- function(vcf_file) {
      # Import VCF file
      vcf <- read.vcfR(vcf_file, verbose=F)
      
      # Extract genotype information
      gt <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
      
      # Extract positions and filter passing variants
      pos <- getPOS(vcf)
      pass_filter <- which(grepl("PASS", vcf@fix[,7]))
      het_filter <- which(rowSums(gt == 1, na.rm = TRUE) > 0)  # heterozygous calls are represented by 1
      
      # Get positions of passing heterozygous variants
      het_pos <- pos[intersect(pass_filter, het_filter)]
      
      # Extract chromosome lengths from VCF header
      contigs <- vcf@meta[grepl("##contig", vcf@meta)]
      chrom_lengths <- sapply(contigs, function(contig) {
        contig <- gsub("##contig=<ID=|>", "", contig)
        parts <- strsplit(contig, ",length=")[[1]]
        chrom <- parts[1]
        length <- as.numeric(parts[2])
        return(c(chrom, length))
      })
      
      chrom_lengths <- as.data.frame(t(chrom_lengths))
      colnames(chrom_lengths) <- c("chrom", "length")
      chrom_lengths$length <- as.numeric(chrom_lengths$length)
      
      # Define 10kbp windows across the genome
      windows <- GRanges(seqnames = rep(chrom_lengths$chrom, chrom_lengths$length %/% 10000 + 1),
                         ranges = IRanges(start = unlist(lapply(chrom_lengths$length, function(x) seq(1, x, by = 10000))),
                                          width = 10000))
      
      # Convert heterozygous positions to GRanges
      het_gr <- GRanges(seqnames = vcf@fix[intersect(pass_filter, het_filter), 1],
                        ranges = IRanges(start = het_pos, width = 1))
      
      # Count heterozygous variants in each window
      counts <- countOverlaps(windows, het_gr)
      
      # Find 95th percentile of heterozygous SNP counts
      threshold <- quantile(counts, 0.95)
      
      # Filter windows with counts > 95th percentile
      high_het_windows <- windows[counts > threshold]
      
      return(high_het_windows)
    }

    # get SV breakpoints from VCF file and create GRanges object with a window around these regions
    get_sv_breakends <- function(svvcf_file, window = 1000) {
  
      # Read the VCF file
      vcf <- read.vcfR(svvcf_file,verbose = F)
      
      # Extract chromosome lengths from VCF header
      contigs <- vcf@meta[grepl("##contig", vcf@meta)]
      chrom_lengths <- sapply(contigs, function(contig) {
            contig <- gsub("##contig=<ID=|>", "", contig)
            parts <- strsplit(contig, ",length=")[[1]]
            chrom <- parts[1]
            length <- as.numeric(parts[2])
            return(c(chrom, length))
      })
          
      # Create Seqinfo object
      seqinfo_obj <- Seqinfo(seqnames = chrom_lengths[1,], seqlengths = as.numeric(chrom_lengths[2,]))
      
      # Extract INFO field with SV types
      svtypes <- extract.info(vcf, element = "SVTYPE")
      
      # Get structural variants of interest
      sv_indices <- which(svtypes %in% c("BND", "DEL", "DUP", "INS", "INV"))
      sv_variants <- vcf[sv_indices, ]
      
      # Get start and end coordinates
      chrom <- getCHROM(sv_variants)
      start_coords <- getPOS(sv_variants)
      end_coords <- sapply(strsplit(extract.info(sv_variants, "END"), ";"), function(x) as.numeric(x[1]))
      
      # Construct data frame
      df <- data.frame(seqnames = chrom,
                      start = start_coords,
                      end = ifelse(is.na(end_coords), start_coords, end_coords),
                      svtype = svtypes[sv_indices])
      
      # Split data frame by SV type
      df_split <- split(df, df$svtype)
      
      # Create GRanges object
      gr <- c(GRanges(df_split$BND$seqnames, IRanges(df_split$BND$start, df_split$BND$end),seqinfo = seqinfo_obj),
              GRanges(df_split$INS$seqnames, IRanges(df_split$INS$start, df_split$INS$end),seqinfo = seqinfo_obj),
              GRanges(df_split$DEL$seqnames, IRanges(df_split$DEL$start, df_split$DEL$start),seqinfo = seqinfo_obj),
              GRanges(df_split$DEL$seqnames, IRanges(df_split$DEL$end, df_split$DEL$end),seqinfo = seqinfo_obj),
              GRanges(df_split$DUP$seqnames, IRanges(df_split$DUP$start, df_split$DUP$start),seqinfo = seqinfo_obj),
              GRanges(df_split$DUP$seqnames, IRanges(df_split$DUP$end, df_split$DUP$end),seqinfo = seqinfo_obj),
              GRanges(df_split$INV$seqnames, IRanges(df_split$INV$start, df_split$INV$start),seqinfo = seqinfo_obj),
              GRanges(df_split$INV$seqnames, IRanges(df_split$INV$end, df_split$INV$end),seqinfo = seqinfo_obj))
      
      # Assign Seqinfo to GRanges object
      seqinfo(gr) <- seqinfo_obj
      
      # Filter chromosomes
      gr <- gr[which(seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y")))]

      # Add flanks and trim. 
      gr <- reduce(trim((flank(gr, width = window, both = TRUE))))
      
      return(gr)
    }

    # remove sex chromosomes
    if (filter.sex.chromosomes){
        if ("chrX" %in% seqnames(gr) || "chrY" %in% seqnames(gr)){
            gr <- gr[-which(seqnames(gr) %in% c("chrX","chrY"))]
        }
    }

    # run fisher test and filter if requested
    if (run.fisher.test){
        diffmeth <- bsseqDiffRegion(aml.bs[,c(namehap1,namehap2)],gr,group1 = namehap1, group2 = namehap2)
        gr <- gr[which(abs(diffmeth$meth.diff)>min.diff)]
    }

    # get divergent regions from VCF file and add to excludeRegions
    if(!is.null(vcffile)){
      excludeRegions <- c(excludeRegions,getDivergentRegions(vcffile))
    }

    # get SV regions from VCF file and add to excludeRegions
    if(!is.null(svvcffile)){
      excludeRegions <- c(excludeRegions,get_sv_breakends(svvcffile, window = sv.breakend.window))
    }

    # remove dmrs in imprinted or highly divergent regions or any other regions to exclude
    # note we could add the Vcf as an input and obtain the highly divergent regions on the fly
    if (length(excludeRegions)>0){
        gr <- subsetByOverlaps(gr,excludeRegions,invert = T)
    }

    # Retain only DMRs where hap1 and hap2 are significantly methylated and unmethylated (or vice versa) via one-sided binomial test
    # vs. either <0.4 or > 0.6 (nominal p-value < 0.05)
    hap1.counts <- cbind(getCoverage(aml.bs[,namehap1],region=gr,type="M",what="perRegionTotal"),getCoverage(aml.bs[,namehap1],region=gr,type="Cov",what="perRegionTotal"))
    hap2.counts <- cbind(getCoverage(aml.bs[,namehap2],region=gr,type="M",what="perRegionTotal"),getCoverage(aml.bs[,namehap2],region=gr,type="Cov",what="perRegionTotal"))
    
    binom_probs <- cbind(apply(hap1.counts, 1, function(row) {
                                k <- row[1]   # Number of successes (first column)
                                n <- row[2]   # Number of trials (second column)
                                binom.test(k, n, p = purity.tolerance,alternative = "l")$p.value }),
                         apply(hap1.counts, 1, function(row) {
                                k <- row[1]   # Number of successes (first column)
                                n <- row[2]   # Number of trials (second column)
                                binom.test(k, n, p = 1-purity.tolerance,alternative = "g")$p.value }),
                        apply(hap2.counts, 1, function(row) {
                                k <- row[1]   # Number of successes (first column)
                                n <- row[2]   # Number of trials (second column)
                                binom.test(k, n, p = purity.tolerance,alternative = "l")$p.value }),
                        apply(hap2.counts, 1, function(row) {
                                k <- row[1]   # Number of successes (first column)
                                n <- row[2]   # Number of trials (second column)
                                binom.test(k, n, p = 1-purity.tolerance,alternative = "g")$p.value }))
    
    gr <- gr[which((binom_probs[,1]<purity.tolerance.p.value & binom_probs[,4]<purity.tolerance.p.value) | (binom_probs[,2]<purity.tolerance.p.value & binom_probs[,3]<purity.tolerance.p.value))]
    
    # create mcol DataFrame with methylation values
    methDf <- getMeth(aml.bs,regions=gr,type="raw",what="perRegion")
    if (!is.null(normal.bs)){
        methDf <- cbind(getMeth(normal.bs,regions=gr,type="raw",what="perRegion"),methDf)
    }

    mcols(gr) <- methDf
    
    # only keep dmrs where one haplotype is >0.4 different from the normal.
    if (!is.null(normal.bs)){
      normal.counts <- cbind(getCoverage(normal.bs,region=gr,type="M",what="perRegionTotal"),getCoverage(normal.bs,region=gr,type="Cov",what="perRegionTotal"))
      normal_binom_probs <- cbind(apply(normal.counts, 1, function(row) {
                                  k <- row[1]   # Number of successes (first column)
                                  n <- row[2]   # Number of trials (second column)
                                  binom.test(k, n, p = purity.tolerance,alternative = "l")$p.value }),
                           apply(normal.counts, 1, function(row) {
                                  k <- row[1]   # Number of successes (first column)
                                  n <- row[2]   # Number of trials (second column)
                                  binom.test(k, n, p = 1-purity.tolerance,alternative = "g")$p.value }))
                           
      gr <- gr[which(rowMins(normal_binom_probs)<purity.tolerance.p.value)]

      mcols(gr) <- data.frame(mcols(gr),check.names = F) %>%
        mutate(type = case_when(
            abs(!!sym(normal.sample.name) - !!sym(namehap1)) > abs(!!sym(normal.sample.name) - !!sym(namehap2)) & !!sym(namehap1) < !!sym(normal.sample.name) ~ "hap1_hypo",
            abs(!!sym(normal.sample.name) - !!sym(namehap1)) > abs(!!sym(normal.sample.name) - !!sym(namehap2)) & !!sym(namehap1) > !!sym(normal.sample.name) ~ "hap1_hyper",
            abs(!!sym(normal.sample.name) - !!sym(namehap2)) > abs(!!sym(normal.sample.name) - !!sym(namehap1)) & !!sym(namehap2) < !!sym(normal.sample.name) ~ "hap2_hypo",
            abs(!!sym(normal.sample.name) - !!sym(namehap2)) > abs(!!sym(normal.sample.name) - !!sym(namehap1)) & !!sym(namehap2) > !!sym(normal.sample.name) ~ "hap2_hyper",
            TRUE ~ NA))
    }   
    return(gr)
}



#' filterDmrs
#' Filter Differentially Methylated Regions
#'
#' This function filters differentially methylated regions (DMRs) based on various criteria,
#' comparing methylation levels between a sample and optionally a normal sample.
#'
#' @param gr A GRanges object containing the DMRs to be filtered.
#' @param bs A BSseq object containing methylation data for the sample.
#' @param normal.bs A BSseq object containing methylation data for the normal sample (optional).
#' @param min.diff Numeric value specifying the minimum methylation difference to consider (default: 0.4).
#' @param min.cpgs Integer specifying the minimum number of CpGs required in a region (default: 5).
#' @param min.length Integer specifying the minimum length of a region in base pairs (default: 100).
#' @param excludeRegions A GRanges object specifying regions to exclude from the analysis.
#' @param filter.sex.chromosomes Logical indicating whether to filter out sex chromosomes (default: TRUE).
#'
#' @return A filtered GRanges object containing the DMRs that meet the specified criteria.
#'
#' @details
#' The function performs the following main steps:
#' 1. Checks input validity and prepares data.
#' 2. Optionally removes sex chromosomes.
#' 3. Excludes specified regions from the analysis.
#' 4. Filters regions based on various criteria (number of CpGs, length, coverage, methylation difference).
#' 5. Classifies regions as hypomethylated or hypermethylated compared to the normal sample.
#'
#' @examples
#' # Basic usage:
#' filtered_dmrs <- filterDmrs(gr, bs, normal.bs, excludeRegions = excluded_regions)
#'
#' # With custom parameters:
#' filtered_dmrs <- filterDmrs(gr, bs, normal.bs, min.diff = 0.3, min.cpgs = 10,
#'                             excludeRegions = excluded_regions,
#'                             filter.sex.chromosomes = FALSE)
#'
#' @importFrom methylKit getMeth
#' @importFrom bsseq sampleNames getCoverage
#' @importFrom GenomicRanges subsetByOverlaps
#' @importFrom dplyr mutate case_when
#'
#' @export
filterDmrs <- function(gr,bs,normal.bs=NULL,min.diff=.4,min.cpgs=5,min.length=100,excludeRegions,filter.sex.chromosomes=TRUE){
    require(methylKit)
    require(bsseq)
    
    sample.name <- NULL
    normal.sample.name <- NULL
    # if nlbmbs is not null check number of samples. If >1 then collapse
    if (!is.null(normal.bs) && ncol(normal.bs)>1){
        # exit with error that says 'nlbmbs must be a single sample'
        stop("normal.bs must be a single sample")
    } else if(!is.null(normal.bs) && ncol(normal.bs)==1){
        normal.sample.name <- sampleNames(normal.bs)[1]
    }

    # if bs is not null check number of samples. If >1 then collapse
    if (!is.null(bs) && ncol(bs)>1){
        # exit with error that says 'bs must be a single sample'
        stop("bs must be a single sample")
    } else if(!is.null(bs) && ncol(bs)==1){
        sample.name <- sampleNames(bs)[1]
    }

    # remove sex chromosomes
    if (filter.sex.chromosomes){
        if ("chrX" %in% seqnames(gr) || "chrY" %in% seqnames(gr)){
            gr <- gr[-which(seqnames(gr) %in% c("chrX","chrY"))]
        }
    }
  
    # remove dmrs in imprinted or highly divergent regions or any other regions to exclude
    # note we could add the Vcf as an input and obtain the highly divergent regions on the fly
    if (!is.null(excludeRegions)){
        gr <- subsetByOverlaps(gr,excludeRegions,invert = T)
    }

    # create mcol DataFrame with methylation values
    methDf <- getMeth(bs,regions=gr,type="raw",what="perRegion")
    methDf <- cbind(getMeth(normal.bs,regions=gr,type="raw",what="perRegion"),methDf)
    mcols(gr) <- methDf

    # filter based on minimum number of observed CpGs, length, total coverage, and minimum difference in methylation
    gr <- gr[countOverlaps(gr,bs)>=min.cpgs]
    gr <- gr[which(width(gr)>=min.length)]
    gr <- gr[which(abs(rowDiffs(as.matrix(as.data.frame(mcols(gr)[,c(sample.name,normal.sample.name)]))))>min.diff)]
    
    # only keep dmrs where one haplotype is >0.4 different from the normal.
    mcols(gr) <- data.frame(mcols(gr),check.names = F) %>%
        mutate(type = case_when(
            !!sym(sample.name) < !!sym(normal.sample.name) ~ "hypo",
            !!sym(sample.name) > !!sym(normal.sample.name) ~ "hyper",
            TRUE ~ NA))
            
    return(gr)
}

#' inferGRanges
#' Infer GRanges Object from Various Inputs
#'
#' This function infers a GRanges object from a file or data frame, providing flexibility
#' in how genomic regions are specified.
#'
#' @param input Either a file path (.bed, .txt, .tsv, .Rda) or a data frame.
#' @param chr_col Character string specifying the column name for chromosome information in the data frame (optional).
#' @param start_col Character string specifying the column name for start position information in the data frame (optional).
#' @param end_col Character string specifying the column name for end position information in the data frame (optional).
#'
#' @return A GRanges object representing the genomic regions.
#'
#' @details
#' The function can handle various input types:
#' 1. Data frame: Directly converts to GRanges, using specified or default column names.
#' 2. BED file: Reads and converts to GRanges, assuming standard BED format.
#' 3. Text or TSV file: Reads as a table and converts to GRanges.
#' 4. RDA file: Loads an R object named 'dmrs' and converts to GRanges.
#'
#' If column names are not specified, the function attempts to use standard names
#' ("chr", "start", "end"). For BED files, it handles additional standard BED columns if present.
#'
#' @examples
#' # From a data frame
#' df <- data.frame(chr = c("chr1", "chr2"), start = c(100, 200), end = c(200, 300))
#' gr <- inferGRanges(df)
#'
#' # From a BED file
#' gr <- inferGRanges("path/to/file.bed")
#'
#' # From a text file with custom column names
#' gr <- inferGRanges("path/to/file.txt", chr_col = "chromosome", start_col = "begin", end_col = "stop")
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom data.table fread
#' @importFrom tools file_ext
#'
#' @export
inferGRanges <- function(input, chr_col = NULL, start_col = NULL, end_col = NULL) {
  
  # Helper function to check and create GRanges
  create_granges <- function(df, chr_col, start_col, end_col) {
    if (is.null(chr_col) | is.null(start_col) | is.null(end_col)) {
      stop("chr_col, start_col, and end_col must be specified if column names are not standard.")
    }
    
    # Check if specified columns exist in the dataframe
    if (!all(c(chr_col, start_col, end_col) %in% colnames(df))) {
      warning("Specified columns not found. Using first three columns instead.")
      chr_col <- colnames(df)[1]
      start_col <- colnames(df)[2]
      end_col <- colnames(df)[3]
    }
    
    # Ensure that mcols have the same number of rows as the ranges
    if (is.data.table(df)) {
      mcols <- df[, !(colnames(df) %in% c(chr_col, start_col, end_col)), with = FALSE]
    } else {
      mcols <- df[, !(colnames(df) %in% c(chr_col, start_col, end_col))]
    }
    
    gr <- GRanges(seqnames = df[[chr_col]],
                  ranges = IRanges(start = df[[start_col]], end = df[[end_col]]))
    
    if (nrow(mcols) == length(gr)) {
      mcols(gr) <- mcols
    } else {
      warning("Metadata columns do not match the number of ranges. Metadata columns will be ignored.")
    }
    
    return(gr)
  }
  
  # Check if input is a data frame
  if (is.data.frame(input)) {
    # If columns are not specified, use default column names
    if (is.null(chr_col)) chr_col <- "chr"
    if (is.null(start_col)) start_col <- "start"
    if (is.null(end_col)) end_col <- "end"
    # Create and return GRanges object
    return(create_granges(input, chr_col, start_col, end_col))
  }
  
  # Check if input is a file
  if (is.character(input) && file.exists(input)) {
    # Determine file extension
    file_ext <- tools::file_ext(input)
    df <- NULL
    
    if (file_ext == "bed") {
      df <- fread(input, header = FALSE)
      colnames(df)[1:3] <- c("chr", "start", "end")
      if (ncol(df) > 3) {
        additional_cols <- c("name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
        colnames(df)[4:ncol(df)] <- additional_cols[1:(ncol(df) - 3)]
      }
    } else if (file_ext %in% c("txt", "tsv")) {
      df <- fread(input, header = TRUE)
    } else if (file_ext == "rda") {
      load(input)
      if (!exists("dmrs")) {
        stop("The .rda file does not contain an object named 'dmrs'.")
      }
      df <- dmrs
    } else {
      stop("Unsupported file format.")
    }
    
    # If columns are not specified, use default column names
    if (is.null(chr_col)) chr_col <- "chr"
    if (is.null(start_col)) start_col <- "start"
    if (is.null(end_col)) end_col <- "end"
    
    # Create and return GRanges object
    return(create_granges(df, chr_col, start_col, end_col))
  }
  
  stop("Input must be either a data frame or a valid file path.")
}


#' Parse bamstats output file
#'
#' @param statsfile path to bamstats output file
#'
#' @return a list with basic stats and histogram of read lengths
#'
#' @importFrom dplyr mutate filter pull select
#' @importFrom tibble tibble
#' @importFrom readr read.table
#'
parseBamStats <- function(statsfile){
  stats <- list()
  dat <- read.table(statsfile,sep = "\t",comment.char = "#",fill=T) %>% mutate(V2=gsub(":","",V2),V3=as.numeric(V3))
  
  # basic stats
  stats$reads <- dat %>% filter(V2=="raw total sequences") %>% pull(V3)
  stats$mapped_reads <- dat %>% filter(V2=="reads mapped") %>% pull(V3)
  stats$total_bases <- dat %>% filter(V2=="total length") %>% pull(V3)
  stats$mapped_bases <- dat %>% filter(V2=="bases mapped (cigar)") %>% pull(V3)
  stats$error_rate <- dat %>% filter(V2=="error rate") %>% pull(V3)
  stats$mean_read_length <- dat %>% filter(V2=="average length") %>% pull(V3)
  stats$mean_coverage <- stats$mapped_bases/3143893127
  
  # read length histogram
  stats$read_lengths <- dat %>% filter(V1=="FRL") %>% tibble(length=as.numeric(V2),count=as.numeric(V3))
  
  # n50
  stats$n50 <- stats$read_lengths %>% mutate(bases=as.numeric(length)*as.numeric(count)) %>% mutate(cumbases=cumsum(bases)) %>% mutate(n50=ifelse(cumbases<stats$total_bases/2,1,0)) %>% filter(n50==0) %>% dplyr::slice(1) %>% pull(length)
  
  # coverage distn
  stats$coverage_distn <- dat %>% filter(V1=="COV") %>% mutate(coverage_bin=V2,coverage=as.numeric(V3),nt=as.numeric(V4),percent=as.numeric(V4)/3143893127) %>% select(coverage_bin,coverage,nt,percent)
  
  return(stats)
}

#' Annotate DMRs with genome annotations
#' 
#' Annotate DMRs with genome annotations, including custom annotations, and methylation values from user-provided BSseq objects.
#' 
#' @param gr character or GRanges object. The DMRs to annotate. If a character, assumed to be a path to a .rda file with a GRanges object named 'filtered_dmrs'.
#' @param annots_beds character vector of file paths to custom annotation BED files.
#' @param annots_names character vector of names for custom annotations. If NULL, the names will be the base names of the BED files.
#' @param bsseq_paths character vector of file paths to BSseq objects.
#' @param hg38_features character vector of genome annotation features to include. Default: c("hg38_genes_promoters", "hg38_genes_intergenic", "hg38_cpg_islands", "hg38_basicgenes").
#' @return GRanges object with annotations.
#' @export
annotateDmrs <- function(gr, annots_beds = NULL, annots_names = NULL, 
                         bsseq_paths = NULL, 
                         hg38_features = c("hg38_genes_promoters", "hg38_genes_intergenic", 
                                           "hg38_cpg_islands", "hg38_basicgenes")) {
  # Load DMRs, if gr is a file path to a .rda file. Else assume gr is a GRanges object
  if (is.character(gr) && grepl("\\.rda$", gr)) {
    if (!file.exists(gr)) {
      stop("DMRs file does not exist.")
    }
    load(gr)
    dmrs <- filtered_dmrs
  } else {
    dmrs <- gr
  }

  # Process BSseq objects, ensure sample names are unique
  bsseq_obj_list <- list()
  all_sample_names <- character(0)
  
  if (!is.null(bsseq_paths)) {
    for (i in seq_along(bsseq_paths)) {
      # Load the BSseq object
      bs <- HDF5Array::loadHDF5SummarizedExperiment(dir = bsseq_paths[i])
      bs_samples <- sampleNames(bs)
      
      # Ensure unique sample names
      for (j in seq_along(bs_samples)) {
        original_name <- bs_samples[j]
        new_name <- original_name
        suffix <- 1
        while (new_name %in% all_sample_names) {
          new_name <- paste0(original_name, "_", suffix)
          suffix <- suffix + 1
        }
        bs_samples[j] <- new_name
      }
      
      # Update sample names in the BSseq object if they were changed
      if (!identical(sampleNames(bs), bs_samples)) {
        sampleNames(bs) <- bs_samples
      }
      
      # Extract methylation data for each sample and add it to the DMRs object
      for (j in seq_along(bs_samples)) {
        sample_name <- bs_samples[j]
        sample_mtx <- getMeth(bs, regions = dmrs, type = "raw", what = "perRegion")
        mcols(dmrs)[[sample_name]] <- sample_mtx
      }
      
      bsseq_obj_list[[i]] <- bs
      all_sample_names <- c(all_sample_names, bs_samples)
    }
  }
  
  # Handle custom annotations
  custom_annotations_granges <- list()
  
  # Set default annotation names based on file basenames if not provided
  if (is.null(annots_names) && !is.null(annots_beds)) {
    annots_names <- basename(annots_beds)
    annots_names <- sub("\\..*", "", annots_names)
  }
  
  # Ensure the number of annotation names matches the number of custom annotation files
  if (!is.null(annots_names) && !is.null(annots_beds)) {
    if (length(annots_names) != length(annots_beds)) {
      stop("The number of annotations names must match the number of custom annotation files.")
    }
  }
  
  # Read custom annotations into GRanges objects
  for (i in seq_along(annots_beds)) {
    custom_annotations_granges[[i]] <- read_annotations(con = annots_beds[i], genome = 'hg38', 
                                                        name = annots_names[i], format = 'bed')
  }
  
  # Prepare custom annotation names
  custom_annotations_names <- paste0("hg38_custom_", annots_names)
  
  # Combine default hg38 features with custom annotations
  annotation_features <- c(hg38_features, custom_annotations_names)
  print(paste("Annotating DMRs with:", paste(annotation_features, collapse = ", ")))
  # Build annotation object
  ann <- build_annotations(genome = 'hg38', annotations = annotation_features)
  
  # Annotate DMRs with genome and custom annotations
  dmrs_annotated <- annotate_regions(regions = dmrs, annotations = ann, ignore.strand = TRUE, 
                                     quiet = FALSE)
  
  return(dmrs_annotated)
}
