% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{filterHapDmrs}
\alias{filterHapDmrs}
\title{filterHapDmrs
Filter Haplotype-Specific Differentially Methylated Regions}
\usage{
filterHapDmrs(
  gr,
  aml.bs,
  normal.bs = NULL,
  normal.sample.name = NULL,
  vcffile = NULL,
  svvcffile = NULL,
  sv.breakend.window = 1000,
  min.diff = 0.4,
  purity.tolerance = 0.4,
  purity.tolerance.p.value = 0.01,
  excludeRegions = GRanges(),
  namehap1 = grep("hap1", sampleNames(aml.bs), value = T)[1],
  namehap2 = grep("hap2", sampleNames(aml.bs), value = T)[1],
  filter.sex.chromosomes = TRUE,
  run.fisher.test = FALSE
)
}
\arguments{
\item{gr}{A GRanges object containing the regions to be filtered.}

\item{aml.bs}{A BSseq object containing methylation data for the AML sample.}

\item{normal.bs}{A BSseq object containing methylation data for the normal sample (optional).}

\item{normal.sample.name}{Character string specifying the name of the normal sample (optional).}

\item{vcffile}{Path to a VCF file for identifying divergent regions (optional).}

\item{svvcffile}{Path to SV VCF file for filtering out DMRs near SVs (optional).}

\item{sv.breakend.window}{bp to flank SV breakends for filtering DMRs (default: 1000).}

\item{min.diff}{Numeric value specifying the minimum methylation difference to consider (default: 0.4).}

\item{purity.tolerance}{Numeric value specifying the maximum allowable purity deviation (default: 0.4).}

\item{purity.tolerance.p.value}{One-sided P-value cutoff for binomial test of methylation counts against purity.tolerance or 1-purity.tolerance (default: 0.01).}

\item{excludeRegions}{A GRanges object specifying regions to exclude from the analysis (default: empty GRanges).}

\item{namehap1}{Character string specifying the name of haplotype 1 in the aml.bs object (default: first sample name containing "hap1").}

\item{namehap2}{Character string specifying the name of haplotype 2 in the aml.bs object (default: first sample name containing "hap2").}

\item{filter.sex.chromosomes}{Logical indicating whether to filter out sex chromosomes (default: TRUE).}

\item{run.fisher.test}{Logical indicating whether to run Fisher's exact test for differential methylation (default: FALSE).}
}
\value{
A filtered GRanges object containing the haplotype-specific DMRs that meet the specified criteria.
}
\description{
This function filters haplotype-specific differentially methylated regions (DMRs) based on various criteria.
It can compare methylation levels between two haplotypes and optionally against a normal sample.
}
\details{
The function performs the following main steps:
1. Checks input validity and prepares data.
2. Optionally removes sex chromosomes.
3. Optionally runs Fisher's exact test for differential methylation.
4. Identifies and excludes divergent regions if a VCF file is provided.
5. Filters regions based on various criteria (number of CpGs, length, coverage, methylation difference).
6. If a normal sample is provided, further filters and classifies regions based on comparison to the normal sample.
}
\examples{
# Basic usage:
filtered_dmrs <- filterHapDmrs(gr, aml.bs)

# With normal sample and VCF file:
filtered_dmrs <- filterHapDmrs(gr, aml.bs, normal.bs, normal.sample.name = "normal",
                               vcffile = "path/to/vcf.file", min.diff = 0.3)

}
