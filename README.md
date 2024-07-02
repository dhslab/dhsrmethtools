# dhsrmethtools

R package with utility functions for converting methylation data to `bsseq` objects, performing differential methylation analysis, and filtering DMRs.

## Installation

```R
devtools::install_github("dhslab/dhsrmethtools")
```

## Usage

library(dhsrmethtools)

## Custom functions

`bed2bsseq()` - Converts simplified bed files with methylation data to `bsseq` objects.

`modkit2bsseq()` - Converts ONT bedmethyl files to `bsseq` objects.

`pbcpg2bsseq()` - Converts pb-CpG tools bedmethyl files to `bsseq` objects.

`hapBSseq()` - Creates `bsseq` objects from haplotype-specific and combined bedmethyl files based on platform type.

`callHapDmrs()` - Performs differential methylation analysis between two haplotypes within the same sample.

`callDmrs()` - Conducts differential methylation analysis between a sample and a normal reference, identifying DMRs.

`filterHapDmrs()` - Filters haplotype-specific DMRs based on various criteria, optionally comparing against a normal sample.

`filterDmrs()` - Filters DMRs based on coverage, methylation differences, and other criteria.

`inferGRanges()` - Infers `GRanges` objects from files or data frames, supporting various input formats.


## To-Do List

- [x] Update callHapDmrs function to include cpg and length parameters
- [ ] Update filterHapDmrs function
- [ ] Update documentation
