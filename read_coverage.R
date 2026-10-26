# Read coverage per sample from an allele table
#
# Usage:  Rscript read_coverage.R <allele_table> [output_file]
#   or set allele.table.path (and output.path) below and run the script in R or RStudio.
#

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)

# Allele table (tsv/txt or csv), one row per allele. Needs a sample ID column
# (sample_name, sampleID, sample_id or NIDA), reads, and target_name or locus. Run is optional
allele.table.path = if (length(args) >= 1) args[1] else "/path/allele_data.txt"

# Output table; by default next to the allele table (e.g. allele_data_read_coverage.txt)
output.path = if (length(args) >= 2) args[2] else sub("(\\.[^./]+)?$", "_read_coverage.txt", allele.table.path)

# Remove sequencing suffixes from sample names (_S12, _L001), as in estCNV.qmd,
# so that the lanes of a sample are counted together
standardise.sample.name = TRUE

# TRUE: loci of the run with no reads in a sample count as 0 reads (as in estCNV_PostAnalysis.qmd).
# FALSE: only the loci with reads in the sample are used
count.missing.loci.as.zero = TRUE


# Read as text so sample IDs keep their format (e.g. 2076817.1)
allele.table <- if (grepl("\\.csv$", allele.table.path, ignore.case = TRUE)) {
  read.csv(allele.table.path, colClasses = "character", check.names = FALSE)
} else {
  read.delim(allele.table.path, colClasses = "character", check.names = FALSE)
}

# Sample ID columns that are recognised (case-insensitive), in order of preference
id.aliases = c("sample_name", "sampleid", "sample_id", "nida")
id.column <- colnames(allele.table)[tolower(colnames(allele.table)) %in% id.aliases]
if (length(id.column) == 0) {
  stop("No sample ID column. Expected one of: sample_name, sampleID, sample_id, NIDA")
}
id.column <- id.column[order(match(tolower(id.column), id.aliases))][1]
locus.column <- intersect(c("target_name", "locus"), colnames(allele.table))[1]
if (is.na(locus.column)) stop("No target_name or locus column")
if (!"reads" %in% colnames(allele.table)) stop("No reads column")
has.run <- "Run" %in% colnames(allele.table)
if (!has.run) allele.table$Run <- NA_character_
cat("Sample ID column:", id.column, "| locus column:", locus.column,
    "| Run column:", if (has.run) "yes" else "no (one run)", "\n")

alleles <- allele.table %>%
  transmute(SampleID = trimws(.data[[id.column]]),
            Run,
            locus = .data[[locus.column]],
            reads = suppressWarnings(as.numeric(reads)))
if (any(!is.finite(alleles$reads) | alleles$reads < 0)) stop("Invalid read counts in the reads column")
if (standardise.sample.name) {
  alleles$SampleID <- sub("_L00\\d+$", "", sub("_S\\d+(_.*)?$", "", alleles$SampleID))
}

# Reads of each locus in each sample (and run), summed across alleles
locus.reads <- alleles %>%
  group_by(SampleID, Run, locus) %>%
  summarise(locus_reads = sum(reads), .groups = "drop")

if (count.missing.loci.as.zero) {
  # Loci of each run: those with reads in at least one sample of the run
  run.loci <- locus.reads %>%
    filter(locus_reads > 0) %>%
    distinct(Run, locus)
  locus.reads <- locus.reads %>%
    distinct(SampleID, Run) %>%
    inner_join(run.loci, by = "Run", relationship = "many-to-many") %>%
    left_join(locus.reads, by = c("SampleID", "Run", "locus")) %>%
    mutate(locus_reads = replace_na(locus_reads, 0))
} else {
  locus.reads <- filter(locus.reads, locus_reads > 0)
}

# One row per sample (and run). Min and max are the reads of the locus with the fewest
# and with the most reads in the sample
coverage <- locus.reads %>%
  group_by(SampleID, Run) %>%
  summarise(ReadCoverage = round(mean(locus_reads), 2),
            MinReads = min(locus_reads),
            MaxReads = max(locus_reads),
            TotalReads = sum(locus_reads),
            NLoci = n(),
            NLociWithReads = sum(locus_reads > 0),
            .groups = "drop") %>%
  arrange(Run, SampleID)
if (!has.run) coverage <- select(coverage, -Run)

cat("Samples:", n_distinct(coverage$SampleID), "| sample-runs:", nrow(coverage), "\n")
if (has.run) {
  cat("Loci per run:\n")
  locus.reads %>% distinct(Run, locus) %>% count(Run, name = "NLoci") %>% print(n = Inf)
} else {
  cat("Loci:", n_distinct(locus.reads$locus), "\n")
}
print(coverage)

old.options <- options(scipen = 999)
write.table(coverage, output.path, sep = "\t", quote = FALSE, row.names = FALSE)
options(old.options)
cat("Read coverage per sample saved as:", output.path, "\n")
