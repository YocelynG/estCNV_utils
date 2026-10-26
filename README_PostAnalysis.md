# CNV estimation: post-analysis

`estCNV_PostAnalysis.qmd` takes the per-gene fold changes produced by [`estCNV.qmd`] and links them to sample metadata and parasitemia. For each sample it keeps one sequencing run and adds the read depth
The code generates plots fold changes per gene and, if you choose, compares them between groups (e.g. regions). Finally, it lists the samples that pass the deletion and duplication thresholds

## Prerequisites

1. **Run `estCNV.qmd` first.** The post-analysis reads its final output, `fold_changes_final_output.csv`, and the same allele tables it used.
3. **R packages**:
   - `tidyverse`, `readxl`, `knitr`, `rmarkdown` (always needed)
   - `rstatix`, `ggpubr` (only needed when `location = TRUE`)

```r
install.packages(c("tidyverse", "readxl", "knitr", "rmarkdown", "rstatix", "ggpubr"))
```

## Inputs

### 1. Fold changes (`fold.changes.path`)

The `fold_changes_final_output.csv` file written by `estCNV.qmd`. It has one row per sample, run and gene:

| Column           | Description                               |
|------------------|-------------------------------------------|
| `sample_name`    | Sample identifier                         |
| `Run`            | Sequencing run                            |
| `gene`           | Gene (e.g. `CRT`, `HRP2`, `HRP3`, `MDR1`, `PM`) |
| `FoldChangeGene` | Per-gene fold change                      |

### 2. Allele tables (`allele.data.path`)

The folder that holds the allele tables used in `estCNV.qmd`. The script searches it recursively for files whose names start with `allele_data` or `allele_table`. Each file needs these columns:

- a sample ID column
- `reads`
- `target_name` or `locus`

If a file has no `Run` column, the name of its parent folder is used as the run. Sample names are cleaned (suffixes such as `_S12` and `_L001` are removed)

The allele tables are used for two things:

- **Run selection**: when a sample passed QC in more than one run, the run with more total reads is kept
- **Read depth**: the depth columns and plots described in [Workflow overview](#workflow-overview)

### 3. Metadata (`metadata.path`)

A CSV, TSV/TXT or Excel (`.xlsx`) file with **one row per sample**. It must have a sample ID column and contain parasitemia and location information. If `location = TRUE`, it must also have the grouping column named in `location.variable` 

### 4. Panel information (`panel.info.dir`)

The `panel_information/` folder bundled with this repository:

```
panel_information/
├── validated_target_filter.tsv   # Validated loci, used for the all-loci read depth
├── targets_of_interest.tsv       # target_name and group, used for the per-gene read depth
├── D1.1/
│   └── D1.1_amplicon_info.tsv    # Maps amplicon coordinates (locus) to target_name
└── ...
```

In `targets_of_interest.tsv`, groups are collapsed to their gene by dropping everything after the first `_`. For example, `HRP2_P4`, `HRP2_DD2` and `HRP2_P34` all become `HRP2`.

## Sample ID matching

The script looks for the sample ID column in each table. Column names are matched case-insensitively, and the first match in the order below is used:

| Table                          | Recognised ID columns |
|--------------------------------|-----------------------|
| Fold changes, allele tables    | `sample_name`, `sampleID`, `sample_id`, `NIDA` |
| Metadata, parasitemia          | `NIDA`, `sample_name`, `sampleID`, `sample_id` |

Before matching, NIDA-style IDs are converted to a common form:

The report prints how many sample-runs matched the metadata and lists those that did not. Unmatched samples are left out of the analysis. Controls are recognised by `control.pattern` and are always left out of plots, tests and exports.

## Configuration

Edit the **Prepare Inputs** chunk at the top of `estCNV_PostAnalysis.qmd`:

```r
fold.changes.path = "/path/to/fold_changes_final_output.csv"
allele.data.path = "/path/to/allele_tables/"
metadata.path = "/path/to/metadata.txt"
parasitemia.path = "/path/to/parasitemia_db.txt"
output.dir = "/path/to/output/"           # Created if it does not exist
```

| Parameter | Example | Description |
|-----------|---------|-------------|
| `location` | `TRUE` | Run the **Plots by Region or other** section (`TRUE`) or skip it (`FALSE`). |
| `location.variable` | `"region"` | Metadata column used to group samples, matched case-insensitively (e.g. `"region"`, `"district"`, `"us"`). Its values are treated as categories. |
| `location.order` | `c("North", "Centre", "South")` | Optional order of the groups on the plots. Values not listed follow in alphabetical order. |
| `parasitemia.decimal.mark` | `","` | Decimal separator used in `parasitemia.path`. |
| `parasitemia.extreme` | `1e9` | In the per-gene parasitemia plot, points above this value are drawn in navy. They are still used in the regression. |
| `control.pattern` | `"3D7\|Dd2\|HB3\|PM2"` | Regular expression (case-insensitive) that identifies control samples by name. |
| `expected.cnv` | see the file | Expected CNV per gene, drawn as a dotted line in the group boxplots. `duplication` is drawn at fold change 2, `deletion` at 0. |
| `gene.colors` | see the file | One color per gene. Every gene in the fold-change file needs a color, otherwise the script stops. |


## Outputs

All files are written to `output.dir`.

| File | Description |
|------|-------------|
| `fold_changes_final_metadata.txt` | One row per sample and gene for the selected run, with metadata, parasitemia and read-depth columns. Controls are excluded. |
| `CNV_ReadDepth_PerGene.txt` | Per-gene locus depth and all-loci depth for each sample and run, with the number of amplicons used and the number that had reads. |
| `CNV_Boxplots_FoldChange_PerGene.pdf` | Fold change per gene, all samples pooled. |
| `CNV_FoldChange_Parasitemia_PerGene.pdf` | Fold change against parasitemia (log10), per gene. |
| `CNV_FoldChange_AllLociReadDepth_PerGene.pdf` | Fold change against mean reads across all validated loci, per gene. |
| `CNV_FoldChange_ReadDepth_PerGene.pdf` | Fold change against mean reads across the gene's targets of interest. |
| `CNV_Outlier_Summary.txt` | Samples and genes that pass the deletion or duplication thresholds. Columns: `SampleID`, `FoldChange`, `Gene`, `GeneReads`, `SampleTotalReads`, `Parasitemia`, `Type`, then the original metadata columns. |

When `location = TRUE`, three more files are written. Their names end with `location.variable` (e.g. `_region`):

| File | Description |
|------|-------------|
| `CNV_Tests_<variable>.txt` | Dunn's test for each pair of groups within each gene, with the group medians and the gene's Kruskal–Wallis result. |
| `CNV_Boxplots_FoldChange_<variable>.pdf` | Fold change per gene and group, with Dunn p-values. |
| `CNV_Regression_FoldChange_Parasitemia_<variable>.pdf` | Fold change against log10(parasitemia), per gene and group. |
