# Tumor / cfDNA Concordance Toolkit

A collection of standalone Python scripts for comparing paired tumor-tissue
and cell-free DNA (liquid biopsy) whole-genome sequencing results: variant
concordance, substitution spectra, copy-number profiling, oncoprint-style
mutation matrices, and gene-set / functional-enrichment overlap analysis.

Each script is a self-contained command-line tool with its own inputs,
outputs, and options. They share no runtime state, so you can use any
subset of them independently, and each can be read on its own.

## Contents

| Script | Purpose |
|---|---|
| `snv_substitution_spectrum.py` | Single-base substitution spectrum (SBS6-style) from VCFs, overall / per-sample / per-group. |
| `variant_concordance_analysis.py` | Shared vs. unique variants and concordance (overlap or Jaccard) between paired VCFs. |
| `oncoplot_from_annovar.py` | Gene x sample oncoprint from ANNOVAR-style exonic variant annotations. |
| `cnv_profiling_pipeline.py` | Genome-wide CNV profiles, GISTIC-style segment plots, recurrent-region detection, and summary statistics. |
| `gene_set_overlap_analysis.py` | Union, pairwise, and exhaustive N-way overlap of gene lists across two sample groups. |
| `functional_enrichment_overlap.py` | Shared terms between two functional-enrichment result sets, per ontology/category. |
| `mutation_annotation_parser.py` | Parse free-text mutation annotations (e.g. `"C>T(3)"`) into a structured gene x category matrix. |
| `upset_intersection_plot.py` | Two-group UpSet-style intersection plot from a precomputed intersection table. |
| `vcf_utils.py` | Shared VCF-reading helpers used by the two VCF-based scripts. |

## Requirements

- Python 3.9+
- See `requirements.txt` (pandas, numpy, matplotlib, openpyxl)

```bash
python -m venv .venv
source .venv/bin/activate   # on Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## General usage pattern

Every script follows the same shape:

```bash
python <script>.py --input-... <path> [--output-dir results/<name>] [options...]
```

Run `python <script>.py --help` for the full list of options for that
script. Output directories are created automatically and default to a
subfolder under `results/`; none of the scripts write outside the
directory you point them at.

None of the scripts use randomness (no seeds required for reproducibility);
given the same inputs and options, outputs are deterministic.

## Script details

### `snv_substitution_spectrum.py`

Counts single-base substitutions in a set of VCFs and collapses them onto
the 6 canonical, strand-symmetric substitution classes (C>A, C>G, C>T,
T>A, T>C, T>G — the same convention used by COSMIC's SBS signatures).

**Input:** a sample sheet CSV with columns `sample_id, group, vcf_path`.
**Output:** a class x sample count table, plus bar charts for the overall
spectrum, the per-sample contribution, and the per-group comparison.

### `variant_concordance_analysis.py`

For paired samples (e.g. a tumor and a matched liquid-biopsy draw from the
same case), compares the two raw variant call sets and reports shared /
unique counts and concordance, using either the fraction of the reference
sample's variants (`--metric overlap`, default) or the Jaccard index
(`--metric jaccard`). Both metrics are always written to the output table.

**Input:** a pairing sheet CSV with columns `pair_id, reference_vcf, other_vcf`.
**Output:** a concordance table (CSV) plus a bar chart, a stacked
shared/unique bar chart, and a dot plot.

### `oncoplot_from_annovar.py`

Builds a gene x sample oncoprint from ANNOVAR-style
`exonic_variant_function` annotation files: filters to a configurable set
of mutation categories, picks the most frequently mutated genes, and
draws each cell as a proportional stack of the mutation categories
observed for that gene in that sample, with marginal bar charts showing
per-sample and per-gene mutation burden.

**Input:** a directory of `*.exonic_variant_function` files (one per
sample), or a sample sheet CSV (`sample_id, file_name`) if your file
names don't already match the sample IDs you want to use.
**Output:** the gene x sample count matrix (CSV), the long-format
gene/sample/category counts (CSV), and the oncoprint figure.

### `cnv_profiling_pipeline.py`

The most involved script: reads per-sample copy-number ratio tracks and
discrete CNV segment calls (Control-FREEC-style `*_ratio.txt` / `*_CNVs`
output), lays samples out on a shared genome-position axis, and produces
per-sample genome-wide and GISTIC-style plots, a combined multi-sample
overview figure, recurrent-region tables (exact-coordinate and
reciprocal-overlap-merged), summary statistics, and a multi-sheet Excel
workbook.

**Input:** a sample sheet CSV (`sample_id, group, sample_type, ratio_file,
cnv_file`). If your reference genome uses accession-style contig names
(e.g. RefSeq `NC_XXXXXX.X`) rather than plain chromosome numbers, supply
`--chrom-map` (a 2-column TSV of `accession<TAB>label`) to translate them.
Plot order defaults to a natural sort of the chromosome labels found in
the data; pass `--chrom-order` (one label per line) to control it explicitly.
**Output:** `results/cnv/genome_wide_profiles/`, `results/cnv/gistic_style_plots/`,
`combined_overview.png`, and `cnv_analysis_summary.xlsx`.

Note: the per-sample profile plots use `log2(1.25)` / `log2(0.75)` as the
default gain/loss guide lines, while the combined overview figure uses a
flat ±0.25 by default (both are configurable via `--gain-threshold` /
`--loss-threshold` / `--overview-guide-threshold`) -- set them equal if you
want the two figure types to use the same visual threshold.

### `gene_set_overlap_analysis.py`

Starting from one gene-list CSV per sample, builds a combined gene matrix
and characterizes overlap within and between two sample groups: a
frequency table across all samples, group-vs-group shared/unique genes,
and exhaustive pairwise/N-way overlap (inclusive, in the CSV output, and
exclusive-only, in a printed/saved text report) within each group.

**Input:** a directory of per-sample gene-list CSVs (gene symbols in the
first column); group membership is inferred from a configurable
case-insensitive prefix match on the file name.
**Output:** `gene_matrix.csv`, `gene_frequency.csv`, `group_overlap.csv`,
`group_unique_genes.csv`, `group_a_specific_genes.csv`,
`pairwise_multiway_overlap.csv`, `exclusive_nway_report.txt`.

### `functional_enrichment_overlap.py`

For each ontology/category (default: GO Biological Process, Cellular
Component, Molecular Function, KEGG, Reactome), finds the enrichment terms
shared between two group-specific result files and writes the matching
rows from both.

**Input:** an input directory containing `{category}_{group_a_suffix}.csv`
and `{category}_{group_b_suffix}.csv` per category, each with a `Term`
column.
**Output:** `{category}_shared_terms.csv` per category.

### `mutation_annotation_parser.py`

Parses free-text mutation annotations (e.g. `"C>T(3)"`, or an indel count
folded into a "Mutation Type" column) into structured records, aggregates
them into a gene x category matrix, cross-checks the parsed total against
a reported total if available, and plots mutation frequency by gene.

**Input:** a CSV with `Gene name`, `Sample`, `Mutation Type`, and a
mutation-detail column (name configurable; falls back to the 5th column if
not found, since spreadsheet-to-CSV exports often leave this column
unnamed).
**Output:** `mutation_records.csv`, `mutation_frequency_by_gene.csv`, and a
bar chart.

### `upset_intersection_plot.py`

Draws a two-group UpSet-style intersection plot (paired bar chart + dot
matrix + set-size bars) from a precomputed intersection table. This
script visualizes intersection sizes you supply; it does not compute
them, so it can be reused for any kind of paired-group set comparison
(variants, genes, or anything else).

**Input:** a CSV with per-combination size columns for each group and
0/1 membership columns per sample (see the script's `--help` for the
exact column-naming convention).
**Output:** the plot (PNG).

## Customization

- All file paths accepted by these scripts are relative to either the
  current working directory or (for sample/pairing sheets) the sheet's own
  location -- there are no hard-coded absolute paths anywhere in this
  repository.
- Mutation-category lists, overlap thresholds, top-N cutoffs, and plot
  thresholds are all exposed as CLI options with the defaults noted above.
- None of the scripts bundle example data. They are written against the
  documented input formats above; use your own sequencing/annotation
  outputs.

## Reproducibility

- No randomness is used anywhere in this toolkit, so outputs are fully
  deterministic given the same inputs and options.
- Dependency versions are pinned to minimums in `requirements.txt`; for
  strict reproducibility, pin exact versions with `pip freeze` after
  installing.
- All scripts write their outputs under a single `--output-dir` (or
  `--output`), so a full run's results are self-contained and easy to
  diff or archive.

## Known limitations

- `oncoplot_from_annovar.py` assumes the standard 12-column ANNOVAR
  `exonic_variant_function` layout; files with a different column count
  are padded/truncated to 12 columns, which may misalign extra or missing
  fields for non-standard exports.
- `cnv_profiling_pipeline.py`'s reciprocal-overlap region merging is a
  simple greedy clustering (first-seen anchor absorbs any sufficiently
  overlapping event); it is not a formal interval-clustering algorithm and
  results can depend on row order within ties.
- `gene_set_overlap_analysis.py` and `oncoplot_from_annovar.py` compare
  gene symbols as case-insensitive strings; they do not resolve gene
  synonyms or aliases.

## License

MIT -- see `LICENSE`.
