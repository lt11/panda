# panda

## Overview

The `panda` pipeline is devised for building and inspecting a feature-level pangenome from pairwise genome alignments and genome annotations, with a current focus on *Saccharomyces cerevisiae* gene annotations.

The core workflow lives in [mrk/panda.Rmd]. It stages GFF annotations, converts them to BED, runs `impg` against a PAF alignment, reduces overlapping feature intervals into nonredundant sub-blocks, computes presence statistics, and checks how systematic yeast gene names are represented in the resulting pangenome.

The different chunks of code can be also found, e.g. [mrk/c2.R].

The markdown [mrk/panda.Rmd] reports a detailed description of the chunks and extensive comments.

## Installation

This repository is not yet distributed as a packaged software. Instead, the workflow is executed in a local environment where all required tools and R packages are pre-installed.

All the following step are documented below:

- clone the repository
```
git clone https://github.com/lt11/panda.git
cd panda
```

- install dependencies (e.g. `impg`) and make sure they are accessible in your `$PATH`

- install R packages, e.g.:
```
Rscript -e 'install.packages(c(
  "data.table",
  "rmarkdown"
))'
```

- configure external inputs

- check `impg` version compatibility

## Repository Layout

Top-level directories:

- [`ids/`] input lists of strain-haplotype identifiers used to define the analysis set and output column order.
- [`aln/`] input alignments in PAF format.
- [`anno/`] staged annotation files produced during the workflow.
- [`mrk/`] the R Markdown analysis and helper scripts.
- [`png/`] primary pangenome output tables produced by the notebook.
- [`sts/`] downstream summary statistics, currently focused on gene sub-blocks.
- [`check-sys/`] validation outputs comparing systematic gene names in the pangenome against source GFF annotations.

Key source files:

- [mrk/panda.Rmd]: end-to-end notebook.
- [mrk/c1.sh]: copies and normalises input GFF files.
- [mrk/c2.R]: converts GFF annotations to BED.
- [mrk/c3.R]: computes the feature pangenome from `impg` output.
- [mrk/c4.R]: checks systematic genes between annotations and the pangenome.
- [mrk/c5.R]: additional systematic-gene validation/reporting step.
- [mrk/run-profvis.R]: helper to profile an R script with `profvis`.

## Inputs

The workflow expects three main inputs:

1. A single PAF alignment in [`aln/`]. E.g. you can use a PAF file generated from [`PGGB`](https://github.com/pangenome/pggb) or directly from [`wfmash`](https://github.com/waveygang/wfmash).
2. A list of strain-haplotype identifiers in [`ids/ids-ps.txt`], one per line, such as `SGDref-0` or `CMF-1`.
3. A directory of source GFF annotations outside this repository, configured in [mrk/c1.sh] and in the corresponding chunk inside [mrk/panda.Rmd].

Naming conventions matter:

- Contigs in the PAF use PanSN-style names with `#`, for example `SGDref#0#chrI`.
- Annotation files are named with `-`, for example `SGDref-0-features.gff`.
- `ids/ids-ps.txt` defines both the haplotypes (e.g. for example `SGDref-0` to analyse and the order of the output columns.

## Workflow

The notebook runs the analysis in four main stages.

1. Stage GFF files.
   `mrk/c1.sh` copies the requested annotation files into `anno/gff/`, optionally merging mitochondrial annotations.
2. Convert GFF to BED.
   `mrk/c2.R` filters feature classes, normalises coordinates, and writes BED files to `anno/bed/`.
3. Build the pangenome.
   `mrk/c3.R` runs `impg` for each haplotype BED file against the single PAF alignment, merges overlapping hits into nonredundant feature sub-blocks, and writes the main pangenome tables to `png/`.
4. Summarise and validate gene sub-blocks.
   `mrk/c4.R` derives gene-level summary tables in `sts/gene/` and compares systematic yeast gene names in the pangenome against the original GFF annotations, writing reports to `check-sys/`.

## Dependencies

This repository is not packaged as an R package. The analysis assumes a local environment with:

- `Rscript`
- `rmarkdown` to render the notebook
- `impg`

R packages used in the scripts include:

- `data.table`
- `this.path`
- `scriptName`
- `GenomicRanges`
- `tictoc`
- `profvis`

The notebook notes an `impg` version difference:

- `impg 0.2.0`: `impg -I -p file.paf -b file.bed`
- `impg 0.2.3`: `impg query -I -p file.paf -b file.bed`

Adjust the command in [mrk/c3.R] or the matching chunk in [mrk/panda.Rmd] to match the installed version.

## Running an Analysis

There is no executable or standalone CLI entrypoint; the repository is driven from the notebook or by running a single chunk, as described below.

Before rendering:

- place exactly one input `.paf` file in [`aln/`]
- update [`ids/ids-ps.txt`] for the genomes/haplotypes to analyse
- set the external annotations directory in the first notebook chunk
- decide whether mitochondrial annotations should be included with `with_mito`
- confirm the `impg` command variant for the installed version
- review notebook chunks marked with `eval = FALSE`

Then render:

```bash
Rscript -e 'rmarkdown::render("mrk/panda.Rmd")'
```

If you prefer to run the stages individually, the scripts map to the notebook sections:

```bash
bash mrk/c1.sh
Rscript mrk/c2.R
Rscript mrk/c3.R
Rscript mrk/c4.R
Rscript mrk/c5.R
```

This can be useful if you want to profile the scripts on your machine.

Profiling helper:

```bash
Rscript mrk/run-profvis.R c3.R
```

Run that command from `mrk/` if you want `run-profvis.R` to resolve the target script as written.

## Main Outputs

Primary outputs are written to [`png/`]:

- [`png/pan-features.txt`]: the feature pangenome table. The first columns are `Class_id` and `Features_id`; remaining columns are haplotypes such as `SGDref#0` and `S288C#0`, containing genomic intervals or `MA` (i.e. missing annotation) when absent.
- [`png/generators.txt`]: the `class:feature#strand` entries that generated each output row.
- `pan-features.RData`: an R object with the same table.

Gene-level summaries are written to [`sts/gene/`]:

- [`sts/gene/sts-by-haplos.txt`]: haplotype-based presence table with `Ν_pres`, `F_pres`, `N_feats`, `N_feats_sys`, and `N_feats_rid`.
- [`sts/gene/sts-by-genomes.txt`]: genome-collapsed version of the same summary.
- `n-sblocks-regs.txt` and `n-rid-private.txt`: additional sub-block count summaries.

Systematic-gene checks are written to [`check-sys/`]:

- [`check-sys/sys-summary.txt`]: per-genome counts of systematic genes seen in the pangenome versus the source GFF, plus gains/losses.
- `gene-lost-*.txt` and `gene-gained-*.txt`: per-genome discrepancy lists.

## Notes

- The analysis currently focuses on `gene` and `pseudogene` annotations in the BED conversion step, although the scripts contain commented alternatives for additional feature classes.
- `mrk/panda.Rmd` documents a known warning from `impg 0.2.0` related to asymmetric `wfmash` alignments.

## License

This repository is distributed under the [MIT License].
