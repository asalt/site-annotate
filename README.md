# site-annotate

A command-line tool for 15mer site rollup and annotation.


this is designed to be modular and easy to modify how filtering and aggregation are performed

# Input file formats

## metadata

a tab separated file with sample information of the form:

| name (must be unique) | recno | runno | searchno | label | group |
| ---- | ----- | ----- | -------- | ----- | ----- |
| ctrl-1 | 12345 | 1 | 1 | 126 | ctrl |
| ctrl-2 | 12345 | 1 | 1 | 127N | ctrl |
| treat-1 | 12345 | 1 | 1 | 126 | treatment |
| treat-2 | 12345 | 1 | 1 | 127N | treatment |

The last `group` column is optional but is a default standard variable
name used internally for faceting / grouping.
This can be overwritten by specifying variables in the config file.
Multiple additional columns of metadata can be specified and referenced
within the config file as well.

## config

a base toml config file can be generated with:
`site-annotate get-config`

# Running

## site level summarization

note this can be skipped entirely using fragpipe output (described more later)

starting with `psm.tsv` file from fragpipe output
and reference database.
the reference database format is currently limited to a custom
style (to be updated and elaborated on later)

```
# install
pip install -e .

# run (writes outputs next to the PSM; use --output-dir to override)
site-annotate run \
  --fasta path/to/proteome.fa \
  --psms  path/to/psms.tsv \
  --cores 8 \
  --refresh-uniprot \
  --output-dir out
```

Outputs include per‑mod reduced files (e.g., `sty_79_9663_site_annotation_reduced.tsv`) and PSP‑mapped files (e.g., `sty_79_9663_site_annotation_reduced_mapped.tsv`). Reduced outputs carry a `uniprot_id` column.

UniProt mapping sources, in order:

- PSM `Protein.Ids` (parses sp|ACC|…/tr|ACC|…)
- ENSP→UniProt (cache/online)
- Fallback via Entrez `geneid`, `symbol` (human), `ENSG`, `ENST`

Use `--refresh-uniprot` to enable online lookups via mygene.info. Otherwise only the local cache is used.

## Cache management

Warm or inspect the local mapping cache:

```
site-annotate cache ls --limit 20 --filter ENSP0000
site-annotate cache rm ENSP00000486824 ENSP00000350398
site-annotate cache clear --yes
site-annotate cache import --file mappings.tsv   # TSV with columns: query, uniprot.Swiss-Prot
```

## PSP data

Place PhosphoSitePlus files under `data/phosphositeplus/`:

- `Phosphosite_seq.fasta` (required for UniProt sequence checks)
- `Phosphorylation_site_dataset`, `Acetylation_site_dataset`, `Regulatory_sites` (optional but recommended)

Mapped files match by UniProt + 15‑mer. When UniProt is missing, mapping falls back to symbol + 15‑mer.

## Doctor summary

Each run writes `site_annotation_run_summary.txt` under the output directory with:

- Aligned PSM rows
- Overall UniProt coverage (PSM)
- Per‑mod rows with non‑null UniProt (reduced)
- Per‑mod rows with PSP `acc_id`/UniProt (mapped)



## Analysis and Reporting


### Merging data
Combining fragpipe site-level expression data with meta data into a single gct file
The resulting gct file is used for creating reports.

`site-annotate merge-meta metadata.tsv /data-dir/`

The `data-dir` should have rec_run_search prefixed site-level data
corresponding to the recno | runno | searchno specified in the metadata file.

so `data-dir` may look like:
`
  data-dir/
    - 12345_1_1_abundance_single-site_MD.tsv
`
the gct file will be saved within the same data directory.

TODO: describe normalization options

### Reporting


`site-annotate report -m metadata.tsv -c config.toml --gct /path/to/gct/file`


## Utilities

- Dry-run preview: validate metadata and matching files without running the full merge.
  `site-annotate dry-run -m metadata.tsv -d /path/to/data_dir [-o outdir] [--latex]`
  Writes a concise text summary (and an optional PDF if LaTeX is available).

- Site vs Protein comparison: compute site/protein ratios from GCT files.
  `site-annotate compare-protein --site-gct site.gct --protein-gct protein.gct [-j auto] [--metric log2ratio] [--write-gct/--no-write-gct] [--write-excel] [-o outdir] [--latex]`
  - Auto-selects the best join key across rdesc columns if `-j auto` (default).
  - Produces ratio and log2ratio TSVs plus a summary, and optionally a GCT/Excel of the chosen metric for downstream analysis.
- Optional per-gene heatmaps: `--gene ABC1 --gene-file genes.txt` creates ComplexHeatmap PDFs (metric-matched) under `genes/`.
  - Generates 4 clustering variants per gene (rcT/F_ccT/F) and, by default, also row-wise z-scored versions (suffix `_z`).
  - Control with `--zscore/--no-zscore`. Discover existing genes from output folders with `--discover-genes`.


# Developer notes

Please keep a running JSONL worklog in `dev_log.jsonl` whenever you modify
the repository.  Each entry should be a single line JSON object containing at
least a timestamp (ISO-8601 preferred), an `event` label, and a short
`details` description of the noteworthy action you performed.  Append new
entries to the file rather than rewriting existing history.  This log is used
by downstream agents to understand recent context, so try to update it after
commits or other significant steps.
## LIMMA Top Tables (optional)

- Render a compact PDF and text summary from one or multiple LIMMA top-table TSVs.
  `site-annotate limma-report --limma-tsv results.tsv --latex --ollama --ollama-model llama3.2:3b`
  Generates a summary, optional LaTeX PDF (if lualatex/pdflatex available), and an optional local LLM summary if `ollama` is present.
