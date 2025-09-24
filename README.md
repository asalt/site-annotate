# site-annot

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

`site-annotate run --uniprot-check -p psm.tsv -f ref.fa`



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

