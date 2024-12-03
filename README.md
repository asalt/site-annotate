# site-annot

A command-line tool for 15mer site rollup and annotation.


this is designed to be modular and easy to modify how filtering and aggregation are performed

# Input file formats

## metadata

a tab separated file with sample information of the form:

| name (must be unique) | recno | runno | searchno | label | group |
| ---- | ----- | ----- | -------- | ----- | ----- |
| ctrl-1 | 123 | 1 | 1 | 126 | ctrl |
| ctrl-2 | 123 | 1 | 1 | 127N | ctrl |
| treat-1 | 123 | 1 | 1 | 126 | treatment |
| treat-2 | 123 | 1 | 1 | 127N | treatment |

The last `group` column is optional but is a default standard variable
name used internally for faceting / grouping.
This can be overwritten by specifying variables in the config file.
Multiple additional columns of metadata can be specified and referenced
within the config file as well.

## config

a base config file can be generated with:
`site-annotate run`

# Running

## site level summarization

note this can be skipped entirely using fragpipe output (described more later)

starting with `psm.tsv` file from fragpipe output
and reference database.
the reference database format is currently limited to a custom
style (to be updated and elaborated on later)

`site-annotate run --uniprot-check -p psm.tsv -f ref.fa`



## report




