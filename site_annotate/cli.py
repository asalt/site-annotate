import sys
import os
import re
import glob
from functools import partial
from pprint import pprint
import subprocess
import functools
import logging
import pathlib
from pathlib import Path
import click
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import pandas as pd
from itertools import product

from rapidfuzz import process, fuzz

import pyfaidx
from pyfaidx import Fasta

import janitor

from . import log
from .io import io
from . import io_external
from .io.io import get_reader, validate_metadata, load_and_validate_files
from .io.io import read_gct
from .io.io import write_gct as write_gct_file, write_excel as write_excel_file
from .compare import compute_site_vs_protein, generate_gene_heatmaps
from .constants import get_all_columns
from . import misc
from . import modisite
from .utils import data_generator

# from .constants import VALID_MODI_COLS, get_all_columns
from .runner import run_pipeline
from . import mapper
from . import reduce
from . import tasks
from .tasks import maybe_compile_latex, compile_latex_from_template, maybe_ollama_summarize
from .spectra_plot import plot_top_psms_from_fragpipe

logger = log.get_logger(__file__)

TEMPLATE_PATH = pathlib.Path(__file__).parent.parent / "scripts"  # 2 levels up
BASE_CONFIG = pathlib.Path(__file__).parent.parent / "config" / "base.toml"


@click.group(chain=True)
def main():
    pass


@main.command()
@click.option("-n", "--name", default="site-annotate.toml")
def get_config(name):
    import shutil

    new_file = pathlib.Path(".").absolute() / name
    if new_file.exists():
        logger.error(f"File {new_file} already exists")
        raise FileExistsError(f"File {new_file} already exists")
    print(f"Writing {new_file} ")
    shutil.copy(BASE_CONFIG, new_file)


def get_templates(TEMPLATE_PATH):
    """
    premade Rmd templates
    """

    if not TEMPLATE_PATH.exists():
        logger.error(f"Template path not found: {TEMPLATE_PATH}")
        raise
    REPORT_TEMPLATES = TEMPLATE_PATH.glob("*Rmd")
    REPORT_TEMPLATES = {x.stem: x for x in REPORT_TEMPLATES}
    # logger.info(REPORT_TEMPLATES)
    return REPORT_TEMPLATES


REPORT_TEMPLATES = get_templates(TEMPLATE_PATH)


# Deprecated flag support: --uniprot-check -> --refresh-uniprot
def _deprecated_uniprot_check(ctx, param, value):
    if value:
        logger.warning("--uniprot-check is deprecated; use --refresh-uniprot")
        # ensure the new flag reflects this
        ctx.params["refresh_uniprot"] = True
    return None


def prepare_params(
    template, config, data_dir, output_dir, metadata, gct, root_dir, save_env
):
    params_dict = {}

    # Resolve paths
    if config:
        params_dict["config_file"] = str(pathlib.Path(config).absolute())
    if gct:
        if not os.path.exists(gct):
            raise FileNotFoundError(f"{gct} does not exist")
        params_dict["gct_file"] = str(pathlib.Path(gct).absolute())
    if data_dir:
        params_dict["data_dir"] = str(pathlib.Path(data_dir).absolute())
    if output_dir is None:
        output_dir = data_dir or "."
    output_dir = pathlib.Path(output_dir).absolute()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    params_dict["output_dir"] = str(output_dir)
    params_dict["root_dir"] = str(root_dir)
    params_dict["save_env"] = str(save_env).lower()

    # Handle metadata
    if metadata and not gct:
        metadata_path = pathlib.Path(metadata).absolute()
        validated_meta = validate_meta(metadata_path, data_dir)
        validated_meta_path = metadata_path.parent / (
            metadata_path.stem + "_validated.tsv"
        )
        validated_meta.to_csv(validated_meta_path, sep="\t", index=False)
        params_dict["metadata"] = str(validated_meta_path)
        logger.info(f"Wrote validated metadata: {validated_meta_path}")

    return params_dict


def _summarize_dry_run(meta_df: pd.DataFrame, found: dict, mappings: dict) -> str:
    lines = []
    lines.append("DRY-RUN SUMMARY")
    lines.append("")
    lines.append(f"metadata rows: {len(meta_df)}")
    lines.append(f"distinct rec_run_search: {meta_df['rec_run_search'].nunique()}")
    lines.append("")

    all_rrs = sorted(meta_df["rec_run_search"].unique())
    missing_rrs = [rrs for rrs in all_rrs if not found.get(rrs)]
    if missing_rrs:
        lines.append(f"missing data files for: {', '.join(missing_rrs)}")
    else:
        lines.append("all rec_run_search have at least one matching data file")

    lines.append("")
    for rrs in all_rrs:
        files = found.get(rrs) or []
        lines.append(f"[ {rrs} ] found files: {len(files)}")
        for f in files:
            m = mappings.get((rrs, f))
            if not m:
                lines.append(f"  - {os.path.basename(f)} : no label mapping (no labels provided?)")
                continue
            expected = len(m)
            matched = sum(1 for v in m.values() if v)
            missing = [k for k, v in m.items() if not v]
            lines.append(
                f"  - {os.path.basename(f)} : mapped {matched}/{expected} labels"
            )
            if missing:
                lines.append(f"      missing: {', '.join(missing)}")
        lines.append("")

    return "\n".join(lines)


@main.command()
@click.option(
    "-e",
    "--extended",
    is_flag=True,
    default=False,
    show_default=True,
    help="print extended info about templates",
)
def show_templates(extended):
    for k, v in REPORT_TEMPLATES.items():
        print(k, v)
    # end


def common_options(f):
    f = click.option(
        "-c",
        "--config",
        type=click.Path(exists=True, dir_okay=False),
        help=".toml file with additional parameters for report",
    )(f)
    f = click.option(
        "-m",
        "--metadata",
        type=click.Path(exists=True, dir_okay=False),
        help="tsv file with rec run search info used to associate with data",
    )(f)
    f = click.option(
        "-o",
        "--output-dir",
        type=click.Path(exists=False, file_okay=False, dir_okay=True),
        required=False,
        default=None,
        show_default=True,
        help="Output directory, defaults to data_dir if not explicitly set",
    )(f)
    f = click.option(
        "-d",
        "--data-dir",
        type=click.Path(exists=True, file_okay=False, dir_okay=True),
        required=True,
        default=pathlib.Path(".").absolute(),
        show_default=True,
    )(f)
    f = click.option(
        "-g",
        "--gct",
        type=click.Path(exists=True, file_okay=True, dir_okay=True),
        required=False,
        default=None,
        show_default=True,
    )(f)
    f = click.option(
        "-r",
        "--root-dir",
        type=click.Path(exists=True, file_okay=False, dir_okay=True),
        required=True,
        default=pathlib.Path(".").absolute(),
        show_default=True,
        help="Root directory, default is location of command invocation",
    )(f)
    f = click.option(
        "--save-env",
        type=click.BOOL,
        required=False,
        default=False,
        show_default=True,
        help="save r environ to output directory",
    )(f)
    return f


# @click.option(
#     "-d",
#     "--data-dir",
#     type=click.Path(exists=True, file_okay=False, dir_okay=True),
#     required=True,
#     default=pathlib.Path(".").absolute(),
#     show_default=True,
# )
# @click.option(
#     "-o",
#     "--
#     type=click.Path(exists=True, file_okay=False, dir_okay=True),
#     required=False,
#     default=None,
#     show_default=True,
#     help="Output directory, defaults to data_dir if not explictely set",
# )
# @click.option("-m", "--metadata", type=click.Path(exists=True, dir_okay=False))


def conf_to_dataframe(conf_file, set_index=False):
    # this should be easy to test

    # Initialize config parser
    import configparser

    config = configparser.ConfigParser()
    config.read(conf_file)

    # Prepare dictionary for storing data
    data = {}

    # Loop through each section to create rows
    for section in config.sections():
        section_data = config[section]
        row = {key: section_data[key] for key in section_data}
        row["name"] = section  # Include section name as 'name' column
        data[section] = row

    # Convert dictionary to DataFrame
    df = pd.DataFrame.from_dict(data, orient="index")

    # Optional: Set 'name' as index
    if set_index:
        df.set_index("name", inplace=True)

    return df


def validate_meta(
    metadata_file: pathlib.Path, data_dir: pathlib.Path, ensure_data=True
) -> dict:
    """
    top level function to read metadata file and validate it
    """

    reader = io.get_reader(str(metadata_file))
    meta_df = reader(metadata_file, dtype={"label": "str"})
    meta_df = io.validate_metadata(meta_df)  # or .pipe

    rec_run_searches = meta_df.rec_run_search.unique()
    expr_data = io.find_expr_files(rec_run_searches, data_dir)
    expr_data = {k: v for (k, v) in expr_data.items() if v}
    if len(expr_data) == 0 and not ensure_data:  # no data
        return

    meta_df_final_collection = io.validate_expr_files(expr_data, meta_df)

    logger.info(f"successfully validated {metadata_file}")

    return meta_df_final_collection


@main.command(name="plot-spectra")
@click.option(
    "--psm-tsv",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="FragPipe psm.tsv file to select PSMs from.",
)
@click.option(
    "--raw-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="Directory containing mzML files (e.g. FragPipe raw/IMAC).",
)
@click.option(
    "-o",
    "--out-dir",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default="spectra",
    show_default=True,
    help="Output directory for spectrum PNGs.",
)
@click.option(
    "-n",
    "--n-top",
    type=int,
    default=10,
    show_default=True,
    help="Number of top PSMs to plot, sorted by Hyperscore.",
)
@click.option(
    "--qvalue-cutoff",
    type=float,
    default=0.01,
    show_default=True,
    help="Maximum Qvalue for PSM selection.",
)
@click.option(
    "--min-localization",
    type=float,
    default=0.0,
    show_default=True,
    help="Minimum localization probability for phospho sites.",
)
@click.option(
    "--phospho-loc-col",
    type=str,
    default="STY:79.9663 Best Localization",
    show_default=True,
    help="Column name used for phospho localization score.",
)
@click.option(
    "--annotate/--no-annotate",
    default=False,
    show_default=True,
    help="Use spectrum_utils styling for spectra (no fragment labels yet).",
)
@click.option(
    "--mz-min",
    type=float,
    default=None,
    show_default=True,
    help="Minimum m/z to display (applied to all plots).",
)
@click.option(
    "--mz-max",
    type=float,
    default=None,
    show_default=True,
    help="Maximum m/z to display (applied to all plots).",
)
@click.option(
    "--max-labels-plain",
    "--plain-max-labels",
    type=str,
    default="auto",
    show_default=True,
    help="For plain plots, maximum number of peaks to label with m/z; use 'auto' for heuristic labeling.",
)
@click.option(
    "--max-aa-length",
    type=int,
    default=None,
    show_default=True,
    help="Only consider PSMs with Peptide Length <= this value when selecting top PSMs.",
)
@click.option(
    "--html/--no-html",
    default=False,
    show_default=True,
    help="Also write interactive HTML spectra (annotated only).",
)
@click.option(
    "--verbose/--no-verbose",
    default=False,
    show_default=True,
    help="Print progress information while reading mzML and plotting spectra.",
)
def plot_spectra(
    psm_tsv,
    raw_dir,
    out_dir,
    n_top,
    qvalue_cutoff,
    min_localization,
    phospho_loc_col,
    annotate,
    mz_min,
    mz_max,
    max_labels_plain,
    max_aa_length,
    html,
    verbose,
):
    """
    Plot top-scoring phospho MS/MS spectra from a FragPipe psm.tsv file.

    Examples:

      site-annotate plot-spectra \\
        --psm-tsv path/to/psm.tsv \\
        --raw-dir path/to/raw/IMAC \\
        --out-dir spectra \\
        --n-top 4
    """
    psm_path = pathlib.Path(psm_tsv).absolute()
    raw_dir_path = pathlib.Path(raw_dir).absolute()
    out_dir_path = pathlib.Path(out_dir).absolute()
    out_dir_path.mkdir(parents=True, exist_ok=True)

    # interpret max-labels-plain
    plain_max_labels = None
    if isinstance(max_labels_plain, str):
        if max_labels_plain.lower() != "auto":
            try:
                plain_max_labels = int(max_labels_plain)
            except ValueError:
                raise click.BadParameter(
                    f"--max-labels-plain must be an integer or 'auto', got {max_labels_plain!r}"
                )
    else:
        plain_max_labels = max_labels_plain

    written = plot_top_psms_from_fragpipe(
        psm_path=psm_path,
        raw_dir=raw_dir_path,
        out_dir=out_dir_path,
        n_top=n_top,
        qvalue_cutoff=qvalue_cutoff,
        min_localization=min_localization,
        phospho_loc_col=phospho_loc_col,
        annotated=annotate,
        mz_min=mz_min,
        mz_max=mz_max,
        max_labels=plain_max_labels,
        max_length=max_aa_length,
        write_html=html,
        verbose=verbose,
    )
    logger.info(f"Wrote {len(written)} spectra to {out_dir_path}")


@main.command()
@click.argument("metadata", type=click.Path(exists=True, dir_okay=False))
@click.argument("data_dir", type=click.Path(exists=True, file_okay=False))
def check_meta(metadata, data_dir):
    """
    used to check if metadata tabular file matches with tmt-integrator output file
    used internally in `merge_meta`
    looks for rec_run_search string in metadtaa file (e.g. 12345_1_1) and matching
    files (e.g. 12345_1_1_abundance_single-site_MD.tsv)
    a searchno and runno are not strictly necessary in the metadata, and will take on default values(1, 7) if not provided.
    but critically, the `rec_runno_searchno` pattern string needs to match
    TODO: allow partial matches and let runno and searchno take on default values if not provided
    """
    meta_validated_collection = validate_meta(metadata, data_dir)
    for meta_validated in meta_validated_collection:
        print(meta_validated.to_string())
    print("Success")


@main.command()
@click.option("--cores", default=1, show_default=True)
@click.option("--taxon", type=str, default=None, show_default=True)
@click.option(
    "--make-nr/--no-make-nr", type=bool, default=True, is_flag=True, show_default=True
)
@click.option(
    "--result-dir", type=click.Path(exists=False, file_okay=False, dir_okay=True)
)
@click.argument("metadata", type=click.Path(exists=True, dir_okay=False))
@click.argument(
    "data_dir",
    type=click.Path(exists=True, file_okay=False),
    default=click.Path("."),
    # show_default=True,
)
def merge_meta(cores, taxon, make_nr, result_dir, metadata, data_dir):
    r"""
    [DATA_DIR] is the directory where the data files are located

    default is current directory

    Merge expression matrix from fragpipe output:
    starting with this file:
        `tmt-report/abundance_single-site_MD.tsv`
        and a config file with recno 12345
    add the recno prefix and run:
        `tmt-report/12345_abundance_single-site_MD.tsv`
    and run:
        `site-annotate merge-meta ./path/to/meta.tsv ./siteinfo-dir`
    which will produce a gct file with metadata

    Note only works with a single "plex". (not true anymore)
    recno information is only available from the prefix (12345) appended to the front of the emat

    """
    data_dir = pathlib.Path(data_dir).absolute()

    # if result_dir is None:
    result_dir = result_dir or data_dir.parent
    result_dir = pathlib.Path(result_dir).absolute()

    metadata = pathlib.Path(metadata).absolute()
    # can add try/except here to provide better error msg
    meta_validated_collection = validate_meta(metadata, data_dir)

    for name, meta_validated in meta_validated_collection.items():
        meta_validated_fname = metadata.parent / (
            metadata.stem + f"_{name}_validated.tsv"
        )
        logger.info(f"writing {meta_validated_fname}")
        meta_validated.to_csv(meta_validated_fname, sep="\t", index=False)

        taxon_str, nr_str = "", ""
        if taxon is not None:
            taxon_str = f"_{taxon}"
        if make_nr:
            nr_str = "_nr"

        outname = result_dir / (
            metadata.stem + f"_{name}_siteinfo_combined{taxon_str}{nr_str}"
        )
        # ===
        logger.info(f"preparing {outname}")

        res = io.merge_metadata(
            meta_validated, outname, taxon=taxon, make_nr=make_nr, cores=cores
        )  # writes gct file to output


def match_columns(sample_names: list, df: pd.DataFrame, threshold=97):
    """
    Matches `sample` column values from annotation_df to df.columns
    using fuzzy string matching with preprocessing.

    :param annotation_df: DataFrame containing a "sample" column.
    :param df: DataFrame whose columns need to be matched.
    :param threshold: Similarity threshold for a confident match (default=97%).
    :return: A dictionary mapping `sample` values to best matches in `df.columns`.
    """
    # Extract and preprocess column names and samples to match

    def preprocess_string(s):
        """
        Preprocesses a string by removing spaces, ...,  and normalizing case.
        """
        return re.sub(r"[\s]+", "", s).lower()

    samples = [preprocess_string(sample) for sample in sample_names]
    columns = [preprocess_string(col) for col in df.columns]

    # Dictionary to store matches
    matches = {}

    # Perform fuzzy matching
    for original_sample, processed_sample in zip(sample_names, samples):
        # Get best match using token_sort_ratio for better handling of word order
        match, score, _ = process.extractOne(
            processed_sample, columns, scorer=fuzz.token_sort_ratio
        )

        if (
            score / 100 >= threshold / 100
        ):  # Normalize score to 0-1 range for comparison
            # Find the original column name corresponding to the match
            original_match = df.columns[columns.index(match)]
            matches[original_sample] = original_match
        else:
            logger.warning(
                f"No high-confidence match for '{original_sample}' (Best: '{match}' with {score}%)"
            )
            matches[original_sample] = None  # No confident match

    return matches


@main.command()
@click.option(
    "-e",
    "--experiment-annotation",
    type=click.Path(exists=True, file_okay=True),
    required=True,
    help="Tabular file with columns: plex, channel, sample.",
)
@click.argument(
    "target",
    type=click.Path(exists=True, file_okay=True),
    default=".",
)
def add_experiment_annotation(experiment_annotation, target):
    """
    (experimental)
    For FragPipe TMT-integrator output.
    Merges FragPipe TMT-integrator output with an experiment annotation file.

    [experiment-annotation] is a tabular file generated from FragPipe.

    [TARGET] target file or directory

    :experiment-annotation:
    expected columns are : `plex`, `channel` and `sample`
    other columns such as `condition` and `replicate`, which come from FragPipe, are not used
    information about the "recno" should either be extractable from the `plex` column
    or a separate `recno` column can be manually added

    :target:
    example tmt-report/
    -> reads in abundance_xx_yy.tsv
        xx is a "level" such as gene, signle_site, protein.
            we are primarily interested in gene and signle_site
        we compare all of this to the config file, experiment_annotation)
        and attempt to merge.



    """
    # Validate inputs
    try:
        reader = get_reader(str(experiment_annotation))
        annotation_df = reader(experiment_annotation)
    except Exception as e:
        raise click.ClickException(f"Error reading annotation file: {e}")

    required_columns = {"plex", "channel", "sample"}
    if not required_columns.issubset(annotation_df.columns):
        raise click.ClickException(
            f"Annotation file missing required columns: {required_columns - set(annotation_df.columns)}"
        )

    click.echo(f"Processing target: {target}")
    click.echo(f"Loaded annotation file with {len(annotation_df)} rows.")

    # Placeholder for additional processing
    click.echo("Processing...")

    def search_files(base_dir, middle_part="gene"):
        """
        Searches for all files in the base_dir that match the middle_part pattern.

        :param base_dir: Directory to search in.
        :param middle_part: Middle part of the filename to match.
        :return: List of matching file paths.
        """
        # Define patterns
        prefixes = ["abundance", "ratio"]
        suffixes = ["GN.tsv", "MD.tsv", "None.tsv"]

        # Generate all possible filenames using product
        filenames = [
            f"{prefix}_{middle_part}_{suffix}"
            for prefix, suffix in product(prefixes, suffixes)
        ]

        # Check for file existence
        matching_files = [
            str(Path(base_dir) / filename)
            for filename in filenames
            if (Path(base_dir) / filename).exists()
        ]

        return matching_files

    files = search_files(target)

    for file in files:
        print(file)
        df = get_reader(file)(file)
        cid_cols = set(annotation_df["sample"]) & set(df.columns)
        # these should match
        # now do the merge and save the results
        sample_mapping = match_columns(annotation_df["sample"], df)
        # normalize names
        df = df.rename(columns={v: k for k, v in sample_mapping.items()})

        from .io import write_gct

        write_gct(
            df,  # this is not done yet
            cdesc=annotation_df,
        )


@main.command(name="dry-run")
@click.option("-m", "--metadata", type=click.Path(exists=True, dir_okay=False), required=True)
@click.option("-d", "--data-dir", type=click.Path(exists=True, file_okay=False, dir_okay=True), required=True)
@click.option("-o", "--output-dir", type=click.Path(exists=False, file_okay=False, dir_okay=True), default=None, show_default=True)
@click.option("-t", "--threshold", type=int, default=97, show_default=True, help="Fuzzy match confidence threshold")
@click.option("--latex/--no-latex", is_flag=True, default=False, show_default=True, help="Attempt to write a PDF summary if LaTeX is available")
def dry_run(metadata, data_dir, output_dir, threshold, latex):
    """
    Preview metadata and expression file alignment, and summarize label-to-column mapping.
    Writes a text summary; optionally compiles a PDF if lualatex/pdflatex is available.
    """
    data_dir = pathlib.Path(data_dir).absolute()
    metadata_path = pathlib.Path(metadata).absolute()
    if output_dir is None:
        output_dir = data_dir
    outdir = pathlib.Path(output_dir).absolute()
    outdir.mkdir(parents=True, exist_ok=True)

    # Load and standardize metadata
    reader = io.get_reader(str(metadata_path))
    meta_df = reader(metadata_path, dtype={"label": "str"})
    meta_df = io.validate_metadata(meta_df)

    rec_run_searches = meta_df.rec_run_search.unique().tolist()
    found = io.find_expr_files(rec_run_searches, str(data_dir))

    # For each found file, attempt label mapping
    mappings = {}
    for rrs, files in found.items():
        if not files:
            continue
        meta_rrs = meta_df[meta_df.rec_run_search == rrs]
        if "label" not in meta_rrs.columns or meta_rrs["label"].isna().all():
            for f in files:
                mappings[(rrs, f)] = None
            continue

        labels = meta_rrs["label"].dropna().astype(str).tolist()
        for f in files:
            try:
                df_head = get_reader(f)(f, nrows=5)
            except Exception:
                mappings[(rrs, f)] = {lab: None for lab in labels}
                continue
            # use fuzzy matching helper
            m = match_columns(labels, df_head, threshold=threshold)
            mappings[(rrs, f)] = m

    summary = _summarize_dry_run(meta_df, found, mappings)
    txt_path = outdir / "dry_run_summary.txt"
    txt_path.write_text(summary)
    click.echo(f"Wrote summary: {txt_path}")

    if latex:
        pdf_path = maybe_compile_latex(summary, str(outdir), basename="dry_run_summary", title="site-annot dry-run")
        if pdf_path:
            click.echo(f"Wrote PDF: {pdf_path}")
        else:
            click.echo("LaTeX engine not available or failed; skipped PDF.")


@main.command()
@common_options
@click.option(
    "-t",
    "--template",
    type=click.Choice(REPORT_TEMPLATES.keys()),
    required=False,
    default="analyze_modi",
    show_default=True,
)
@click.option("--interactive", default=False, show_default=True, is_flag=True)
@click.option("--save-env", default=False, show_default=True)
def report(
    template,
    config,
    data_dir,
    output_dir,
    metadata,
    gct,
    root_dir,
    interactive,
    save_env,
    **kwargs,
):

    params_dict = prepare_params(
        template, config, data_dir, output_dir, metadata, gct, root_dir, save_env
    )
    tasks.run_r_code_with_params(params_dict, interactive=interactive)
    return

    print(template)

    template_file = REPORT_TEMPLATES[template]

    params_dict = dict()

    if config is not None:
        config = pathlib.Path(config).absolute()
        params_dict["config"] = str(config)

    if gct is not None:
        if not os.path.exists(gct):
            logger.error(f"{gct} does not exist")
            raise FileNotFoundError(f"{gct} does not exist")
        params_dict["gct"] = pathlib.Path(gct).absolute()

    if data_dir is not None:
        data_dir = pathlib.Path(data_dir).absolute()
        params_dict["data_dir"] = str(data_dir)

    if metadata is not None and gct is None:
        metadata = pathlib.Path(metadata).absolute()
        # now we check if metadata is aligned with available experimental data
        meta_validated = validate_meta(metadata, data_dir)
        meta_validated_fname = metadata.parent / (metadata.stem + "_validated.tsv")
        meta_validated.to_csv(meta_validated_fname, sep="\t", index=False)

        # ===
        logger.info(f"wrote {meta_validated_fname}")
        params_dict["metadata"] = str(meta_validated_fname)
    if output_dir is None:
        output_dir = data_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_dir = str(pathlib.Path(output_dir).absolute())
    #

    # Create a dictionary with all parameters
    for k, v in kwargs.items():
        params_dict[k] = v
    params_dict.update(
        {
            # "config": str(config),
            # "data_dir": str(data_dir),
            # "metadata": str(meta_validated_fname),
            "output_dir": str(output_dir),
            "root_dir": str(root_dir),
            "save_env": str.lower(str(save_env)),
        }
    )

    if (
        template == "merge_modis" or template == "heatmap"
    ):  # eventually will move this out
        params = (
            "list(" + ", ".join(f"'{k}' = '{v}'" for k, v in params_dict.items()) + ")"
        )
        cmd = [
            f"Rscript",
            "-e",
            f"""library(rmarkdown)
            rmarkdown::render("{template_file}",
            output_dir="{output_dir}",
            params={params},
            )""",
        ]
        logger.info(cmd)
        subprocess.run(cmd)
        return

    import rpy2.robjects as robjects

    for k, v in params_dict.items():
        robjects.r.assign(k, str(v))

    runfile = pathlib.Path(__file__).parent.parent / "R/run.R"
    robjects.r.assign("here_dir", str(runfile.parent))
    robjects.r("setwd(here_dir)")
    robjects.r.source(str(runfile))

    os.chdir(str(runfile.parent))

    robjects.r(
        """
        data_dir <- ifelse(exists("data_dir"), data_dir, '.')
        output_dir <- ifelse(exists("output_dir"), output_dir, '.')
        config_file <- ifelse(exists("config"), config, NA) # have to make it NA not null
        gct_file <- ifelse(exists("gct"), gct, NA) # have to make it NA not null
        save_env <- ifelse(exists("save_env"), save_env, FALSE)

        print(paste0('output dir is: ', output_dir))
        run(
            data_dir = data_dir,
            output_dir = output_dir,
            config_file = config_file,
            gct_file = gct_file,
            save_env = save_env,
               )
               """
    )
    return


# we are not using this
# can probably remove
@main.command()
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=False,
    default=None,
    show_default=True,
    help="Output directory, defaults to data_dir if not explicitly set",
)
@click.option(
    "-d",
    "--data-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    default=pathlib.Path(".").absolute(),
    show_default=True,
)
@click.option("--modi-abbrev", default="p", show_default=True)
def reduce_sites(output_dir, data_dir, modi_abbrev, **kwargs):
    raise ValueError("depreciated")

    template_file = REPORT_TEMPLATES["site_annotate_post_not_reduced"]

    data_dir = pathlib.Path(data_dir).absolute()
    if output_dir is None:
        output_dir = data_dir

    params_dict = {
        "data_dir": str(data_dir),
        "modi_abbrev": modi_abbrev,
    }
    params = "list(" + ", ".join(f"'{k}' = '{v}'" for k, v in params_dict.items()) + ")"
    cmd = [
        f"Rscript",
        "-e",
        f"""library(rmarkdown)
        rmarkdown::render("{template_file}",
        output_dir="{output_dir}",
        params={params},
        )""",
    ]

    logger.info(cmd)
    subprocess.run(cmd)


@main.command()
@click.option("--cores", default=1, show_default=True)
@click.option(
    "-p",
    "--psms",
    type=click.Path(exists=True, dir_okay=False),
    multiple=True,
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(exists=False, file_okay=False),
    default=None,
    show_default=True,
)
@click.option(
    "--prefix",
    type=str,
    default=None,
    show_default=True,
    help="prefix to append to beginning of files",
)
@click.option(
    "--refresh-uniprot/--no-refresh-uniprot",
    default=False,
    show_default=True,
    is_flag=True,
    help="Perform an additional UniProt refresh pass (cache/online if available). Initial mapping from PSM/FASTA/cache always runs.",
)
@click.option(
    "--uniprot-check/--no-uniprot-check",
    is_flag=True,
    default=False,
    expose_value=False,
    callback=_deprecated_uniprot_check,
    help="DEPRECATED: use --refresh-uniprot",
)
@click.option(
    "-f", "--fasta", type=click.Path(exists=True, dir_okay=False), help="fasta file"
)
def run(cores, psms, output_dir, prefix, refresh_uniprot, fasta):
    """
    Build site-level tables from a PSM file and a FASTA.

    - Parses site probabilities into per-mod DataFrames (e.g., sty_79_9663, m_15_9949).
    - Aligns proteins to FASTA headers; keeps uniprot_id in reduced outputs.
    - UniProt mapping sources (in order): Protein.Ids, ENSP (cache/online),
      geneid (Entrez), symbol (human), ENSG, ENST. Use --refresh-uniprot to enable
      online lookups; otherwise cached/offline fills only.
    - Writes both reduced and PSP-mapped files; mapped uses UniProt+15mer or
      symbol+15mer when UniProt is missing.
    """
    if output_dir is not None:
        outdir_path = pathlib.Path(output_dir)
        outdir_path.mkdir(parents=True, exist_ok=True)
    else:
        outdir_path = None
    if not psms:
        logger.warning("Warning: No PSM files provided.")
        return
    if len(psms) > 1:
        logger.error("more than 1 psms not supported yet")
        raise NotImplementedError("more than 1 psms not supported yet")

    logger.info("Loading and validating inputs…")
    df, fasta_data, fa_psp_ref = load_and_validate_files(psms[0], fasta, refresh_uniprot)
    logger.info("Aligned PSM rows: %d", len(df))
    logger.info("Running site extraction across proteins (cores=%s)…", cores)
    fullres = run_pipeline(df, fasta_data, fa_psp_ref, cores)
    if not fullres:
        # Produce a summary and a placeholder output instead of exiting silently
        logger.error(
            "No site-level results produced. Check that PSMs contain modification columns (e.g., sty_79_9663) and that protein IDs align to FASTA keys."
        )
        # Write a brief run summary to the chosen outdir or PSM folder
        base_dir = outdir_path if outdir_path is not None else pathlib.Path(psms[0]).parent
        base_dir.mkdir(parents=True, exist_ok=True)
        summary = base_dir / f"{prefix or ''}site_annotation_run_summary.txt"
        mods_present = sorted(get_all_columns(df.columns))
        with open(summary, "w") as h:
            h.write("site-annotate run summary\n")
            h.write(f"psm: {psms[0]}\n")
            h.write(f"fasta: {fasta}\n")
            h.write(f"aligned_rows: {len(df)}\n")
            h.write(f"mods_in_psm_columns: {', '.join(mods_present) if mods_present else 'none'}\n")
            if "uniprot_id" in df.columns:
                h.write(f"uniprot_non_null: {int(df['uniprot_id'].notna().sum())}\n")
        # Fallback: write basic per-mod site tables without requiring FASTA alignment
        mods = [m for m in mods_present if m in df.columns and df[m].notna().any()]
        if mods:
            logger.info("Writing basic per-mod outputs without FASTA alignment")
            basic = {}
            for col in mods:
                res = misc.extract_and_transform_data(df, col)
                if res is None or res.empty:
                    continue
                merged = res.merge(df, left_on="original_index", right_index=True)
                # Minimal ordering: peptide, protein, AA, position_relative, prob + useful IDs
                keep = [x for x in ("peptide","protein","AA","position_relative","prob","modified_peptide","spectrum","intensity","uniprot_id","ENSP","ENST","ENSG","geneid","symbol") if x in merged.columns]
                basic[col] = merged[keep]
            if basic:
                save_results(basic, psms[0], name="site_annotation_basic_no_fasta", prefix=prefix, outdir=base_dir)
                logger.info("Wrote basic per-mod outputs (no FASTA alignment)")
        logger.info(f"Wrote summary: {summary}")
        return
    # type fullres is list of dicts
    # fullres[0].keys()
    # dict_keys(['sty_79_9663'])
    # fullres[0]['sty_79_9663'].columns
    # Index(['position_relative', 'AA', 'prob', 'original_index', 'spectrum', 'spectrum_file', 'peptide', 'modified_peptide',
    #     'extended_peptide', 'prev_aa', 'next_aa', 'peptide_length', 'charge', 'retention', 'observed_mass', 'calibrated_observed_mass',
    #     'observed_m_z', 'calibrated_observed_m_z', 'calculated_peptide_mass', 'calculated_m_z', 'delta_mass', 'spectralsim', 'rtscore',
    #     'expectation', 'hyperscore', 'nextscore', 'peptideprophet_probability', 'number_of_enzymatic_termini',
    #     'number_of_missed_cleavages', 'protein_start', 'protein_end', 'intensity', 'assigned_modifications', 'observed_modifications',
    #     'm_15_9949', 'm_15_9949_best_localization', 'sty_79_9663', 'sty_79_9663_best_localization', 'purity', 'is_unique', 'protein',
    #     'protein_id', 'entry_name', 'gene', 'protein_description', 'mapped_genes', 'mapped_proteins', 'TMT_126', 'TMT_127_N',
    #     'TMT_127_C', 'TMT_128_N', 'TMT_128_C', 'TMT_129_N', 'TMT_129_C', 'TMT_130_N', 'TMT_130_C', 'TMT_131_N', 'TMT_131_C',
    #     'TMT_132_N', 'TMT_132_C', 'TMT_133_N', 'TMT_133_C', 'TMT_134_N', 'position_absolut', 'fifteenmer', 'protein_length', 'tmt_sum',
    #     'TMT_126_ratio', 'TMT_127_N_ratio', 'TMT_127_C_ratio', 'TMT_128_N_ratio', 'TMT_128_C_ratio', 'TMT_129_N_ratio',
    #     'TMT_129_C_ratio', 'TMT_130_N_ratio', 'TMT_130_C_ratio', 'TMT_131_N_ratio', 'TMT_131_C_ratio', 'TMT_132_N_ratio',
    #     'TMT_132_C_ratio', 'TMT_133_N_ratio', 'TMT_133_C_ratio', 'TMT_134_N_ratio', 'TMT_126_intensity', 'TMT_127_N_intensity',
    #     'TMT_127_C_intensity', 'TMT_128_N_intensity', 'TMT_128_C_intensity', 'TMT_129_N_intensity', 'TMT_129_C_intensity',
    #     'TMT_130_N_intensity', 'TMT_130_C_intensity', 'TMT_131_N_intensity', 'TMT_131_C_intensity', 'TMT_132_N_intensity',
    #     'TMT_132_C_intensity', 'TMT_133_N_intensity', 'TMT_133_C_intensity', 'TMT_134_N_intensity', 'highest_prob'],
    #     dtype='object')
    finalres = process_results(fullres, decoy_label="rev_")
    if not finalres:
        logger.error("No results after postprocessing; nothing to write")
        base_dir = outdir_path if outdir_path is not None else pathlib.Path(psms[0]).parent
        base_dir.mkdir(parents=True, exist_ok=True)
        summary = base_dir / f"{prefix or ''}site_annotation_run_summary.txt"
        mods_present = sorted(get_all_columns(df.columns))
        with open(summary, "w") as h:
            h.write("site-annotate run summary\n")
            h.write(f"psm: {psms[0]}\n")
            h.write(f"fasta: {fasta}\n")
            h.write(f"aligned_rows: {len(df)}\n")
            h.write(f"mods_in_psm_columns: {', '.join(mods_present) if mods_present else 'none'}\n")
            if "uniprot_id" in df.columns:
                h.write(f"uniprot_non_null: {int(df['uniprot_id'].notna().sum())}\n")
        placeholder = base_dir / f"{prefix or ''}site_annotation_NO_RESULTS.txt"
        with open(placeholder, "w") as h:
            h.write("No site-level results after postprocessing (e.g., only decoys). See run summary for details.\n")
        logger.info(f"Wrote summary: {summary}")
        logger.info(f"Wrote placeholder: {placeholder}")
        return
    logger.info("Mods with results: %s", ", ".join(sorted(finalres.keys())))
    # finalres is a dict with keys modi ( sty_79_9663, k_42_0106, ... )
    # and values the concatenated dataframe of all sites, along with tmt quant if applicable
    # this "not reduced" file can be huge ( 15g + )
    # save_results(finalres, psms[0])

    site_reduced = reduce.reduce_sites(
        finalres,
        decoy_label="rev_",
    )

    save_results(site_reduced, psms[0], name="site_annotation_reduced", prefix=prefix, outdir=outdir_path)

    site_reduced_mapped = mapper.add_annotations(site_reduced)
    save_results(site_reduced_mapped, psms[0], name="site_annotation_reduced_mapped", prefix=prefix, outdir=outdir_path)

    # Doctor summary: write a concise run report and log it
    try:
        summary_path = _write_doctor_summary(df, site_reduced, site_reduced_mapped, outdir_path or pathlib.Path(psms[0]).parent, fasta)
        logger.info("Wrote run summary: %s", summary_path)
    except Exception as exc:
        logger.warning("Failed to write run summary: %s", exc)

    return


def _write_doctor_summary(df, site_reduced, site_reduced_mapped, outdir: pathlib.Path, fasta_path: str) -> pathlib.Path:
    outdir.mkdir(parents=True, exist_ok=True)
    lines = []
    lines.append("site-annotate doctor summary")
    lines.append("")
    lines.append(f"aligned_rows: {len(df)}")
    uni_nonnull = int(df['uniprot_id'].notna().sum()) if 'uniprot_id' in df.columns else 0
    lines.append(f"uniprot_non_null_psm: {uni_nonnull}")
    # Per-mod stats for reduced
    lines.append("")
    lines.append("per-mod reduced stats:")
    for mod, rdf in sorted(site_reduced.items()):
        total = len(rdf)
        nonnull = int(rdf['uniprot_id'].notna().sum()) if 'uniprot_id' in rdf.columns else 0
        lines.append(f"  {mod}: rows={total}, uniprot_non_null={nonnull}")
    # Per-mod stats for mapped
    lines.append("")
    lines.append("per-mod mapped stats:")
    for mod, mdf in sorted(site_reduced_mapped.items()):
        total = len(mdf)
        uni_non = int(mdf['uniprot_id'].notna().sum()) if 'uniprot_id' in mdf.columns else 0
        acc_non = int(mdf['acc_id'].notna().sum()) if 'acc_id' in mdf.columns else 0
        lines.append(f"  {mod}: rows={total}, uniprot_non_null={uni_non}, psp_acc_non_null={acc_non}")
    # PSP data presence
    lines.append("")
    lines.append(f"fasta: {fast_a_path}")
    lines.append(f"psp_dir: {io_external.data_dir}")
    # Write file
    summary = outdir / "site_annotation_run_summary.txt"
    summary.write_text("\n".join(lines))
    return summary


def process_results(fullres, decoy_label="rev_"):
    """
    fullres is a list of dicts
    """
    finalres = {}
    all_frames = defaultdict(list)
    for result_dict in fullres:
        for modi_id, df in result_dict.items():
            all_frames[modi_id].append(df)
        # for col in VALID_MODI_COLS:
        #     frames = [items.get(col) for items in fullres if items.get(col) is not None]
        # if frames:
        #     finalres[col] = pd.concat(frames, ignore_index=True)
        # else:
        #     pass  # not all need to be present
        # logger.warning(f"No results returned for {col}")
    for k, v in all_frames.items():
        _df = pd.concat(v)  # filter decoys here
        _df = _df[~_df.protein.str.startswith(decoy_label)]
        finalres[k] = _df
    return finalres


def save_results(finalres, input_file, name="site_annotation_notreduced", prefix=None, outdir: pathlib.Path | None = None):
    if prefix is None:
        prefix = ""
    infile = pathlib.Path(input_file)
    for key, val in finalres.items():
        base_dir = outdir if outdir is not None else infile.parent
        outfile = base_dir / f"{prefix}{key}_{name}.tsv"
        logger.info(f"Writing {outfile}")
        val.to_csv(outfile, sep="\t", index=False)


@main.command(name="compare-protein")
@click.option("--site-gct", type=click.Path(exists=True, dir_okay=False), required=True, help="Site-level GCT file")
@click.option("--protein-gct", type=click.Path(exists=True, dir_okay=False), required=True, help="Protein/Gene-level GCT file")
@click.option("-j", "--join-on", type=str, default="auto", show_default=True, help="Join key or 'auto' to select best match (e.g., symbol, ENSG, geneid)")
@click.option("-o", "--output-dir", type=click.Path(exists=False, file_okay=False, dir_okay=True), default=None, show_default=True)
@click.option("--agg", type=click.Choice(["mean", "median", "sum"]), default="mean", show_default=True, help="Aggregation across multiple protein rows per join key")
@click.option("--metric", type=click.Choice(["ratio", "log2ratio"]), default="log2ratio", show_default=True, help="Which matrix to emit as GCT/Excel when writing files")
@click.option("--write-gct/--no-write-gct", "write_gct_flag", is_flag=True, default=True, show_default=True, help="Write output matrix to GCT as well as TSV")
@click.option("--write-excel/--no-write-excel", "write_excel_flag", is_flag=True, default=False, show_default=True, help="Write output matrix to Excel with row metadata")
@click.option("--drop-unmatched/--keep-unmatched", is_flag=True, default=False, show_default=True, help="Drop unmatched site rows from outputs")
@click.option("--eps", type=float, default=1e-6, show_default=True, help="Small constant to avoid division by zero")
@click.option("--gene", "genes", multiple=True, help="Gene/geneID identifiers to plot heatmaps for")
@click.option("--genes-file", type=click.Path(exists=True, dir_okay=False), help="File containing gene IDs (one per line) to plot")
@click.option("--gene-col", type=str, default=None, help="Explicit rdesc column to match genes (case-insensitive)")
@click.option("--gene-debug/--no-gene-debug", is_flag=True, default=False, show_default=True, help="Print per-gene match diagnostics")
@click.option("--force/--no-force", is_flag=True, default=False, show_default=True, help="Overwrite existing outputs (both data and plots)")
@click.option("--force-plots/--no-force-plots", "force_plots", is_flag=True, default=False, show_default=True, help="Overwrite existing plot PDFs")
@click.option("--force-data/--no-force-data", "force_data", is_flag=True, default=False, show_default=True, help="Overwrite existing data files (TSVs, GCT/Excel)")
@click.option("--zscore/--no-zscore", is_flag=True, default=True, show_default=True, help="Also render row-wise z-scored heatmaps")
@click.option("--discover-genes/--no-discover-genes", is_flag=True, default=False, show_default=True, help="Discover gene names from existing subfolders under outdir/genes")
@click.option("--latex/--no-latex", is_flag=True, default=False, show_default=True, help="Attempt to write a PDF summary if LaTeX is available")
def compare_protein(site_gct, protein_gct, join_on, output_dir, agg, metric, write_gct_flag, write_excel_flag: bool, drop_unmatched, eps, genes, genes_file, gene_col, gene_debug, force, force_plots, force_data, zscore, discover_genes, latex):
    """
    Compare site-level intensities to protein/gene-level abundances by computing per-sample ratios and log2 ratios.
    """
    site_gct = pathlib.Path(site_gct).absolute()
    protein_gct = pathlib.Path(protein_gct).absolute()
    outdir = pathlib.Path(output_dir or site_gct.parent).absolute()
    outdir.mkdir(parents=True, exist_ok=True)

    gene_list = [g for g in genes if g]
    if genes_file:
        with open(genes_file, "r", encoding="utf-8") as handle:
            extra_genes = [line.strip() for line in handle if line.strip()]
        gene_list.extend(extra_genes)
    if gene_list:
        gene_list = list(dict.fromkeys(gene_list))
    elif discover_genes:
        # Discover gene names from outdir/genes/* subdirectories
        genes_dir = outdir / "genes"
        if genes_dir.exists() and genes_dir.is_dir():
            subdirs = [p.name for p in genes_dir.iterdir() if p.is_dir() and not p.name.startswith(".")]
            if subdirs:
                gene_list = sorted(subdirs)

    # Resolve force behavior
    eff_force_plots = bool(force or force_plots)
    eff_force_data = bool(force or force_data)

    emat_s, cdesc_s, rdesc_s = read_gct(str(site_gct))
    emat_p, cdesc_p, rdesc_p = read_gct(str(protein_gct))

    try:
        ratio_df, log2ratio_df, rdesc_out, samples, site_join_key, prot_join_key, matched_rows, total_rows = compute_site_vs_protein(
            emat_s, rdesc_s, emat_p, rdesc_p, join_on=join_on, agg=agg, eps=eps, drop_unmatched=drop_unmatched
        )
    except ValueError as e:
        raise click.ClickException(str(e))

    # Write TSVs
    base = outdir / f"site_vs_protein_{len(ratio_df)}x{len(samples)}"
    ratio_tsv = str(base) + "_ratio.tsv"
    log2ratio_tsv = str(base) + "_log2ratio.tsv"
    ratio_out = rdesc_out.join(ratio_df)
    log2ratio_out = rdesc_out.join(log2ratio_df)
    if not Path(ratio_tsv).exists() or eff_force_data:
        ratio_out.to_csv(ratio_tsv, sep="\t", index=True, header=True)
    else:
        click.echo(f"Exists, skipping: {ratio_tsv}")
    if not Path(log2ratio_tsv).exists() or eff_force_data:
        log2ratio_out.to_csv(log2ratio_tsv, sep="\t", index=True, header=True)
    else:
        click.echo(f"Exists, skipping: {log2ratio_tsv}")

    # Text summary
    def _median_ignore_na(a):
        try:
            return float(np.nanmedian(a))
        except Exception:
            return float("nan")

    medians = {s: _median_ignore_na(ratio_df[s].to_numpy()) for s in samples}
    summary_lines = [
        "SITE vs PROTEIN COMPARISON",
        f"site GCT: {site_gct}",
        f"protein GCT: {protein_gct}",
        f"join_on: site='{site_join_key}', protein='{prot_join_key}'; agg: {agg}; eps: {eps}",
        f"samples (intersection): {len(samples)}",
        f"matched site rows: {matched_rows}/{total_rows}",
        "",
        "median(site/protein) by sample:",
    ]
    for s in samples:
        summary_lines.append(f"  - {s}: {medians[s]}")
    summary_text = "\n".join(summary_lines)
    summary_path = outdir / "compare_protein_summary.txt"
    if not summary_path.exists() or eff_force_data:
        summary_path.write_text(summary_text)
    else:
        click.echo(f"Exists, skipping: {summary_path}")

    click.echo(f"Wrote: {ratio_tsv}")
    click.echo(f"Wrote: {log2ratio_tsv}")
    click.echo(f"Wrote summary: {summary_path}")

    if latex:
        pdf_path = maybe_compile_latex(summary_text, str(outdir), basename="compare_protein_summary", title="site-annot site vs protein")
        if pdf_path:
            click.echo(f"Wrote PDF: {pdf_path}")
        else:
            click.echo("LaTeX engine not available or failed; skipped PDF.")

    # Optional GCT/Excel outputs of the chosen metric
    if write_gct_flag or write_excel_flag:
        emat_out = ratio_df if metric == "ratio" else log2ratio_df
        # Ensure sample metadata aligns
        cdesc_out = cdesc_s.loc[samples]
        gct_base = outdir / f"site_vs_protein_{metric}"
        rows, cols = emat_out.shape
        if write_gct_flag:
            gct_out = Path(f"{gct_base}_{rows}x{cols}.gct")
            if not gct_out.exists() or eff_force_data:
                write_gct_file(emat_out, cdesc=cdesc_out, rdesc=rdesc_out, filename=str(gct_base))
            else:
                click.echo(f"Exists, skipping: {gct_out}")
        if write_excel_flag:
            xlsx_out = Path(f"{gct_base}_{rows}x{cols}.xlsx")
            if not xlsx_out.exists() or eff_force_data:
                combined = rdesc_out.join(emat_out)
                write_excel_file(combined, filename=str(gct_base), shape=emat_out.shape, column_metadata=cdesc_out)
            else:
                click.echo(f"Exists, skipping: {xlsx_out}")

    if gene_list:
        emat_for_plot = ratio_df if metric == "ratio" else log2ratio_df
        try:
            generated, missing, dbg = generate_gene_heatmaps(
                emat_for_plot,
                rdesc_s,
                samples,
                outdir,
                gene_list,
                metric,
                force_data=eff_force_data,
                force_plots=eff_force_plots,
                gene_col=gene_col,
                debug=gene_debug,
                generate_zscore=zscore,
            )
            click.echo(f"Generated heatmaps for {generated}/{len(gene_list)} genes.")
            if missing:
                click.echo("Missing genes: " + ", ".join(missing))
            if gene_debug and dbg:
                for line in dbg:
                    click.echo(line)
        except FileNotFoundError as exc:
            click.echo(f"Skipping gene heatmaps: {exc}")
        except subprocess.CalledProcessError as exc:
            click.echo(f"Gene heatmap generation failed: {exc}")


@main.command(name="limma-report")
@click.option("--limma-tsv", type=click.Path(exists=True, dir_okay=False), multiple=True, required=True, help="Limma toptable TSV(s) with columns like id, logFC, adj.P.Val, P.Value [and optional symbol, contrast]")
@click.option("-n", "--top-n", type=int, default=25, show_default=True)
@click.option("-o", "--output-dir", type=click.Path(exists=False, file_okay=False, dir_okay=True), default=None, show_default=True)
@click.option("--latex/--no-latex", is_flag=True, default=True, show_default=True)
@click.option("--ollama/--no-ollama", is_flag=True, default=False, show_default=True)
@click.option("--ollama-model", type=str, default="llama3.2:3b", show_default=True)
def limma_report(limma_tsv, top_n, output_dir, latex, ollama, ollama_model):
    """Render a small LaTeX/PDF report of top limma hits and optionally ask a local LLM (Ollama) to summarize."""
    outdir = pathlib.Path(output_dir or pathlib.Path(limma_tsv[0]).parent).absolute()
    outdir.mkdir(parents=True, exist_ok=True)

    # Collect sections
    import pandas as pd
    sections = []
    text_summary_lines = []
    for tsv in limma_tsv:
        df = pd.read_csv(tsv, sep="\t")
        label = pathlib.Path(tsv).stem
        if "contrast" in df.columns and df["contrast"].nunique() > 1:
            for cname, sub in df.groupby("contrast"):
                sub2 = sub.sort_values(["adj.P.Val", "P.Value"], ascending=[True, True]).head(top_n)
                headers = [h for h in ["id", "symbol", "logFC", "adj.P.Val", "P.Value"] if h in sub2.columns]
                rows = [[str(x) for x in row] for row in sub2[headers].itertuples(index=False, name=None)]
                sections.append({"title": f"{label} ({cname})", "headers": headers, "rows": rows})
                text_summary_lines.append(f"{label} ({cname}) top {len(rows)} hits: best adj.P.Val = {sub2['adj.P.Val'].min()}")
        else:
            df2 = df.sort_values(["adj.P.Val", "P.Value"], ascending=[True, True]).head(top_n)
            headers = [h for h in ["id", "symbol", "logFC", "adj.P.Val", "P.Value"] if h in df2.columns]
            rows = [[str(x) for x in row] for row in df2[headers].itertuples(index=False, name=None)]
            sections.append({"title": label, "headers": headers, "rows": rows})
            text_summary_lines.append(f"{label} top {len(rows)} hits: best adj.P.Val = {df2['adj.P.Val'].min() if 'adj.P.Val' in df2 else 'NA'}")

    # Save text summary
    summary_text = "\n".join(text_summary_lines)
    (outdir / "limma_summary.txt").write_text(summary_text)

    # Optional LaTeX PDF via Jinja2
    if latex:
        context = {"title": "LIMMA Top Hits", "sections": sections, "columns_spec": "l" + "c" * (len(sections[0]['headers']) - 1) if sections and sections[0]['headers'] else "lc"}
        pdf_path = compile_latex_from_template(context, str(outdir), basename="limma_top_hits")
        if pdf_path:
            click.echo(f"Wrote PDF: {pdf_path}")
        else:
            click.echo("LaTeX engine not available or Jinja2 missing; skipped PDF.")

    # Optional Ollama summary
    if ollama:
        ollama_out = maybe_ollama_summarize(summary_text, str(outdir), basename="limma_ollama", model=ollama_model)
        if ollama_out:
            click.echo(f"Wrote Ollama summary: {ollama_out}")
        else:
            click.echo("Ollama not available or failed; skipped model summary.")


@main.command(name="cache", help="Inspect and manage the local UniProt mapping cache")
@click.argument("action", type=click.Choice(["ls", "rm", "clear", "import"]))
@click.argument("keys", nargs=-1)
@click.option("--limit", type=int, default=20, show_default=True, help="Maximum keys to display (for ls)")
@click.option("--filter", "filter_pat", type=str, default=None, help="Substring to filter keys (for ls)")
@click.option("--yes", is_flag=True, help="Confirm deletion for clear")
@click.option("--file", "import_file", type=click.Path(exists=True, dir_okay=False), help="TSV file with two columns: query and uniprot.Swiss-Prot (for import)")
def cache_cmd(action, keys, limit, filter_pat, yes, import_file):
    if action == "ls":
        with mapper.get_db() as db:
            total = len(db)
            shown = 0
            click.echo(f"Cache file: {mapper.sqlitedict_filename}")
            click.echo(f"Total keys: {total}")
            for k in db.keys():
                ks = str(k)
                if filter_pat and filter_pat not in ks:
                    continue
                click.echo(ks)
                shown += 1
                if shown >= limit:
                    break
            if shown == 0:
                click.echo("No keys match filter" if filter_pat else "Cache is empty or unreadable")
        return

    if action == "rm":
        if not keys:
            raise click.UsageError("Provide at least one ENSP key to remove")
        removed = 0
        with mapper.get_db() as db:
            for key in keys:
                if key in db:
                    del db[key]
                    removed += 1
            db.commit()
        click.echo(f"Removed {removed} keys from cache")
        return

    if action == "clear":
        path = mapper.sqlitedict_filename
        if not yes:
            click.echo(f"This will delete {path}. Re-run with --yes to confirm.")
            return
        try:
            if os.path.exists(path):
                os.remove(path)
                click.echo(f"Deleted {path}")
            else:
                click.echo(f"Cache file does not exist: {path}")
        except Exception as exc:
            click.echo(f"Failed to delete {path}: {exc}")
        return

    if action == "import":
        if not import_file:
            raise click.UsageError("Provide --file path to TSV with columns: query, uniprot.Swiss-Prot")
        import pandas as pd
        df = pd.read_csv(import_file, sep="\t")
        required = {"query", "uniprot.Swiss-Prot"}
        if not required.issubset(df.columns):
            raise click.UsageError(f"Missing required columns: {required - set(df.columns)}")
        rows = [
            {"query": str(row["query"]), "uniprot": {"Swiss-Prot": str(row["uniprot.Swiss-Prot"])}}
            for _, row in df.iterrows()
        ]
        with mapper.get_db() as db:
            mapper.update_db(db, rows)
        click.echo(f"Imported {len(rows)} rows into cache")
