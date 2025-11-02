import pandas as pd
import numpy as np
from typing import Tuple, List, Optional
from pathlib import Path
import re
import subprocess


def auto_select_join_keys(rdesc_site: pd.DataFrame, rdesc_prot: pd.DataFrame, requested: str = "auto") -> Tuple[str, str]:
    """
    Return (site_key, prot_key) to use for joining site rows to protein rows.
    - If requested != 'auto', validate presence and return (requested, requested)
    - Otherwise choose a case-insensitive matching pair that maximizes matched site rows
    """
    if requested != "auto":
        if requested not in rdesc_site.columns:
            raise ValueError(f"Join key '{requested}' not found in site rdesc")
        if requested not in rdesc_prot.columns:
            raise ValueError(f"Join key '{requested}' not found in protein rdesc")
        return requested, requested

    site_lut = {c.lower(): c for c in rdesc_site.columns}
    prot_lut = {c.lower(): c for c in rdesc_prot.columns}
    common_lc = sorted(set(site_lut) & set(prot_lut))

    scores = []
    if common_lc:
        for lc in common_lc:
            s_key = site_lut[lc]
            p_key = prot_lut[lc]
            p_set = set(rdesc_prot[p_key].dropna().astype(str))
            s_vals = rdesc_site[s_key].dropna().astype(str)
            matched = s_vals.isin(p_set).sum()
            scores.append((matched, s_key, p_key))
    else:
        exact = sorted(set(rdesc_site.columns) & set(rdesc_prot.columns))
        if not exact:
            raise ValueError("No overlapping join keys between site and protein rdesc")
        for c in exact:
            p_set = set(rdesc_prot[c].dropna().astype(str))
            s_vals = rdesc_site[c].dropna().astype(str)
            matched = s_vals.isin(p_set).sum()
            scores.append((matched, c, c))

    scores.sort(key=lambda x: x[0], reverse=True)
    best = scores[0]
    return best[1], best[2]


def aggregate_protein_matrix(emat_p: pd.DataFrame, rdesc_p: pd.DataFrame, prot_key: str, agg: str) -> pd.DataFrame:
    """Aggregate protein matrix by a row descriptor column using mean/median/sum."""
    rdesc_key = rdesc_p[[prot_key]].copy()
    rdesc_key.index = emat_p.index
    rdesc_key = rdesc_key.dropna()
    emat_p_keyed = emat_p.loc[rdesc_key.index]

    if agg == "mean":
        prot_mat = emat_p_keyed.groupby(rdesc_key[prot_key]).mean()
    elif agg == "median":
        prot_mat = emat_p_keyed.groupby(rdesc_key[prot_key]).median()
    else:
        prot_mat = emat_p_keyed.groupby(rdesc_key[prot_key]).sum()
    return prot_mat


def compute_site_vs_protein(
    emat_s: pd.DataFrame,
    rdesc_s: pd.DataFrame,
    emat_p: pd.DataFrame,
    rdesc_p: pd.DataFrame,
    join_on: str = "auto",
    agg: str = "mean",
    eps: float = 1e-6,
    drop_unmatched: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str], str, str, int, int]:
    """
    Compute site/protein ratio and log2ratio matrices.

    Returns: (ratio_df, log2ratio_df, rdesc_out, samples, site_key, prot_key, matched_rows, total_rows)
    """
    # Convert to numeric to avoid object arrays and pd.NA issues
    emat_s = emat_s.apply(pd.to_numeric, errors="coerce")
    emat_p = emat_p.apply(pd.to_numeric, errors="coerce")

    # Align samples
    samples = [s for s in emat_s.columns if s in emat_p.columns]
    if not samples:
        raise ValueError("No overlapping samples between site and protein matrices")
    emat_s = emat_s[samples]
    emat_p = emat_p[samples]

    # Determine join keys and aggregate protein matrix
    site_key, prot_key = auto_select_join_keys(rdesc_s, rdesc_p, join_on)
    prot_mat = aggregate_protein_matrix(emat_p, rdesc_p, prot_key, agg)

    # Build mapping and compute ratios
    site_keys = rdesc_s[[site_key]].copy()
    site_keys.index = emat_s.index
    matched_mask = site_keys[site_key].isin(prot_mat.index)
    matched_rows = int(matched_mask.sum())
    total_rows = int(len(site_keys))

    # Convert to float arrays for numeric ops
    prot_mat_aligned = prot_mat.reindex(site_keys[site_key]).to_numpy(dtype=float, copy=False)
    site_mat = emat_s.to_numpy(dtype=float, copy=False)

    denom = prot_mat_aligned + float(eps)
    ratio_mat = site_mat / denom
    log2ratio_mat = np.log2(site_mat + float(eps)) - np.log2(denom)

    ratio_df = pd.DataFrame(ratio_mat, index=emat_s.index, columns=samples)
    log2ratio_df = pd.DataFrame(log2ratio_mat, index=emat_s.index, columns=samples)

    if drop_unmatched:
        keep_idx = site_keys.index[matched_mask]
        ratio_df = ratio_df.loc[keep_idx]
        log2ratio_df = log2ratio_df.loc[keep_idx]
        rdesc_out = rdesc_s.loc[keep_idx]
    else:
        rdesc_out = rdesc_s

    return ratio_df, log2ratio_df, rdesc_out, samples, site_key, prot_key, matched_rows, total_rows


GENE_COLUMN_CANDIDATES = [
    "geneid",
    "gene_id",
    "gene",
    "symbol",
    "gene name",
    "gene_name",
    "ensg",
    "enst",
    "ensp",
    "uniprot",
    "uniprot_id",
    "protein",
    "proteinid",
    "protein_id",
    "entry_name",
]


def _sanitize_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return cleaned or "unnamed"


def _candidate_gene_columns(columns: List[str]) -> List[str]:
    lc = [c.lower() for c in columns]
    base = set(GENE_COLUMN_CANDIDATES)
    # add any columns containing 'symbol' or 'gene'
    for c, orig in zip(lc, columns):
        if (
            "symbol" in c
            or c == "gene"
            or c.startswith("gene_")
            or c.endswith("_gene")
            or "geneid" in c
            or c == "gene id"
            or c in ("ensg", "enst", "ensp")
            or "uniprot" in c
            or c.startswith("protein")
        ):
            base.add(c)
    # Map back to original names preserving first occurrence
    out = []
    seen = set()
    for name_lc, orig in zip(lc, columns):
        if name_lc in base and name_lc not in seen:
            out.append(orig)
            seen.add(name_lc)
    return out


def _tokenize(val: str) -> List[str]:
    return [t for t in re.split(r"[\s,;|]+", val) if t]


def _match_gene_rows(rdesc: pd.DataFrame, gene: str, prefer_col: Optional[str] = None) -> pd.Index:
    if gene is None:
        return pd.Index([])
    gene = str(gene).strip()
    if not gene:
        return pd.Index([])

    gene_lower = gene.lower()
    mask = pd.Series(False, index=rdesc.index)

    candidate_cols = _candidate_gene_columns(list(rdesc.columns))
    if prefer_col:
        # find matching column in a case-insensitive way
        lut = {c.lower(): c for c in rdesc.columns}
        if prefer_col.lower() in lut:
            candidate_cols = [lut[prefer_col.lower()]] + [c for c in candidate_cols if c != lut[prefer_col.lower()]]
    for col in candidate_cols:
        col_lower = col.lower()
        series = rdesc[col]
        series = series.dropna().astype(str)
        if series.empty:
            continue
        # exact match (case-sensitive and insensitive)
        exact_idx = series.index[series == gene]
        if not exact_idx.empty:
            mask.loc[exact_idx] = True
        lower_idx = series.index[series.str.lower() == gene_lower]
        if not lower_idx.empty:
            mask.loc[lower_idx] = True
        # tokenized match (split multi-valued cells)
        tok_idx = series.index[
            series.apply(lambda s: gene_lower in [t.lower() for t in _tokenize(s)])
        ]
        if not tok_idx.empty:
            mask.loc[tok_idx] = True

    return mask[mask].index


def _call_heatmap_script(
    matrix_path: Path,
    output_pdf: Path,
    title: str,
    metric: str,
    cluster_rows: bool,
    cluster_cols: bool,
    zscore: bool,
    scale_min: Optional[float] = None,
    scale_max: Optional[float] = None,
) -> None:
    script_path = Path(__file__).parent.parent / "R" / "plot_gene_heatmap.R"
    if not script_path.exists():
        raise FileNotFoundError(f"Heatmap script not found: {script_path}")

    cmd = [
        "Rscript",
        str(script_path),
        "--matrix",
        str(matrix_path),
        "--output",
        str(output_pdf),
        "--title",
        title,
        "--metric",
        metric,
        "--cluster-rows",
        str(bool(cluster_rows)).upper(),
        "--cluster-cols",
        str(bool(cluster_cols)).upper(),
        "--zscore",
        str(bool(zscore)).upper(),
    ]
    if scale_min is not None and scale_max is not None:
        cmd.extend(["--scale-min", str(float(scale_min))])
        cmd.extend(["--scale-max", str(float(scale_max))])
    res = subprocess.run(cmd, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if res.returncode != 0:
        out = res.stdout.decode("utf-8", errors="replace")
        raise RuntimeError(f"R heatmap failed (exit {res.returncode}):\n{out}")


def generate_gene_heatmaps(
    emat: pd.DataFrame,
    rdesc: pd.DataFrame,
    samples: List[str],
    output_dir: Path,
    genes: List[str],
    metric: str,
    force_data: bool = False,
    force_plots: bool = False,
    gene_col: Optional[str] = None,
    debug: bool = False,
    generate_zscore: bool = True,
) -> Tuple[int, List[str], List[str]]:
    """
    Generate per-gene heatmaps using ComplexHeatmap via an R script.

    Returns (count_generated, missing_genes)
    """
    if not genes:
        return 0, [], []

    gene_root = output_dir / "genes"
    gene_root.mkdir(parents=True, exist_ok=True)

    generated = 0
    missing = []
    debug_lines: List[str] = []

    # Pre-compute global color scales
    emat_values = emat.to_numpy(dtype=float)
    # raw_min = np.nanmin(emat_values) if emat_values.size else None
    # raw_max = np.nanmax(emat_values) if emat_values.size else None
    raw_min = np.nanpercentile(emat_values, 0.5) if emat_values.size else None
    raw_max = np.nanpercentile(emat_values, 95) if emat_values.size else None


    z_min = z_max = None
    if generate_zscore and emat_values.size:
        # Row-wise zscore across entire dataset
        row_means = np.nanmean(emat_values, axis=1)
        row_sds = np.nanstd(emat_values, axis=1, ddof=0)
        row_sds = np.where(np.isnan(row_sds) | (row_sds == 0), 1.0, row_sds)
        z = (emat_values - row_means[:, None]) / row_sds[:, None]
        z_min = np.nanmin(z)
        z_max = np.nanmax(z)

    for gene in genes:
        if gene is None:
            continue
        gene_str = str(gene).strip()
        if not gene_str:
            continue

        idx = _match_gene_rows(rdesc, gene_str, prefer_col=gene_col)
        # ensure indices exist in emat
        idx = idx.intersection(emat.index)
        if idx.empty:
            missing.append(gene_str)
            if debug:
                debug_lines.append(
                    f"{gene_str}: matched 0 rows; columns searched: {', '.join(_candidate_gene_columns(list(rdesc.columns)))}"
                )
            continue

        sub_mat = emat.loc[idx]
        sub_mat = sub_mat.dropna(how="all")
        if sub_mat.empty:
            missing.append(f"{gene_str} (all NaN)")
            continue

        gene_dir = gene_root / _sanitize_name(gene_str)
        gene_dir.mkdir(parents=True, exist_ok=True)

        matrix_path = gene_dir / f"{metric}_matrix.tsv"

        if debug:
            debug_lines.append(
                f"{gene_str}: matched {len(idx)} rows; columns searched: {', '.join(_candidate_gene_columns(list(rdesc.columns)))}"
            )

        # write inputs only if forcing or missing
        if force_data or not matrix_path.exists():
            sub_mat.to_csv(matrix_path, sep="\t")
        rdesc_path = gene_dir / "rdesc.tsv"
        if force_data or not rdesc_path.exists():
            rdesc.loc[idx].to_csv(rdesc_path, sep="\t")

        # Render combinations of clustering for raw and optional z-score
        variants = [(metric, False)]
        if generate_zscore:
            variants.append((f"{metric}_z", True))

        for label, zflag in variants:
            for cr in (True, False):
                for cc in (True, False):
                    suffix = f"rc{'T' if cr else 'F'}_cc{'T' if cc else 'F'}"
                    pdf_path = gene_dir / f"{label}_heatmap_{suffix}.pdf"
                    if pdf_path.exists() and not force_plots:
                        generated += 1
                        continue
                    try:
                        smin = raw_min if not zflag else z_min
                        smax = raw_max if not zflag else z_max
                        _call_heatmap_script(
                            matrix_path,
                            pdf_path,
                            title=gene_str,
                            metric=metric,
                            cluster_rows=cr,
                            cluster_cols=cc,
                            zscore=zflag,
                            scale_min=smin,
                            scale_max=smax,
                        )
                        generated += 1
                    except Exception as exc:
                        missing.append(f"{gene_str} (heatmap failed {label}_{suffix}: {exc}")

    return generated, missing, debug_lines
