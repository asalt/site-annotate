import logging
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from pyteomics import mzml  # noqa: E402
from pyteomics import mass as mass_mod  # noqa: E402


logger = logging.getLogger(__name__)


ANNOTATION_SETTINGS = {
    "fragment_tol_mass": 0.05,
    "fragment_tol_mode": "Da",
    "ion_types": "by",
    "max_ion_charge": 2,
    "neutral_losses": {
        "NH3": -17.026549,
        "H2O": -18.010565,
        "H3PO4": -97.976896,
    },
}


def parse_spectrum_id(spectrum: str) -> Tuple[str, int, int]:
    """
    Parse FragPipe Spectrum identifier of the form:
        basename.scan.scan.charge

    Returns (basename, scan_number, charge).
    """
    parts = spectrum.split(".")
    if len(parts) < 4:
        raise ValueError(f"Unexpected spectrum id format: {spectrum}")

    basename = ".".join(parts[:-3]) or parts[0]
    scan_number = int(parts[-3])
    charge = int(parts[-1])
    return basename, scan_number, charge


def find_mzml_for_basename(raw_dir: Path, basename: str) -> Path:
    """
    Prefer calibrated mzML, fall back to uncalibrated.
    """
    candidates = [
        raw_dir / f"{basename}_calibrated.mzML",
        raw_dir / f"{basename}.mzML",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"No mzML found for {basename} in {raw_dir} "
        f"(tried: {', '.join(str(c) for c in candidates)})"
    )


def build_ms2_index(mzml_path: Path, target_scans: Iterable[int]) -> Dict[int, dict]:
    """
    Build a mapping from scan number -> mzML spectrum dict for MS2 scans.

    Only scans in `target_scans` are stored to keep memory usage small.
    """
    target_scans = set(target_scans)
    index: Dict[int, dict] = {}

    logger.info(f"Indexing MS2 spectra from {mzml_path}")
    with mzml.MzML(str(mzml_path)) as reader:
        for spectrum in reader:
            if spectrum.get("ms level") != 2:
                continue

            scan_id = spectrum.get("id")
            scan_number = None
            if isinstance(scan_id, str) and "scan=" in scan_id:
                try:
                    scan_number = int(scan_id.split("scan=")[1])
                except Exception:
                    scan_number = None

            if scan_number is None:
                scan_number = spectrum["index"] + 1

            if scan_number in target_scans and scan_number not in index:
                index[scan_number] = spectrum
                if len(index) == len(target_scans):
                    break

    missing = target_scans - set(index.keys())
    if missing:
        logger.warning(
            f"{len(missing)} target scans not found in {mzml_path}: "
            f"{sorted(list(missing))[:10]}{'...' if len(missing) > 10 else ''}"
        )
    return index


def select_top_phospho_psms(
    psm_path: Path,
    n_top: int = 4,
    qvalue_cutoff: float = 0.01,
    min_localization: float = 0.0,
    phospho_loc_col: str = "STY:79.9663 Best Localization",
    max_length: Optional[int] = None,
) -> pd.DataFrame:
    """
    Read FragPipe psm.tsv and select top phospho PSMs by Hyperscore.
    """
    logger.info(f"Reading PSMs from {psm_path}")
    psm_df = pd.read_table(psm_path)

    if "Is Decoy" not in psm_df.columns or "Qvalue" not in psm_df.columns:
        raise ValueError("Expected columns 'Is Decoy' and 'Qvalue' not found.")

    if phospho_loc_col not in psm_df.columns:
        raise ValueError(f"Expected column {phospho_loc_col!r} not found.")

    mask = (psm_df["Is Decoy"] == False) & (psm_df["Qvalue"] <= qvalue_cutoff)  # noqa: E712
    phospho_scores = psm_df[phospho_loc_col].fillna(0)
    mask &= phospho_scores > min_localization

    if max_length is not None and "Peptide Length" in psm_df.columns:
        mask &= psm_df["Peptide Length"] <= max_length

    filtered = psm_df[mask]
    if filtered.empty:
        raise ValueError("No PSMs passed the filtering criteria.")

    if "Hyperscore" not in filtered.columns:
        raise ValueError("Expected column 'Hyperscore' not found.")

    top_psms = filtered.sort_values("Hyperscore", ascending=False).head(n_top).copy()
    logger.info(
        f"Selected {len(top_psms)} PSMs (n_top={n_top}, "
        f"qvalue<={qvalue_cutoff}, localization>{min_localization})"
    )
    return top_psms


def _format_title(
    modified_peptide: str,
    charge: int,
    scan_number: int,
    hyperscore: float,
    qvalue: float,
) -> str:
    """
    Build a compact multi-line title including modified sequence.
    """
    line1 = f"{modified_peptide}  z={charge}  scan={scan_number}"
    line2 = f"Hyperscore={hyperscore:.1f}  Q={qvalue:.3g}"
    return f"{line1}\n{line2}"


def _format_ms_title(
    charge: int,
    scan_number: int,
    hyperscore: float,
    qvalue: float,
) -> str:
    """
    Simpler title for MS2 panel, without the (long) modified sequence.
    """
    line1 = f"scan={scan_number}  z={charge}"
    line2 = f"Hyperscore={hyperscore:.1f}  Q={qvalue:.3g}"
    return f"{line1}\n{line2}"


def _import_spectrum_utils():
    """
    Import spectrum_utils with a numba cache workaround.

    The environment in which this project runs can trigger a numba
    caching error when spectrum_utils is imported. To avoid that,
    we monkeypatch numba's FunctionCache before importing
    spectrum_utils.
    """
    try:
        from numba.core import dispatcher
    except Exception as exc:  # pragma: no cover - defensive
        raise ImportError("numba is required for spectrum_utils plotting") from exc

    class DummyCache:
        def __init__(self, *args, **kwargs):
            self._memo = {}

        def load_overload(self, *args, **kwargs):
            return None

        def save_overload(self, *args, **kwargs):
            return

    dispatcher.FunctionCache = DummyCache

    import spectrum_utils.spectrum as sus  # type: ignore
    import spectrum_utils.plot as sup  # type: ignore
    from spectrum_utils import fragment_annotation as fa  # type: ignore

    return sus, sup, fa


def _compute_sequence_coverage(sequence: str, annotations) -> Dict[str, List[bool]]:
    """
    Compute per-residue coverage for b and y ions based on annotations.

    annotations: iterable where each element is either an Ion annotation
    object or a string representation (e.g. 'b5', 'y7-H2O').
    """
    L = len(sequence)
    cov_b = [False] * L
    cov_y = [False] * L

    ion_re = re.compile(r"([abyxyz])(\d+)")

    for ann in annotations:
        if ann is None:
            continue
        s = str(ann)
        m = ion_re.search(s)
        if not m:
            continue
        ion_type = m.group(1)
        n = int(m.group(2))
        if ion_type == "b":
            # prefix of length n -> residues 0..n-1
            for i in range(min(n, L)):
                cov_b[i] = True
        elif ion_type == "y":
            # suffix of length n -> residues L-n..L-1
            start = max(0, L - n)
            for i in range(start, L):
                cov_y[i] = True

    return {"b": cov_b, "y": cov_y}


def _draw_sequence_ribbon(
    ax,
    sequence: str,
    coverage: Dict[str, List[bool]],
    mod_masses: Optional[Dict[int, List[float]]] = None,
):
    """
    Draw a sequence coverage ribbon on the given axis.
    """
    L = len(sequence)
    ax.set_xlim(0, L)
    ax.set_ylim(0, 1)
    ax.axis("off")

    box_height = 0.6
    y_center = 0.5

    mod_masses = mod_masses or {}
    mod_positions_set = set(mod_masses.keys())

    for i, aa in enumerate(sequence):
        x = i + 0.5
        has_b = coverage["b"][i]
        has_y = coverage["y"][i]

        if has_b and has_y:
            facecolor = "#ccccff"  # both
        elif has_b:
            facecolor = "#a6cee3"  # b-only
        elif has_y:
            facecolor = "#fb9a99"  # y-only
        else:
            facecolor = "#f0f0f0"  # no coverage

        pos1 = i + 1
        is_mod = pos1 in mod_positions_set

        rect = plt.Rectangle(
            (i, y_center - box_height / 2),
            1.0,
            box_height,
            facecolor=facecolor,
            edgecolor="black" if not is_mod else "red",
            linewidth=0.5 if not is_mod else 1.2,
        )
        ax.add_patch(rect)

        ax.text(
            x,
            y_center,
            aa,
            ha="center",
            va="center",
            fontsize=10,
        )

    # 1-based residue indices below ribbon (forward, for b-ions)
    for i in range(L):
        x = i + 0.5
        ax.text(
            x,
            y_center - box_height / 2 - 0.1,
            str(i + 1),
            ha="center",
            va="top",
            fontsize=7,
        )

    # reverse indices above ribbon (L..1, for y-ions)
    for i in range(L):
        x = i + 0.5
        ax.text(
            x,
            y_center + box_height / 2 + 0.1,
            str(L - i),
            ha="center",
            va="bottom",
            fontsize=7,
        )

    # modification mass deltas (if any) just below the residue letter
    for pos1, masses in mod_masses.items():
        if not masses:
            continue
        i = pos1 - 1
        if not (0 <= i < L):
            continue
        x = i + 0.5
        mass = masses[0]
        ax.text(
            x,
            y_center - box_height / 4,
            f"[+{mass:.2f}]",
            ha="center",
            va="top",
            fontsize=7,
        )

    # sequence string above ribbon (unmodified peptide)
    ax.text(0, 1.02, sequence, fontsize=10, va="bottom", ha="left")


def _label_annotated_peaks(ax, msms, max_labels: int = 60, label_intensity_frac: float = 0.05):
    """
    Add labels for annotated fragments, including ion type, charge, and m/z.
    """
    if msms.annotation is None or not len(msms.annotation):
        return

    mz = msms.mz
    inten = msms.intensity
    ann = msms.annotation

    max_int = float(np.max(inten)) if len(inten) else 0.0
    if max_int <= 0:
        return
    thr = label_intensity_frac * max_int

    # all annotated indices
    all_idx = [i for i, a in enumerate(ann) if a is not None and str(a) != "?"]
    if not all_idx:
        return

    # determine charge state for each annotated fragment
    charge_re = re.compile(r"\^(\d+)\+")

    def get_z(i: int) -> int:
        s = str(ann[i])
        m = charge_re.search(s)
        return int(m.group(1)) if m else 1

    # high-intensity annotated peaks
    base_idx = [i for i in all_idx if inten[i] >= thr]
    base_one = [i for i in base_idx if get_z(i) == 1]
    base_multi = [i for i in base_idx if get_z(i) > 1]

    indices = set(base_one) | set(base_multi)
    if not indices:
        # fall back to all annotated peaks, 1+ preferred
        one_plus = [i for i in all_idx if get_z(i) == 1]
        if one_plus:
            indices = set(one_plus)
        else:
            indices = set(all_idx)

    # enforce label budget, favor 1+ first then multi-charged
    if len(indices) > max_labels:
        labels: List[int] = []

        if base_one:
            base_one_arr = np.array(base_one)
            order = np.argsort(inten[base_one_arr])[::-1]
            base_one_arr = base_one_arr[order]
            take_one = min(max_labels, len(base_one_arr))
            labels.extend(base_one_arr[:take_one].tolist())
        remaining = max_labels - len(labels)

        if remaining > 0 and base_multi:
            base_multi_arr = np.array(base_multi)
            order = np.argsort(inten[base_multi_arr])[::-1]
            base_multi_arr = base_multi_arr[order]
            labels.extend(base_multi_arr[:remaining].tolist())

        indices = sorted(set(labels))
    else:
        indices = sorted(indices)

    if not indices:
        return

    y_min, y_max = ax.get_ylim()
    for i in indices:
        ion_str = str(ann[i])
        label = f"{ion_str} {mz[i]:.1f}"
        frac = inten[i] / max_int if max_int > 0 else 0.0
        if frac >= 0.5:
            fs = 6
        elif frac >= 0.2:
            fs = 5
        else:
            fs = 4

        # place label close to peak but keep inside axis
        y_pos = min(inten[i] * 1.05, y_max * 0.9)

        # thin connector line from peak top to label
        ax.plot(
            [mz[i], mz[i]],
            [inten[i], y_pos],
            color="0.5",
            linewidth=0.3,
        )
        ax.text(
            mz[i],
            y_pos,
            label,
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=fs,
        )


def _build_proforma_from_fragpipe(peptide: str, assigned_mods: str) -> Optional[str]:
    """
    Build a simple ProForma string from FragPipe 'Assigned Modifications'.

    Example Assigned Modifications:
        N-term(304.2071),9S(79.9663)
    """
    if not isinstance(assigned_mods, str) or not assigned_mods.strip():
        return None

    mods_by_pos: Dict[int, List[float]] = {}
    parts = [p.strip() for p in assigned_mods.split(",") if p.strip()]

    for part in parts:
        m_term = re.match(r"^(N-term|C-term)\(([-0-9.]+)\)$", part)
        if m_term:
            term = m_term.group(1)
            mass = float(m_term.group(2))
            if term == "N-term":
                pos = 0
            else:
                pos = len(peptide) + 1
            mods_by_pos.setdefault(pos, []).append(mass)
            continue

        m_pos = re.match(r"^(\d+)([A-Z])\(([-0-9.]+)\)$", part)
        if m_pos:
            pos = int(m_pos.group(1))
            aa = m_pos.group(2)
            mass = float(m_pos.group(3))
            if 1 <= pos <= len(peptide):
                if peptide[pos - 1] != aa:
                    logger.warning(
                        f"AA mismatch for modification {part}: "
                        f"peptide has {peptide[pos - 1]} at position {pos}"
                    )
                mods_by_pos.setdefault(pos, []).append(mass)
            else:
                logger.warning(f"Position {pos} out of range for peptide of length {len(peptide)}")
            continue

        logger.warning(f"Could not parse Assigned Modifications token: {part}")

    proforma: List[str] = []

    # N-term modifications at position 0: use ProForma n[...] syntax
    for mass in mods_by_pos.get(0, []):
        proforma.append(f"n[+{mass:.4f}]")

    for i, aa in enumerate(peptide, start=1):
        proforma.append(aa)
        for mass in mods_by_pos.get(i, []):
            proforma.append(f"[+{mass:.4f}]")

    # C-term modifications at position len(peptide) + 1: use c[...] syntax
    for mass in mods_by_pos.get(len(peptide) + 1, []):
        proforma.append(f"c[+{mass:.4f}]")

    return "".join(proforma)


def _get_mod_masses_from_assigned(peptide: str, assigned_mods: str) -> Dict[int, List[float]]:
    """
    Extract 1-based residue positions and mass deltas for positional
    modifications from a FragPipe Assigned Modifications string.

    N-term and C-term modifications are ignored for the ribbon.
    """
    masses: Dict[int, List[float]] = {}
    if not isinstance(assigned_mods, str) or not assigned_mods.strip():
        return masses

    tokens = [t.strip() for t in assigned_mods.split(",") if t.strip()]
    for tok in tokens:
        m_pos = re.match(r"^(\d+)([A-Z])\(([-0-9.]+)\)$", tok)
        if not m_pos:
            continue
        pos = int(m_pos.group(1))
        mass = float(m_pos.group(3))
        if 1 <= pos <= len(peptide):
            masses.setdefault(pos, []).append(mass)
    return masses


def _parse_mod_masses_full(peptide: str, assigned_mods: str) -> Dict[int, List[float]]:
    """
    Parse Assigned Modifications into positional mass deltas.

    Positions:
      0           -> N-term
      1..len(seq) -> residues
      len(seq)+1  -> C-term
    """
    mods: Dict[int, List[float]] = {}
    if not isinstance(assigned_mods, str) or not assigned_mods.strip():
        return mods

    L = len(peptide)
    tokens = [t.strip() for t in assigned_mods.split(",") if t.strip()]
    for tok in tokens:
        if tok.startswith("N-term(") and tok.endswith(")"):
            try:
                m = float(tok[len("N-term(") : -1])
                mods.setdefault(0, []).append(m)
            except Exception:
                continue
            continue
        if tok.startswith("C-term(") and tok.endswith(")"):
            try:
                m = float(tok[len("C-term(") : -1])
                mods.setdefault(L + 1, []).append(m)
            except Exception:
                continue
            continue
        # positional AA(mass)
        try:
            head, mass_part = tok.split("(", 1)
            mass = float(mass_part.rstrip(")"))
        except Exception:
            continue
        pos_str = "".join(ch for ch in head if ch.isdigit())
        if not pos_str:
            continue
        pos = int(pos_str)
        if 1 <= pos <= L:
            mods.setdefault(pos, []).append(mass)
    return mods


def compute_mass_ladder(peptide: str, assigned_mods: str) -> pd.DataFrame:
    """
    Compute theoretical b and y ladders (1+ and 2+) for a peptide
    with positional modifications encoded as mass deltas.
    """
    L = len(peptide)
    mods = _parse_mod_masses_full(peptide, assigned_mods)

    rows = []

    def prefix_delta(n: int) -> float:
        # N-term (0) plus positions 1..n
        return sum(
            d for pos, ds in mods.items() if pos <= n and pos >= 0 for d in ds
        )

    def suffix_delta(n: int) -> float:
        # positions L-n+1..L plus C-term
        start = L - n + 1
        return sum(
            d
            for pos, ds in mods.items()
            if (pos == L + 1 or (start <= pos <= L) and pos > 0)
            for d in ds
        )

    # b-ions
    for n in range(1, L):
        seq = peptide[:n]
        for z in (1, 2):
            base = mass_mod.fast_mass(seq, ion_type="b", charge=z)
            mz = base + prefix_delta(n) / z
            rows.append({"ion": "b", "n": n, "z": z, "mz": mz})

    # y-ions
    for n in range(1, L):
        seq = peptide[L - n :]
        for z in (1, 2):
            base = mass_mod.fast_mass(seq, ion_type="y", charge=z)
            mz = base + suffix_delta(n) / z
            rows.append({"ion": "y", "n": n, "z": z, "mz": mz})

    df = pd.DataFrame(rows).sort_values(["ion", "n", "z"])
    return df


def format_mass_ladder_text(ladder: pd.DataFrame) -> str:
    """
    Format mass ladder as aligned 4-column text:
      n_b, b1+, b2+, n_y, y1+, y2+

    where n_y runs in reverse (L-1..1) so large y-ions are easy to spot.
    """
    if ladder is None or ladder.empty:
        return ""
    max_n = int(ladder["n"].max())
    lines = []
    lines.append(" n_b      b1+        b2+     n_y      y1+        y2+")
    for n in range(1, max_n + 1):
        sub_b = ladder[ladder["ion"] == "b"]
        sub_y = ladder[ladder["ion"] == "y"]

        ny = max_n - n + 1  # reverse index for y-ions

        def get_mz(ion, z):
            vals = ladder[(ladder["ion"] == ion) & (ladder["n"] == (n if ion == "b" else ny)) & (ladder["z"] == z)][
                "mz"
            ].values
            return f"{vals[0]:.4f}" if len(vals) else "-"

        s_b1 = get_mz("b", 1)
        s_b2 = get_mz("b", 2)
        s_y1 = get_mz("y", 1)
        s_y2 = get_mz("y", 2)
        lines.append(
            f"{n:4d}  {s_b1:>10}  {s_b2:>10}  {ny:5d}  {s_y1:>10}  {s_y2:>10}"
        )
    return "\n".join(lines)


def _plot_ms2_interactive_html(
    mz_array,
    intensity_array,
    charge: int,
    scan_number: int,
    precursor_mz: Optional[float],
    peptide: Optional[str],
    assigned_mods: Optional[str],
    hyperscore: float,
    qvalue: float,
    out_path: Path,
) -> None:
    """
    Write a standalone interactive HTML spectrum using spectrum_utils.iplot.
    """
    sus, _, _ = _import_spectrum_utils()
    import spectrum_utils.iplot as supi  # type: ignore

    pmz = float(precursor_mz) if precursor_mz is not None else 0.0

    msms = sus.MsmsSpectrum(
        identifier=f"scan={scan_number}",
        precursor_mz=pmz,
        precursor_charge=charge,
        mz=mz_array,
        intensity=intensity_array,
        retention_time=0.0,
    )

    if peptide and assigned_mods:
        proforma = _build_proforma_from_fragpipe(peptide, assigned_mods)
        if proforma:
            try:
                msms.annotate_proforma(proforma, **ANNOTATION_SETTINGS)
            except Exception as exc:  # pragma: no cover
                logger.warning(
                    f"Failed to annotate spectrum for HTML (scan {scan_number}): {exc}"
                )

    chart = supi.spectrum(msms)
    chart = chart.properties(
        title=f"scan={scan_number} z={charge}  Hyperscore={hyperscore:.1f}  Q={qvalue:.3g}"
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    chart.save(str(out_path))


def plot_ms2_plain(
    mz_array,
    intensity_array,
    modified_peptide: str,
    charge: int,
    scan_number: int,
    hyperscore: float,
    qvalue: float,
    out_path: Path,
    mz_min: Optional[float] = None,
    mz_max: Optional[float] = None,
    max_labels: Optional[int] = None,
    label_intensity_frac: float = 0.05,
    info_text: Optional[str] = None,
) -> None:
    """
    Simple stick spectrum plot using matplotlib, with sequence/mods in the title.

    Highlights peaks above a relative intensity threshold (up to max_labels)
    with their m/z values.
    """
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.vlines(mz_array, 0, intensity_array, color="black", linewidth=0.5)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")

    if mz_min is not None or mz_max is not None:
        ax.set_xlim(mz_min, mz_max)

    # dynamic selection of peaks to label
    labeled_idx = set()
    # auto mode: choose a generous cap based on spectrum size
    if max_labels is None:
        max_labels_eff = min(len(mz_array), 120)
    else:
        max_labels_eff = max_labels

    if max_labels_eff and len(mz_array):
        max_int = float(np.max(intensity_array))
        if max_int > 0:
            threshold = label_intensity_frac * max_int
            cand_idx = np.where(intensity_array >= threshold)[0]
            if len(cand_idx) > max_labels_eff:
                # take top max_labels among candidates
                top_idx = cand_idx[np.argsort(intensity_array[cand_idx])[-max_labels_eff:]]
                base_idx = np.sort(top_idx)
            else:
                base_idx = np.sort(cand_idx)

            for i in base_idx:
                labeled_idx.add(int(i))

            # ensure some labels at high m/z even if intensities are low
            remaining = max_labels_eff - len(labeled_idx)
            if remaining > 0:
                try:
                    high_thresh_mz = max(
                        1300.0, float(np.percentile(mz_array, 75))
                    )
                except Exception:
                    high_thresh_mz = 1300.0
                high_idx = np.where(mz_array >= high_thresh_mz)[0]
                if len(high_idx):
                    # pick highest-intensity peaks in high region that are not yet labeled
                    high_idx = [i for i in high_idx if i not in labeled_idx]
                    if high_idx:
                        high_idx_arr = np.array(high_idx)
                        take = min(remaining, len(high_idx_arr))
                        top_high = high_idx_arr[
                            np.argsort(intensity_array[high_idx_arr])[-take:]
                        ]
                        for i in top_high:
                            labeled_idx.add(int(i))

            for i in sorted(labeled_idx):
                mz_val = mz_array[i]
                inten = intensity_array[i]
                ax.text(
                    mz_val,
                    inten,
                    f"{mz_val:.1f}",
                    rotation=90,
                    va="bottom",
                    ha="center",
                    fontsize=6,
                )

    # add optional info text block
    # info text will be placed in figure coordinates after layout

    # cosmetics: remove top/right spines and leave some headroom
    for spine in ("top", "right"):
        if spine in ax.spines:
            ax.spines[spine].set_visible(False)
    if len(intensity_array):
        ax.set_ylim(0, float(np.max(intensity_array)) * 1.2)

    ax.set_title(_format_title(modified_peptide, charge, scan_number, hyperscore, qvalue))
    fig.tight_layout()

    # Place info text just below the x-axis in figure coordinates
    if info_text:
        fig.text(
            0.99,
            0.01,
            info_text,
            ha="right",
            va="bottom",
            fontsize=8,
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    logger.info(f"Wrote {out_path}")


def plot_ms2_annotated(
    mz_array,
    intensity_array,
    modified_peptide: str,
    charge: int,
    scan_number: int,
    hyperscore: float,
    qvalue: float,
    out_path: Path,
    precursor_mz: Optional[float] = None,
    mz_min: Optional[float] = None,
    mz_max: Optional[float] = None,
    peptide: Optional[str] = None,
    assigned_mods: Optional[str] = None,
    info_text: Optional[str] = None,
) -> None:
    """
    Plot MS2 spectrum using spectrum_utils styling, with a sequence coverage ribbon.

    Uses FragPipe peptide + Assigned Modifications to construct a simple
    ProForma string for fragment annotation when possible.
    """
    sus, sup, fa = _import_spectrum_utils()

    pmz = float(precursor_mz) if precursor_mz is not None else 0.0

    msms = sus.MsmsSpectrum(
        identifier=f"scan={scan_number}",
        precursor_mz=pmz,
        precursor_charge=charge,
        mz=mz_array,
        intensity=intensity_array,
        retention_time=0.0,
    )

    if peptide and assigned_mods:
        proforma = _build_proforma_from_fragpipe(peptide, assigned_mods)
        if proforma:
            try:
                msms.annotate_proforma(proforma, **ANNOTATION_SETTINGS)
            except Exception as exc:  # pragma: no cover - defensive
                logger.warning(f"Failed to annotate spectrum for scan {scan_number}: {exc}")

    # figure with sequence ribbon (top) + MS2 spectrum (bottom) + ladder panel (right)
    width = 16
    if peptide:
        L = len(peptide)
        if L >= 40:
            width = 18
        if L >= 60:
            width = 20
    fig = plt.figure(figsize=(width, 5))
    gs = fig.add_gridspec(
        nrows=2, ncols=2, height_ratios=[1, 4], width_ratios=[3, 1], hspace=0.05, wspace=0.3
    )
    ax_seq = fig.add_subplot(gs[0, 0])
    ax_ms = fig.add_subplot(gs[1, 0])
    ax_ladder = fig.add_subplot(gs[:, 1])
    ax_ladder.axis("off")

    # sequence coverage ribbon if peptide available and annotations present
    if peptide and getattr(msms, "annotation", None) is not None:
        coverage = _compute_sequence_coverage(peptide, msms.annotation)
        mod_masses = _get_mod_masses_from_assigned(peptide, assigned_mods or "")
        _draw_sequence_ribbon(ax_seq, peptide, coverage, mod_masses=mod_masses)
    else:
        ax_seq.axis("off")

    # spectrum plot
    sup.spectrum(msms, ax=ax_ms, grid=False)

    # ensure x-axis covers annotated fragments; extend user-specified mz_max if needed
    if mz_min is not None or mz_max is not None:
        x_min, x_max = ax_ms.get_xlim()
        if mz_min is not None:
            x_min = mz_min
        if mz_max is not None:
            x_max = mz_max
        if getattr(msms, "annotation", None) is not None:
            ann_mz = [
                msms.mz[i]
                for i, a in enumerate(msms.annotation)
                if a is not None and str(a) != "?"
            ]
            if ann_mz:
                max_ann = max(ann_mz)
                x_max = max(x_max, max_ann * 1.02)
        ax_ms.set_xlim(x_min, x_max)

    # label annotated peaks with ion type, charge and m/z
    _label_annotated_peaks(ax_ms, msms)

    # cosmetics similar to plain plotting (no manual y-limit: spectrum_utils uses % scale)
    for spine in ("top", "right"):
        if spine in ax_ms.spines:
            ax_ms.spines[spine].set_visible(False)

    # No title on annotated spectra to avoid crowding the ribbon;
    # scan / score info is available in the filename and info text.
    # mass ladder panel for this peptide / modification state
    if peptide and assigned_mods:
        try:
            ladder = compute_mass_ladder(peptide, assigned_mods or "")
        except Exception as exc:  # pragma: no cover
            logger.warning(f"Failed to compute mass ladder for scan {scan_number}: {exc}")
            ladder = None
        if ladder is not None and not ladder.empty:
            text = format_mass_ladder_text(ladder)
            ax_ladder.text(
                0.0,
                1.0,
                text,
                ha="left",
                va="top",
                fontsize=7,
                family="monospace",
            )

    fig.tight_layout()

    # Place info text just below the x-axis in figure coordinates
    if info_text:
        fig.text(
            0.99,
            0.01,
            info_text,
            ha="right",
            va="bottom",
            fontsize=8,
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    logger.info(f"Wrote {out_path}")


def plot_top_psms_from_fragpipe(
    psm_path: Path,
    raw_dir: Path,
    out_dir: Path,
    n_top: int = 4,
    qvalue_cutoff: float = 0.01,
    min_localization: float = 0.0,
    phospho_loc_col: str = "STY:79.9663 Best Localization",
    annotated: bool = False,
    mz_min: Optional[float] = None,
    mz_max: Optional[float] = None,
    max_labels: Optional[int] = 10,
    max_length: Optional[int] = None,
    write_html: bool = False,
    verbose: bool = False,
) -> List[Path]:
    """
    High-level helper: select top PSMs from a FragPipe psm.tsv and plot them.

    Returns list of written image paths.
    """
    top_psms = select_top_phospho_psms(
        psm_path=psm_path,
        n_top=n_top,
        qvalue_cutoff=qvalue_cutoff,
        min_localization=min_localization,
        phospho_loc_col=phospho_loc_col,
        max_length=max_length,
    )

    if verbose:
        logger.info(
            f"Selected {len(top_psms)} PSMs (n_top={n_top}, max_aa_length={max_length})"
        )

    written: List[Path] = []

    groups: Dict[str, List[Tuple[int, int, int]]] = {}
    for idx, row in top_psms.iterrows():
        basename, scan_number, charge = parse_spectrum_id(row["Spectrum"])
        groups.setdefault(basename, []).append((idx, scan_number, charge))

    for basename, entries in groups.items():
        mzml_path = find_mzml_for_basename(raw_dir, basename)
        target_scans = [scan_number for _, scan_number, _ in entries]
        if verbose:
            logger.info(
                f"Indexing MS2 for {basename} ({len(entries)} PSMs, {len(target_scans)} target scans)"
            )
        ms2_index = build_ms2_index(mzml_path, target_scans)
        if verbose:
            logger.info(
                f"Indexed {len(ms2_index)} MS2 spectra for {basename}"
            )

        for idx, scan_number, charge in entries:
            row = top_psms.loc[idx]
            spectrum = ms2_index.get(scan_number)
            if spectrum is None:
                logger.warning(
                    f"Missing spectrum for scan {scan_number} in {mzml_path}, skipping."
                )
                continue

            mz_array = spectrum["m/z array"]
            intensity_array = spectrum["intensity array"]

            modified_peptide = row.get("Modified Peptide", row["Peptide"])
            hyperscore = float(row["Hyperscore"])
            qvalue = float(row["Qvalue"])

            precursor_mz = row.get("Calibrated Observed M/Z", None)
            if pd.isna(precursor_mz):
                precursor_mz = row.get("Observed M/Z", None)

            peptide = row.get("Peptide", None)
            assigned_mods = row.get("Assigned Modifications", None)

            # Build a small info text block with precursor / MS1 details
            info_parts: List[str] = []
            if precursor_mz is not None and not pd.isna(precursor_mz):
                info_parts.append(f"prec={float(precursor_mz):.4f} m/z")
            if charge:
                info_parts.append(f"z={charge}")
                if precursor_mz is not None and not pd.isna(precursor_mz):
                    mass = float(precursor_mz) * charge - charge * 1.007276
                    info_parts.append(f"M={mass:.1f} Da")

            parent_scan = row.get("Parent Scan Number", None)
            if parent_scan is not None and not pd.isna(parent_scan):
                try:
                    info_parts.append(f"MS1 scan={int(parent_scan)}")
                except Exception:
                    pass

            apex_rt = row.get("Apex Retention Time", None)
            if apex_rt is not None and not pd.isna(apex_rt):
                try:
                    rt_min = float(apex_rt) / 60.0
                    info_parts.append(f"Apex RT={rt_min:.2f} min")
                except Exception:
                    pass

            info_text = "  ".join(info_parts) if info_parts else None

            suffix = "annot" if annotated else "plain"
            out_name = (
                f"{basename}_scan_{scan_number}_z{charge}_"
                f"hyperscore_{hyperscore:.1f}_{suffix}.pdf"
            )
            out_path = out_dir / out_name

            if verbose:
                pep_len = len(peptide) if isinstance(peptide, str) else "NA"
                logger.info(
                    f"Plotting scan={scan_number} z={charge} len={pep_len} "
                    f"Hyperscore={hyperscore:.1f} -> {out_path.name}"
                )

            if annotated:
                plot_ms2_annotated(
                    mz_array=mz_array,
                    intensity_array=intensity_array,
                    modified_peptide=modified_peptide,
                    charge=charge,
                    scan_number=scan_number,
                    hyperscore=hyperscore,
                    qvalue=qvalue,
                    out_path=out_path,
                    precursor_mz=precursor_mz,
                    mz_min=mz_min,
                    mz_max=mz_max,
                    peptide=peptide,
                    assigned_mods=assigned_mods,
                    info_text=info_text,
                )

                if write_html:
                    html_name = out_name.replace(".pdf", ".html")
                    html_path = out_dir / html_name
                    try:
                        _plot_ms2_interactive_html(
                            mz_array=mz_array,
                            intensity_array=intensity_array,
                            charge=charge,
                            scan_number=scan_number,
                            precursor_mz=precursor_mz,
                            peptide=peptide,
                            assigned_mods=assigned_mods,
                            hyperscore=hyperscore,
                            qvalue=qvalue,
                            out_path=html_path,
                        )
                        if verbose:
                            logger.info(f"Wrote interactive HTML: {html_path.name}")
                    except Exception as exc:  # pragma: no cover
                        logger.warning(
                            f"Failed to write HTML for scan {scan_number}: {exc}"
                        )
            else:
                plot_ms2_plain(
                    mz_array=mz_array,
                    intensity_array=intensity_array,
                    modified_peptide=modified_peptide,
                    charge=charge,
                    scan_number=scan_number,
                    hyperscore=hyperscore,
                    qvalue=qvalue,
                    out_path=out_path,
                    mz_min=mz_min,
                    mz_max=mz_max,
                    max_labels=max_labels,
                    info_text=info_text,
                )

            written.append(out_path)

    return written
