import re


RENAME = {
    "sample_01": "TMT_126",
    "sample_02": "TMT_127_N",
    "sample_03": "TMT_127_C",
    "sample_04": "TMT_128_N",
    "sample_05": "TMT_128_C",
    "sample_06": "TMT_129_N",
    "sample_07": "TMT_129_C",
    "sample_08": "TMT_130_N",
    "sample_09": "TMT_130_C",
    "sample_10": "TMT_131_N",
    "sample_11": "TMT_131_C",
    "sample_12": "TMT_132_N",
    "sample_13": "TMT_132_C",
    "sample_14": "TMT_133_N",
    "sample_15": "TMT_133_C",
    "sample_16": "TMT_134_N",
    "sample_17": "TMT_134_C",
    "sample_18": "TMT_135_N",
    "ion_126_128": "TMT_126",
    "ion_127_125": "TMT_127_N",
    "ion_127_131": "TMT_127_C",
    "ion_128_128": "TMT_128_N",
    "ion_128_134": "TMT_128_C",
    "ion_129_131": "TMT_129_N",
    "ion_129_138": "TMT_129_C",
    "ion_130_135": "TMT_123_N",
    "ion_130_141": "TMT_130_C",
    "ion_131_138": "TMT_131_N",
}
RENAME_SHORT = {
    "sample_01": "TMT_126",
    "sample_02": "TMT_127_N",
    "sample_03": "TMT_128_C",
    "sample_04": "TMT_129_N",
    "sample_05": "TMT_130_C",
    "sample_06": "TMT_131_N",
    "TMT_127": "TMT_127_N",
    "TMT_128": "TMT_128_C",
    "TMT_129": "TMT_129_N",
    "TMT_130": "TMT_130_C",
    "TMT_131": "TMT_131_N",
    "126": "TMT_126",
    "127": "TMT_127_N",
    "128": "TMT_128_C",
    "129": "TMT_129_N",
    "130": "TMT_130_C",
    "131": "TMT_131_N",
}
RENAME.update({y:y for y in RENAME.values()})
RENAME.update({x.lstrip("TMT_") : y for x,y in RENAME.items()})
RENAME.update({x.replace("_", "-"): y for x, y in RENAME.items()})
RENAME.update({x.replace("_", ""): y for x, y in RENAME.items()})

RENAME.update({x + "_intensity": y for x, y in RENAME.items()})
RENAME.update({v + "_intensity": v for v in RENAME.values()})


RENAME_SHORT.update({y:y for y in RENAME.values()})
RENAME_SHORT.update({x.lstrip("TMT_") : y for x,y in RENAME.items()})
RENAME_SHORT.update({x.replace("_", "-"): y for x, y in RENAME_SHORT.items()})
RENAME_SHORT.update({x.replace("_", ""): y for x, y in RENAME.items()})
RENAME_SHORT.update({x + "_intensity": y for x, y in RENAME.items()})
RENAME_SHORT.update(
    {v + "_intensity": v for v in RENAME_SHORT.values() if not v.endswith("_intensity")}
)

VALID_MODI_COLS = [
    "sty_79_9663",
    "k_42_0106",
    "k_43_0058",
    "k_43_0058",
    "k_114_0429",
]


PHOSPHOSITEPLUS_ANNOTATIONS = {
    "k_42_0106": "Acetylation_site_dataset",
    "sty_79_9663": "Phosphorylation_site_dataset",
    "disease": "Disease-associated_sites",
    "regulatory_sites": "Regulatory_sites",
}

# "Acetylation_site_dataset",  "Disease-associated_sites",  "Phosphorylation_site_dataset",  "PhosphoSitePlugin.jar",  "Regulatory_sites",


# TODO simply list the masses and calculate / look for matches without considering the aas
# to allow for modi mass on other aas not listed
MODI_ABBREVS = {
    "k_42_0106": "ac",
    "sty_79_9663": "p",
    "m_15_9949": "ox",
    "kr_14.0156": "me1",
    "kr_28_0313": "me2",
    "k_14_0156": "me1",  # this is redundant with above
    "k_28.0313": "me2",
    "k_42_0470": "me3",
    "acdefghiklmnpqrstvwyk_114_0429": "gg",  # this could be present
    "k_114_0429": "gg",  # this could be present
}


pattern = re.compile(r"[a-z]+_\d+_\d+$")


def get_all_columns(iterable):
    matches = [pattern.match(it) for it in iterable]
    matches = [x.group() for x in filter(None, matches)]
    return matches


# Allowed residues by modification column prefix
ALLOWED_AA_PREFIX = {
    "sty": set("STY"),
    "k": set("K"),
    "m": set("M"),
    "kr": set("KR"),
}


def allowed_residues_for_col(col: str) -> set | None:
    """Return allowed amino-acid letters for a given modification column.

    Falls back to K for known GG/ubiquitin remnants (114.0429) regardless of prefix.
    Returns None if no restriction should be applied.
    """
    if not isinstance(col, str):
        return None
    # Ubiquitin GG remnant on Lys
    if col.endswith("_114_0429"):
        return set("K")
    prefix = col.split("_", 1)[0]
    return ALLOWED_AA_PREFIX.get(prefix)


def get_all_possible_id_cols():

    return [
        "gene",
        "acc_id",
        "hu_chr_loc",
        "mod_rsd",
        "site_grp_id",
        "organism",
        "mw_kd",
        "domain",
        "site_+_7_aa",
        "lt_lit",
        "ms_lit",
        "ms_cst",
        "cst_cat#",
        "ambiguous_site",
        "psp_position",
        "fifteenmer",
        "sitename",
        "uniprot_id",
        "ENSP",
        "ENST",
        "ENSG",
        "geneid",
        "Gene",
        "ProteinID",
        "taxon",
        "symbol",
        "protein_start",
        "protein_end",
        "protein_start_psp",
        "position_relative",
        "position_absolut",
        "position_absolut_psp",
        "sitename2",
        # "all_possible_positions",
        # "all_possible_positions_psp",
        "AA",
        "site_id",
        # "hyperscore_best",
        # "rtscore_best",
        # "delta_mass_best",
        # "highest_prob_best",
    ]


POSSIBLE_SITE_ID_COLS = get_all_possible_id_cols()
