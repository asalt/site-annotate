import re

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
