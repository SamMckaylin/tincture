"""Registry of public data files used by tincture case studies.

Checksums: figshare publishes MD5s (copied from the figshare API). GEO and Sanger do not
publish checksums, so their SHA-256 values were computed on the first complete download
and are recorded with the date in ``docs/data_status.md``.
"""

from __future__ import annotations

from tincture.datasets.cache import RemoteFile

# ---------------------------------------------------------------------------------------
# D1 sci-Plex (Srivatsan et al., Science 2020), GEO GSE139944, sample GSM4150378 (sciPlex3)
# ---------------------------------------------------------------------------------------
_GEO = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4150nnn/GSM4150378/suppl/GSM4150378_"
_GEO_LICENCE = "GEO public data; cite Srivatsan et al., Science 2020"


def _geo(name: str, size: int, sha256: str) -> RemoteFile:
    return RemoteFile(
        "sciplex", f"GSM4150378_{name}", _GEO + name, sha256=sha256, size=size, licence=_GEO_LICENCE
    )


SCIPLEX3: dict[str, RemoteFile] = {
    "umi_matrix": _geo(
        "sciPlex3_A549_MCF7_K562_screen_UMI.count.matrix.gz",
        3_096_224_606,
        "7d632716aa6ed0fc1780996003d0abc440ff78340609bceaaa4c8ade9345d00a",
    ),
    "cell_annotations": _geo(
        "sciPlex3_A549_MCF7_K562_screen_cell.annotations.txt.gz",
        2_770_264,
        "5db11591747a49f07ceb90a070971d60816591911f1022171adba7d8d032f379",
    ),
    "gene_annotations": _geo(
        "sciPlex3_A549_MCF7_K562_screen_gene.annotations.txt.gz",
        806_168,
        "fbe43028cfb75dc5ebf383cbbbb100da6b51e7be8e2e74f6354f7c59b9fda7c2",
    ),
    "pdata": _geo(
        "sciPlex3_pData.txt.gz",
        84_673_357,
        "860ba5c21846e2cfe97c39c805895982ef9db4b108cdf8b841a4a91be65024c2",
    ),
    "hash_metadata": _geo(
        "sciPlex3_A549_MCF7_K562_hashTable_metadata.txt.gz",
        64_715,
        "98a65191d343d149c5742cf6e7b734d403e853d165da800fa2daea41b978a2ff",
    ),
}

# ---------------------------------------------------------------------------------------
# D2 PRISM Repurposing 20Q2 (Corsello et al., Nat Cancer 2020), figshare 20564034 v1
# ---------------------------------------------------------------------------------------
_FIG = "https://ndownloader.figshare.com/files/"
_PRISM_LICENCE = "CC BY 4.0 (figshare 10.6084/m9.figshare.20564034.v1)"


def _fig(file_id: int, name: str, size: int, md5: str) -> RemoteFile:
    return RemoteFile("prism", name, f"{_FIG}{file_id}", md5=md5, size=size, licence=_PRISM_LICENCE)


PRISM_SECONDARY: dict[str, RemoteFile] = {
    "readme": _fig(
        36794613,
        "prism-repurposing-20q2-secondary-screen-readme",
        8046,
        "47419770d1edfcc36e0b11465f7d1fd0",
    ),
    "cell_line_info": _fig(
        36794589,
        "prism-repurposing-20q2-secondary-screen-cell-line-info.csv",
        84730,
        "4c3a5fbee769e8384064a8fbe4b4ed7d",
    ),
    "pooling_info": _fig(
        36794610,
        "prism-repurposing-20q2-secondary-screen-pooling-info.csv",
        331675,
        "8a484f5f4b704abd06cd3201e1d7692b",
    ),
    "curve_params": _fig(
        36794595,
        "prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv",
        290170269,
        "eb0c7ee3ccb9e480148c71fcdc97312e",
    ),
    "treatment_info": _fig(
        36794622,
        "prism-repurposing-20q2-secondary-screen-replicate-treatment-info.csv",
        16127430,
        "610fa36ba5b8c12d735917a7d555f65d",
    ),
    "logfold_change": _fig(
        36794607,
        "prism-repurposing-20q2-secondary-screen-logfold-change.csv",
        308481930,
        "51402a49933c652c7ba261db4fd076b6",
    ),
    "collapsed_treatment_info": _fig(
        36794619,
        "prism-repurposing-20q2-secondary-screen-replicate-collapsed-treatment-info.csv",
        4284816,
        "162da467eb97a67abdbdbee7af101091",
    ),
}


# ---------------------------------------------------------------------------------------
# D3 GDSC2 raw monotherapy data, release 8.5 (27 Oct 2023), Sanger
# ---------------------------------------------------------------------------------------
GDSC2_RAW = RemoteFile(
    "gdsc",
    "GDSC2_public_raw_data_27Oct23.csv.zip",
    "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_public_raw_data_27Oct23.csv.zip",
    sha256="7bc484eb0e859d58bfca30f72121cd353b26d35b28d18549464d496f3fc1bfc6",
    size=122_676_850,
    licence="GDSC terms: research use, no redistribution (could not be re-checked; see "
    "docs/data_status.md)",
)


# --- D9: LINCS L1000 phase II (GEO GSE70138), checksums from GSE70138_SHA512SUMS.txt.gz ---
_LINCS = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/"
_LINCS_LICENCE = "GEO public data; cite Subramanian et al., Cell 2017 (CMap / LINCS L1000)"


def _lincs(name: str, sha512: str) -> RemoteFile:
    return RemoteFile("lincs", name, _LINCS + name, sha512=sha512, licence=_LINCS_LICENCE)


LINCS_GSE70138 = {
    "level5": _lincs(
        "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz",
        "9d078903d3d028f37b0bae584c1afe48f5bab0082e61f6e542de471674c78946cff710132b8e95beb997a39dc2dd6442fdda0fbb1cc9b2399ef21cf636975c9f",
    ),
    "sig_info": _lincs(
        "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
        "8d16dfc3a8fdfba572797f2519f9bc96b687c3666c78d386b53fc997164dd2a8ab5953c0b323795ec7f62de5de4ddbd0285e79f8be48487156479dc8ee0128c9",
    ),
    "sig_metrics": _lincs(
        "GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz",
        "b9a79bc8835aa19bb46a824e1458d15ecc6d5ceff16afee1f02e759f495a7c2e8b4c99761e1d00ca8db10fc4789a6ba65958f9aabf3daa557681fd542a5ab89a",
    ),
    "pert_info": _lincs(
        "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz",
        "b0df8aaa5653c581bd69a3dddd6eaa64e6296a70f40f50b1b80a2854abbee17a9d66e874dc00657c045be738df521f2e673115c1aae9f5972e820e88f0032ddc",
    ),
    "gene_info": _lincs(
        "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz",
        "b2df105ebf7a69799b23130d03b7da769bfd6ab83fd69cbb419c91c78ae069ae3e5d474e09d36e3ac4271f3b7fc56b24caaf1ac5ece9554ea7738b15535dfefd",
    ),
}


# --- D7 / D8: kinase panels. Publishers block scripted downloads (docs/data_status.md), so
# these are placed by hand; checksums pin the exact files the loaders were written against.
NANOBRET_VASTA2018 = RemoteFile(
    "nanobret",
    "vasta2018_mmc2.xlsx",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC5814754/",
    sha256="b6348a34e7c62946eb5701d0adac35d15b73be0354de143c4aca0eecec677ac3",
    size=26112,
    licence="Table S1 of Vasta et al., Cell Chem Biol 2018 (NIH author manuscript)",
    manual=True,
)
DAVIS2011_TABLE4 = RemoteFile(
    "davis",
    "davis2011_table4.xls",
    "https://www.nature.com/articles/nbt.1990",
    sha256="e47821b2e959ff22aa6d1613b236a63fb530f0da19cfa8d45408598cf27498d7",
    size=304640,
    licence="Supplementary Table 4 of Davis et al., Nat Biotechnol 2011",
    manual=True,
)


# --- D4: GDSC anchored combinations (Jaaks, Coker et al., Nature 2022), figshare. The
# figshare licence field says CC BY 4.0, but every item's description restricts use to
# internal research and education (no resale, no commercial services): treated as the
# restrictive GDSC terms, so outputs stay local (docs/data_status.md).
_JAAKS_LICENCE = (
    "Jaaks et al. 2022 (figshare): internal research and educational use only; "
    "no resale or commercial services"
)


def _jaaks(file_id: int, name: str, size: int, md5: str) -> RemoteFile:
    return RemoteFile(
        "jaaks2022", name, f"{_FIG}{file_id}", md5=md5, size=size, licence=_JAAKS_LICENCE
    )


JAAKS2022: dict[str, RemoteFile] = {
    "glossary": _jaaks(34006715, "Glossary.xlsx", 12193, "8706ea6fbab6dce6210697e8849c6dbf"),
    "original_fitted": _jaaks(
        34006655,
        "Original screen_All tissues_fitted data.csv",
        198893607,
        "786589e5eb556670cc7d24356487f8fd",
    ),
    "validation_fitted": _jaaks(
        34007258,
        "Validation screen_All tissues_fitted data.csv",
        11774397,
        "7cbc398261624ff3d9965669c1645a51",
    ),
}

# --- D10: MIX-seq (McFarland, Paolella et al., Nat Commun 2020), figshare 10298696 v3
_MIX_LICENCE = "CC BY 4.0 (figshare 10.6084/m9.figshare.10298696.v3)"


def _mix(file_id: int, name: str, size: int, md5: str) -> RemoteFile:
    return RemoteFile("mixseq", name, f"{_FIG}{file_id}", md5=md5, size=size, licence=_MIX_LICENCE)


MIXSEQ: dict[str, RemoteFile] = {
    "readme": _mix(24848579, "README.txt", 5515, "681292e2dc103b37f8ebb5d3e3b16e34"),
    "supp_tables": _mix(
        18716510, "Supplementary Tables.xlsx", 87927, "42b8df97cfb407c81619c877a03f2e47"
    ),
    "trametinib_tc": _mix(
        20000609, "trametinib_tc_expt5.zip", 291873754, "97b3712228cc1bf474d5b0fa2777c48a"
    ),
}

# --- D5: GDSC full 7x7 matrix combinations (Bashi et al., Cancer Discov 2024), figshare
# project 163378 (named in the paper's Data Availability statement). CC BY 4.0 on figshare;
# GDSC data are still treated as internal-use (outputs stay local, docs/data_status.md).
_BASHI_LICENCE = "CC BY 4.0 (figshare project 163378; Bashi et al., Cancer Discov 2024)"


def _bashi(file_id: int, name: str, size: int, md5: str) -> RemoteFile:
    return RemoteFile(
        "bashi2024", name, f"{_FIG}{file_id}", md5=md5, size=size, licence=_BASHI_LICENCE
    )


BASHI2024: dict[str, RemoteFile] = {
    "well_stats": _bashi(
        45451363,
        "matrix_well_stats_GDSC_007_010_plus_Spiper1.zip",
        241088863,
        "ff3dab95ad59a1c956e3b11c397dcbd8",
    ),
    "gdsc007": _bashi(
        39806602, "gdsc-007_matrix_results.csv", 14947830, "0223b9ad0921ac4974bd4479aea8853b"
    ),
    "gdsc008": _bashi(
        39806650, "gdsc-008_matrix_results.csv", 16897258, "ad402a31c2749652632e2dd94e57c6ee"
    ),
    "gdsc009": _bashi(
        39806656, "gdsc-009_matrix_results.csv", 15492311, "9119449ae1d9567a876ab8b87ca4005a"
    ),
    "gdsc010": _bashi(
        44311796, "gdsc-010_matrix_results.csv", 8026089, "a50cab826d6fda7853a54bad97f23e15"
    ),
    "sandpiper01": _bashi(
        39806584, "sandpiper-01_matrix_results.csv", 1293674, "3d3bdc4765989c0ab0f7abf781394cb4"
    ),
}

# --- D6: GDSC pan-cancer anchored (Vis et al., Cell Rep Med 2024), open-access supplement
# (Data S1 holds pancan_combi_fitted_data_jan21.tsv.gz). Europe PMC serves the supplement as
# one zip and publishes no checksum: SHA-256 computed on the first download (2026-11-18).
VIS2024_SUPPLEMENT = RemoteFile(
    "vis2024",
    "PMC11384948_supplementary.zip",
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11384948/supplementaryFiles",
    sha256="5f19097023fbc446269b1c6f0d317d3aefa53de87e1d36fdcbf1650722e1f358",
    size=53745336,
    licence="Vis et al., Cell Rep Med 2024 (open access, CC BY); GDSC data: internal use",
)

# --- Cell Model Passports model annotations (MSI status, tissue) for GDSC lines. The URL is
# "latest", so the checksum pins the copy the case studies were run with (2026-11-18).
CMP_MODEL_LIST = RemoteFile(
    "cmp",
    "model_list_latest.csv.gz",
    "https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz",
    sha256="333affaa690c9f7e517467a6f78b00ffb418cb8449e2eb1cc5d93d85f7b0862e",
    size=145029,
    licence="Cell Model Passports: non-commercial use",
)
