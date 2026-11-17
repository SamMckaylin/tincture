# API reference

_Generated from docstrings by `docs/make_api.py` (tincture 0.0.1.dev0)._ Do not edit by hand: run `python docs/make_api.py`.

```python
import tincture as tc
```

Top-level names: `tc.PlotResult`, `tc.TinctureDataError`, `tc.TinctureError`, `tc.TinctureWarning`, `tc.__version__`, `tc.compass`, `tc.demo`, `tc.dose_paths`, `tc.effect_ladder`, `tc.export`, `tc.fingerprint`, `tc.plate_atlas`, `tc.render_mpl`, `tc.retrievability`, `tc.theme`.

## Plots

Every plot returns a `PlotResult` (figure, caption, methods, warnings, classic fallback). Call them as `tc.<name>(...)`.

### `tincture.plots.retrievability`

Retrievability: which compounds have a real, reproducible phenotype? (HANDOVER §5.5)

#### `tc.retrievability(profiles: PerturbationSet, group_by: Sequence[str] = ('compound_id', 'dose'), focus: Sequence[str] = (), negatives: str = 'controls', null_size: int = 10000, alpha: float = 0.05, seed: int = 0, title: str = 'Retrievability', footer: Sequence[str] = (), null: str = 'permutation') -> PlotResult`

Retrievability gate for replicate-level profiles.

``focus`` holds group keys (``"<compound_id>|<dose>"``) to highlight; the first one
also gets the CAT panel and is used for the classic replicate scatter.

#### `compute(profiles: PerturbationSet, group_by: Sequence[str] = ('compound_id', 'dose'), negatives: str = 'controls', null_size: int = 10000, alpha: float = 0.05, seed: int = 0, null: str = 'permutation') -> MapResult`

Retrievability statistics for replicate-level profiles (controls = dose 0).

#### `gate_mask(profiles: PerturbationSet, result: MapResult, group_by: Sequence[str] = ('compound_id', 'dose')) -> pd.Series`

Boolean per sample: True if its condition is retrievable (vehicle always True).

Also stored in ``profiles.adata.uns['tincture_retrievable']`` for other plots' ``gate``.

### `tincture.plots.plate_atlas`

Plate Atlas: can we trust these plates? (HANDOVER §5.6; QC / supplementary figure)

#### `tc.plate_atlas(pt: PlateTable, neg: Selector | None = None, pos: Selector | None = None, edge: Selector | str = 'auto', plate: str | None = None, title: str = 'Plate Atlas', footer: tuple[str, ...] = ()) -> PlotResult`

Plate QC figure. See :func:`tincture.stats.zprime.plate_qc` for the selectors.

### `tincture.plots.compass`

Compass: what known mechanism does this compound resemble, and how strongly? (HANDOVER §5.1)

#### `tc.compass(signatures: Mapping[str, Signature], library: ReferenceLibrary, focus: str | None = None, gate: set[str] | None = None, spoke_labels: Mapping[str, str] | None = None, max_spokes: int = 8, label_top: int = 4, n_boot: int = 200, title: str = 'Compass', footer: Sequence[str] = (), gate_mode: str = 'drop') -> PlotResult`

Compass plot.

``signatures``: compound id → :class:`Signature` (one per compound, e.g. top dose).
``gate``: compound ids that passed retrievability; others are drawn hollow grey.

#### `select_spokes(lib: ReferenceLibrary, max_spokes: int = 8) -> tuple[ReferenceLibrary, list[str]]`

Keep the ``max_spokes`` most *coherent* references, then order them by similarity.

Coherence = mean pairwise cosine between a class's member signatures (do its members
agree?) × √(members). Magnitude alone would favour classes with one or two very strong
compounds and drop well-defined classes with moderate effects (e.g. Aurora kinase
inhibitors in sci-Plex), which would then wrongly look "novel". Classes whose members do
not agree are never used: mean cosine ≥ ``MIN_COHERENCE`` with ≥ 3 members, or
≥ ``MIN_COHERENCE_PAIR`` for a two-member class.
Returns the ordered library and the names left out.

#### `class_coherence(lib: ReferenceLibrary, name: str) -> tuple[float, int]`

(mean pairwise cosine between members, number of members).

#### `short_name(name: str) -> str`

Drop salt forms and code synonyms in parentheses: 'Abexinostat (PCI-24781)' →
'Abexinostat'. Names that are only a code are kept.

### `tincture.plots.dose_paths`

Dose Paths: does a compound's mechanism change with dose? (HANDOVER §5.2)

#### `tc.dose_paths(sigset: SignatureSet, compounds: Sequence[str], library: ReferenceLibrary, stress_genes: Collection[str] | None = None, focus: str | None = None, y_axis: str = 'stress', cell_loss: Mapping[tuple[str, float], float] | None = None, gate: Collection[str] | None = None, spoke_labels: Mapping[str, str] | None = None, reference: str | None = None, n_boot: int = 500, title: str = 'Dose Paths', footer: Sequence[str] = ()) -> PlotResult`

Dose Paths for ``compounds`` from one stratum's :class:`SignatureSet`.

``stress_genes``: gene symbols for the stress axis (see
:func:`tincture.stats.paths.stress_gene_set`). ``y_axis="cell_loss"`` uses
``cell_loss[(compound, dose)]`` instead. ``gate``: compound ids that passed the
retrievability gate; others are dropped and the caption says so. ``reference``: use
this reference for x in every panel (default: the focus compound's top reference for one
panel, each compound's own top reference in small multiples).

#### `compute(sigset: SignatureSet, compounds: Sequence[str], library: ReferenceLibrary, stress_mask: np.ndarray | None, stratum: _Stratum, cell_loss: Mapping[tuple[str, float], float] | None = None, shared_reference: str | None = None, n_boot: int = 500) -> tuple[dict[str, DosePath], dict[str, str], dict[str, PathAxes]]`

Dose path per compound, its top reference and the axes used.

### `tincture.plots.fingerprint`

Fingerprint Strip: what did each compound do, at pathway level? (HANDOVER §5.3)

#### `tc.fingerprint(data: SignatureSet | Sequence[Signature] | Mapping[str, Signature], gene_sets: pd.DataFrame | None = None, focus: str | None = None, references: Sequence[str] = (), dose: float | str | None = None, symbols: np.ndarray | None = None, alpha: float = 0.05, clip: float | None = None, order: str = 'cluster', drop_empty: bool = True, max_jaccard: float = 0.5, gate: Collection[str] | None = None, source: str = 'MSigDB Hallmark', title: str = 'Fingerprint Strip', footer: Sequence[str] = ()) -> PlotResult`

Fingerprint Strip for a set of signatures.

``data``: a :class:`SignatureSet` (top dose per compound unless ``dose`` is given) or
signatures. ``gene_sets``: long ``source, target`` table; default MSigDB Hallmark from
the cache. ``references``: compound ids drawn first, above a gap, in the given order.
``order``: ``"cluster"`` (average-linkage on score correlation, optimal leaf order) or
``"given"``. ``gate``: compound ids that passed retrievability; others are dropped.

#### `nice_clip(scores: np.ndarray) -> float`

Symmetric colour limit: the 95th percentile of |score|, rounded up to a nice value.

#### `compute(sigs: list[Signature], gene_sets: pd.DataFrame, symbols: np.ndarray | None, max_jaccard: float = 0.5, min_size: int = 5) -> tuple[pd.DataFrame, pd.DataFrame, list[tuple[str, str, float]]]`

ULM scores, BH q-values (rows × pathways) and the dedup decisions.

### `tincture.plots.effect_ladder`

Effect Ladder: which features moved, and how confident are we? (HANDOVER §5.4)

#### `tc.effect_ladder(data: Signature | SignatureSet, compound: str | None = None, dose: float | None = None, gene_sets: pd.DataFrame | None = None, pathway: str | None = None, top_n: int = 30, alpha: float = 0.05, name: str | None = None, title: str = 'Effect Ladder', footer: Sequence[str] = ()) -> PlotResult`

Effect Ladder for one signature (a :class:`Signature`, or ``compound``/``dose`` of a
:class:`SignatureSet`). ``gene_sets``: long ``source, target`` table (e.g.
``load_geneset("hallmark").table``) for grouping and ``pathway=`` drill-down.

#### `compute(sig: Signature, pvalues: np.ndarray | None, symbols: np.ndarray | None, gene_sets: pd.DataFrame | None, top_n: int = 30, pathway: str | None = None, alpha: float = 0.05) -> tuple[pd.DataFrame, dict[str, int]]`

Table of the chosen features (in drawing order) and transcriptome-wide counts.

#### `volcano_spec(sig: Signature, pvalues: np.ndarray | None, symbols: np.ndarray | None, alpha: float = 0.05) -> PlotSpec`

Classic volcano: shrunk log2 FC vs −log10 p, q < alpha coloured.

## Results and export

### `tincture.result`

PlotResult: what every public plot function returns.

#### class `PlotResult(spec: PlotSpec, classic: PlotSpec | None = None)`

A built plot: spec, caption, methods, warnings and an optional classic fallback.

The figure is rendered lazily (and cached per preset/height) with
:func:`tincture.render.mpl.render_mpl`.

- `caption`
- `methods`
- `warnings`
- `meta`
- `figure(preset: str | JournalPreset = 'nature_double', height_mm: float | None = None) -> Figure`
- `classic_figure(preset: str | JournalPreset = 'nature_double', height_mm: float | None = None) -> Figure | None`
- `save(path: str | Path, preset: str | JournalPreset = 'nature_double', height_mm: float | None = None, formats: tuple[str, ...] = ('pdf', 'png'), classic: bool = True, check: bool = True) -> dict[str, Path]`: Export figure, caption/methods sidecar and (optionally) the classic fallback.

### `tincture.export`

Export figures at journal sizes, with caption/methods sidecars.

#### `check_figure(fig: Figure, preset: str | JournalPreset) -> None`

Raise if the figure breaks the accessibility rules for ``preset``.

#### `export(obj: PlotResult | PlotSpec | Figure, path: str | Path, preset: str | JournalPreset = 'nature_double', height_mm: float | None = None, formats: tuple[str, ...] = ('pdf', 'png'), classic: bool = True, check: bool = True) -> dict[str, Path]`

Export a plot at a journal preset.

Parameters
----------
obj
    A :class:`PlotResult`, a :class:`PlotSpec` or any matplotlib ``Figure`` (resized
    to the preset width; no sidecar).
path
    Output path *without* extension (an extension, if given, is replaced).
preset, height_mm
    Size. See :mod:`tincture.theme.presets`.
formats
    File formats to write.
classic
    Also export the classic fallback when the plot has one.
check
    Enforce the minimum font size and the colour-blind check.

Returns
-------
Mapping of written file kind (``"pdf"``, ``"md"``, ``"classic_pdf"`` ...) to path.

## Data loaders

### `tincture.io.sciplex`

Loader for sci-Plex3 (Srivatsan et al., Science 2020; GEO GSE139944, GSM4150378).

#### `clean_smiles(s: object) -> str | None`

Strip a ChemAxon CXSMILES extension (``|...|``) to plain SMILES.

#### `read_cell_metadata(root: Path | None = None, download: bool = True) -> pd.DataFrame`

Per-cell metadata in tincture's vocabulary, indexed by the 0-based matrix column.

Columns: ``cell``, ``well``, ``plate``, ``cell_line``, ``compound_id``,
``compound_name``, ``dose`` (M), ``time`` (h), ``replicate``, ``target``, ``pathway``,
``pathway_level_1``, ``pathway_level_2``, ``smiles``, ``n_umi``, ``assigned``.

#### `read_genes(root: Path | None = None, download: bool = True) -> pd.DataFrame`

Gene table indexed by 0-based matrix row, with a ``human`` flag.

#### `build_pseudobulk(root: Path | None = None, out: Path | None = None, chunksize: int = 20000000, progress: Callable[[str], None] | None = <built-in function print>) -> Path`

Stream the triplet matrix once and sum counts per well (human genes only).

Writes an ``.h5ad`` with raw summed counts in ``X`` (CSR, float32) and returns its path.

#### `load_sciplex3(root: Path | None = None, cache: Path | None = None, cell_lines: list[str] | None = None, time: float | list[float] | None = 24.0, compounds: list[str] | None = None, build: bool = False) -> PerturbationSet`

sci-Plex3 as a per-well pseudobulk :class:`PerturbationSet` (raw summed UMI counts).

Parameters
----------
cell_lines, time, compounds
    Optional filters. ``time=24`` (default) keeps the 24 h screen; A549 also has 72 h.
    ``compounds`` are catalog ids or names; vehicle is always kept.
build
    Build the pseudobulk cache if missing (streams the 3.1 GB matrix; tens of minutes).

### `tincture.io.prism`

Loader for the PRISM Repurposing 20Q2 secondary screen (Corsello et al., Nat Cancer 2020).

#### `compound_core(broad_id: pd.Series) -> pd.Series`

'BRD-K12345678-001-02-3' -> 'BRD-K12345678' (compound, not batch).

#### `choose_screens(t: pd.DataFrame, screen: str) -> pd.Series`

Map compound_id -> screen_id to use.

#### `load_prism_secondary(compounds: list[str] | None = None, cell_lines: list[str] | None = None, screen: str = 'auto', include_failed_str: bool = False, root: Path | None = None, download: bool = True) -> DoseResponseSet`

PRISM secondary viability as a :class:`DoseResponseSet`.

Parameters
----------
compounds
    Compound names (e.g. ``"dabrafenib"``) or Broad compound ids (``"BRD-K…"``).
    ``None`` loads everything (≈1,500 compounds; several GB in memory).
cell_lines
    DepMap ids (``"ACH-…"``) or CCLE names to keep.
screen
    ``"auto"`` (per-compound preference MTS010 > HTS002 > MTS006 > MTS005) or one id.

### `tincture.io.gdsc`

Loader for GDSC2 raw monotherapy screening data (release 8.5, 27 Oct 2023, Sanger).

#### `tag_to_well_type(tags: pd.Series) -> pd.Series`

Map GDSC ``TAG`` values to tincture well types; raise on unknown tags.

#### `collapse_duplicate_wells(raw: pd.DataFrame) -> pd.DataFrame`

Collapse wells listed twice (DMSO solvent row + drug row) into one row.

#### `read_raw(path: Path | None = None, plates: list[int] | None = None, cell_lines: list[str] | None = None, max_plates: int | None = None, chunksize: int = 2000000) -> pd.DataFrame`

Read (a filtered subset of) the raw CSV inside the release zip.

#### `load_gdsc2_plates(plates: list[int] | None = None, cell_lines: list[str] | None = None, max_plates: int | None = None, path: Path | None = None) -> PlateTable`

GDSC2 raw wells as a :class:`PlateTable` (value = raw luminescence intensity).

Plates are ordered by scan date (``plate_order``), which is what control-drift
analysis needs. ``tag`` keeps the original GDSC label (e.g. to choose NC-1 vs B for
Z′). Doses are molar.

## Datasets and caches

### `tincture.datasets.cache`

Download-and-cache helpers with checksums.

#### class `DatasetUnavailable`

A registered file is not cached and could not be downloaded.


#### class `RemoteFile(dataset: str, name: str, url: str, md5: str | None = None, sha256: str | None = None, size: int | None = None, licence: str = '')`

One downloadable file.

Attributes
----------
dataset
    Cache folder name (e.g. ``"sciplex"``).
name
    File name in the cache.
url
    Source URL.
md5, sha256
    Expected checksum(s); at least one must be given.
size
    Expected size in bytes (checked before hashing; cheap early failure).
licence
    Short licence note (printed in reports and errors).


#### `data_home() -> Path`

Root of the local data cache.

#### `verify(f: RemoteFile, path: Path) -> bool`

True if ``path`` has the registered size and checksum(s).

#### `is_cached(f: RemoteFile, root: Path | None = None, check: bool = False) -> bool`

Whether the file is present (and, with ``check=True``, verified).

#### `fetch(f: RemoteFile, root: Path | None = None, download: bool = True, timeout: float = 600.0) -> Path`

Return the verified local path of ``f``, downloading it if needed.

Raises
------
DatasetUnavailable
    If the file is missing and cannot be downloaded, or the checksum fails.

### `tincture.datasets.genesets`

Gene-set libraries (MSigDB Hallmark, PROGENy) fetched via decoupler/OmniPath and cached.

#### `load_geneset(name: str, root: Path | None = None, download: bool = True) -> GeneSetLibrary`

Load (fetching once if needed) a human gene-set library: ``hallmark`` or ``progeny``.

## Data model

### `tincture.model.perturbation`

PerturbationSet: samples × features with perturbation metadata, wrapping AnnData.

#### class `PerturbationSet(adata: ad.AnnData, source: str = '')`

Samples (cells, wells or pseudobulk groups) × features with perturbation metadata.

``obs`` holds one row per sample with at least ``compound_id``, ``dose`` (molar; 0 =
vehicle) and ``replicate``. ``var`` holds features (genes, pathways, morphology
features). ``X`` may be dense or sparse. Additional AnnData layers are preserved.

Signatures (``.signatures()``) are added in Milestone 2 (``tincture.stats``).

- `obs`
- `var`
- `obs_names`
- `var_names`: Feature names (e.g. Ensembl gene ids).
- `X`
- `shape`
- `compounds`: Treated compound ids (vehicle excluded), sorted.
- `is_vehicle`
- `doses(compound_id: str) -> np.ndarray`
- `subset(obs: np.ndarray | pd.Series | None = None, features: Sequence[str] | np.ndarray | None = None, **equals: Any) -> PerturbationSet`: Subset samples and/or features.
- `signatures(**kwargs: Any) -> list[Any]`: Signatures vs same-stratum vehicle.
- `pseudobulk(by: Sequence[str], min_samples: int = 1) -> PerturbationSet`: Sum ``X`` over samples sharing the values of ``by`` (e.g. well or replicate).
- `from_long(cls, df: pd.DataFrame, source: str = '') -> PerturbationSet`: Build from the canonical long table (one row per sample × feature).
- `to_long() -> pd.DataFrame`: Inverse of :meth:`from_long` (dense; intended for small sets and tests).

### `tincture.model.signature`

Signature: one perturbation's effect vector with uncertainty.

#### class `Signature(features: np.ndarray, effect: np.ndarray, se: np.ndarray | None = None, boot: np.ndarray | None = None, compound_id: str = '', dose: float = nan, time: float = nan, cell_line: str = '', n_replicates: int = 0, method: str = '', meta: dict[str, Any] = <factory>)`

Effect of one perturbation (compound × dose × time × cell line) on each feature.

Bootstrap draws are stored (not just point estimates) so downstream confidence
intervals are honest. ``effect`` is usually a (shrunk) log fold change vs vehicle.

- `n_features`
- `magnitude`: Euclidean norm of the effect vector (NaNs ignored).
- `ci(level: float = 0.95) -> tuple[np.ndarray, np.ndarray]`: Per-feature CI: bootstrap percentiles if draws exist, else normal from ``se``.
- `align(features: np.ndarray | list[str]) -> Signature`: Reorder/restrict to ``features``; missing features become NaN.

### `tincture.model.reference`

ReferenceLibrary: named reference signatures with mechanism annotations.

#### `library_from_signatures(signatures: Mapping[str, np.ndarray], features: np.ndarray, classes: Mapping[str, Sequence[str]], version: str, source: str = '', citation: str = '', min_members: int = 2) -> ReferenceLibrary`

Build a library from per-compound signatures and class membership.

``signatures`` maps compound id → effect vector (aligned to ``features``); classes with
fewer than ``min_members`` available members are skipped (recorded in ``meta``).

### `tincture.model.dose_response`

Dose-response containers.

#### class `FourPL(bottom: float, top: float, log10_ec50: float, hill: float)`

y = bottom + (top - bottom) / (1 + (x / ec50) ** hill), x in molar.

``top`` is the response at zero dose and ``bottom`` at infinite dose (for an
inhibitor with ``hill > 0``). This matches drc's LL.4 ``d`` (upper) and ``c`` (lower).

- `ec50`
- `inverse(y: float) -> float`: Concentration giving response ``y``; NaN if the curve never reaches it.

#### class `DoseResponse(compound_id: str, dose: np.ndarray, response: np.ndarray, replicate: np.ndarray | None = None, assay: str = '', context: dict[str, Any] = <factory>, fit: FourPL | None = None, r2: float = nan, converged: bool | None = None, ec50_ci: tuple[float, float] | None = None, fit_source: str = '', response_label: str = 'viability (fraction of vehicle)')`

Raw points and (optionally) a fitted curve for one compound in one assay context.

- `n_doses`
- `dose_range`
- `quality_flags() -> list[str]`: Reasons this curve's fit should not be trusted (empty list = usable).
- `fit_ok`

#### class `DoseResponseSet(points: pd.DataFrame, fits: pd.DataFrame | None = None, contexts: pd.DataFrame | None = None, compounds: pd.DataFrame | None = None, assay: str = '', fit_source: str = '', source: str = '')`

Many dose-response curves as two tables.

``points``: one row per measurement with :data:`POINT_COLS` (+ any metadata).
``fits``: one row per curve with :data:`FIT_COLS` (+ provider extras, e.g. auc, ic50).
``contexts``: one row per ``context_id`` (e.g. cell line) with descriptive columns.

- `keys() -> pd.DataFrame`
- `get(compound_id: str, context_id: str) -> DoseResponse`

### `tincture.model.plate`

PlateTable: well-level values with plate geometry.

#### `row_label(i: int) -> str`

0 -> 'A', 25 -> 'Z', 26 -> 'AA', 31 -> 'AF' (1536-well convention).

#### `position_to_rowcol(position: np.ndarray, n_wells: int, one_based: bool = True, order: str = 'row-major') -> tuple[np.ndarray, np.ndarray]`

Convert linear well positions to 0-based (row, col).

#### class `PlateTable(data: pd.DataFrame, n_wells: int, source: str = '', value_label: str = 'intensity')`

Well-level measurements on one or more plates of a single format.

``data`` columns: ``plate``, ``row``, ``col`` (0-based), ``well_type`` (one of
:data:`tincture.model.measurements.WELL_TYPES`), ``value``; optional ``compound_id``,
``dose`` (molar), ``cell_line``, ``tag`` (source-specific label), ``plate_order``.

Large loaders (GDSC: 13.5 M wells) store text columns as pandas categoricals to save
memory; use :meth:`relative_to` rather than ``.map`` on them for per-plate maths.

- `shape`
- `plates`
- `is_edge() -> np.ndarray`: Boolean mask of outermost-ring wells.
- `plate_matrix(plate: str, well_type: str | None = None) -> np.ndarray`: Values of one plate as a (rows × cols) array; NaN where absent/filtered.
- `relative_to(mask: np.ndarray | pd.Series) -> pd.Series`: Each well's value divided by its plate's median over the wells in ``mask`` (e.g. negative controls). Works with categorical plate columns.
- `well_name() -> pd.Series`
- `subset_plates(plates: list[str]) -> PlateTable`

### `tincture.model.measurements`

The canonical long measurement table (HANDOVER §4.1) and its validator.

#### `validate_long(df: pd.DataFrame, require_vehicle: bool = False) -> pd.DataFrame`

Validate and normalise a long measurement table.

Returns a copy with canonical dtypes (``dose``/``value``/``time``/``se`` float, ids as
strings). Rows whose ``value`` is NaN are kept (missing is information); rows with a
missing ``dose`` are rejected because dose is a first-class axis.

Raises
------
TinctureDataError
    With a fix suggestion, for missing columns, negative or non-numeric doses,
    unknown ``well_type`` values, or (if ``require_vehicle``) no vehicle rows.

#### `round_dose(dose: np.ndarray | pd.Series, sig: int = 6) -> np.ndarray`

Round molar doses to ``sig`` significant figures.

Unit conversion (e.g. µM / 1e6) leaves binary noise (9.769999999999999e-09), which
breaks grouping by dose; rounding restores the intended values (9.77e-09).

#### `add_constant_dose(df: pd.DataFrame, dose: float = nan) -> pd.DataFrame`

Explicitly add a ``dose`` column for data without dose information.

With the default (NaN) the table will still fail validation: this is deliberate so
that dose-free data can only enter dose-aware plots after the user states a dose.

## Statistics

### `tincture.stats.signatures`

Perturbation signatures from replicate-level count profiles.

#### `logcpm(counts: np.ndarray, prior: float = 0.5) -> np.ndarray`

log2 counts per million (samples in rows).

#### `shrink_normal(effect: np.ndarray, se: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]`

Normal-normal empirical-Bayes shrinkage of effects towards 0.

Prior effect ~ N(0, τ²), τ² = max(mean(effect² − se²), floor). Returns posterior mean,
posterior SD and τ².

#### class `SignatureSet(signatures: list[Signature], stratum: dict[str, Any], features: np.ndarray, df_prior: float, df_residual: float, n_vehicle: int, method: str, pvalues: dict[tuple[str, float], np.ndarray] = <factory>, feature_symbols: np.ndarray | None = None)`

All signatures of one stratum plus the shared statistics behind them.

- `get(compound_id: str, dose: float) -> Signature`
- `matrix(shrunk: bool = True) -> pd.DataFrame`: Signatures × features effect matrix (shrunk by default).

#### `compute_signatures(ps: PerturbationSet, strata: Sequence[str] = ('cell_line', 'time'), min_cpm: float = 5.0, min_frac: float = 0.25, trend: bool = True, shrink: bool = True, min_replicates: int = 1) -> list[SignatureSet]`

Compute one :class:`SignatureSet` per stratum (see module docstring).

Raises
------
TinctureDataError
    No vehicle samples in a stratum, or no feature passes the expression filter.

#### `replicate_profiles(ps: PerturbationSet, strata: Sequence[str] = ('cell_line', 'time'), min_cpm: float = 5.0, min_frac: float = 0.25, n_features: int | None = 2000) -> PerturbationSet`

Per-sample profiles for retrieval: logCPM minus same-stratum vehicle mean, scaled
by the vehicle SD per feature (robust z vs vehicle).

Feature selection is **label-free** (highest variance across all samples of the
stratum after the expression filter), so retrieval is not inflated by choosing
features that separate the very conditions being retrieved.

### `tincture.stats.ebayes`

Empirical-Bayes moderation of gene-wise variances (Smyth 2004, limma ``squeezeVar``).

#### `trigamma_inverse(x: float) -> float`

Solve trigamma(y) = x for y > 0 (Newton iteration as in limma::trigammaInverse).

#### class `FDistFit(s0_sq: np.ndarray, d0: float)`

Prior for the variances: ``s0²`` (scale, per gene if trended) and ``d0`` (df).


#### `fit_f_dist(s2: np.ndarray, df: float | np.ndarray, covariate: np.ndarray | None = None) -> FDistFit`

Moment estimation of the scaled-F prior (limma ``fitFDist``).

#### `squeeze_var(s2: np.ndarray, df: float | np.ndarray, covariate: np.ndarray | None = None) -> SqueezedVar`

Posterior (moderated) variances. ``df_total = df + d0`` (capped as in limma).

### `tincture.stats.pathways`

Pathway activity scores: univariate linear model (ULM), as in decoupler.

#### class `PathwayScores(score: pd.DataFrame, pvalue: pd.DataFrame, n_features: pd.Series, source: str = '')`

Scores and p-values, signatures × pathways.


#### `network_matrix(net: pd.DataFrame, features: np.ndarray, min_size: int = 5) -> tuple[np.ndarray, list[str], np.ndarray]`

Features × pathways weight matrix from a long ``source, target[, weight]`` table.

#### `ulm(Y: pd.DataFrame, net: pd.DataFrame, min_size: int = 5, source: str = '') -> PathwayScores`

ULM scores for each row of ``Y`` (signatures × features) and each gene set.

#### `dedup_gene_sets(net: pd.DataFrame, max_jaccard: float = 0.5) -> tuple[pd.DataFrame, list[tuple[str, str, float]]]`

Drop the smaller of any pair of sets with Jaccard > ``max_jaccard``.

Returns the filtered table and the list of (dropped, kept, jaccard) decisions.

#### `symbols_for(features: np.ndarray, symbols: np.ndarray | None) -> np.ndarray`

Map feature ids to gene symbols (first occurrence wins; unmapped keep their id).

### `tincture.stats.map_retrieval`

Replicate retrievability: mean average precision (mAP) with a permutation null.

#### `ap_from_relevance(rel: np.ndarray) -> np.ndarray`

AP for each row of a 0/1 relevance matrix already sorted by decreasing score.

#### `random_ap(n_draws: int, n_pos: int, total: int, rng: np.random.Generator) -> np.ndarray`

Null APs: random rankings with ``n_pos`` positives among ``total`` items.

#### `average_precision(X: np.ndarray, groups: Sequence[str] | np.ndarray, is_control: np.ndarray, negatives: str = 'controls') -> pd.DataFrame`

Per-profile AP (one row per non-control profile that has ≥ 1 positive).

#### class `MapResult(table: pd.DataFrame, ap: pd.DataFrame, null_size: int, alpha: float, negatives: str, excluded: list[str], null: str = 'copied', min_p: np.ndarray | None = None, design_can_pass: bool = True)`

Per-group retrievability.


#### `permutation_null(S: np.ndarray, group_idx: np.ndarray, control_idx: np.ndarray, null_size: int, rng: np.random.Generator) -> tuple[np.ndarray, bool]`

Null mAPs by relabelling group members among group ∪ controls.

Returns (null values, exact): ``exact`` means all C(pool, n_rep) relabelings were
enumerated (including the observed labelling).

#### `min_achievable_p(n_rep: int, n_controls: int, null: str, null_size: int) -> float`

Smallest p-value the design can produce (perfect retrieval).

#### `mean_average_precision(ap: pd.DataFrame, null_size: int = 10000, alpha: float = 0.05, seed: int = 0, all_groups: Sequence[str] | None = None, negatives: str = 'controls') -> MapResult`

mAP per group with the copied permutation null and BH q-values.

#### `retrievability(X: np.ndarray, groups: Sequence[str] | np.ndarray, is_control: np.ndarray, negatives: str = 'controls', null_size: int = 10000, alpha: float = 0.05, seed: int = 0, null: str = 'permutation') -> MapResult`

AP per profile + mAP, p, q per group in one call (see module docstring for nulls).

#### `cat_curve(a: np.ndarray, b: np.ndarray, ks: Sequence[int] | None = None, level: float = 0.95) -> CatCurve`

Correspondence at the top: fraction of the top-k features (by |value|) shared
between two replicate profiles, with an exact hypergeometric null band (overlap of
two random k-subsets of n features).

### `tincture.stats.zprime`

Plate quality: Z′, edge effect and control drift.

#### `zprime(pos: np.ndarray, neg: np.ndarray, robust: bool = False) -> float`

Z′ factor; NaN when either control group has fewer than 2 finite values.

#### class `Selector(well_types: Sequence[str] = (), tags: Sequence[str] = ())`

Selects wells by ``well_type`` values and/or source ``tag`` values (OR).

- `mask(d: pd.DataFrame) -> np.ndarray`
- `describe() -> str`

#### `plate_qc(pt: PlateTable, neg: Selector | None = None, pos: Selector | None = None, edge: Selector | str = 'auto', min_controls: int = 3, threshold: float = 0.5, seeded_frac: float = 0.1) -> PlateQC`

Per-plate Z′, edge ratio and control level, plus drift across plate order.

``edge="auto"`` uses negative controls when every plate has ≥ ``min_controls`` of them
both on the edge and inside; otherwise all drug-free cell wells (neg_ctrl + empty).

### `tincture.stats.projection`

Projection of signatures onto reference mechanisms (Compass statistics).

#### `nnls_project(y: np.ndarray, R: np.ndarray, se: np.ndarray | None = None) -> Projection`

Non-negative least squares of ``y`` on the rows of ``R`` (refs × features).

With ``se`` (per-feature standard errors of ``y``), magnitude and R² are corrected for
noise (see module docstring).

#### `spoke_angles(n: int, start: float = 1.5707963267948966) -> np.ndarray`

Evenly spaced spoke angles, the first pointing up, going clockwise.

#### `compass_xy(p: Projection, angles: np.ndarray) -> tuple[float, float, float]`

(x, y, concentration) for a projection: direction from weights, radius = magnitude.

#### `bootstrap_ellipse(effect: np.ndarray, sd: np.ndarray, R: np.ndarray, angles: np.ndarray, n_boot: int = 200, level: float = 0.95, seed: int = 0) -> tuple[float, float, float, float, float]`

95% ellipse (cx, cy, width, height, angle_deg) of the compass position under
posterior draws ``effect + sd · z``.

#### `cov_ellipse(pts: np.ndarray, level: float = 0.95) -> tuple[float, float, float, float, float]`

Normal-theory ellipse (cx, cy, width, height, angle_deg) covering ``level`` of 2-D
points (n × 2).

#### `seriate(sim: np.ndarray) -> np.ndarray`

Order references so similar ones are adjacent (optimal leaf ordering).

### `tincture.stats.paths`

Dose-path statistics: 2-D positions per dose and the bend test.

#### `stress_gene_set(hallmark_table, sets: Sequence[str] = ('HALLMARK_APOPTOSIS', 'HALLMARK_P53_PATHWAY', 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE')) -> set[str]`

Union of the stress Hallmark sets (``source``/``target`` table, names with or
without the ``HALLMARK_`` prefix).

#### class `PathAxes(ref_unit: np.ndarray, stress_mask: np.ndarray | None, scale: tuple[float, float] = (1.0, 1.0))`

How to turn an effect vector into (x, y).

- `x(eff: np.ndarray) -> np.ndarray`
- `y_weights`: y as a linear map: y = eff @ w (stress mean minus overall mean).
- `y(eff: np.ndarray) -> np.ndarray | float`

#### `signed_turns(xy: np.ndarray) -> np.ndarray`

Signed turning angles at interior points of a path that starts at the origin.

``xy``: (..., n_points, 2) for the doses; the origin (vehicle) is prepended.
Returns (..., n_points - 1).

#### `dose_path(sigs: Sequence[Signature], axes: PathAxes, y_values: np.ndarray | None = None, n_boot: int = 500, seed: int = 0) -> DosePath`

Positions and bend test for one compound's signatures (any order; sorted by dose).

``y_values`` overrides the y axis (cell-loss axis): one value per signature, no draws on
y (cell counts carry their own, smaller uncertainty, not modelled here).

## Theme

### `tincture.theme.presets`

Journal export presets.

#### class `JournalPreset(name: str, width_mm: float, max_height_mm: float, min_font_pt: float = 6.0, base_font_pt: float = 7.0, dpi: int = 600, note: str = '')`

Physical figure size and typography for one journal format.

- `width_in`
- `max_height_in`
- `figsize(height_mm: float | None = None, aspect: float = 0.62) -> tuple[float, float]`: Figure size in inches.

#### `get_preset(preset: str | JournalPreset) -> JournalPreset`

Look up a preset by name (or pass one through).

### `tincture.theme.colorblind`

Colour-vision-deficiency check for figures.

#### `simulate(rgb: np.ndarray, kind: str) -> np.ndarray`

Simulate a colour-vision deficiency on sRGB colours in [0, 1], shape (..., 3).

#### `srgb_to_lab(rgb: np.ndarray) -> np.ndarray`

sRGB in [0, 1] -> CIE L*a*b* (D65).

#### `delta_e2000(lab1: np.ndarray, lab2: np.ndarray) -> float`

CIEDE2000 colour difference (Sharma, Wu & Dalal 2005).

#### class `ColorblindReport(colours: list[str], threshold: float, failures: list[tuple[str, str, str, float]] = <factory>, min_delta: dict[str, float] = <factory>)`

Result of :func:`check_colorblind`.

- `passed`

#### `check_colours(colours: list[str], threshold: float = 10.0, kinds: tuple[str, ...] = ('normal', 'deuteranopia', 'protanopia')) -> ColorblindReport`

Check a list of colours pairwise under each simulated vision type.

#### `check_colorblind(fig: Figure, threshold: float = 10.0, raise_on_fail: bool = True) -> ColorblindReport`

Check that a figure's categorical colours stay distinguishable for colour-blind readers.

Parameters
----------
fig
    A matplotlib figure (tincture-rendered or not).
threshold
    Minimum CIEDE2000 difference between any two categorical colours.
raise_on_fail
    Raise :class:`TinctureColorblindError` if any pair collides (default), otherwise
    return the report.

## Errors

### `tincture.errors`

Errors and warnings raised by tincture.

#### class `TinctureError`

Base class for tincture errors.


#### class `TinctureDataError(message: str, fix: str | None = None)`

Input data cannot support the requested plot or statistic.

Parameters
----------
message
    What is wrong with the data.
fix
    A concrete suggestion for how to fix it (always shown to the user).


#### class `TinctureSpecError`

A PlotSpec is internally inconsistent (a bug in a plot builder or user spec).


#### class `TinctureStyleError`

A styling request violates the theme rules (e.g. a rainbow colormap).


#### class `TinctureColorblindError`

Raised by ``check_colorblind(..., raise_on_fail=True)`` when colours collide.


#### class `TinctureWarning`

A result is weak, reduced or degraded; the plot still renders but says so.
