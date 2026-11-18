# API reference

_Generated from docstrings by `docs/make_api.py` (tincture 0.0.1.dev0)._ Do not edit by hand: run `python docs/make_api.py`.

```python
import tincture as tc
```

Top-level names: `tc.ComboMatrix`, `tc.PlotResult`, `tc.TinctureDataError`, `tc.TinctureError`, `tc.TinctureWarning`, `tc.__version__`, `tc.combination_lift`, `tc.combination_surface`, `tc.compass`, `tc.demo`, `tc.dose_paths`, `tc.effect_ladder`, `tc.engagement_ladder`, `tc.export`, `tc.fingerprint`, `tc.neighbourhood`, `tc.onset_timeline`, `tc.plate_atlas`, `tc.render_mpl`, `tc.retrievability`, `tc.selectivity`, `tc.theme`, `tc.window`.

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

### `tincture.plots.window`

Window: is there a dose that works without the toxicity? (HANDOVER §5.7)

#### `tc.window(on: DoseResponse, tox: DoseResponse | Sequence[DoseResponse], cmax: tuple[float, float] | None = None, on_frac: float = 0.9, tox_drop: float = 0.1, n_boot: int = 200, compound: str | None = None, title: str = 'Window', footer: Sequence[str] = ()) -> PlotResult`

Window plot. ``cmax``: free Cmax range in molar (low, high), drawn as a band; omit if
no exposure data are available (never pass a guess).

#### class `CurveCall(dr: DoseResponse, hill: HillFit, ec: float, ec_ci: tuple[float, float], problems: list[str] = <factory>, censored: bool = False)`

Fit, effective concentration and verdict for one curve.

- `usable`

### `tincture.plots.neighbourhood`

Neighbourhood Card: which annotated compounds does this one resemble? (HANDOVER §5.9)

#### `tc.neighbourhood(query: Signature, library: ReferenceLibrary, smiles: str | None = None, k: int = 10, alpha: float = 0.05, query_replicates: np.ndarray | None = None, name: str | None = None, null: np.ndarray | None = None, null_kind: str = 'all', title: str = 'Neighbourhood Card', footer: Sequence[str] = ()) -> PlotResult`

Neighbourhood Card for ``query`` against an annotated ``library``.

``query_replicates``: (n, features) signatures of the same compound (e.g. other cell
lines) for bootstrap intervals; otherwise the query's SEs are used if present.
``null``: precomputed library-pair similarities (reuse across many queries).
``null_kind``: ``"all"`` library pairs (HANDOVER default) or ``"unrelated"`` pairs (pairs
sharing an MoA annotation removed; see :mod:`tincture.stats.neighbours`).

#### `structure_image(smiles: str | None, size: int = 360) -> np.ndarray | None`

RGBA image of the molecule, or None if no SMILES / RDKit unavailable / unparsable.

### `tincture.plots.selectivity`

Selectivity Spectrum: how specific is this compound? (HANDOVER §5.8)

#### `tc.selectivity(table: PotencyTable, compound: str, intended: str | Sequence[str], cmax: tuple[float, float] | None = None, allow_mixed: bool = False, panel_size: int | None = None, panel_name: str = 'targets', top_labels: int = 8, title: str = 'Selectivity Spectrum', footer: Sequence[str] = ()) -> PlotResult`

Selectivity Spectrum for one compound. ``cmax``: free Cmax range (molar); never pass
a guess. ``panel_size``: size of the full family (e.g. ~518 human kinases) for the
coverage footer.

### `tincture.plots.engagement`

Engagement Ladder: does potency hold up from protein to cell to phenotype, and where does
it break? (HANDOVER §5.10)

#### `tc.engagement_ladder(table: PotencyTable, pairs: Sequence[tuple[str, str]], km_atp: Mapping[str, KmATP] | None = None, atp: tuple[float, float] = (0.001, 0.005), title: str = 'Engagement Ladder', footer: Sequence[str] = ()) -> PlotResult`

Engagement Ladder for (compound, target) ``pairs``.

``table`` needs a ``rung`` column with values from :data:`RUNGS`; optional columns
``px_lo``/``px_hi`` (interval), ``cell_system`` and ``label``. ``km_atp``: ATP Km
(molar) per target, one value or a (low, high) range, for the Cheng–Prusoff wedge;
never pass a guess.

#### `atp_shift(km_atp: KmATP, atp: tuple[float, float] = (0.001, 0.005)) -> tuple[float, float]`

Expected Cheng–Prusoff shift (log10) for ATP-competitive inhibitors.

``km_atp`` is one Km or a (low, high) range (e.g. unphosphorylated vs phosphorylated
enzyme); the wedge then spans every combination, log10(1 + ATP_low / Km_high) to
log10(1 + ATP_high / Km_low), so the shift flag fires only beyond the widest expectation.

### `tincture.plots.combination`

Combination Surface and Combination Lift (HANDOVER §5.11).

#### `tc.combination_surface(cm: ComboMatrix, model: Literal['bliss', 'zip', 'hsa'] = 'bliss', single_agent: Literal['fitted', 'observed'] = 'fitted', noise_sd: float | Literal['matrix'] | None = None, cmax_a: tuple[float, float] | None = None, cmax_b: tuple[float, float] | None = None, min_excess: float = 0.1, alpha: float = 0.05, n_boot: int = 200, title: str = 'Combination Surface', footer: Sequence[str] = ()) -> PlotResult`

Combination Surface for one full concentration matrix. ``cmax_a``/``cmax_b``: free
Cmax ranges (molar) of each drug; never pass a guess. ``noise_sd``: well-level SD of
viability, used only when there is a single replicate matrix; ``"matrix"`` estimates it
from the matrix itself (``stats.synergy.matrix_noise_sd``).

#### `lift_table(d: pd.DataFrame, group: str | None = None, facet: str | None = None) -> pd.DataFrame`

Per cell line (and facet): mean viabilities over replicates, ΔEmax with a t-interval
when there are ≥ 2 replicates. ``fit_ok`` (optional, bool) marks usable fits.

#### `tc.combination_lift(table: pd.DataFrame, drug_a: str, drug_b: str, group: str | None = None, facet: str | None = None, threshold: float = 0.2, model_note: str = '', title: str = 'Combination Lift', footer: Sequence[str] = ()) -> PlotResult`

Combination Lift across cell lines. ``table``: one row per cell line × replicate (×
facet, e.g. anchor concentration) with viabilities ``single_a``, ``single_b``, ``combo``
at Emax and optional ``fit_ok``, group and facet columns.

### `tincture.plots.onset`

Onset Timeline: what happens first, and what follows? (HANDOVER §5.12)

#### `tc.onset_timeline(table: pd.DataFrame, feature: str = 'feature', threshold: float = 0.0, alpha: float = 0.05, time_unit: str = 'h', score_label: str = 'effect score', style: Literal['gantt', 'horizon'] = 'gantt', title: str = 'Onset Timeline', footer: Sequence[str] = ()) -> PlotResult`

Onset Timeline from a long table (``feature``, ``time``, ``score``, ``q`` or ``p``).

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

### `tincture.io.lincs`

Loader for LINCS L1000 phase II, Level 5 signatures (GEO GSE70138).

#### `parse_dose(s: object) -> float`

``"10.0 um"`` → 1e-5 (molar); unknown or control doses (``-666``) → NaN.

#### `build_landmark_cache(root: Path | None = None, out: Path | None = None, download: bool = True, block: int = 4000) -> Path`

Download Level 5, keep the 978 landmark genes, write the h5ad cache.

#### `load_lincs(cell_lines: list[str] | None = None, pert_type: str | None = 'trt_cp', time: float | None = 24.0, min_tas: float | None = None, cache: Path | None = None, build: bool = False) -> ad.AnnData`

LINCS GSE70138 Level 5 landmark signatures (signatures × 978 genes) with metadata
(``pert_id``, ``pert_iname``, ``cell_id``, ``dose`` in molar, ``time`` in hours, ``tas``,
``canonical_smiles``, ``inchi_key``).

### `tincture.io.davis`

Loader for Davis et al. 2011 kinase binding constants (Nat Biotechnol 29:1046; D8).

#### `classify(construct: str) -> str`

``"mutant"``, ``"nonhuman"`` or ``"wildtype"`` for a Davis assay construct name.

#### `load_davis2011(compounds: list[str] | None = None, per_gene: bool = True, root: Path | None = None) -> PotencyTable`

Davis 2011 Kd values as a :class:`PotencyTable` (assay type ``Kd``).

### `tincture.io.nanobret`

Loader for NanoBRET live-cell kinase occupancy (Vasta et al., Cell Chem Biol 2018; D7).

#### `read_table_s1(root: Path | None = None) -> pd.DataFrame`

Table S1 as a long table: kinase, compound, replicate, occupancy (%, NaN = not
determined), no_binding flag.

#### `occupancy(root: Path | None = None) -> pd.DataFrame`

Mean occupancy (fraction) per compound × kinase over determined replicates (kinases
listed twice keep the row with determined values).

#### `apparent_ic50(f: np.ndarray, dose: float) -> np.ndarray`

C·(1−f)/f with f clipped to [F_MIN, F_MAX] (1:1 binding, Hill 1).

#### `load_nanobret(root: Path | None = None, davis_names: bool = True) -> PotencyTable`

Apparent cellular IC50s as a :class:`PotencyTable` (assay type ``IC50_cell``).

``davis_names``: rename kinases to Davis 2011 gene symbols (curated synonyms in
:mod:`tincture.io.kinase_names`) so the two panels join; the original name is kept in
``nanobret_name``.

### `tincture.io.kinase_names`

Curated kinase name synonyms between panels (HANDOVER §6.2: map by curated synonyms
when no shared identifier exists; unmatched names are listed in docs/data_status.md).

### `tincture.io.gdsc_combo`

GDSC combination screens: full 7×7 matrices (D5, Bashi et al. 2024) and anchored 2×7
screens (D4, Jaaks et al. 2022; D6, Vis et al. 2024). File facts in docs/data_status.md.

#### `bashi_results() -> pd.DataFrame`

All per-screen matrix summaries, with a ``screen`` column (e.g. 'gdsc-008').

#### `load_bashi_matrices(drug_1: str, drug_2: str, cell_lines: Sequence[str] | None = None) -> dict[str, ComboMatrix]`

7×7 matrices for one drug pair (orientation as screened: drug A = lib1).

#### `load_bashi_replicated(screens: Sequence[str] = ('gdsc-008', 'gdsc-010'), min_replicates: int = 3) -> list[ComboMatrix]`

Every cell line × drug pair screened as ≥ ``min_replicates`` replicate 7×7 matrices
(used to calibrate single-replicate noise estimates).

#### `load_jaaks2022(screen: str = 'original') -> pd.DataFrame`

Jaaks et al. 2022 fitted anchored data (one row per anchor conc × library × line ×
replicate plate).

#### `load_vis2024() -> pd.DataFrame`

Vis et al. 2024 pan-cancer anchored fitted data (Data S1).

#### `cell_model_annotations() -> pd.DataFrame`

Cell Model Passports model list (model_id = SIDM, msi_status, tissue …).

#### `anchored_lift_table(d: pd.DataFrame, anchor: str, library: str, max_rmse: float | None = None) -> pd.DataFrame`

Lift input from a GDSC anchored table for one anchor × library pair: viabilities at
the highest library concentration (anchor alone, library Emax, combination Emax) per
cell line × anchor concentration × replicate. ``fit_ok`` is False when any value is
missing or (if given) a curve RMSE exceeds ``max_rmse``.

### `tincture.io.mixseq`

MIX-seq (McFarland, Paolella et al., Nat Commun 2020; figshare 10298696 v3). Facts in
docs/data_status.md (D10).

#### `pseudobulk_trametinib_tc(min_cells: int = 10) -> tuple[pd.DataFrame, pd.DataFrame]`

Returns ``(counts, groups)``: counts = genes (symbol index, Ensembl in column
``ensembl``… dropped) × groups; groups = one row per column with ``depmap_id``,
``condition`` (DMSO / trametinib / untreated), ``time_h`` and ``n_cells``. Groups with
fewer than ``min_cells`` cells are dropped.

## Datasets and caches

### `tincture.datasets.cache`

Download-and-cache helpers with checksums.

#### class `DatasetUnavailable`

A registered file is not cached and could not be downloaded.


#### class `RemoteFile(dataset: str, name: str, url: str, md5: str | None = None, sha256: str | None = None, size: int | None = None, licence: str = '', sha512: str | None = None, manual: bool = False)`

One downloadable file.

Attributes
----------
dataset
    Cache folder name (e.g. ``"sciplex"``).
name
    File name in the cache.
url
    Source URL.
md5, sha256, sha512
    Expected checksum(s); at least one must be given (GEO publishes SHA-512).
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

Gene-set libraries (MSigDB Hallmark, PROGENy) fetched from OmniPath and cached.

#### `fetch_hallmark(url: str = 'https://static.omnipathdb.org/tables/msigdb-hallmark.tsv.gz', timeout: float = 120.0) -> pd.DataFrame`

Hallmark gene sets from OmniPath's static table (no decoupler needed).

Same processing as ``decoupler.op.hallmark(organism="human")``: drop the ``HALLMARK_``
prefix, split complexes (``COMPLEX:A_B``) into member genes, de-duplicate.

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

#### `library_from_compounds(signatures: Mapping[str, np.ndarray], features: np.ndarray, annotations: Mapping[str, Mapping[str, str]], version: str, source: str = '', citation: str = '') -> ReferenceLibrary`

One reference per compound (each its own single member), with annotations such as
MoA. Leave-one-out (:meth:`ReferenceLibrary.without`) then removes a query compound from
its own library, as the Neighbourhood Card requires.

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

### `tincture.model.potency`

PotencyTable: compound × target (or rung) potency with assay type and context.

#### class `PotencyTable(data: pd.DataFrame, source: str = '')`

Long table with columns ``compound_id, target, px, assay_type`` (+ optional
``censored`` ('', '>' = weaker than the limit, '<' = more potent than the limit),
``context`` and any metadata).

- `from_molar(cls, df: pd.DataFrame, value_col: str = 'value_molar', source: str = '') -> PotencyTable`
- `compound(compound_id: str) -> pd.DataFrame`
- `assay_types`

### `tincture.model.combo`

ComboMatrix: drug A × drug B concentration grid of viability (HANDOVER §4.2).

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

### `tincture.stats.hill`

4-parameter log-logistic (Hill) fits with bootstrap uncertainty.

#### `fit_4pl(dose: np.ndarray, response: np.ndarray, start: FourPL | None = None) -> HillFit`

Fit a 4PL to treated points (dose > 0, molar). ``start`` warm-starts a single
optimisation (used by the bootstrap); otherwise six starting points are tried.

#### `ec_relative(fit: FourPL, frac: float) -> float`

Concentration giving ``frac`` of the curve's span (EC90: frac = 0.9).

#### `ec_drop(fit: FourPL, drop: float) -> float`

Concentration at which the response is ``drop`` below the top (inf if never).

#### `median_of_units(dose: np.ndarray, response: np.ndarray, cluster: np.ndarray) -> tuple[np.ndarray, np.ndarray]`

Per-dose median over units of each unit's mean response (dose > 0).

#### `bootstrap_4pl(dose: np.ndarray, response: np.ndarray, n_boot: int = 200, seed: int = 0, cluster: np.ndarray | None = None) -> HillBootstrap`

Refit ``n_boot`` times on points resampled with replacement within each dose, or,
with ``cluster``, on the median-of-units curve of resampled units.

### `tincture.stats.neighbours`

Nearest annotated neighbours of a query signature, with an empirical null.

#### `library_null(lib: ReferenceLibrary, kind: str = 'all', key: str = 'moa', max_pairs: int = 2000000, seed: int = 0) -> np.ndarray`

Cosine similarities of distinct library pairs (subsampled above ``max_pairs``).

``kind="unrelated"`` keeps only pairs whose ``annotations[key]`` differ (pairs with a
missing annotation are kept).

#### `neighbours(query: np.ndarray, lib: ReferenceLibrary, query_replicates: np.ndarray | None = None, query_se: np.ndarray | None = None, null: np.ndarray | None = None, n_boot: int = 200, seed: int = 0, null_kind: str = 'all') -> NeighbourResult`

Rank library references by cosine similarity to ``query`` (aligned to lib.features).

### `tincture.stats.synergy`

Combination effect models: Bliss (default), HSA and ZIP, with honest uncertainty.

#### `grid_noise_sd(grid: np.ndarray, lo: float = 0.1, hi: float = 0.9, min_triplets: int = 10) -> float | None`

Well-level SD from one concentration grid: second differences y₋₁ − 2y₀ + y₊₁ along
rows and columns (removes linear dose trends), over triplets whose three viabilities all
lie in (lo, hi) (saturated wells understate noise); robust scale
σ = 1.4826 · median|Δ²| / √6. None if fewer than ``min_triplets`` triplets qualify.

#### `matrix_noise_sd(cm: ComboMatrix) -> float | None`

Median over replicate matrices of :func:`grid_noise_sd` (single-well SD).

#### `synergy(cm: ComboMatrix, model: Model = 'bliss', single_agent: Literal['fitted', 'observed'] = 'fitted', noise_sd: float | Literal['matrix'] | None = None, n_boot: int = 200, alpha: float = 0.05, min_excess: float = 0.1, ceiling: float = 0.9, seed: int = 0) -> SynergyResult`

Excess over a reference model on a full concentration matrix.

### `tincture.stats.onset`

Onset, peak and resolution of a response over sparse time points (HANDOVER §5.12).

#### `onset_table(d: pd.DataFrame, feature: str = 'feature', threshold: float = 0.0, alpha: float = 0.05) -> tuple[list[Onset], pd.DataFrame]`

Onsets for every feature. ``d`` has columns ``feature``, ``time``, ``score`` and ``q``
(or ``p``).

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
