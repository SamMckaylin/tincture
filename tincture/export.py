"""Export figures at journal sizes, with caption/methods sidecars.

``export(result, "figs/fig2", preset="nature_double")`` writes::

    figs/fig2.pdf, figs/fig2.png        # the figure, exactly preset width
    figs/fig2.md                        # caption, methods, footer, warnings
    figs/fig2_classic.pdf, ..._classic.png  # classic fallback, when one exists

Before writing, the figure is checked: every text must be ≥ the preset's minimum font
size (6 pt by default) and categorical colours must pass the colour-blind check.
"""

from __future__ import annotations

from pathlib import Path

from matplotlib.figure import Figure

from tincture.errors import TinctureStyleError
from tincture.render.mpl import iter_texts, min_font_size
from tincture.result import PlotResult
from tincture.spec import PlotSpec
from tincture.theme.colorblind import check_colorblind
from tincture.theme.presets import MM_PER_INCH, JournalPreset, get_preset
from tincture.theme.style import theme_context

ALLOWED_FORMATS = ("pdf", "png", "svg", "tiff", "eps")


def check_figure(fig: Figure, preset: str | JournalPreset) -> None:
    """Raise if the figure breaks the accessibility rules for ``preset``."""
    p = get_preset(preset)
    smallest = min_font_size(fig)
    if smallest < p.min_font_pt - 1e-6:
        offenders = sorted(
            {
                f"{t.get_text()[:30]!r} ({float(t.get_fontsize()):.1f} pt)"
                for t in iter_texts(fig)
                if float(t.get_fontsize()) < p.min_font_pt - 1e-6
            }
        )
        raise TinctureStyleError(
            f"Text below the {p.min_font_pt:g} pt minimum for {p.name}: {', '.join(offenders[:5])}"
        )
    check_colorblind(fig, raise_on_fail=True)


def _save_fig(
    fig: Figure, stem: Path, formats: tuple[str, ...], p: JournalPreset
) -> dict[str, Path]:
    out: dict[str, Path] = {}
    with theme_context(p):
        for fmt in formats:
            if fmt not in ALLOWED_FORMATS:
                raise TinctureStyleError(
                    f"Unsupported format {fmt!r}; use one of {ALLOWED_FORMATS}."
                )
            path = stem.with_suffix(f".{fmt}")
            # bbox="standard": never crop — the exported width must equal the preset width.
            fig.savefig(path, format=fmt, dpi=p.dpi, bbox_inches=None)
            out[fmt] = path
    return out


def sidecar_markdown(spec: PlotSpec, preset: JournalPreset, classic_name: str | None) -> str:
    lines = [f"# {spec.title or spec.plot}", ""]
    lines += ["## Caption", "", spec.caption or "_(no caption)_", ""]
    if spec.warnings:
        lines += ["## Warnings", ""] + [f"- {w}" for w in spec.warnings] + [""]
    lines += ["## Methods", "", spec.methods or "_(no methods paragraph)_", ""]
    if spec.footer:
        lines += ["## Footer", ""] + [f"- {f}" for f in spec.footer] + [""]
    lines += [
        "## Export",
        "",
        f"- Preset: `{preset.name}` ({preset.width_mm:g} mm wide, {preset.dpi} dpi)",
        f"- Classic fallback: {classic_name or 'none'}",
        "",
    ]
    return "\n".join(lines)


def export(
    obj: PlotResult | PlotSpec | Figure,
    path: str | Path,
    preset: str | JournalPreset = "nature_double",
    height_mm: float | None = None,
    formats: tuple[str, ...] = ("pdf", "png"),
    classic: bool = True,
    check: bool = True,
) -> dict[str, Path]:
    """Export a plot at a journal preset.

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
    """
    p = get_preset(preset)
    stem = Path(path)
    if stem.suffix in {f".{f}" for f in ALLOWED_FORMATS} | {".md"}:
        stem = stem.with_suffix("")
    stem.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(obj, Figure):
        fig = obj
        w_in, h_in = fig.get_size_inches()
        height = height_mm if height_mm is not None else p.width_mm * h_in / w_in
        fig.set_size_inches(p.width_in, min(height, p.max_height_mm) / MM_PER_INCH)
        if check:
            check_figure(fig, p)
        return _save_fig(fig, stem, formats, p)

    result = obj if isinstance(obj, PlotResult) else PlotResult(spec=obj)
    fig = result.figure(p, height_mm)
    if check:
        check_figure(fig, p)
    written = _save_fig(fig, stem, formats, p)

    from tincture.report.classic import classic_name

    cname = classic_name(result.spec.plot) if result.classic is not None else None
    md_path = stem.with_suffix(".md")
    md_path.write_text(sidecar_markdown(result.spec, p, cname), encoding="utf-8")
    written["md"] = md_path

    if classic and result.classic is not None:
        cfig = result.classic_figure(p, height_mm)
        assert cfig is not None
        if check:
            check_figure(cfig, p)
        cstem = stem.with_name(stem.name + "_classic")
        for fmt, fp in _save_fig(cfig, cstem, formats, p).items():
            written[f"classic_{fmt}"] = fp
    return written
