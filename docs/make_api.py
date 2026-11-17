"""Generate docs/api.md from the package's public docstrings (no extra dependencies).

    python docs/make_api.py            # rewrite docs/api.md
    python docs/make_api.py --check    # exit 1 if docs/api.md is stale (used by the tests)

Public = listed in ``tincture.__all__`` for the top level; in submodules, every function or
class defined in that module whose name has no leading underscore.
"""

from __future__ import annotations

import importlib
import inspect
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent / "api.md"

SECTIONS: list[tuple[str, str, list[str]]] = [
    (
        "Plots",
        "Every plot returns a `PlotResult` (figure, caption, methods, warnings, classic "
        "fallback). Call them as `tc.<name>(...)`.",
        [
            "tincture.plots.retrievability",
            "tincture.plots.plate_atlas",
            "tincture.plots.compass",
            "tincture.plots.dose_paths",
            "tincture.plots.fingerprint",
            "tincture.plots.effect_ladder",
        ],
    ),
    ("Results and export", "", ["tincture.result", "tincture.export"]),
    ("Data loaders", "", ["tincture.io.sciplex", "tincture.io.prism", "tincture.io.gdsc"]),
    ("Datasets and caches", "", ["tincture.datasets.cache", "tincture.datasets.genesets"]),
    (
        "Data model",
        "",
        [
            "tincture.model.perturbation",
            "tincture.model.signature",
            "tincture.model.reference",
            "tincture.model.dose_response",
            "tincture.model.plate",
            "tincture.model.measurements",
        ],
    ),
    (
        "Statistics",
        "",
        [
            "tincture.stats.signatures",
            "tincture.stats.ebayes",
            "tincture.stats.pathways",
            "tincture.stats.map_retrieval",
            "tincture.stats.zprime",
            "tincture.stats.projection",
            "tincture.stats.paths",
        ],
    ),
    ("Theme", "", ["tincture.theme.presets", "tincture.theme.colorblind"]),
    ("Errors", "", ["tincture.errors"]),
]

#: Spec builders and plumbing used by the plot functions; not part of the user API.
INTERNAL = {
    "build_spec",
    "build_classic_spec",
    "sidecar_markdown",
    "compute_points",
    "place_labels",
}

TOP_LEVEL = {
    "tincture.plots.retrievability": "retrievability",
    "tincture.plots.plate_atlas": "plate_atlas",
    "tincture.plots.compass": "compass",
    "tincture.plots.dose_paths": "dose_paths",
    "tincture.plots.fingerprint": "fingerprint",
    "tincture.plots.effect_ladder": "effect_ladder",
}


def _first_paragraph(doc: str | None) -> str:
    if not doc:
        return "_No docstring._"
    doc = inspect.cleandoc(doc)
    return doc


def _signature(obj: object) -> str:
    """Signature with annotations as written in the source (they are strings under
    ``from __future__ import annotations``) and defaults as Python literals."""
    try:
        sig = inspect.signature(obj)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return "(...)"
    parts, star_done = [], False
    for p in sig.parameters.values():
        if p.name == "self" or p.name.startswith("_"):
            continue
        kind = p.kind
        if kind is p.KEYWORD_ONLY and not star_done:
            parts.append("*")
            star_done = True
        prefix = {p.VAR_POSITIONAL: "*", p.VAR_KEYWORD: "**"}.get(kind, "")
        if kind is p.VAR_POSITIONAL:
            star_done = True
        txt = prefix + p.name
        if p.annotation is not p.empty:
            txt += f": {p.annotation if isinstance(p.annotation, str) else p.annotation.__name__}"
        if p.default is not p.empty:
            txt += f" = {p.default!r}"
        parts.append(txt)
    out = f"({', '.join(parts)})"
    ret = sig.return_annotation
    if ret is not sig.empty:
        out += f" -> {ret if isinstance(ret, str) else getattr(ret, '__name__', ret)}"
    return out


def _members(mod: object) -> list[tuple[str, object]]:
    out = []
    for name, obj in vars(mod).items():
        if name.startswith("_"):
            continue
        if not (inspect.isfunction(obj) or inspect.isclass(obj)):
            continue
        if getattr(obj, "__module__", None) != mod.__name__:  # type: ignore[attr-defined]
            continue
        if name in INTERNAL or not _own_doc(obj):  # spec builders / undocumented helpers
            continue
        out.append((name, obj))
    # the plot function first, then the rest in source order
    return sorted(out, key=lambda kv: kv[0] != "plot")


def _own_doc(obj: object) -> str | None:
    """The object's docstring, ignoring dataclasses' auto-generated ``Name(field: ...)``."""
    d = inspect.getdoc(obj)
    if d and inspect.isclass(obj) and d.startswith(f"{obj.__name__}("):
        return None
    return d


def _module_block(modname: str) -> list[str]:
    mod = importlib.import_module(modname)
    lines = [f"### `{modname}`", ""]
    doc = inspect.getdoc(mod)
    if doc:
        lines += [doc.split("\n\n")[0], ""]
    for name, obj in _members(mod):
        if modname in TOP_LEVEL and name == "plot":
            title = f"`tc.{TOP_LEVEL[modname]}{_signature(obj)}`"
        elif inspect.isclass(obj):
            # constructor signature only where the class defines one (stable across Pythons)
            own_init = "__init__" in vars(obj) or hasattr(obj, "__dataclass_fields__")
            sig = _signature(obj).split(" -> ")[0] if own_init else ""
            title = f"class `{name}{sig}`"
        else:
            title = f"`{name}{_signature(obj)}`"
        lines += [f"#### {title}", "", _first_paragraph(_own_doc(obj)), ""]
        if inspect.isclass(obj):
            for mname, m in vars(obj).items():
                if mname.startswith("_") or not (
                    inspect.isfunction(m) or isinstance(m, property | classmethod | staticmethod)
                ):
                    continue
                target = m.fget if isinstance(m, property) else getattr(m, "__func__", m)
                sig = "" if isinstance(m, property) else _signature(target)
                d = inspect.getdoc(target)
                summary = d.split("\n\n")[0].replace("\n", " ") if d else ""
                lines.append(f"- `{mname}{sig}`" + (f": {summary}" if summary else ""))
            lines.append("")
    return lines


def render() -> str:
    import tincture

    lines = [
        "# API reference",
        "",
        f"_Generated from docstrings by `docs/make_api.py` (tincture {tincture.__version__})._ "
        "Do not edit by hand: run `python docs/make_api.py`.",
        "",
        "```python",
        "import tincture as tc",
        "```",
        "",
        "Top-level names: " + ", ".join(f"`tc.{n}`" for n in tincture.__all__) + ".",
        "",
    ]
    for title, intro, mods in SECTIONS:
        lines += [f"## {title}", ""]
        if intro:
            lines += [intro, ""]
        for m in mods:
            lines += _module_block(m)
    return "\n".join(lines).rstrip() + "\n"


def main(argv: list[str]) -> int:
    text = render()
    if "--check" in argv:
        if not OUT.exists() or OUT.read_text(encoding="utf-8") != text:
            print("docs/api.md is stale: run python docs/make_api.py", file=sys.stderr)
            return 1
        return 0
    OUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUT} ({len(text.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
