"""Tiny helper for case-study reports (figures + caption + methods + checks in markdown)."""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import tincture as tc

OUT = Path(__file__).resolve().parent / "out"


class Report:
    def __init__(self, name: str, title: str):
        self.dir = OUT / name
        self.dir.mkdir(parents=True, exist_ok=True)
        self.lines = [
            f"# {title}",
            "",
            f"_Generated {dt.date.today().isoformat()} with tincture {tc.__version__}._",
            "",
        ]

    def text(self, *lines: str) -> None:
        self.lines += [*lines, ""]

    def figure(self, name: str, result: tc.PlotResult, preset: str = "nature_double") -> None:
        written = tc.export(result, self.dir / name, preset=preset, formats=("png", "pdf"))
        self.lines += [
            f"![{name}]({written['png'].name})",
            "",
            f"**Caption.** {result.caption}",
            "",
        ]
        if result.classic is not None:
            self.lines += [
                f"Classic fallback: [{written['classic_png'].name}]({written['classic_png'].name})",
                "",
            ]
        self.lines += [f"<details><summary>Methods</summary>\n\n{result.methods}\n\n</details>", ""]

    def check(self, name: str, passed: bool, detail: str) -> bool:
        self.lines += [f"- {'✅' if passed else '❌'} **{name}**: {detail}"]
        return passed

    def write(self) -> Path:
        p = self.dir / "report.md"
        p.write_text("\n".join(self.lines) + "\n", encoding="utf-8")
        return p


def md_table(df, floatfmt: str = ".2f") -> str:
    """Markdown table without extra dependencies."""
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = [format(v, floatfmt) if isinstance(v, float) else str(v) for v in r.tolist()]
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)
