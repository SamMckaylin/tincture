"""Captions, methods paragraphs and classic fallbacks."""

from tincture.report.captions import TEMPLATES, Caption, fmt_conc, fmt_p
from tincture.report.classic import build_classic, register_classic
from tincture.report.methods import Methods, register_block

__all__ = [
    "TEMPLATES",
    "Caption",
    "Methods",
    "build_classic",
    "fmt_conc",
    "fmt_p",
    "register_block",
    "register_classic",
]
