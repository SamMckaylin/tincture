"""Gene-set libraries (MSigDB Hallmark, PROGENy) fetched from OmniPath and cached.

Hallmark is read directly from OmniPath's static table (the same file and the same
processing as ``decoupler.op.hallmark``), so the core install needs no ``[omics]`` extra.
PROGENy uses decoupler (``pip install 'tincture[omics]'``).

OmniPath tables can change between fetches, so the first fetch is cached under
``<data_home>/genesets/`` with a JSON sidecar recording the SHA-256, fetch date, source
and decoupler version. Later loads use the cache and verify it; the sidecar's version
string is printed on every figure that uses the gene sets.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from tincture.datasets.cache import DatasetUnavailable, data_home
from tincture.errors import TinctureDataError

HALLMARK_URL = "https://static.omnipathdb.org/tables/msigdb-hallmark.tsv.gz"

SOURCES = {
    "hallmark": "MSigDB Hallmark (Liberzon et al., Cell Syst 2015) via OmniPath",
    "progeny": "PROGENy top-500 (Schubert et al., Nat Commun 2018) via OmniPath/decoupler",
}


@dataclass(frozen=True)
class GeneSetLibrary:
    name: str
    table: pd.DataFrame  # source, target[, weight]
    version: str

    @property
    def n_sets(self) -> int:
        return int(self.table["source"].nunique())


def _paths(name: str, root: Path | None) -> tuple[Path, Path]:
    base = (root or data_home()) / "genesets"
    return base / f"{name}_human.tsv", base / f"{name}_human.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch_hallmark(url: str = HALLMARK_URL, timeout: float = 120.0) -> pd.DataFrame:
    """Hallmark gene sets from OmniPath's static table (no decoupler needed).

    Same processing as ``decoupler.op.hallmark(organism="human")``: drop the ``HALLMARK_``
    prefix, split complexes (``COMPLEX:A_B``) into member genes, de-duplicate.
    """
    import gzip
    import io
    import urllib.request

    req = urllib.request.Request(url, headers={"User-Agent": "tincture-data-fetch"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read()
    hm = pd.read_csv(io.BytesIO(gzip.decompress(raw)), sep="\t")[["geneset", "genesymbol"]]
    hm["geneset"] = hm["geneset"].str.replace("HALLMARK_", "", regex=False)
    hm["genesymbol"] = (
        hm["genesymbol"].astype(str).str.replace("COMPLEX:", "", regex=False).str.split("_")
    )
    hm = hm.explode("genesymbol").rename(columns={"geneset": "source", "genesymbol": "target"})
    return hm.drop_duplicates(["source", "target"]).reset_index(drop=True)


def load_geneset(name: str, root: Path | None = None, download: bool = True) -> GeneSetLibrary:
    """Load (fetching once if needed) a human gene-set library: ``hallmark`` or ``progeny``."""
    if name not in SOURCES:
        raise TinctureDataError(
            f"Unknown gene-set library {name!r}.", fix=f"Use one of {list(SOURCES)}."
        )
    tsv, meta = _paths(name, root)
    if not tsv.exists():
        if not download:
            raise DatasetUnavailable(f"Gene-set library {name!r} not cached at {tsv}.")
        fetched_with = "urllib (OmniPath static table)"
        try:
            if name == "hallmark":
                table = fetch_hallmark()
            else:
                try:
                    import decoupler as dc
                except ImportError as exc:
                    raise TinctureDataError(
                        "Fetching PROGENy needs decoupler.", fix="pip install 'tincture[omics]'"
                    ) from exc
                table = dc.op.progeny(organism="human", top=500)[["source", "target", "weight"]]
                fetched_with = f"decoupler {getattr(dc, '__version__', '?')}"
        except TinctureDataError:
            raise
        except Exception as exc:  # network / proxy failures surface here
            raise DatasetUnavailable(
                f"Could not fetch {name!r} from OmniPath ({type(exc).__name__}: {exc}). "
                "Allow static.omnipathdb.org and omnipathdb.org, or place a TSV with columns "
                f"source,target[,weight] at {tsv} and a sidecar JSON at {meta}."
            ) from exc
        tsv.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(tsv, sep="\t", index=False)
        meta.write_text(
            json.dumps(
                {
                    "name": name,
                    "source": SOURCES[name],
                    "fetched": dt.date.today().isoformat(),
                    "fetched_with": fetched_with,
                    "sha256": _sha(tsv),
                    "n_sets": int(table["source"].nunique()),
                    "n_rows": len(table),
                },
                indent=2,
            )
        )
    info = json.loads(meta.read_text()) if meta.exists() else {}
    if info.get("sha256") and info["sha256"] != _sha(tsv):
        raise DatasetUnavailable(
            f"{tsv} does not match the checksum recorded at fetch time; delete both files "
            "and re-fetch."
        )
    table = pd.read_csv(tsv, sep="\t")
    version = (
        f"{info.get('source', SOURCES[name])}; fetched {info.get('fetched', '?')}; "
        f"sha256 {str(info.get('sha256', '?'))[:12]}"
    )
    return GeneSetLibrary(name=name, table=table, version=version)
