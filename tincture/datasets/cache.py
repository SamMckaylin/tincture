"""Download-and-cache helpers with checksums.

Data never lives in the repository. Files are cached under ``$TINCTURE_DATA`` (default:
``<repo>/data/raw`` when run from a checkout, else ``~/.cache/tincture``), one folder per
dataset. Every registered file has a pinned checksum (MD5 when the source publishes
one, e.g. figshare; otherwise SHA-256 recorded at first verified download). A file
whose checksum does not match is never used.

If a download fails, :func:`fetch` raises :class:`DatasetUnavailable` with the URL and
the exact place to put the file by hand; it never substitutes other data.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from tincture.errors import TinctureError


class DatasetUnavailable(TinctureError):
    """A registered file is not cached and could not be downloaded."""


@dataclass(frozen=True)
class RemoteFile:
    """One downloadable file.

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
    """

    dataset: str
    name: str
    url: str
    md5: str | None = None
    sha256: str | None = None
    size: int | None = None
    licence: str = ""
    sha512: str | None = None

    def __post_init__(self) -> None:
        if not (self.md5 or self.sha256 or self.sha512):
            raise ValueError(f"RemoteFile {self.name!r} needs an md5, sha256 or sha512 checksum.")


def data_home() -> Path:
    """Root of the local data cache."""
    env = os.environ.get("TINCTURE_DATA")
    if env:
        return Path(env).expanduser()
    repo = Path(__file__).resolve().parents[2]
    if (repo / "pyproject.toml").exists() and (repo / "tincture").is_dir():
        return repo / "data" / "raw"
    return Path.home() / ".cache" / "tincture"


def local_path(f: RemoteFile, root: Path | None = None) -> Path:
    return (root or data_home()) / f.dataset / f.name


def file_digest(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(f: RemoteFile, path: Path) -> bool:
    """True if ``path`` has the registered size and checksum(s)."""
    if not path.is_file():
        return False
    if f.size is not None and path.stat().st_size != f.size:
        return False
    if f.md5 and file_digest(path, "md5") != f.md5:
        return False
    if f.sha512 and file_digest(path, "sha512") != f.sha512:
        return False
    return not (f.sha256 and file_digest(path, "sha256") != f.sha256)


def is_cached(f: RemoteFile, root: Path | None = None, check: bool = False) -> bool:
    """Whether the file is present (and, with ``check=True``, verified)."""
    p = local_path(f, root)
    return verify(f, p) if check else p.is_file()


def fetch(
    f: RemoteFile,
    root: Path | None = None,
    download: bool = True,
    timeout: float = 600.0,
) -> Path:
    """Return the verified local path of ``f``, downloading it if needed.

    Raises
    ------
    DatasetUnavailable
        If the file is missing and cannot be downloaded, or the checksum fails.
    """
    dest = local_path(f, root)
    if dest.is_file():
        if verify(f, dest):
            return dest
        raise DatasetUnavailable(
            f"{dest} exists but its checksum does not match the registry "
            f"(expected md5={f.md5}, sha256={f.sha256}, sha512={f.sha512}). "
            "It may be truncated or a "
            f"different release. Delete it and re-download from {f.url}."
        )
    if not download:
        raise DatasetUnavailable(
            f"{f.name} is not cached. Download it from {f.url} and place it at {dest}."
        )
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{f.name}.", dir=dest.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        req = urllib.request.Request(f.url, headers={"User-Agent": "tincture-data-fetch"})
        with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as out:
            shutil.copyfileobj(resp, out, length=1 << 20)
    except OSError as exc:  # URLError, HTTPError, timeouts are OSError subclasses
        tmp.unlink(missing_ok=True)
        raise DatasetUnavailable(
            f"Could not download {f.url} ({exc}). Licence: {f.licence or 'see source'}. "
            f"Download it manually and place it at {dest}, then re-run."
        ) from exc
    if not verify(f, tmp):
        tmp.unlink(missing_ok=True)
        raise DatasetUnavailable(
            f"Downloaded {f.url} but the checksum does not match the registry. The source "
            "may have changed; record this in docs/data_status.md before updating the pin."
        )
    tmp.replace(dest)
    return dest
