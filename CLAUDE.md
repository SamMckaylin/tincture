# Working on tincture

The full brief is in `HANDOVER.md`. Read it before doing anything. Key rules (HANDOVER §0):

- Work milestone by milestone (§8). Commit at the end of each, write `docs/milestone_<n>.md`, then stop and summarise.
- Never invent data. Synthetic data lives only in `tests/fixtures/synthetic_*`. If a download fails, stop and record it in `docs/data_status.md`.
- Never commit raw data (`data/` is git-ignored).
- Verify "verify" dataset facts against the source before writing loaders; record the results in `docs/data_status.md`.
- Ask before adding heavy dependencies (R, Java, CUDA, more than ~50 MB) and before settling §9 decisions.
- Every red-team guard in §5 needs its own test.

Dev loop: `source .venv/bin/activate && ruff check . && ruff format --check . && mypy && pytest`.
Architecture: stats -> `PlotSpec` (`tincture/spec`) -> `render_mpl` (`tincture/render/mpl.py`). No matplotlib objects outside `render/`, `theme/` and `export.py`.
