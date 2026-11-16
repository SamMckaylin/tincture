# Data status

Ground rules: no data is committed (`data/` is git-ignored), no data is invented, and every "verify" fact in HANDOVER §6.1 is checked against its source before a loader is written. Anything that differs from the handover is recorded here.

## Environment

| Date | Finding |
|---|---|
| 2026-11-16 (M0) | The build environment's network allowlist blocked nature.com, cell.com and pubs.acs.org (DNS failure through the proxy). Dataset hosts (GEO/NCBI, figshare, cancerrxgene / Sanger, DepMap) have **not been tested yet**; that happens at the start of M1. If they are blocked, the M1 loaders will be written and tested on synthetic fixtures only, and the human will be asked to place the files in `data/raw/<dataset>/` or allow the hosts. |

## Datasets

| ID | Dataset | Status | Verified facts | Differences from handover |
|---|---|---|---|---|
| D1 | sci-Plex (GSE139944) | not started (M1) | | |
| D2 | PRISM Repurposing secondary | not started (M1) | | |
| D3 | GDSC monotherapy raw | not started (M1) | | |
| D4–D12 | | not started (M5–M7) | | |

## Unmatched compound names

None yet.
