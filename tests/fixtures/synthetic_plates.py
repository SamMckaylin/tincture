"""SYNTHETIC plates for Plate Atlas tests only (ground rule 2)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from tincture.model import GEOMETRY, PlateTable


def synthetic_plates(
    n_plates: int = 6,
    n_wells: int = 384,
    edge_factor: float = 1.0,
    unseeded_outer: bool = False,
    drift_per_plate: float = 0.0,
    bad_plates: tuple[int, ...] = (),
    cell_lines: tuple[str, ...] = ("L1",),
    line_brightness: tuple[float, ...] = (1.0,),
    drop_wells: int = 0,
    seed: int = 0,
) -> PlateTable:
    rng = np.random.default_rng(seed)
    nr, nc = GEOMETRY[n_wells]
    rows = []
    for p in range(n_plates):
        line = cell_lines[p % len(cell_lines)]
        bright = line_brightness[p % len(line_brightness)]
        level = 1000 * bright * (1 + drift_per_plate * p)
        noise = 0.04 if p not in bad_plates else 0.25
        for r in range(nr):
            for c in range(nc):
                depth = min(r, c, nr - 1 - r, nc - 1 - c)
                if c == 1 and 2 <= r < nr - 2:
                    wt, mu = "blank", 0.02 * level
                elif c == nc // 2 and 2 <= r < nr - 2:
                    wt, mu = "neg_ctrl", level
                elif c == nc // 2 + 1 and 2 <= r < nr - 2:
                    wt, mu = "pos_ctrl", 0.1 * level
                elif depth == 0:
                    wt, mu = "empty", (0.02 if unseeded_outer else edge_factor) * level
                elif unseeded_outer and depth == 1:
                    wt, mu = "empty", edge_factor * level
                elif depth <= 1:
                    wt, mu = "empty", level
                else:
                    wt, mu = "compound", level * rng.uniform(0.3, 1.0)
                rows.append((f"P{p:02d}", r, c, wt, mu * (1 + rng.normal(0, noise)), p, line))
    d = pd.DataFrame(
        rows, columns=["plate", "row", "col", "well_type", "value", "plate_order", "cell_line"]
    )
    if drop_wells:
        d = d.drop(index=rng.choice(d.index, drop_wells, replace=False)).reset_index(drop=True)
    return PlateTable(d, n_wells, source="synthetic")
