"""Select a pre-specified, diverse Phase-5 validation set from a Pareto union.

The selection is deterministic and is performed before any high-fidelity
solver run. It preserves the three single-objective extremes and the TOPSIS
compromise, then fills the remaining slots by maximin coverage in normalized
objective space. Every selected case must be re-evaluated and reported; this
tool must not be used for post-hoc selection of favorable validation results.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


OBJECTIVES = [
    "dose_pred_5min_mSv",
    "massa_total_kg",
    "custo_total_unidade",
]
REQUIRED = OBJECTIVES + [
    "pareto_id",
    "run_id",
    "seed",
    "ballistic_constraint_violation_mm",
    "ballistic_feasible_95PI",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def topsis_iqr(F: np.ndarray) -> int:
    q1 = np.percentile(F, 25, axis=0)
    q3 = np.percentile(F, 75, axis=0)
    scale = np.where((q3 - q1) > 1e-15, q3 - q1, 1.0)
    Z = (F - q1[None, :]) / scale[None, :]
    ideal = np.min(Z, axis=0)
    anti_ideal = np.max(Z, axis=0)
    distance_ideal = np.linalg.norm(Z - ideal[None, :], axis=1)
    distance_anti = np.linalg.norm(Z - anti_ideal[None, :], axis=1)
    score = distance_anti / (distance_ideal + distance_anti + 1e-15)
    return int(np.argmax(score))


def select_candidates(df: pd.DataFrame, n_select: int) -> pd.DataFrame:
    missing = [column for column in REQUIRED if column not in df.columns]
    if missing:
        raise ValueError(f"Pareto source is missing columns: {missing}")
    if n_select < 5:
        raise ValueError("Use at least five candidates to retain all pre-specified anchors.")
    if n_select > len(df):
        raise ValueError(f"Requested {n_select} candidates, but Pareto union has {len(df)} rows.")
    if not df[REQUIRED].notna().all().all():
        raise ValueError("Pareto source has missing values in required columns.")
    if not np.allclose(df["ballistic_constraint_violation_mm"].astype(float), 0.0):
        raise ValueError("Pareto source contains a ballistic-infeasible solution.")
    if not df["ballistic_feasible_95PI"].astype(int).eq(1).all():
        raise ValueError("Pareto source contains a non-feasible conformal prediction.")

    out = df.copy().reset_index(drop=True)
    F = out[OBJECTIVES].to_numpy(dtype=float)
    selected: dict[int, list[str]] = {}

    def add(index: int, reason: str) -> None:
        selected.setdefault(int(index), []).append(reason)

    add(int(np.argmin(F[:, 0])), "MIN_DOSE")
    add(int(np.argmin(F[:, 1])), "MIN_MASS")
    add(int(np.argmin(F[:, 2])), "MIN_COST")
    add(topsis_iqr(F), "TOPSIS_IQR")

    lo = np.min(F, axis=0)
    hi = np.max(F, axis=0)
    span = np.where((hi - lo) > 1e-15, hi - lo, 1.0)
    Z = (F - lo[None, :]) / span[None, :]

    while len(selected) < n_select:
        selected_idx = np.asarray(sorted(selected), dtype=int)
        remaining = np.asarray(
            [i for i in range(len(out)) if i not in selected], dtype=int
        )
        distances = np.linalg.norm(
            Z[remaining, None, :] - Z[selected_idx][None, :, :], axis=2
        )
        add(int(remaining[int(np.argmax(np.min(distances, axis=1)))]), "MAXIMIN_PARETO_COVERAGE")

    rows = []
    for index in sorted(selected):
        row = out.iloc[index].copy()
        row["validation_selection_reason"] = " | ".join(selected[index])
        rows.append(row)
    selected_df = pd.DataFrame(rows).reset_index(drop=True)
    selected_df.insert(0, "validation_candidate_id", np.arange(1, len(selected_df) + 1))
    return selected_df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pareto", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--n-candidates", type=int, default=30)
    args = parser.parse_args()

    if not args.pareto.is_file():
        raise FileNotFoundError(args.pareto)
    if args.output_csv.exists():
        raise FileExistsError(f"Refusing to overwrite: {args.output_csv}")

    selected = select_candidates(pd.read_csv(args.pareto), args.n_candidates)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output_csv, index=False)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection_is_pre_specified": True,
        "post_hoc_favorable_selection_prohibited": True,
        "candidate_source": str(args.pareto),
        "candidate_source_sha256": sha256_file(args.pareto),
        "n_candidates": int(len(selected)),
        "selection_rule": (
            "MIN_DOSE, MIN_MASS, MIN_COST and TOPSIS_IQR anchors, followed "
            "by deterministic maximin coverage in normalized three-objective space"
        ),
        "all_selected_cases_must_be_reported": True,
    }
    manifest_path = args.output_csv.with_name("manifesto_selecao_validacao_fase5.json")
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(args.output_csv)
    print(manifest_path)


if __name__ == "__main__":
    main()
