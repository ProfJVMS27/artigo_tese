#!/usr/bin/env python3
"""Measure batch inference latency for every optimized Phase-3 candidate.

This is a reporting utility only.  It reconstructs each candidate with its
saved optimized hyperparameters, refits on its original training partition,
and times prediction on unlabelled training-design inputs.  Calibration and
test responses are never read or used.  Latency is hardware-dependent and is
reported as a median batch time divided by the number of designs.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
PHASE3 = BASE / "Fase_3_Metamodelagem_Tese_Final.py"
OUT = BASE / "analise_desempenho_modelos (Fase 3) - MCNP63_DEDF_AP"


def load_phase3():
    spec = importlib.util.spec_from_file_location("phase3_latency", PHASE3)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {PHASE3}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=300)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.repeats < 5 or args.batch_size < 10:
        raise SystemExit("Use at least 5 repetitions and batch size of at least 10.")
    summary_path = OUT / "benchmark_cv_summary.csv"
    if not summary_path.exists():
        raise SystemExit(f"Final Phase-3 summary not found: {summary_path}")
    p3 = load_phase3()
    summary = pd.read_csv(summary_path)
    records = []
    for target in p3.TARGETS:
        expected = p3.EXPECTED_CAMPAIGN_CASES if target == p3.TARGET_DOSE else p3.EXPECTED_LSDYNA_VALID_CASES
        path = p3.DOSE_DATASET_PATH if target == p3.TARGET_DOSE else p3.PENETRATION_DATASET_PATH
        data, _ = p3.load_consolidated_target_dataset(path, target, expected)
        train = p3.subset_by_ids(data, p3.make_splits(data)["train"])
        probe = train[p3.RAW_INPUT_COLUMNS].head(min(args.batch_size, len(train))).copy()
        y_train = train[target].to_numpy(dtype=float)
        candidates = summary[summary["target"] == target]
        for row in candidates.itertuples(index=False):
            params = json.loads(row.best_params)
            pipeline = p3.make_pipeline(row.feature_set, p3.model_factories()[row.model]())
            pipeline.set_params(**params)
            pipeline.fit(train[p3.RAW_INPUT_COLUMNS], y_train)
            pipeline.predict(probe)  # warm-up; excluded from timing
            elapsed = []
            for _ in range(args.repeats):
                t0 = time.perf_counter()
                pipeline.predict(probe)
                elapsed.append(time.perf_counter() - t0)
            median_batch = float(np.median(elapsed))
            records.append({
                "target": target,
                "feature_set": row.feature_set,
                "model": row.model,
                "probe_batch_size": len(probe),
                "repeats": args.repeats,
                "median_batch_inference_seconds": median_batch,
                "median_inference_seconds_per_case": median_batch / len(probe),
                "p95_batch_inference_seconds": float(np.quantile(elapsed, 0.95)),
                "timing_scope": "Prediction only; post-fit; unlabelled training-design batch",
            })
            print(f"{target} | {row.feature_set} | {row.model}: {1e3 * median_batch / len(probe):.6f} ms/case", flush=True)
    table_dir = OUT / "tabelas_artigo"
    table_dir.mkdir(parents=True, exist_ok=True)
    output = table_dir / "supplementary_inference_all_optimized_candidates.csv"
    pd.DataFrame(records).to_csv(output, index=False)
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "repeats": args.repeats,
        "batch_size_requested": args.batch_size,
        "calibration_or_test_responses_used": False,
    }
    (table_dir / "inference_benchmark_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
