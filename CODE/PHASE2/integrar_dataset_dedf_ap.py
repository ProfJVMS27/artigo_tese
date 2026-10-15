#!/usr/bin/env python3
"""Create a versioned MCNP DE/DF-AP + LS-DYNA dataset without overwriting sources."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

EXPECTED_MCNP = set(range(1, 2001))
EXPECTED_LSDYNA_ML = 1996


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def require_ids(frame: pd.DataFrame, label: str, expected: set[int], id_col: str = "sim_id") -> None:
    if id_col not in frame.columns:
        raise ValueError(f"{label}: missing {id_col}")
    ids = set(pd.to_numeric(frame[id_col], errors="raise").astype(int))
    if ids != expected or frame[id_col].duplicated().any():
        raise ValueError(f"{label}: IDs are not exactly the required set or are duplicated")


def update_mcnp_columns(frame: pd.DataFrame, audit: pd.DataFrame, source_result_csv: str) -> pd.DataFrame:
    out = frame.copy()
    id_col = "sim_id" if "sim_id" in out.columns else "id_simulacao"
    merged = out.merge(
        audit[["sim_id", "dose_efetiva_mSv_seg", "erro_relativo"]].rename(columns={"dose_efetiva_mSv_seg": "dedf_rate", "erro_relativo": "dedf_rel"}),
        left_on=id_col, right_on="sim_id", how="left", validate="one_to_one",
    )
    rate = pd.to_numeric(merged["dedf_rate"], errors="raise")
    rel = pd.to_numeric(merged["dedf_rel"], errors="raise")
    if rate.isna().any() or rel.isna().any():
        raise ValueError("MCNP DE/DF merge produced missing dose values")
    # Dataset-specific names are intentionally updated only if present.
    for col in ("dose_efetiva_mSv_seg", "mcnp_dose_efetiva_mSv_seg"):
        if col in merged.columns:
            merged[col] = rate
    for col in ("dose_efetiva_5min_mSv", "mcnp_dose_efetiva_5min_mSv"):
        if col in merged.columns:
            merged[col] = 300.0 * rate
    for col in ("erro_relativo_mcnp", "mcnp_erro_relativo"):
        if col in merged.columns:
            merged[col] = rel
    for col in ("mcnp_valid",):
        if col in merged.columns:
            merged[col] = True
    if "mcnp_dose_5min_method" in merged.columns:
        merged["mcnp_dose_5min_method"] = "DE4_DF4_ICRP116_AP_FM4_7400_mSv_s_x_300_s"
    if "mcnp_source_result_csv" in merged.columns:
        merged["mcnp_source_result_csv"] = source_result_csv
    if "mcnp_design_precision_note" in merged.columns:
        merged["mcnp_design_precision_note"] = (
            "MCNP6.3 campaign with frozen normalized template; DE4/DF4 ICRP-116 AP response; FM4=7400 mSv/s."
        )
    drops = ["dedf_rate", "dedf_rel"]
    if id_col != "sim_id":
        drops.append("sim_id")
    return merged.drop(columns=drops)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dataset-dir", required=True, type=Path)
    parser.add_argument("--dedf-results", required=True, type=Path)
    parser.add_argument("--dedf-audit-json", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    src, out = args.source_dataset_dir.resolve(), args.output_dir.resolve()
    if out.exists():
        raise SystemExit(f"Refusing to overwrite existing output directory: {out}")
    audit_report = json.loads(args.dedf_audit_json.read_text(encoding="utf-8"))
    if audit_report.get("result") != "PASS" or audit_report.get("valid_cases") != 2000:
        raise SystemExit("DE/DF audit is not a complete PASS; integration aborted.")
    audit = read_csv(args.dedf_results)
    require_ids(audit, "DE/DF audit results", EXPECTED_MCNP)
    if set(audit["status"].astype(str)) != {"PASS"}:
        raise SystemExit("DE/DF audit contains non-PASS cases; integration aborted.")
    for col in ("input_sha256_match", "dedf_contract", "output_exists", "nps_100m", "tally4_present", "tfc_bin_error_pass"):
        if not audit[col].astype(str).str.lower().eq("true").all():
            raise SystemExit(f"DE/DF audit contract failed in column {col}; integration aborted.")

    names = {
        "mcnp": "MCNP63_RESULTADOS_2000_FINAL.csv",
        "dose": "DATASET_ML_DOSE_MCNP63_2000_FINAL.csv",
        "master": "DATASET_MESTRE_CONJUNTO_MCNP63_LSDYNA_2000_FINAL.csv",
        "joint": "DATASET_ML_CONJUNTO_MCNP63_LSDYNA_1996_FINAL.csv",
        "penetration": "DATASET_ML_PENETRACAO_LSDYNA_1996_FINAL.csv",
    }
    source_paths = {key: src / value for key, value in names.items()}
    missing = [str(path) for path in source_paths.values() if not path.exists()]
    if missing:
        raise SystemExit(f"Missing source datasets: {missing}")
    frames = {key: read_csv(path) for key, path in source_paths.items()}
    require_ids(frames["mcnp"], "mcnp", EXPECTED_MCNP, id_col="id_simulacao")
    for key in ("dose", "master"):
        require_ids(frames[key], key, EXPECTED_MCNP)
    for key in ("joint", "penetration"):
        if len(frames[key]) != EXPECTED_LSDYNA_ML or frames[key]["sim_id"].duplicated().any():
            raise SystemExit(f"{key}: expected {EXPECTED_LSDYNA_ML} non-duplicated rows")

    out.mkdir(parents=True)
    audit_rate = audit[["sim_id", "dose_efetiva_mSv_seg", "erro_relativo"]].copy()
    audit_rate["dose_efetiva_mSv_seg"] = pd.to_numeric(audit_rate["dose_efetiva_mSv_seg"], errors="raise")
    audit_rate["erro_relativo"] = pd.to_numeric(audit_rate["erro_relativo"], errors="raise")
    source_label = "MCNP63_DEDF_AP_RESULTADOS_2000_FINAL.csv"
    updated = {key: update_mcnp_columns(frame, audit_rate, source_label) for key, frame in frames.items()}

    # Preserve original MCNP status columns while replacing only the scientific response.
    updated["mcnp"]["dose_efetiva_mSv_seg"] = audit_rate.set_index("sim_id").loc[updated["mcnp"]["id_simulacao"].astype(int), "dose_efetiva_mSv_seg"].to_numpy()
    updated["mcnp"]["erro_relativo"] = audit_rate.set_index("sim_id").loc[updated["mcnp"]["id_simulacao"].astype(int), "erro_relativo"].to_numpy()
    updated["mcnp"]["arquivo_output"] = [
        audit.loc[audit["sim_id"].astype(int).eq(int(sim)), "output_relpath"].iloc[0]
        for sim in updated["mcnp"]["id_simulacao"]
    ]

    output_names = {
        "mcnp": "MCNP63_DEDF_AP_RESULTADOS_2000_FINAL.csv",
        "dose": "DATASET_ML_DOSE_MCNP63_DEDF_AP_2000_FINAL.csv",
        "master": "DATASET_MESTRE_CONJUNTO_MCNP63_DEDF_AP_LSDYNA_2000_FINAL.csv",
        "joint": "DATASET_ML_CONJUNTO_MCNP63_DEDF_AP_LSDYNA_1996_FINAL.csv",
        "penetration": "DATASET_ML_PENETRACAO_LSDYNA_1996_FINAL.csv",
    }
    for key, filename in output_names.items():
        updated[key].to_csv(out / filename, index=False)
    audit.copy().to_csv(out / "AUDITORIA_CAMPANHA_DEDF_AP_2000.csv", index=False)

    manifest = {
        "created_utc": utc_now(), "integration_status": "PASS",
        "response_contract": audit_report["response_contract"],
        "mcnp_cases": 2000, "lsdyna_ml_cases": EXPECTED_LSDYNA_ML,
        "source_dataset_dir": str(src), "dedf_results": str(args.dedf_results.resolve()),
        "dedf_audit_json": str(args.dedf_audit_json.resolve()),
        "source_file_sha256": {key: sha256(path) for key, path in source_paths.items()},
        "output_file_sha256": {key: sha256(out / filename) for key, filename in output_names.items()},
        "invariants": {
            "mcnp_ids_exact_1_to_2000": True, "mcnp_valid_cases": 2000,
            "lsdyna_ml_rows_preserved": EXPECTED_LSDYNA_ML,
            "lsdyna_columns_unchanged_except_joined_mcnp_response": True,
            "source_datasets_not_overwritten": True,
        },
    }
    (out / "MANIFESTO_INTEGRACAO_DEDF_AP.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
