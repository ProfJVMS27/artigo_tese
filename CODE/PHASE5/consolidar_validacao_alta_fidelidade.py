"""Consolidate pre-specified MCNP6.3 / LS-DYNA Phase-5 validations."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-root", required=True, type=Path)
    parser.add_argument(
        "--phase3-model-dir", type=Path,
        default=Path(
            r"C:\Users\Admin\Desktop\Scripts Tese 1\Fase 3"
            r"\modelos_treinados (Fase 3) - MCNP63_DEDF_AP"
        ),
    )
    args = parser.parse_args()
    root = args.study_root
    prepared = pd.read_csv(root / "preparation_manifest.csv")
    mcnp = pd.read_csv(root / "mcnp_run_manifest.csv")
    lsdyna = pd.read_csv(root / "lsdyna" / "study_validation_summary.csv")
    n_cases = len(prepared)
    if n_cases < 1 or len(mcnp) != n_cases or len(lsdyna) != n_cases:
        raise ValueError(
            "Os manifestos de preparação, MCNP e LS-DYNA devem conter o mesmo "
            "número não nulo de casos."
        )
    out = prepared.merge(mcnp, on=["validation_candidate_id", "validation_sim_id"], validate="one_to_one")
    out = out.merge(lsdyna, left_on="validation_sim_id", right_on="sim_id", validate="one_to_one")
    required = (
        out["status_execucao_mcnp"].eq("COMPLETED")
        & out["status_extracao"].eq("COMPLETED")
        & out["qa_tally"].eq("PASS")
        & out["preliminary_decision"].eq("PASS_PRELIMINARY")
    )
    if not required.all():
        bad = out.loc[~required, ["validation_candidate_id", "status_execucao_mcnp", "status_extracao", "qa_tally", "preliminary_decision"]]
        raise ValueError("Validação incompleta/reprovada:\n" + bad.to_string(index=False))
    dose_meta = json.loads((args.phase3_model_dir / "surrogate_dose_metadata.json").read_text(encoding="utf-8"))
    pen_meta = json.loads((args.phase3_model_dir / "surrogate_penetration_metadata.json").read_text(encoding="utf-8"))
    dose_q = float(dose_meta["conformal_q_abs"])
    pen_q = float(pen_meta["conformal_q_abs"])
    out["hf_dose_5min_mSv"] = out["dose_efetiva_mSv_seg"].astype(float) * 300.0
    out["dose_error_mSv"] = out["hf_dose_5min_mSv"] - out["pred_dose_5min_mSv"]
    out["penetration_error_mm"] = out["max_penetration_mm"].astype(float) - out["pred_penetration_mm"]
    out["dose_within_conformal_95"] = out["dose_error_mSv"].abs() <= dose_q
    out["penetration_within_conformal_95"] = out["penetration_error_mm"].abs() <= pen_q
    out["abs_dose_error_mSv"] = out["dose_error_mSv"].abs()
    out["abs_penetration_error_mm"] = out["penetration_error_mm"].abs()
    dest = root / "validacao_alta_fidelidade_consolidada.csv"
    out.to_csv(dest, index=False)
    report = {
        "n_cases": int(n_cases), "dose_conformal_q95_mSv": dose_q, "penetration_conformal_q95_mm": pen_q,
        "dose_mae_mSv": float(out["abs_dose_error_mSv"].mean()),
        "penetration_mae_mm": float(out["abs_penetration_error_mm"].mean()),
        "dose_within_conformal_95": int(out["dose_within_conformal_95"].sum()),
        "penetration_within_conformal_95": int(out["penetration_within_conformal_95"].sum()),
        "result_csv": str(dest),
    }
    (root / "auditoria_validacao_alta_fidelidade.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(dest)


if __name__ == "__main__":
    main()
