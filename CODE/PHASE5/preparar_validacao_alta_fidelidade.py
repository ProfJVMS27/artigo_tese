"""Prepare pre-specified Phase-4 high-fidelity validation cases without running solvers.

This tool is deliberately generation-only.  It creates new, isolated MCNP
inputs and LS-DYNA keywords from the frozen Phase-2 generators, records hashes,
and never changes the original campaign, its keywords, or either frozen model.
Use the existing ``run_lsdyna_study.py`` only after reviewing this manifest and
obtaining an explicit decision to execute LS-DYNA.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType

import pandas as pd


N_LAYERS = 7
TOTAL_THICKNESS_CM = 6.0
MIN_THICKNESS_CM = 0.2
MAX_THICKNESS_CM = 2.0
MATERIALS = {"SiC", "UHMWPE", "Al", "Fe", "B4C-HDPE", "W-Epoxy", "Bi2O3-HDPE"}
VALIDATION_SIM_START = 9001


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Não foi possível importar {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def material_columns() -> list[str]:
    return [f"material_camada_{i}" for i in range(1, N_LAYERS + 1)]


def thickness_columns() -> list[str]:
    return [f"espessura_camada_{i}" for i in range(1, N_LAYERS + 1)]


def validation_sim_ids(n_candidates: int) -> list[int]:
    if n_candidates < 1:
        raise ValueError("A validação deve conter pelo menos um candidato.")
    return list(range(VALIDATION_SIM_START, VALIDATION_SIM_START + n_candidates))


def validate_candidates(df: pd.DataFrame) -> pd.DataFrame:
    required = ["validation_candidate_id", "pareto_id", "run_id", "seed"]
    required += material_columns() + thickness_columns()
    required += [
        "dose_pred_5min_mSv", "massa_total_kg", "custo_total_unidade",
        "max_penetration_pred_mm", "penetration_pi_upper95_mm",
        "ballistic_constraint_violation_mm", "ballistic_feasible_95PI",
        "validation_selection_reason",
    ]
    missing = [name for name in required if name not in df.columns]
    if missing:
        raise ValueError(f"CSV de candidatos sem colunas obrigatórias: {missing}")
    if df["validation_candidate_id"].duplicated().any():
        raise ValueError("validation_candidate_id duplicado.")
    if df[required].isna().any().any():
        bad = df[required].isna().sum()
        raise ValueError(f"Dados ausentes: {bad[bad.gt(0)].to_dict()}")

    out = df.sort_values("validation_candidate_id").copy().reset_index(drop=True)
    expected_ids = list(range(1, len(out) + 1))
    if out["validation_candidate_id"].astype(int).tolist() != expected_ids:
        raise ValueError(
            "Os candidatos devem estar identificados sequencialmente de 1 a N."
        )
    for row_index, row in out.iterrows():
        mats = [str(row[col]) for col in material_columns()]
        ths = [float(row[col]) for col in thickness_columns()]
        if set(mats) != MATERIALS or len(set(mats)) != N_LAYERS:
            raise ValueError(f"Candidato {row_index + 1}: materiais inválidos.")
        if any(not math.isfinite(value) or value < MIN_THICKNESS_CM - 1e-12 or value > MAX_THICKNESS_CM + 1e-12 for value in ths):
            raise ValueError(f"Candidato {row_index + 1}: espessura fora de [0,2; 2,0] cm.")
        if not math.isclose(sum(ths), TOTAL_THICKNESS_CM, rel_tol=0.0, abs_tol=1e-6):
            raise ValueError(f"Candidato {row_index + 1}: soma das espessuras não é 6,0 cm.")
        if float(row["ballistic_constraint_violation_mm"]) != 0.0 or int(row["ballistic_feasible_95PI"]) != 1:
            raise ValueError(f"Candidato {row_index + 1}: não é viável na restrição balística Fase 4.")
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_mcnp_runner(path: Path, sim_ids: list[int]) -> None:
    """Write the institutional runner; task count defaults to eight threads."""
    path.write_text(
        "@echo off\r\n"
        "setlocal EnableExtensions\r\n"
        "if \"%~1\"==\"\" (\r\n"
        "  echo Uso: run_mcnp_validation.cmd ^\"C:\\caminho\\mcnp6.exe^\" [tasks]\r\n"
        "  exit /b 2\r\n"
        ")\r\n"
        "set \"MCNP_EXE=%~1\"\r\n"
        "set \"MCNP_TASKS=%~2\"\r\n"
        "if \"%MCNP_TASKS%\"==\"\" set \"MCNP_TASKS=8\"\r\n"
        f"for %%I in ({' '.join(str(sim_id) for sim_id in sim_ids)}) do (\r\n"
        "  echo ===== MCNP sim_%%I; tasks=%MCNP_TASKS% =====\r\n"
        "  \"%MCNP_EXE%\" i=\"mcnp\\sim_%%I\\input.i\" o=\"mcnp\\sim_%%I\\output_%%I.o\" r=\"mcnp\\sim_%%I\\runtpe_%%I.r\" tasks %MCNP_TASKS%\r\n"
        "  if errorlevel 1 exit /b 1\r\n"
        ")\r\n"
        "echo Finalizado. Copie de volta todas as pastas mcnp\\sim_XXXX.\r\n",
        encoding="ascii",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--study-root", type=Path, required=True)
    parser.add_argument("--prepare", action="store_true", help="Cria inputs isolados; sem esta flag nada é escrito.")
    parser.add_argument(
        "--refresh-mcnp-inputs", action="store_true",
        help="Regenera somente os inputs MCNP e seus hashes; não toca nos casos LS-DYNA.",
    )
    parser.add_argument(
        "--mcnp-reference", type=Path,
        default=Path(r"C:\Users\Admin\Desktop\Scripts Tese 1\Fase 2\mcnp63_canonical_validation.py"),
    )
    parser.add_argument(
        "--lsdyna-model", type=Path,
        default=Path(r"C:\Users\Admin\Desktop\Scripts Tese 1\Fase 2\LS-DYNA\lsdyna_modelo_mecanico_congelado.py"),
    )
    parser.add_argument(
        "--dedf-response-builder", type=Path,
        default=Path(
            r"C:\Users\Admin\Desktop\Scripts Tese 1\Fase 2\MCNP"
            r"\CAMPANHA_DEDF_AP_MCNP63_20261001\prepare_campaign_dedf_ap.py"
        ),
        help="Gerador congelado do bloco DE4/DF4 AP, FM4=7400.",
    )
    args = parser.parse_args()

    if args.prepare and args.refresh_mcnp_inputs:
        raise ValueError("Use apenas uma de --prepare ou --refresh-mcnp-inputs.")

    for path in (
        args.candidates,
        args.mcnp_reference,
        args.lsdyna_model,
        args.dedf_response_builder,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    candidates = validate_candidates(pd.read_csv(args.candidates))
    if not args.refresh_mcnp_inputs and args.study_root.exists() and any(args.study_root.iterdir()):
        raise FileExistsError(f"Diretório de estudo já existe ou não está vazio: {args.study_root}")

    # These imports and MCNP builds are a dry validation of the frozen sources.
    mcnp = load_module("frozen_mcnp_reference", args.mcnp_reference)
    dedf = load_module("frozen_dedf_response", args.dedf_response_builder)
    lsdyna = load_module("frozen_lsdyna_model", args.lsdyna_model)
    for _, row in candidates.iterrows():
        mats, ths = mcnp.validate_design_row(row)
        canonical = mcnp.build_mcnp_input(
            9000 + int(row["validation_candidate_id"]), mats, ths
        )
        mcnp.validate_input_text(canonical)
        converted = dedf.replace_response_block(canonical)
        if "DE4" not in converted or "DF4" not in converted or "FM4 7400" not in converted:
            raise ValueError("Falha ao aplicar a resposta DE4/DF4 AP ao deck MCNP.")

    if args.refresh_mcnp_inputs:
        source = args.study_root / "preparation_manifest.csv"
        rows = list(csv.DictReader(source.open(encoding="utf-8-sig", newline=""))) if source.is_file() else []
        if len(rows) != len(candidates):
            raise ValueError("O manifesto existente deve conter o mesmo número de candidatos.")
        candidate_by_id = {int(row["validation_candidate_id"]): row for _, row in candidates.iterrows()}
        for record in rows:
            candidate = candidate_by_id[int(record["validation_candidate_id"])]
            mats, ths = mcnp.validate_design_row(candidate)
            target = Path(record["mcnp_input"])
            canonical = mcnp.build_mcnp_input(
                int(record["validation_sim_id"]), mats, ths
            )
            mcnp.validate_input_text(canonical)
            target.write_text(dedf.replace_response_block(canonical), encoding="utf-8")
            record["mcnp_input_sha256"] = sha256_file(target)
            record["status"] = "MCNP_INPUTS_REGENERATED_PENDING_EXECUTION"
        write_csv(source, rows)
        manifest_path = args.study_root / "manifesto_preparacao.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["mcnp_reference_sha256"] = sha256_file(args.mcnp_reference)
        manifest["dedf_response_builder_sha256"] = sha256_file(args.dedf_response_builder)
        manifest["mcnp_inputs_regenerated_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["mcnp_executed"] = False
        manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        write_mcnp_runner(
            args.study_root / "run_mcnp_validation.cmd",
            validation_sim_ids(len(candidates)),
        )
        print(f"Inputs MCNP regenerados: {source}")
        return

    print(f"Candidatos: PASS ({len(candidates)}/{len(candidates)})")
    print("Referência MCNP e modelo LS-DYNA congelados: PASS")
    if not args.prepare:
        print("DRY VALIDATION ONLY — use --prepare to create isolated inputs; no solver is called by this tool.")
        return

    args.study_root.mkdir(parents=True, exist_ok=False)
    mcnp_root = args.study_root / "mcnp"
    lsdyna_root = args.study_root / "lsdyna"
    mcnp_root.mkdir()
    lsdyna_root.mkdir()
    candidate_copy = args.study_root / "candidatos_validacao_alta_fidelidade_fase4.csv"
    candidates.to_csv(candidate_copy, index=False)

    rows: list[dict] = []
    sim_ids = validation_sim_ids(len(candidates))
    for (_, row), validation_sim_id in zip(candidates.iterrows(), sim_ids):
        candidate_id = int(row["validation_candidate_id"])
        mats, ths = mcnp.validate_design_row(row)

        mcnp_case = mcnp_root / f"sim_{validation_sim_id:04d}"
        mcnp_case.mkdir()
        mcnp_input = mcnp_case / "input.i"
        canonical = mcnp.build_mcnp_input(validation_sim_id, mats, ths)
        mcnp.validate_input_text(canonical)
        mcnp_input.write_text(dedf.replace_response_block(canonical), encoding="utf-8")

        lsdyna_case = lsdyna_root / f"sim_{validation_sim_id:04d}"
        lsdyna_case.mkdir()
        lsdyna.configure_case(
            output_dir=lsdyna_case,
            sim_id=validation_sim_id,
            layer_materials=mats,
            layer_thicknesses_mm=[10.0 * value for value in ths],
        )
        lsdyna_summary = lsdyna.write_keyword_file()
        lsdyna_keyword = Path(lsdyna_summary["output_k"])
        if not lsdyna_keyword.is_file():
            raise RuntimeError(f"Keyword LS-DYNA não gerado: {lsdyna_keyword}")

        rows.append({
            "validation_candidate_id": candidate_id,
            "validation_sim_id": validation_sim_id,
            "pareto_id": int(row["pareto_id"]),
            "run_id": int(row["run_id"]),
            "seed": int(row["seed"]),
            "selection_reason": row["validation_selection_reason"],
            "mcnp_input": str(mcnp_input),
            "mcnp_input_sha256": sha256_file(mcnp_input),
            "lsdyna_keyword": str(lsdyna_keyword),
            "lsdyna_keyword_sha256": sha256_file(lsdyna_keyword),
            "lsdyna_precheck": str(lsdyna_case / "PRE_RUN_CHECK.txt"),
            "pred_dose_5min_mSv": float(row["dose_pred_5min_mSv"]),
            "pred_penetration_mm": float(row["max_penetration_pred_mm"]),
            "pred_penetration_upper95_mm": float(row["penetration_pi_upper95_mm"]),
            "status": "INPUTS_PREPARED_NOT_EXECUTED",
        })

    write_csv(args.study_root / "preparation_manifest.csv", rows)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "Fase 4 high-fidelity validation preparation only",
        "candidate_source": str(args.candidates),
        "candidate_source_sha256": sha256_file(args.candidates),
        "mcnp_reference": str(args.mcnp_reference),
        "mcnp_reference_sha256": sha256_file(args.mcnp_reference),
        "dedf_response_builder": str(args.dedf_response_builder),
        "dedf_response_builder_sha256": sha256_file(args.dedf_response_builder),
        "mcnp_response_contract": "DE4/DF4 ICRP-116 AP, DF pSv cm2, FM4=7400 mSv/s",
        "lsdyna_frozen_model": str(args.lsdyna_model),
        "lsdyna_frozen_model_sha256": sha256_file(args.lsdyna_model),
        "n_candidates": len(rows),
        "mcnp_executed": False,
        "lsdyna_executed": False,
        "original_campaign_modified": False,
        "original_keywords_modified": False,
    }
    (args.study_root / "manifesto_preparacao.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    run_cmd = args.study_root / "run_mcnp_validation.cmd"
    write_mcnp_runner(run_cmd, sim_ids)
    print(args.study_root / "preparation_manifest.csv")
    print(args.study_root / "manifesto_preparacao.json")
    print(run_cmd)


if __name__ == "__main__":
    main()
