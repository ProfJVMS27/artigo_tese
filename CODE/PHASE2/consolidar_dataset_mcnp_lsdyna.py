"""Consolida a campanha MCNP e o mestre LS-DYNA sem alterar as fontes.

Fonte canônica de dose: ``dose_efetiva_mSv_seg`` do resultado primário MCNP.
O alvo de cinco minutos é calculado deterministicamente como taxa * 300 s.
As quatro etiquetas balísticas não validadas permanecem ausentes somente nos
arquivos ML; os valores brutos continuam preservados no dataset mestre.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPECTED_CASES = 2000
EXPECTED_LSDYNA_VALID = 1996
PENDING_LSDYNA_IDS = {733, 832, 1277, 1567}
DOSE_SECONDS = 300.0


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV sem cabeçalho: {path}")
        return list(reader), list(reader.fieldnames)


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def parse_int(value: str, label: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: inteiro inválido {value!r}") from exc


def parse_float(value: str, label: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: número inválido {value!r}") from exc


def csv_index(rows: list[dict[str, str]], id_column: str, source_name: str) -> dict[int, dict[str, str]]:
    indexed: dict[int, dict[str, str]] = {}
    for row in rows:
        sim_id = parse_int(row.get(id_column, ""), f"{source_name}.{id_column}")
        if sim_id in indexed:
            raise ValueError(f"{source_name}: sim_id duplicado {sim_id}")
        indexed[sim_id] = row
    return indexed


def required_columns(columns: list[str], names: list[str], source_name: str) -> None:
    missing = [name for name in names if name not in columns]
    if missing:
        raise ValueError(f"{source_name}: colunas ausentes: {missing}")


def expected_ids(indexed: dict[int, dict[str, str]], source_name: str) -> None:
    actual = set(indexed)
    expected = set(range(1, EXPECTED_CASES + 1))
    if actual != expected:
        missing = sorted(expected - actual)[:20]
        extra = sorted(actual - expected)[:20]
        raise ValueError(f"{source_name}: IDs incompatíveis; ausentes={missing}; extras={extra}")


def normalize_status(value: str) -> str:
    return str(value).strip().upper()


def format_number(value: float) -> str:
    return format(value, ".12g")


def build(args: argparse.Namespace) -> dict[str, Any]:
    lsdyna_path = args.lsdyna_master.resolve()
    mcnp_path = args.mcnp_results.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    lsdyna_rows, lsdyna_columns = read_csv(lsdyna_path)
    mcnp_rows, mcnp_columns = read_csv(mcnp_path)
    layer_materials = [f"material_camada_{i}" for i in range(1, 8)]
    layer_thicknesses = [f"espessura_camada_{i}_cm" for i in range(1, 8)]
    required_columns(
        lsdyna_columns,
        ["sim_id", "dataset_role", "working_max_penetration_mm", "working_final_vx_m_s",
         "working_failed_elements_total", "working_energy_ratio_final", "working_ballistic_status"]
        + layer_materials + layer_thicknesses,
        "mestre LS-DYNA",
    )
    required_columns(
        mcnp_columns,
        ["id_simulacao", "status_geracao_input", "status_execucao_mcnp", "status_extracao",
         "dose_efetiva_mSv_seg", "erro_relativo", "materiais", "espessuras"],
        "resultados MCNP",
    )

    lsdyna = csv_index(lsdyna_rows, "sim_id", "mestre LS-DYNA")
    mcnp = csv_index(mcnp_rows, "id_simulacao", "resultados MCNP")
    expected_ids(lsdyna, "mestre LS-DYNA")
    expected_ids(mcnp, "resultados MCNP")

    mcnp_failures: list[int] = []
    material_mismatches: list[int] = []
    max_thickness_diff_cm = 0.0
    dose_min = float("inf")
    dose_max = float("-inf")
    relative_error_min = float("inf")
    relative_error_max = float("-inf")
    master_rows: list[dict[str, Any]] = []
    ml_rows: list[dict[str, Any]] = []
    dose_rows: list[dict[str, Any]] = []
    penetration_rows: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []

    ml_fields = (
        ["sim_id"] + layer_materials + layer_thicknesses
        + ["peso_total_kg_csv", "custo_total_unidade_csv", "massa_blindagem_modelo_final_kg",
           "dose_efetiva_mSv_seg", "dose_efetiva_5min_mSv", "erro_relativo_mcnp", "mcnp_valid",
           "max_penetration_mm", "final_fsp_vx_m_s", "failed_elements_total", "energy_ratio",
           "ballistic_status", "lsdyna_valid", "joint_ml_eligible", "missing_target_reason"]
    )

    for sim_id in range(1, EXPECTED_CASES + 1):
        ls = lsdyna[sim_id]
        mc = mcnp[sim_id]
        statuses_ok = all(
            normalize_status(mc[key]) in {"SUCESSO", "SUCCESS", "COMPLETED"}
            for key in ("status_geracao_input", "status_execucao_mcnp", "status_extracao")
        )
        if not statuses_ok:
            mcnp_failures.append(sim_id)
        rate = parse_float(mc["dose_efetiva_mSv_seg"], f"MCNP dose sim {sim_id}")
        rel_error = parse_float(mc["erro_relativo"], f"MCNP erro relativo sim {sim_id}")
        if rate < 0.0 or rel_error < 0.0:
            raise ValueError(f"MCNP: dose ou erro relativo negativo na simulação {sim_id}")
        dose_5min = rate * DOSE_SECONDS
        dose_min, dose_max = min(dose_min, rate), max(dose_max, rate)
        relative_error_min, relative_error_max = min(relative_error_min, rel_error), max(relative_error_max, rel_error)

        mcnp_materials = [x.strip() for x in mc["materiais"].split(",")]
        mcnp_thicknesses = [parse_float(x.strip(), f"espessura MCNP sim {sim_id}") for x in mc["espessuras"].split(",")]
        if len(mcnp_materials) != 7 or len(mcnp_thicknesses) != 7:
            raise ValueError(f"MCNP: sete camadas esperadas na simulação {sim_id}")
        for index in range(7):
            if mcnp_materials[index] != ls[layer_materials[index]]:
                material_mismatches.append(sim_id)
            difference = abs(mcnp_thicknesses[index] - parse_float(ls[layer_thicknesses[index]], "espessura LS-DYNA"))
            max_thickness_diff_cm = max(max_thickness_diff_cm, difference)
        if max_thickness_diff_cm > 0.0050001:
            raise ValueError("MCNP x LS-DYNA: diferença de espessura superior ao arredondamento documentado")

        lsdyna_valid = ls["dataset_role"].strip().upper() == "VALID_FOR_ANALYSIS"
        if not lsdyna_valid and sim_id not in PENDING_LSDYNA_IDS:
            raise ValueError(f"LS-DYNA: simulação não válida fora da lista documentada: {sim_id}")
        if lsdyna_valid and sim_id in PENDING_LSDYNA_IDS:
            raise ValueError(f"LS-DYNA: simulação pendente marcada como válida: {sim_id}")
        missing_reason = "" if lsdyna_valid else "LS-DYNA: etiqueta balística não validada; não usar como alvo ML"

        master = dict(ls)
        master.update({
            "mcnp_status_geracao_input": mc["status_geracao_input"],
            "mcnp_status_execucao": mc["status_execucao_mcnp"],
            "mcnp_status_extracao": mc["status_extracao"],
            "mcnp_dose_efetiva_mSv_seg": format_number(rate),
            "mcnp_dose_efetiva_5min_mSv": format_number(dose_5min),
            "mcnp_dose_5min_method": "mcnp_dose_efetiva_mSv_seg_x_300_s",
            "mcnp_erro_relativo": format_number(rel_error),
            "mcnp_valid": str(statuses_ok).upper(),
            "mcnp_source_result_csv": str(mcnp_path),
            "mcnp_design_precision_note": "camadas exportadas com 2 casas; projeto canônico preservado do mestre LS-DYNA",
            "lsdyna_target_valid": str(lsdyna_valid).upper(),
            "joint_ml_eligible": str(statuses_ok and lsdyna_valid).upper(),
            "joint_record_status": "VALID_BOTH_TARGETS" if statuses_ok and lsdyna_valid else "VALID_MCNP_LSDYNA_LABEL_PENDING",
            "joint_missing_target_reason": missing_reason,
        })
        master_rows.append(master)

        ml = {field: "" for field in ml_fields}
        ml["sim_id"] = str(sim_id)
        for field in layer_materials + layer_thicknesses + ["peso_total_kg_csv", "custo_total_unidade_csv", "massa_blindagem_modelo_final_kg"]:
            ml[field] = ls.get(field, "")
        ml.update({
            "dose_efetiva_mSv_seg": format_number(rate),
            "dose_efetiva_5min_mSv": format_number(dose_5min),
            "erro_relativo_mcnp": format_number(rel_error),
            "mcnp_valid": str(statuses_ok).upper(),
            "lsdyna_valid": str(lsdyna_valid).upper(),
            "joint_ml_eligible": str(statuses_ok and lsdyna_valid).upper(),
            "missing_target_reason": missing_reason,
        })
        if lsdyna_valid:
            ml.update({
                "max_penetration_mm": ls["working_max_penetration_mm"],
                "final_fsp_vx_m_s": ls["working_final_vx_m_s"],
                "failed_elements_total": ls["working_failed_elements_total"],
                "energy_ratio": ls["working_energy_ratio_final"],
                "ballistic_status": ls["working_ballistic_status"],
            })
        else:
            missing_rows.append({
                "sim_id": str(sim_id),
                "missing_system": "LS-DYNA",
                "missing_fields": "max_penetration_mm; final_fsp_vx_m_s; failed_elements_total; energy_ratio; ballistic_status",
                "reason": missing_reason,
                "action": "Manter fora do treino/validação do alvo balístico; não imputar.",
            })
        ml_rows.append(ml)

        dose_rows.append({key: ml[key] for key in ml_fields if key not in {
            "max_penetration_mm", "final_fsp_vx_m_s", "failed_elements_total", "energy_ratio",
            "ballistic_status", "lsdyna_valid", "joint_ml_eligible", "missing_target_reason"
        }})
        if lsdyna_valid:
            penetration_rows.append(dict(ml))

    if mcnp_failures:
        raise ValueError(f"MCNP: estados não concluídos nos IDs {mcnp_failures[:20]}")
    valid_count = sum(1 for row in ml_rows if row["lsdyna_valid"] == "TRUE")
    if valid_count != EXPECTED_LSDYNA_VALID:
        raise ValueError(f"LS-DYNA: esperados {EXPECTED_LSDYNA_VALID} válidos; encontrado {valid_count}")
    if {parse_int(row["sim_id"], "missing report") for row in missing_rows} != PENDING_LSDYNA_IDS:
        raise ValueError("Relatório de ausências não coincide com os quatro IDs pendentes")

    master_fields = lsdyna_columns + [
        "mcnp_status_geracao_input", "mcnp_status_execucao", "mcnp_status_extracao",
        "mcnp_dose_efetiva_mSv_seg", "mcnp_dose_efetiva_5min_mSv", "mcnp_dose_5min_method",
        "mcnp_erro_relativo", "mcnp_valid", "mcnp_source_result_csv", "mcnp_design_precision_note",
        "lsdyna_target_valid", "joint_ml_eligible", "joint_record_status", "joint_missing_target_reason",
    ]
    dose_fields = [field for field in ml_fields if field not in {
        "max_penetration_mm", "final_fsp_vx_m_s", "failed_elements_total", "energy_ratio",
        "ballistic_status", "lsdyna_valid", "joint_ml_eligible", "missing_target_reason"
    }]
    missing_fields = ["sim_id", "missing_system", "missing_fields", "reason", "action"]
    outputs = {
        "master": output_dir / "DATASET_MESTRE_CONJUNTO_MCNP_LSDYNA_2000_PROVISORIO.csv",
        "ml_joint": output_dir / "DATASET_ML_CONJUNTO_MCNP_LSDYNA_2000_PROVISORIO.csv",
        "ml_dose": output_dir / "DATASET_ML_DOSE_MCNP_2000.csv",
        "ml_penetration": output_dir / "DATASET_ML_PENETRACAO_LSDYNA_1996_PROVISORIO.csv",
        "missing": output_dir / "MISSING_DATA_REPORT_CONJUNTO.csv",
        "manifest": output_dir / "MANIFESTO_CONSOLIDACAO_MCNP_LSDYNA.json",
        "audit": output_dir / "AUDITORIA_CONSOLIDACAO_MCNP_LSDYNA.json",
    }
    write_csv(outputs["master"], master_rows, master_fields)
    write_csv(outputs["ml_joint"], ml_rows, ml_fields)
    write_csv(outputs["ml_dose"], dose_rows, dose_fields)
    write_csv(outputs["ml_penetration"], penetration_rows, ml_fields)
    write_csv(outputs["missing"], missing_rows, missing_fields)

    hashes = {
        "lsdyna_master_source_sha256": sha256(lsdyna_path),
        "mcnp_primary_results_source_sha256": sha256(mcnp_path),
        **{f"{key}_sha256": sha256(path) for key, path in outputs.items() if key not in {"manifest", "audit"}},
    }
    manifest = {
        "created_utc": utc_now(),
        "campaign": "Tese COPPE/UFRJ - Fase 2 - MCNP + LS-DYNA",
        "status": "PROVISORIO_COM_RASTREABILIDADE",
        "sources": {
            "lsdyna_master": str(lsdyna_path),
            "mcnp_primary_results": str(mcnp_path),
        },
        "canonical_design": "Entradas do mestre LS-DYNA; campos de espessura no resultado MCNP são exportados com 2 casas decimais.",
        "mcnp_target": {
            "name": "dose_efetiva_5min_mSv",
            "source_rate": "dose_efetiva_mSv_seg",
            "calculation": "dose_efetiva_mSv_seg * 300.0",
            "unit": "mSv em 5 min",
            "legacy_direct_5min_column": "Não utilizada: CSV enxuto apresentou escala incompatível com a taxa e a metodologia congelada.",
        },
        "lsdyna_target": {
            "name": "max_penetration_mm",
            "valid_labels": EXPECTED_LSDYNA_VALID,
            "pending_ids": sorted(PENDING_LSDYNA_IDS),
            "policy": "Os quatro casos não recebem imputação e não entram em treino/validação do alvo balístico.",
        },
        "outputs": {key: str(path) for key, path in outputs.items()},
        "sha256": hashes,
    }
    with outputs["manifest"].open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    audit = {
        "created_utc": utc_now(),
        "rows_master_joint": len(master_rows),
        "rows_ml_joint": len(ml_rows),
        "rows_ml_dose": len(dose_rows),
        "rows_ml_penetration": len(penetration_rows),
        "master_unique_ids": len({row["sim_id"] for row in master_rows}),
        "ml_joint_unique_ids": len({row["sim_id"] for row in ml_rows}),
        "mcnp_status_failures": mcnp_failures,
        "mcnp_rate_mSv_s_range": [dose_min, dose_max],
        "mcnp_relative_error_range": [relative_error_min, relative_error_max],
        "mcnp_lsdyna_material_position_mismatches": sorted(set(material_mismatches)),
        "mcnp_lsdyna_max_abs_thickness_difference_cm": max_thickness_diff_cm,
        "lsdyna_pending_ids": sorted(PENDING_LSDYNA_IDS),
        "critical_empty_fields": {
            "dose_efetiva_5min_mSv": sum(not row["dose_efetiva_5min_mSv"] for row in ml_rows),
            "max_penetration_mm_among_lsdyna_valid": sum(
                row["lsdyna_valid"] == "TRUE" and not row["max_penetration_mm"] for row in ml_rows
            ),
        },
        "result": "PASS_WITH_DOCUMENTED_LSDYNA_PENDING_LABELS",
        "sha256": hashes,
    }
    with outputs["audit"].open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    return {"outputs": {key: str(path) for key, path in outputs.items()}, "audit": audit}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lsdyna-master", type=Path, required=True)
    parser.add_argument("--mcnp-results", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    result = build(parse_args())
    print(json.dumps(result, ensure_ascii=False, indent=2))
