"""Audit and consolidate the final MCNP 6.3 campaign for the thesis.

This post-processing program never alters MCNP inputs, LS-DYNA files, or
solver outputs.  It reads the final MCNP6.3 campaign, verifies each of the
2,000 output files against its frozen input manifest, extracts the cell-20
tally-4 response, compares it to the preceding MCNP campaign, and writes
traceable final CSV/JSON datasets.

The MCNP response is the tally value already multiplied by ``FM4`` and is
therefore the effective-dose rate in mSv/s.  The five-minute target is
computed deterministically as ``rate * 300 s``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean, median
from typing import Any


EXPECTED_IDS = set(range(1, 2001))
PENDING_LSDYNA_IDS = {733, 832, 1277, 1567}
DOSE_SECONDS = 300.0
CONFORMAL_Q95_MSV = 0.3103037207640513
FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
CELL20_TALLY = re.compile(
    rf"cell\s+20\s*\n\s*({FLOAT})\s+({FLOAT})", re.IGNORECASE
)
NPS_100M = re.compile(
    r"run terminated when\s+100000000\s+particle histories were done", re.IGNORECASE
)


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


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def as_int(value: str, label: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: inteiro inválido {value!r}") from exc


def as_float(value: str, label: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: número inválido {value!r}") from exc


def format_number(value: float) -> str:
    return format(value, ".12g")


def index_rows(rows: list[dict[str, str]], key: str, label: str) -> dict[int, dict[str, str]]:
    indexed: dict[int, dict[str, str]] = {}
    for row in rows:
        sim_id = as_int(row.get(key, ""), f"{label}.{key}")
        if sim_id in indexed:
            raise ValueError(f"{label}: sim_id duplicado: {sim_id}")
        indexed[sim_id] = row
    if set(indexed) != EXPECTED_IDS:
        raise ValueError(
            f"{label}: cobertura inválida; ausentes={sorted(EXPECTED_IDS-set(indexed))[:20]}, "
            f"extras={sorted(set(indexed)-EXPECTED_IDS)[:20]}"
        )
    return indexed


def parse_output(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="latin-1", errors="replace")
    text_lower = text.lower()
    fatal = "fatal error" in text_lower or "bad trouble" in text_lower
    source_failure = "you need a source subroutine" in text_lower
    nps_100m = bool(NPS_100M.search(text))
    tally_4_present = bool(re.search(r"1tally\s+4\s+nps\s*=", text, re.IGNORECASE))
    tally_match = CELL20_TALLY.search(text)
    rate = float(tally_match.group(1)) if tally_match else None
    relative_error = float(tally_match.group(2)) if tally_match else None
    explicit_bin_pass = "passed all bin error check" in text_lower
    explicit_tfc_fail = "tally        4 tfc bin did not pass" in text_lower
    valid = bool(
        not fatal
        and not source_failure
        and nps_100m
        and tally_4_present
        and rate is not None
        and relative_error is not None
        and relative_error < 0.10
        and explicit_bin_pass
        and not explicit_tfc_fail
    )
    flags: list[str] = []
    if fatal:
        flags.append("FATAL_OR_BAD_TROUBLE")
    if source_failure:
        flags.append("SOURCE_FAILURE")
    if not nps_100m:
        flags.append("NPS_100M_NOT_REACHED")
    if not tally_4_present:
        flags.append("TALLY4_ABSENT")
    if tally_match is None:
        flags.append("CELL20_TALLY_VALUE_ABSENT")
    if relative_error is not None and relative_error >= 0.10:
        flags.append("RELATIVE_ERROR_GE_0_10")
    if not explicit_bin_pass:
        flags.append("BIN_ERROR_QA_NOT_PASSED")
    if explicit_tfc_fail:
        flags.append("TFC_QA_FAILED")
    return {
        "valid": valid,
        "rate_mSv_s": rate,
        "relative_error": relative_error,
        "fatal_or_bad_trouble": fatal,
        "source_failure": source_failure,
        "nps_100m": nps_100m,
        "tally_4_present": tally_4_present,
        "tfc_bin_error_pass": explicit_bin_pass and not explicit_tfc_fail,
        "flags": ";".join(flags),
    }


def audit_campaign(campaign_root: Path, output_dir: Path) -> tuple[dict[int, dict[str, str]], list[dict[str, Any]], dict[str, Any]]:
    manifest_path = campaign_root / "input_manifest.csv"
    manifest_rows, manifest_fields = read_csv(manifest_path)
    expected = index_rows(manifest_rows, "sim_id", "input_manifest")
    audit_rows: list[dict[str, Any]] = []
    extracted_rows: list[dict[str, str]] = []
    bad_ids: list[int] = []
    input_hash_mismatch: list[int] = []

    for sim_id in sorted(EXPECTED_IDS):
        spec = expected[sim_id]
        input_path = campaign_root / spec["input_relpath"]
        output_path = input_path.parent / f"output_{sim_id}.o"
        input_exists = input_path.is_file()
        output_exists = output_path.is_file()
        actual_hash = sha256(input_path) if input_exists else ""
        hash_match = input_exists and actual_hash.lower() == spec["input_sha256"].strip().lower()
        parsed = parse_output(output_path) if output_exists else {
            "valid": False, "rate_mSv_s": None, "relative_error": None,
            "fatal_or_bad_trouble": False, "source_failure": False,
            "nps_100m": False, "tally_4_present": False,
            "tfc_bin_error_pass": False, "flags": "OUTPUT_ABSENT",
        }
        if not hash_match:
            input_hash_mismatch.append(sim_id)
        if not parsed["valid"] or not hash_match:
            bad_ids.append(sim_id)
        material_list = ", ".join(spec[f"material_camada_{layer}"] for layer in range(1, 8))
        thickness_list = ", ".join(spec[f"espessura_camada_{layer}"] for layer in range(1, 8))
        status = "sucesso" if parsed["valid"] and hash_match else "falha_auditoria"
        extracted_rows.append({
            "id_simulacao": str(sim_id),
            "status_geracao_input": status,
            "status_execucao_mcnp": status,
            "status_extracao": status,
            "erro_geracao": "" if input_exists else "input.i ausente",
            "erro_execucao": "" if output_exists else "output MCNP ausente",
            "erro_extracao": parsed["flags"],
            "dose_efetiva_mSv_seg": "" if parsed["rate_mSv_s"] is None else format_number(parsed["rate_mSv_s"]),
            "erro_relativo": "" if parsed["relative_error"] is None else format_number(parsed["relative_error"]),
            "arquivo_input": str(input_path),
            "arquivo_output": str(output_path),
            "materiais": material_list,
            "espessuras": thickness_list,
        })
        audit_rows.append({
            "sim_id": sim_id,
            "input_relpath": spec["input_relpath"],
            "output_relpath": str(output_path.relative_to(campaign_root)),
            "input_exists": input_exists,
            "input_sha256_expected": spec["input_sha256"],
            "input_sha256_actual": actual_hash,
            "input_hash_match": hash_match,
            "output_exists": output_exists,
            "dose_efetiva_mSv_seg": "" if parsed["rate_mSv_s"] is None else format_number(parsed["rate_mSv_s"]),
            "erro_relativo": "" if parsed["relative_error"] is None else format_number(parsed["relative_error"]),
            "nps_100m": parsed["nps_100m"],
            "tally_4_present": parsed["tally_4_present"],
            "tfc_bin_error_pass": parsed["tfc_bin_error_pass"],
            "fatal_or_bad_trouble": parsed["fatal_or_bad_trouble"],
            "source_failure": parsed["source_failure"],
            "validation_flags": parsed["flags"],
            "validation_status": "PASS" if parsed["valid"] and hash_match else "FAIL",
        })

    result_csv = output_dir / "MCNP63_RESULTADOS_2000_FINAL.csv"
    result_fields = list(extracted_rows[0])
    write_csv(result_csv, extracted_rows, result_fields)
    audit_csv = output_dir / "AUDITORIA_CAMPANHA_FINAL_MCNP63_2000.csv"
    write_csv(audit_csv, audit_rows, list(audit_rows[0]))
    rates = [as_float(row["dose_efetiva_mSv_seg"], "dose") for row in extracted_rows if row["dose_efetiva_mSv_seg"]]
    errors = [as_float(row["erro_relativo"], "erro") for row in extracted_rows if row["erro_relativo"]]
    audit = {
        "created_utc": utc_now(),
        "campaign_root": str(campaign_root),
        "input_manifest_sha256": sha256(manifest_path),
        "expected_cases": 2000,
        "inputs_present": sum(row["input_exists"] for row in audit_rows),
        "outputs_present": sum(row["output_exists"] for row in audit_rows),
        "input_hash_mismatch_ids": input_hash_mismatch,
        "valid_cases": sum(row["validation_status"] == "PASS" for row in audit_rows),
        "failed_cases": sorted(set(bad_ids)),
        "dose_rate_mSv_s_range": [min(rates), max(rates)] if rates else None,
        "relative_error_range": [min(errors), max(errors)] if errors else None,
        "quality_criteria": {
            "normal_termination": "run terminated when 100000000 particle histories were done",
            "tally": "tally 4 / cell 20 present",
            "relative_error": "< 0.10",
            "statistical_checks": "passed all bin error check and no tally-4 TFC failure",
            "input_integrity": "SHA256 equals input_manifest.csv",
        },
        "outputs": {
            "results_csv": str(result_csv), "audit_csv": str(audit_csv),
            "results_csv_sha256": sha256(result_csv), "audit_csv_sha256": sha256(audit_csv),
        },
        "result": "PASS" if not bad_ids else "FAIL",
    }
    write_json(output_dir / "AUDITORIA_CAMPANHA_FINAL_MCNP63_2000.json", audit)
    return expected, extracted_rows, audit


def compare_legacy(current: list[dict[str, str]], legacy_path: Path, output_dir: Path) -> dict[str, Any]:
    legacy_rows, _ = read_csv(legacy_path)
    legacy = index_rows(legacy_rows, "id_simulacao", "resultados MCNP anteriores")
    comparison: list[dict[str, Any]] = []
    absolute_errors: list[float] = []
    signed_errors: list[float] = []
    percent_errors: list[float] = []
    for row in current:
        sim_id = as_int(row["id_simulacao"], "resultado atual")
        new_rate = as_float(row["dose_efetiva_mSv_seg"], f"novo {sim_id}")
        old_rate = as_float(legacy[sim_id]["dose_efetiva_mSv_seg"], f"anterior {sim_id}")
        signed = new_rate - old_rate
        absolute = abs(signed)
        percent = 100.0 * signed / old_rate if old_rate != 0 else None
        absolute_errors.append(absolute)
        signed_errors.append(signed)
        if percent is not None:
            percent_errors.append(percent)
        comparison.append({
            "sim_id": sim_id,
            "dose_anterior_mSv_seg": format_number(old_rate),
            "dose_mcnp63_final_mSv_seg": format_number(new_rate),
            "diferenca_assinada_mSv_seg": format_number(signed),
            "diferenca_absoluta_mSv_seg": format_number(absolute),
            "diferenca_percentual_vs_anterior": "" if percent is None else format_number(percent),
            "dose_anterior_5min_mSv": format_number(old_rate * DOSE_SECONDS),
            "dose_mcnp63_final_5min_mSv": format_number(new_rate * DOSE_SECONDS),
            "erro_relativo_anterior": legacy[sim_id].get("erro_relativo", ""),
            "erro_relativo_mcnp63_final": row["erro_relativo"],
        })
    comparison.sort(key=lambda row: float(row["diferenca_absoluta_mSv_seg"]), reverse=True)
    comparison_path = output_dir / "COMPARACAO_MCNP63_VS_CAMPANHA_ANTERIOR_2000.csv"
    write_csv(comparison_path, comparison, list(comparison[0]))
    report = {
        "created_utc": utc_now(),
        "legacy_results": str(legacy_path),
        "legacy_results_sha256": sha256(legacy_path),
        "n_cases": len(comparison),
        "difference_rate_mSv_s": {
            "mean_signed": mean(signed_errors), "median_signed": median(signed_errors),
            "mae": mean(absolute_errors), "median_absolute": median(absolute_errors),
            "max_absolute": max(absolute_errors),
        },
        "difference_percent_vs_legacy": {
            "mean_signed": mean(percent_errors), "median_signed": median(percent_errors),
            "max_absolute": max(abs(value) for value in percent_errors),
        },
        "top_10_absolute_differences": comparison[:10],
        "interpretation": (
            "Diferenças refletem a comparação entre duas campanhas Monte Carlo e devem ser "
            "interpretadas junto aos erros relativos de cada tally; este relatório não altera "
            "resultados nem imputa valores."
        ),
        "comparison_csv": str(comparison_path),
        "comparison_csv_sha256": sha256(comparison_path),
    }
    write_json(output_dir / "RELATORIO_DISCREPANCIAS_MCNP63_VS_ANTERIOR.json", report)
    return report


def build_final_datasets(
    provisional_dir: Path, current: list[dict[str, str]], output_dir: Path
) -> dict[str, Any]:
    master_path = provisional_dir / "DATASET_MESTRE_CONJUNTO_MCNP_LSDYNA_2000_PROVISORIO.csv"
    master_rows, master_fields = read_csv(master_path)
    master = index_rows(master_rows, "sim_id", "mestre conjunto provisório")
    current_index = index_rows(current, "id_simulacao", "resultados MCNP63")
    final_master: list[dict[str, str]] = []
    ml_dose: list[dict[str, str]] = []
    ml_joint: list[dict[str, str]] = []
    missing: list[dict[str, str]] = []
    feature_fields = ["sim_id"] + [f"material_camada_{i}" for i in range(1, 8)] + [
        f"espessura_camada_{i}_cm" for i in range(1, 8)
    ] + ["peso_total_kg_csv", "custo_total_unidade_csv", "massa_blindagem_modelo_final_kg"]
    dose_fields = feature_fields + ["dose_efetiva_mSv_seg", "dose_efetiva_5min_mSv", "erro_relativo_mcnp", "mcnp_valid"]
    joint_fields = dose_fields + [
        "max_penetration_mm", "final_fsp_vx_m_s", "failed_elements_total", "energy_ratio",
        "ballistic_status", "lsdyna_valid", "joint_ml_eligible", "missing_target_reason",
    ]
    for sim_id in sorted(EXPECTED_IDS):
        row = dict(master[sim_id])
        result = current_index[sim_id]
        rate = as_float(result["dose_efetiva_mSv_seg"], f"MCNP63 dose {sim_id}")
        rel = as_float(result["erro_relativo"], f"MCNP63 erro {sim_id}")
        row.update({
            "mcnp_status_geracao_input": "sucesso",
            "mcnp_status_execucao": "sucesso",
            "mcnp_status_extracao": "sucesso",
            "mcnp_dose_efetiva_mSv_seg": format_number(rate),
            "mcnp_dose_efetiva_5min_mSv": format_number(rate * DOSE_SECONDS),
            "mcnp_dose_5min_method": "mcnp_dose_efetiva_mSv_seg_x_300_s",
            "mcnp_erro_relativo": format_number(rel),
            "mcnp_valid": "TRUE",
            "mcnp_source_result_csv": str(output_dir / "MCNP63_RESULTADOS_2000_FINAL.csv"),
            "mcnp_design_precision_note": "Projeto canônico do mestre LS-DYNA; campanha final MCNP6.3, template congelado e frações normalizadas sem alteração de proporções.",
        })
        lsdyna_valid = row.get("dataset_role", "").strip().upper() == "VALID_FOR_ANALYSIS"
        missing_reason = "" if lsdyna_valid else "LS-DYNA: etiqueta balística não validada; não usar como alvo ML"
        row.update({
            "lsdyna_target_valid": str(lsdyna_valid).upper(),
            "joint_ml_eligible": str(lsdyna_valid).upper(),
            "joint_record_status": "VALID_BOTH_TARGETS" if lsdyna_valid else "VALID_MCNP_LSDYNA_LABEL_PENDING",
            "joint_missing_target_reason": missing_reason,
        })
        final_master.append(row)
        dose_record = {field: row.get(field, "") for field in feature_fields}
        dose_record.update({
            "dose_efetiva_mSv_seg": row["mcnp_dose_efetiva_mSv_seg"],
            "dose_efetiva_5min_mSv": row["mcnp_dose_efetiva_5min_mSv"],
            "erro_relativo_mcnp": row["mcnp_erro_relativo"],
            "mcnp_valid": "TRUE",
        })
        ml_dose.append(dose_record)
        if lsdyna_valid:
            joint = dict(dose_record)
            joint.update({
                "max_penetration_mm": row["working_max_penetration_mm"],
                "final_fsp_vx_m_s": row["working_final_vx_m_s"],
                "failed_elements_total": row["working_failed_elements_total"],
                "energy_ratio": row["working_energy_ratio_final"],
                "ballistic_status": row["working_ballistic_status"],
                "lsdyna_valid": "TRUE", "joint_ml_eligible": "TRUE", "missing_target_reason": "",
            })
            ml_joint.append(joint)
        else:
            missing.append({
                "sim_id": str(sim_id), "missing_system": "LS-DYNA",
                "missing_fields": "max_penetration_mm; final_fsp_vx_m_s; failed_elements_total; energy_ratio; ballistic_status",
                "reason": missing_reason,
                "action": "Manter fora do treino/validação do alvo balístico; não imputar.",
            })
    invalid_lsdyna = {int(row["sim_id"]) for row in missing}
    if invalid_lsdyna != PENDING_LSDYNA_IDS:
        raise ValueError(f"IDs LS-DYNA pendentes divergentes: {sorted(invalid_lsdyna)}")
    if len(ml_joint) != 1996:
        raise ValueError(f"Esperados 1996 registros com ambos os alvos, recebido {len(ml_joint)}")
    names = {
        "master": "DATASET_MESTRE_CONJUNTO_MCNP63_LSDYNA_2000_FINAL.csv",
        "ml_dose": "DATASET_ML_DOSE_MCNP63_2000_FINAL.csv",
        "ml_joint": "DATASET_ML_CONJUNTO_MCNP63_LSDYNA_1996_FINAL.csv",
        "ml_penetration": "DATASET_ML_PENETRACAO_LSDYNA_1996_FINAL.csv",
        "missing": "MISSING_DATA_REPORT_CONJUNTO_FINAL.csv",
    }
    paths = {key: output_dir / name for key, name in names.items()}
    write_csv(paths["master"], final_master, master_fields)
    write_csv(paths["ml_dose"], ml_dose, dose_fields)
    write_csv(paths["ml_joint"], ml_joint, joint_fields)
    write_csv(paths["ml_penetration"], ml_joint, joint_fields)
    write_csv(paths["missing"], missing, list(missing[0]))
    hashes = {key: sha256(path) for key, path in paths.items()}
    manifest = {
        "created_utc": utc_now(),
        "status": "FINAL_MCNP63__LSDYNA_TARGET_SPECIFIC",
        "scope": "2000 respostas MCNP6.3 validadas; 1996 etiquetas balísticas LS-DYNA validadas; quatro exclusões balísticas documentadas sem imputação.",
        "outputs": {key: str(path) for key, path in paths.items()},
        "sha256": hashes,
        "source_provisional_master": str(master_path),
        "source_provisional_master_sha256": sha256(master_path),
        "mcnp_target": "dose_efetiva_5min_mSv = dose_efetiva_mSv_seg * 300",
        "lsdyna_pending_ids": sorted(PENDING_LSDYNA_IDS),
        "row_counts": {"master": 2000, "ml_dose": 2000, "ml_joint": 1996, "ml_penetration": 1996},
    }
    manifest_path = output_dir / "MANIFESTO_DATASET_CONJUNTO_MCNP63_LSDYNA_FINAL.json"
    write_json(manifest_path, manifest)
    audit = {
        "created_utc": utc_now(), "result": "PASS_WITH_DOCUMENTED_LSDYNA_EXCLUSIONS",
        "master_rows": len(final_master), "master_unique_ids": len({row["sim_id"] for row in final_master}),
        "ml_dose_rows": len(ml_dose), "ml_joint_rows": len(ml_joint),
        "mcnp_critical_empty": sum(not row["dose_efetiva_5min_mSv"] for row in ml_dose),
        "lsdyna_critical_empty_among_valid": sum(not row["max_penetration_mm"] for row in ml_joint),
        "pending_lsdyna_ids": sorted(PENDING_LSDYNA_IDS), "sha256": hashes,
    }
    write_json(output_dir / "AUDITORIA_DATASET_CONJUNTO_MCNP63_LSDYNA_FINAL.json", audit)
    return {"paths": paths, "manifest": manifest, "audit": audit}


def validate_high_fidelity(hf_root: Path, phase5_root: Path) -> dict[str, Any]:
    prepared_path = phase5_root / "preparation_manifest.csv"
    prepared_rows, _ = read_csv(prepared_path)
    selected = [row for row in prepared_rows if 9001 <= as_int(row["validation_sim_id"], "validation") <= 9005]
    if len(selected) != 5:
        raise ValueError("A preparação da Fase 5 não contém exatamente os casos 9001–9005")
    evidence_root = phase5_root / "mcnp" / "EVIDENCIAS_MCNP63_9001_9005"
    validation_rows: list[dict[str, Any]] = []
    for prepared in selected:
        sim_id = as_int(prepared["validation_sim_id"], "validation_sim_id")
        source_dir = hf_root / f"sim_{sim_id}"
        input_path, output_path = source_dir / "input.i", source_dir / f"output_{sim_id}.o"
        if not input_path.is_file() or not output_path.is_file():
            raise FileNotFoundError(f"Evidência MCNP ausente para {sim_id}: {source_dir}")
        parsed = parse_output(output_path)
        if not parsed["valid"]:
            raise ValueError(f"Validação MCNP {sim_id} não passou: {parsed['flags']}")
        target_dir = evidence_root / f"sim_{sim_id}"
        target_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(input_path, target_dir / "input_normalizado.i")
        shutil.copy2(output_path, target_dir / f"output_{sim_id}.o")
        actual = float(parsed["rate_mSv_s"]) * DOSE_SECONDS
        predicted = as_float(prepared["pred_dose_5min_mSv"], f"dose prevista {sim_id}")
        error = actual - predicted
        validation_rows.append({
            "validation_candidate_id": prepared["validation_candidate_id"], "validation_sim_id": sim_id,
            "selection_reason": prepared["selection_reason"],
            "dose_pred_5min_mSv": format_number(predicted),
            "dose_pi_lower95_mSv": "", "dose_pi_upper95_mSv": "",
            "dose_hf_mcnp63_mSv_seg": format_number(float(parsed["rate_mSv_s"])),
            "dose_hf_5min_mSv": format_number(actual),
            "erro_relativo_hf_mcnp": format_number(float(parsed["relative_error"])),
            "dose_error_assinado_mSv": format_number(error),
            "dose_error_absoluto_mSv": format_number(abs(error)),
            "dentro_intervalo_conformal95": str(abs(error) <= CONFORMAL_Q95_MSV).upper(),
            "mcnp_output_valid": "TRUE", "input_normalizado_sha256": sha256(input_path),
            "output_sha256": sha256(output_path), "evidence_output": str(target_dir / f"output_{sim_id}.o"),
        })
    candidate_csv = phase5_root / "candidatos_validacao_alta_fidelidade_fase4.csv"
    candidates, _ = read_csv(candidate_csv)
    candidate_index = {as_int(row["validation_candidate_id"], "candidato"): row for row in candidates}
    for row in validation_rows:
        original = candidate_index[as_int(str(row["validation_candidate_id"]), "candidate")]
        row["dose_pi_lower95_mSv"] = original["dose_pi_lower95_mSv"]
        row["dose_pi_upper95_mSv"] = original["dose_pi_upper95_mSv"]
    dest = phase5_root / "VALIDACAO_MCNP63_ALTA_FIDELIDADE_9001_9005.csv"
    write_csv(dest, validation_rows, list(validation_rows[0]))
    abs_errors = [float(row["dose_error_absoluto_mSv"]) for row in validation_rows]
    report = {
        "created_utc": utc_now(), "scope": "Cinco pontos definidos pelo autor: 9001–9005.",
        "method": "MCNP6.3, 100 milhões de histórias, tally 4/célula 20; dose em 5 min = taxa * 300 s.",
        "all_outputs_valid": True, "n_points": 5,
        "mae_mSv": mean(abs_errors), "median_absolute_error_mSv": median(abs_errors),
        "max_absolute_error_mSv": max(abs_errors),
        "within_conformal_95_count": sum(row["dentro_intervalo_conformal95"] == "TRUE" for row in validation_rows),
        "conformal_q95_mSv": CONFORMAL_Q95_MSV,
        "csv": str(dest), "csv_sha256": sha256(dest), "evidence_root": str(evidence_root),
    }
    write_json(phase5_root / "AUDITORIA_VALIDACAO_MCNP63_9001_9005.json", report)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-root", required=True, type=Path)
    parser.add_argument("--legacy-results", required=True, type=Path)
    parser.add_argument("--provisional-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--hf-root", required=True, type=Path)
    parser.add_argument("--phase5-root", required=True, type=Path)
    parser.add_argument(
        "--release-final",
        action="store_true",
        help=(
            "Autoriza a escrita dos datasets com sufixo FINAL. Use somente após "
            "a decisão explícita do autor sobre o cenário físico MCNP e a validação "
            "da Fase 3/4 correspondente."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.campaign_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    _, current, campaign_audit = audit_campaign(root, output_dir)
    if campaign_audit["result"] != "PASS":
        raise RuntimeError("Auditoria da campanha não aprovada; dataset final não foi consolidado.")
    discrepancy = compare_legacy(current, args.legacy_results.resolve(), output_dir)
    validation = validate_high_fidelity(args.hf_root.resolve(), args.phase5_root.resolve())
    release_gate = {
        "release_requested": args.release_final,
        "high_fidelity_points_within_current_conformal_interval": validation["within_conformal_95_count"],
        "high_fidelity_points": validation["n_points"],
        "legacy_mean_percent_difference": discrepancy["difference_percent_vs_legacy"]["mean_signed"],
        "status": "RELEASED" if args.release_final else "HOLD_FOR_AUTHOR_MODEL_DECISION",
    }
    write_json(output_dir / "GATE_DE_LIBERACAO_MCNP63.json", release_gate)
    if not args.release_final:
        print(json.dumps({
            "campaign_result": campaign_audit["result"],
            "valid_cases": campaign_audit["valid_cases"],
            "release_status": release_gate["status"],
            "reason": "A campanha passou tecnicamente, mas exige decisão explícita sobre a incompatibilidade com o cenário histórico antes de gerar dados finais.",
        }, ensure_ascii=False, indent=2))
        return
    final_dataset = build_final_datasets(args.provisional_dir.resolve(), current, output_dir)
    summary = {
        "created_utc": utc_now(), "campaign_audit": campaign_audit,
        "discrepancy_report": discrepancy, "final_dataset": final_dataset["audit"],
        "high_fidelity_validation": validation, "release_gate": release_gate,
    }
    write_json(output_dir / "RELATORIO_FINALIZACAO_MCNP63.json", summary)
    print(json.dumps({
        "campaign_result": campaign_audit["result"],
        "valid_cases": campaign_audit["valid_cases"],
        "dataset_master": str(final_dataset["paths"]["master"]),
        "phase5_validation": validation["csv"],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
