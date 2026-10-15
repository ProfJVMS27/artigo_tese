# -*- coding: utf-8 -*-
"""
FASE 3 — SURROGATE MODELING PARA A TESE
========================================

Objetivo científico
-------------------
Treinar e selecionar modelos substitutos SOMENTE para respostas computacionalmente
caras:
    1) dose_efetiva_5min_mSv  <- MCNP (F4:P + FM4)
    2) max_penetration_mm     <- LS-DYNA

Massa e custo NÃO são targets de ML. Permanecem cálculos determinísticos.

Desenho experimental da Fase 3
-------------------------------
- Inputs: 7 materiais + 7 espessuras (DoE permutation–Dirichlet restrito).
- Comparação pareada: RAW vs PHYSICS-INFORMED features.
- Split específico por target: 70% treino / 15% calibração / 15% teste final.
- Dose usa DoE + MCNP válido.
- Penetração usa DoE + LS-DYNA válido.
- Seleção de família/modelo: 5-fold CV APENAS no treino.
- Tuning opcional do vencedor: RandomizedSearchCV APENAS no treino.
- Incerteza: split conformal calibration (95%) no conjunto de calibração.
- Teste final: usado uma única vez para desempenho e cobertura do intervalo.

Compatibilidade LS-DYNA
-----------------------
Este script usa os recortes target-specific do dataset conjunto consolidado:
    - dose MCNP: 2.000 rótulos válidos;
    - penetração LS-DYNA: 1.996 rótulos validados.

Uso
---
1) Coloque este arquivo e feature_engineering.py na mesma pasta.
2) Ajuste somente a seção CONFIGURAÇÃO / CAMINHOS se necessário.
3) Para a execução científica final, mantenha RUN_MODE = "FINAL".
4) Para teste de pipeline com dados parciais, use RUN_MODE = "SMOKE_TEST".
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import math
import os
import platform
import sys
import time
import warnings

import joblib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import wilcoxon

from sklearn.base import clone
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import (
    ExtraTreesRegressor,
    GradientBoostingRegressor,
    HistGradientBoostingRegressor,
    RandomForestRegressor,
)
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, RandomizedSearchCV, train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from feature_engineering_final import (
    MATERIALS,
    MATERIAL_COLS,
    THICKNESS_COLS,
    DENSITY_G_CM3,
    ShieldFeatureEngineer,
    physics_feature_names,
    ballistic_strength23_feature_names,
    validate_design_dataframe,
)
from tabelas_artigo_fase3 import generate_article_tables

try:
    import sklearn
except Exception:
    sklearn = None

try:
    import xgboost as xgb
    XGBOOST_AVAILABLE = True
except Exception:
    xgb = None
    XGBOOST_AVAILABLE = False


# =============================================================================
# 1. CONFIGURAÇÃO
# =============================================================================

SCRIPT_VERSION = "FASE3-THESIS-SURROGATES-v1.4-MCNP63-DEDF-AP-TARGET-SPECIFIC"
RANDOM_STATE = 42
N_FOLDS = 5
TEST_FRACTION = 0.15
CALIBRATION_FRACTION = 0.15
CONFORMAL_ALPHA = 0.05               # 95% prediction interval
RUN_MODE = "FINAL"                  # "FINAL" ou "SMOKE_TEST"
EXPECTED_CAMPAIGN_CASES = 2000
EXPECTED_LSDYNA_VALID_CASES = 1996
LSDYNA_PENDING_IDS = {733, 832, 1277, 1567}
MIN_SMOKE_CASES = 50
SMOKE_SAMPLE_SIZE = 160

# Cada candidato é otimizado formalmente no subconjunto de treino antes do
# ranking.  Calibração e teste permanecem inteiramente bloqueados.
TUNE_ALL_CANDIDATES = True
TUNE_WINNER = True  # compatibilidade com o modo smoke-test
TUNING_ITERATIONS = 20

# MCNP: não removemos automaticamente resultados por erro relativo; apenas auditamos.
MCNP_RELATIVE_ERROR_REVIEW_THRESHOLD = 0.05

BASE_DIR = Path(__file__).resolve().parent

PHASE2_DATASET_DIR = (
    BASE_DIR.parent / "Fase 2" / "DATASET_CONJUNTO_MCNP63_DEDF_AP_LSDYNA_FINAL_20261002"
)
DOSE_DATASET_PATH = PHASE2_DATASET_DIR / "DATASET_ML_DOSE_MCNP63_DEDF_AP_2000_FINAL.csv"
PENETRATION_DATASET_PATH = (
    PHASE2_DATASET_DIR / "DATASET_ML_PENETRACAO_LSDYNA_1996_FINAL.csv"
)

OUTPUT_DIR = BASE_DIR / "analise_desempenho_modelos (Fase 3) - MCNP63_DEDF_AP"
MODEL_DIR = BASE_DIR / "modelos_treinados (Fase 3) - MCNP63_DEDF_AP"
FIGURE_DIR = OUTPUT_DIR / "figuras"

for d in (OUTPUT_DIR, MODEL_DIR, FIGURE_DIR):
    d.mkdir(parents=True, exist_ok=True)

TARGET_DOSE = "dose_efetiva_5min_mSv"
TARGET_PEN = "max_penetration_mm"
TARGETS = [TARGET_DOSE, TARGET_PEN]

# Camada 1 = interna/operador; camada 7 = externa/fonte/impacto.
RAW_INPUT_COLUMNS = list(MATERIAL_COLS) + list(THICKNESS_COLS)


# =============================================================================
# 2. UTILITÁRIOS GERAIS
# =============================================================================


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def resolve_first_existing(candidates: list[Path], label: str) -> Path:
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(
        f"{label} não encontrado. Caminhos procurados:\n" +
        "\n".join(f"  - {p}" for p in candidates)
    )


def normalize_id_column(df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    out = df.copy()
    if "id_simulacao" in out.columns:
        pass
    elif "sim_id" in out.columns:
        out = out.rename(columns={"sim_id": "id_simulacao"})
    else:
        raise ValueError(f"{source_name}: coluna id_simulacao/sim_id ausente.")

    out["id_simulacao"] = pd.to_numeric(out["id_simulacao"], errors="raise").astype(int)
    if out["id_simulacao"].duplicated().any():
        dups = out.loc[out["id_simulacao"].duplicated(), "id_simulacao"].tolist()[:10]
        raise ValueError(f"{source_name}: IDs duplicados: {dups}")
    return out


def accepted_status(series: pd.Series, accepted: set[str]) -> pd.Series:
    normalized = series.astype(str).str.strip().str.upper()
    return normalized.isin({x.upper() for x in accepted})


def smape(y_true, y_pred) -> float:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    denom = np.abs(y_true) + np.abs(y_pred)
    mask = denom > 1e-15
    if not np.any(mask):
        return 0.0
    return float(100.0 * np.mean(2.0 * np.abs(y_pred[mask] - y_true[mask]) / denom[mask]))


def regression_metrics(y_true, y_pred) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=float).ravel()
    y_pred = np.asarray(y_pred, dtype=float).ravel()
    rmse = math.sqrt(mean_squared_error(y_true, y_pred))
    std_y = float(np.std(y_true, ddof=0))
    data_range = float(np.max(y_true) - np.min(y_true))
    return {
        "r2": float(r2_score(y_true, y_pred)),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(rmse),
        "nrmse_std": float(rmse / std_y) if std_y > 0 else np.nan,
        "nrmse_range": float(rmse / data_range) if data_range > 0 else np.nan,
        "smape_pct": smape(y_true, y_pred),
    }


def conformal_quantile(abs_residuals: np.ndarray, alpha: float = 0.05) -> float:
    """Quantil split-conformal com correção finita e método 'higher'."""
    r = np.sort(np.asarray(abs_residuals, dtype=float).ravel())
    if len(r) == 0:
        raise ValueError("Conjunto de calibração vazio.")
    level = math.ceil((len(r) + 1) * (1.0 - alpha)) / len(r)
    level = min(1.0, max(0.0, level))
    try:
        return float(np.quantile(r, level, method="higher"))
    except TypeError:  # numpy antigo
        return float(np.quantile(r, level, interpolation="higher"))


def dump_json(obj, path: Path) -> None:
    def converter(v):
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, (np.floating,)):
            return float(v)
        if isinstance(v, np.ndarray):
            return v.tolist()
        if isinstance(v, Path):
            return str(v)
        raise TypeError(f"Tipo não serializável: {type(v)}")

    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False, default=converter)


# =============================================================================
# 3. CARREGAMENTO E AUDITORIA DOS DADOS
# =============================================================================


def load_configs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = normalize_id_column(df, "DoE")
    validate_design_dataframe(df, MATERIAL_COLS, THICKNESS_COLS)

    # IDs únicos e dentro da campanha.
    if (df["id_simulacao"] <= 0).any():
        raise ValueError("DoE: IDs devem ser positivos.")

    # Mantém somente colunas necessárias nesta fase.
    keep = ["id_simulacao"] + RAW_INPUT_COLUMNS
    for optional in ["peso_total_kg", "custo_total_unidade"]:
        if optional in df.columns:
            keep.append(optional)
    return df[keep].copy()


def load_mcnp(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = normalize_id_column(df, "MCNP")

    # Filtragem somente por status explícito quando disponível.
    mask = pd.Series(True, index=df.index)
    if "status_execucao_mcnp" in df.columns:
        mask &= accepted_status(df["status_execucao_mcnp"], {"SUCESSO", "SUCCESS", "COMPLETED"})
    if "status_extracao" in df.columns:
        mask &= accepted_status(df["status_extracao"], {"SUCESSO", "SUCCESS", "COMPLETED"})
    df = df.loc[mask].copy()

    # Target: preferir explicitamente dose em 5 min. Se só houver taxa, multiplicar por 300 s.
    five_min_aliases = [
        "dose_efetiva_5min_mSv",
        "5 min de Dose efetiva (mSv)",
        "5_min_dose_efetiva_mSv",
    ]
    rate_aliases = [
        "dose_efetiva_mSv_seg",
        "dose_efetiva_mSv_s",
        "dose_rate_mSv_s",
    ]

    source_col = next((c for c in five_min_aliases if c in df.columns), None)
    dose_origin = None
    if source_col is not None:
        df[TARGET_DOSE] = pd.to_numeric(df[source_col], errors="coerce")
        dose_origin = f"coluna direta: {source_col}"
    else:
        rate_col = next((c for c in rate_aliases if c in df.columns), None)
        if rate_col is None:
            raise ValueError(
                "MCNP: não foi encontrada dose de 5 min nem taxa em mSv/s. "
                f"Colunas disponíveis: {df.columns.tolist()}"
            )
        rate = pd.to_numeric(df[rate_col], errors="coerce")
        df[TARGET_DOSE] = rate * 300.0
        dose_origin = f"{rate_col} × 300 s"

    if df[TARGET_DOSE].isna().any():
        raise ValueError(f"MCNP: {int(df[TARGET_DOSE].isna().sum())} doses inválidas/NaN.")
    if (df[TARGET_DOSE] < 0).any():
        raise ValueError("MCNP: dose negativa encontrada.")

    # Erro relativo é auditado, nunca usado como feature.
    rel_error_col = "erro_relativo" if "erro_relativo" in df.columns else None
    if rel_error_col:
        df[rel_error_col] = pd.to_numeric(df[rel_error_col], errors="coerce")
        n_review = int((df[rel_error_col] > MCNP_RELATIVE_ERROR_REVIEW_THRESHOLD).sum())
        if n_review:
            warnings.warn(
                f"MCNP: {n_review} casos com erro relativo > "
                f"{MCNP_RELATIVE_ERROR_REVIEW_THRESHOLD:.2%}. Mantidos para auditoria; "
                "não são usados como feature."
            )

    out_cols = ["id_simulacao", TARGET_DOSE]
    if rel_error_col:
        out_cols.append(rel_error_col)
    out = df[out_cols].copy()
    out.attrs["dose_origin"] = dose_origin
    return out


def load_consolidated_target_dataset(
    path: Path,
    target: str,
    expected_rows: int,
) -> tuple[pd.DataFrame, dict]:
    """Carrega um recorte ML consolidado e recupera os nomes canônicos da Fase 3."""
    if not path.exists():
        raise FileNotFoundError(f"Dataset consolidado não encontrado: {path}")
    df = pd.read_csv(path)
    if "sim_id" not in df.columns:
        raise ValueError(f"{path.name}: coluna sim_id ausente.")
    df = df.rename(columns={"sim_id": "id_simulacao"}).copy()

    rename_thickness = {
        f"espessura_camada_{i}_cm": f"espessura_camada_{i}"
        for i in range(1, 8)
    }
    missing = [
        col for col in ["id_simulacao", target] + list(MATERIAL_COLS)
        + list(rename_thickness) if col not in df.columns
    ]
    if missing:
        raise ValueError(f"{path.name}: colunas obrigatórias ausentes: {missing}")
    df = df.rename(columns=rename_thickness)
    df["id_simulacao"] = pd.to_numeric(df["id_simulacao"], errors="raise").astype(int)
    df[target] = pd.to_numeric(df[target], errors="raise")
    if len(df) != expected_rows:
        raise ValueError(
            f"{path.name}: {len(df)} linhas; esperadas {expected_rows} para {target}."
        )
    if df["id_simulacao"].duplicated().any():
        raise ValueError(f"{path.name}: IDs duplicados.")
    if df[target].isna().any() or (df[target] < 0.0).any():
        raise ValueError(f"{path.name}: target {target} inválido.")
    validate_design_dataframe(df, MATERIAL_COLS, THICKNESS_COLS)

    audit = {
        "source_path": str(path),
        "source_sha256": sha256_file(path),
        "rows_total": int(len(df)),
        "rows_valid": int(len(df)),
        "ids": df["id_simulacao"].astype(int).tolist(),
    }
    if target == TARGET_PEN:
        if "lsdyna_valid" not in df.columns or not df["lsdyna_valid"].astype(str).str.upper().eq("TRUE").all():
            raise ValueError(f"{path.name}: contém registro LS-DYNA não validado.")
        audit["ballistic_status_counts"] = {
            str(k): int(v) for k, v in df["ballistic_status"].value_counts(dropna=False).items()
        }
    if target == TARGET_DOSE:
        if "mcnp_valid" not in df.columns or not df["mcnp_valid"].astype(str).str.upper().eq("TRUE").all():
            raise ValueError(f"{path.name}: contém registro MCNP não validado.")
    return df.sort_values("id_simulacao").reset_index(drop=True), audit


def load_lsdyna(path: Path) -> tuple[pd.DataFrame, dict]:
    df_raw = pd.read_csv(path)
    df = normalize_id_column(df_raw, "LS-DYNA")

    required = [
        "campaign_status",
        "normal_termination",
        "parser_status",
        "validation_status",
        "max_penetration_mm",
        "perforation_bool",
        "ballistic_status",
    ]
    required += [f"material_camada_{i}" for i in range(1, 8)]
    required += [f"espessura_camada_{i}_cm" for i in range(1, 8)]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"LS-DYNA: colunas obrigatórias ausentes: {missing}")

    mask = (
        accepted_status(df["campaign_status"], {"COMPLETED"})
        & (pd.to_numeric(df["normal_termination"], errors="coerce") == 1)
        & accepted_status(df["parser_status"], {"OK"})
        & accepted_status(df["validation_status"], {"PASS"})
    )

    valid = df.loc[mask].copy()
    invalid = df.loc[~mask].copy()

    valid[TARGET_PEN] = pd.to_numeric(valid["max_penetration_mm"], errors="coerce")
    valid["perforation_bool"] = pd.to_numeric(valid["perforation_bool"], errors="coerce")

    if valid[TARGET_PEN].isna().any():
        raise ValueError("LS-DYNA: max_penetration_mm NaN em caso marcado válido.")
    if (valid[TARGET_PEN] < 0).any():
        raise ValueError("LS-DYNA: max_penetration_mm negativo.")

    audit = {
        "rows_total": int(len(df)),
        "rows_valid": int(len(valid)),
        "rows_invalid": int(len(invalid)),
        "invalid_ids_first_50": invalid["id_simulacao"].astype(int).tolist()[:50],
        "ballistic_status_counts": {
            str(k): int(v) for k, v in valid["ballistic_status"].value_counts(dropna=False).items()
        },
        "perforation_counts": {
            str(k): int(v) for k, v in valid["perforation_bool"].value_counts(dropna=False).items()
        },
    }

    keep = [
        "id_simulacao",
        TARGET_PEN,
        "perforation_bool",
        "ballistic_status",
        "campaign_status",
        "parser_status",
        "validation_status",
    ]
    keep += [f"material_camada_{i}" for i in range(1, 8)]
    keep += [f"espessura_camada_{i}_cm" for i in range(1, 8)]
    return valid[keep].copy(), audit


def verify_lsdyna_design_matches_doe(configs: pd.DataFrame, lsdyna: pd.DataFrame) -> None:
    joined = configs.merge(lsdyna, on="id_simulacao", how="inner", suffixes=("_doe", "_ls"))
    if joined.empty:
        raise ValueError("Nenhum ID em comum entre DoE e LS-DYNA.")

    for i in range(1, 8):
        m_doe = f"material_camada_{i}_doe"
        m_ls = f"material_camada_{i}_ls"
        # Como ambos os datasets têm nomes iguais, pandas adiciona suffixes.
        if m_doe not in joined.columns or m_ls not in joined.columns:
            raise RuntimeError(f"Falha na auditoria de material da camada {i}.")
        mismatch = joined[m_doe].astype(str) != joined[m_ls].astype(str)
        if mismatch.any():
            ids = joined.loc[mismatch, "id_simulacao"].tolist()[:10]
            raise ValueError(f"DoE x LS-DYNA: material divergente na camada {i}; IDs {ids}")

        t_doe = f"espessura_camada_{i}"
        t_ls = f"espessura_camada_{i}_cm"
        a = pd.to_numeric(joined[t_doe], errors="raise").to_numpy()
        b = pd.to_numeric(joined[t_ls], errors="raise").to_numpy()
        mismatch_t = ~np.isclose(a, b, atol=5e-7, rtol=0.0)
        if mismatch_t.any():
            ids = joined.loc[mismatch_t, "id_simulacao"].tolist()[:10]
            raise ValueError(f"DoE x LS-DYNA: espessura divergente na camada {i}; IDs {ids}")


def build_dose_dataset(
    configs: pd.DataFrame,
    mcnp: pd.DataFrame,
) -> pd.DataFrame:
    """
    Dataset específico do surrogate radiológico.

    Usa:
        DoE canônico + resultados MCNP válidos.

    NÃO depende da disponibilidade de resultados LS-DYNA.
    """
    df = configs.merge(
        mcnp,
        on="id_simulacao",
        how="inner",
        validate="one_to_one",
    )

    validate_design_dataframe(df, MATERIAL_COLS, THICKNESS_COLS)

    if df[TARGET_DOSE].isna().any():
        raise ValueError(
            f"Dataset de dose contém "
            f"{int(df[TARGET_DOSE].isna().sum())} valores NaN no target."
        )

    if df["id_simulacao"].duplicated().any():
        raise ValueError("Dataset de dose contém IDs duplicados.")

    return (
        df.sort_values("id_simulacao")
        .reset_index(drop=True)
    )


def build_penetration_dataset(
    configs: pd.DataFrame,
    lsdyna: pd.DataFrame,
) -> pd.DataFrame:
    """
    Dataset específico do surrogate balístico.

    Usa:
        DoE canônico + resultados LS-DYNA válidos.

    NÃO depende da disponibilidade de resultados MCNP.
    """

    # Verificação fundamental:
    # a geometria executada no LS-DYNA deve corresponder ao DoE.
    verify_lsdyna_design_matches_doe(configs, lsdyna)

    ls_targets = lsdyna[
        [
            "id_simulacao",
            TARGET_PEN,
            "perforation_bool",
            "ballistic_status",
        ]
    ].copy()

    df = configs.merge(
        ls_targets,
        on="id_simulacao",
        how="inner",
        validate="one_to_one",
    )

    validate_design_dataframe(df, MATERIAL_COLS, THICKNESS_COLS)

    if df[TARGET_PEN].isna().any():
        raise ValueError(
            f"Dataset de penetração contém "
            f"{int(df[TARGET_PEN].isna().sum())} valores NaN no target."
        )

    if df["id_simulacao"].duplicated().any():
        raise ValueError("Dataset de penetração contém IDs duplicados.")

    return (
        df.sort_values("id_simulacao")
        .reset_index(drop=True)
    )


def enforce_run_mode(
    dose_dataset: pd.DataFrame,
    penetration_dataset: pd.DataFrame,
    lsdyna_audit: dict,
) -> None:
    """
    Aplica independentemente o critério mínimo de completude
    aos dois targets.
    """

    n_dose = len(dose_dataset)
    n_pen = len(penetration_dataset)

    if RUN_MODE.upper() == "FINAL":

        if n_dose != EXPECTED_CAMPAIGN_CASES:
            raise RuntimeError(
                "RUN_MODE='FINAL': conjunto radiológico incompleto. "
                f"MCNP válidos associados ao DoE: {n_dose}; "
                f"esperados: {EXPECTED_CAMPAIGN_CASES}."
            )

        if n_pen != EXPECTED_LSDYNA_VALID_CASES:
            raise RuntimeError(
                "RUN_MODE='FINAL': conjunto balístico incompatível com a "
                "campanha provisória documentada. "
                f"LS-DYNA válidos associados ao DoE: {n_pen}; "
                f"esperados: {EXPECTED_LSDYNA_VALID_CASES}."
            )

        if lsdyna_audit["rows_total"] != EXPECTED_LSDYNA_VALID_CASES:
            raise RuntimeError(
                "RUN_MODE='FINAL': recorte balístico contém número inesperado de registros."
            )

        missing_ids = set(range(1, EXPECTED_CAMPAIGN_CASES + 1)) - set(lsdyna_audit["ids"])
        if missing_ids != LSDYNA_PENDING_IDS:
            raise RuntimeError(
                "RUN_MODE='FINAL': IDs ausentes do recorte balístico não coincidem "
                "com a lista documental de pendências."
            )

    elif RUN_MODE.upper() == "SMOKE_TEST":

        if n_dose < MIN_SMOKE_CASES:
            raise RuntimeError(
                "SMOKE_TEST: conjunto de dose possui apenas "
                f"{n_dose} casos; mínimo: {MIN_SMOKE_CASES}."
            )

        if n_pen < MIN_SMOKE_CASES:
            raise RuntimeError(
                "SMOKE_TEST: conjunto de penetração possui apenas "
                f"{n_pen} casos; mínimo: {MIN_SMOKE_CASES}."
            )

    else:
        raise ValueError(
            "RUN_MODE deve ser 'FINAL' ou 'SMOKE_TEST'."
        )


# =============================================================================
# 4. SPLIT 70/15/15 FIXO E REPRODUTÍVEL
# =============================================================================


def make_splits(df: pd.DataFrame) -> dict[str, np.ndarray]:
    ids = df["id_simulacao"].to_numpy()

    ids_traincal, ids_test = train_test_split(
        ids,
        test_size=TEST_FRACTION,
        random_state=RANDOM_STATE,
        shuffle=True,
    )

    cal_fraction_within_traincal = CALIBRATION_FRACTION / (1.0 - TEST_FRACTION)
    ids_train, ids_cal = train_test_split(
        ids_traincal,
        test_size=cal_fraction_within_traincal,
        random_state=RANDOM_STATE + 1,
        shuffle=True,
    )

    return {
        "train": np.sort(ids_train),
        "calibration": np.sort(ids_cal),
        "test": np.sort(ids_test),
    }

def save_split_ids(
    splits: dict[str, np.ndarray],
    output_path: Path,
) -> None:
    """
    Persiste os IDs de treino, calibração e teste de um target.
    """
    rows = []

    for split_name, ids in splits.items():
        rows.extend(
            {
                "id_simulacao": int(sim_id),
                "split": split_name,
            }
            for sim_id in ids
        )

    (
        pd.DataFrame(rows)
        .sort_values("id_simulacao")
        .to_csv(output_path, index=False)
    )

def subset_by_ids(df: pd.DataFrame, ids: np.ndarray) -> pd.DataFrame:
    out = df[df["id_simulacao"].isin(set(map(int, ids)))].copy()
    return out.sort_values("id_simulacao").reset_index(drop=True)


# =============================================================================
# 5. PIPELINES E MODELOS
# =============================================================================


def make_onehot():
    kwargs = dict(handle_unknown="ignore", categories=[list(MATERIALS)] * 7)
    try:
        return OneHotEncoder(sparse_output=False, **kwargs)
    except TypeError:
        return OneHotEncoder(sparse=False, **kwargs)


def numeric_columns_for(feature_set: str) -> list[str]:
    base = list(THICKNESS_COLS)
    if feature_set.upper() == "PHYSICS":
        return base + physics_feature_names()
    if feature_set.upper() == "PHYSICS_STRENGTH23":
        return base + physics_feature_names() + ballistic_strength23_feature_names()
    if feature_set.upper() == "RAW":
        return base
    raise ValueError(feature_set)


def make_preprocessor(feature_set: str) -> ColumnTransformer:
    return ColumnTransformer(
        transformers=[
            ("materials", make_onehot(), list(MATERIAL_COLS)),
            ("numeric", StandardScaler(), numeric_columns_for(feature_set)),
        ],
        remainder="drop",
        sparse_threshold=0.0,
    )


def model_factories() -> dict[str, callable]:
    n_cpu = os.cpu_count() or 1
    factories: dict[str, callable] = {
        "Ridge": lambda: Ridge(alpha=1.0),
        "RandomForest": lambda: RandomForestRegressor(
            n_estimators=400,
            max_depth=None,
            min_samples_leaf=2,
            max_features=0.8,
            random_state=RANDOM_STATE,
            n_jobs=n_cpu,
        ),
        "ExtraTrees": lambda: ExtraTreesRegressor(
            n_estimators=400,
            max_depth=None,
            min_samples_leaf=2,
            max_features=0.9,
            random_state=RANDOM_STATE,
            n_jobs=n_cpu,
        ),
        "GradientBoosting": lambda: GradientBoostingRegressor(
            n_estimators=350,
            learning_rate=0.04,
            max_depth=3,
            min_samples_leaf=4,
            subsample=0.85,
            random_state=RANDOM_STATE,
        ),
        "HistGradientBoosting": lambda: HistGradientBoostingRegressor(
            learning_rate=0.05,
            max_iter=350,
            max_leaf_nodes=31,
            min_samples_leaf=20,
            l2_regularization=0.1,
            random_state=RANDOM_STATE,
        ),
        "MLP": lambda: MLPRegressor(
            hidden_layer_sizes=(128, 64, 32),
            activation="relu",
            solver="adam",
            alpha=1e-4,
            learning_rate_init=1e-3,
            max_iter=2000,
            early_stopping=True,
            validation_fraction=0.15,
            n_iter_no_change=40,
            random_state=RANDOM_STATE,
        ),
    }

    if XGBOOST_AVAILABLE:
        factories["XGBoost"] = lambda: xgb.XGBRegressor(
            n_estimators=500,
            max_depth=5,
            learning_rate=0.04,
            subsample=0.85,
            colsample_bytree=0.85,
            reg_alpha=0.0,
            reg_lambda=1.0,
            objective="reg:squarederror",
            random_state=RANDOM_STATE,
            n_jobs=n_cpu,
            verbosity=0,
        )
    else:
        warnings.warn("xgboost não disponível: XGBoost será omitido do benchmark.")

    if RUN_MODE.upper() == "SMOKE_TEST":
        # Teste operacional curto: confirma carregamento, CV, feature engineering,
        # serialização e gráficos sem produzir artefatos científicos finais.
        return {
            "Ridge": factories["Ridge"],
            "ExtraTrees": lambda: ExtraTreesRegressor(
                n_estimators=64,
                max_depth=None,
                min_samples_leaf=2,
                max_features=0.9,
                random_state=RANDOM_STATE,
                n_jobs=n_cpu,
            ),
        }
    return factories


def make_pipeline(feature_set: str, model) -> Pipeline:
    return Pipeline(
        steps=[
            ("feature_engineering", ShieldFeatureEngineer(feature_set=feature_set)),
            ("preprocess", make_preprocessor(feature_set)),
            ("model", model),
        ]
    )


def feature_sets_for_target(target: str) -> tuple[str, ...]:
    """Pre-specified representations available to each target."""
    if target == TARGET_PEN:
        return ("RAW", "PHYSICS", "PHYSICS_STRENGTH23")
    return ("RAW", "PHYSICS")


# =============================================================================
# 6. BENCHMARK 5-FOLD NO TREINO
# =============================================================================


def benchmark_target(
    train_df: pd.DataFrame,
    target: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    X = train_df[RAW_INPUT_COLUMNS].copy()
    y = train_df[target].to_numpy(dtype=float)

    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    splits = list(kf.split(X))
    factories = model_factories()

    summary_rows = []
    fold_rows = []

    # The Strength23 representation was pre-specified from a strictly
    # training-only mechanical-descriptor study.  It is relevant only to
    # penetration and is unavailable to the radiological-dose target.
    for feature_set in feature_sets_for_target(target):
        for model_name, factory in factories.items():
            print(f"\n[{target}] {feature_set} | {model_name}")
            fold_metrics = []
            total_start = time.time()

            for fold, (tr_idx, va_idx) in enumerate(splits, start=1):
                pipe = make_pipeline(feature_set, factory())
                X_tr = X.iloc[tr_idx]
                X_va = X.iloc[va_idx]
                y_tr = y[tr_idx]
                y_va = y[va_idx]

                t0 = time.time()
                pipe.fit(X_tr, y_tr)
                pred = np.maximum(0.0, pipe.predict(X_va))
                fit_predict_s = time.time() - t0

                m = regression_metrics(y_va, pred)
                m.update(
                    {
                        "target": target,
                        "feature_set": feature_set,
                        "model": model_name,
                        "fold": fold,
                        "n_train_fold": len(tr_idx),
                        "n_val_fold": len(va_idx),
                        "fit_predict_s": fit_predict_s,
                    }
                )
                fold_rows.append(m)
                fold_metrics.append(m)
                print(
                    f"  fold {fold}: R²={m['r2']:.4f} | RMSE={m['rmse']:.6g} | "
                    f"MAE={m['mae']:.6g}"
                )

            elapsed = time.time() - total_start
            summary_rows.append(
                {
                    "target": target,
                    "feature_set": feature_set,
                    "model": model_name,
                    "cv_r2_mean": float(np.mean([m["r2"] for m in fold_metrics])),
                    "cv_r2_std": float(np.std([m["r2"] for m in fold_metrics], ddof=1)),
                    "cv_rmse_mean": float(np.mean([m["rmse"] for m in fold_metrics])),
                    "cv_rmse_std": float(np.std([m["rmse"] for m in fold_metrics], ddof=1)),
                    "cv_mae_mean": float(np.mean([m["mae"] for m in fold_metrics])),
                    "cv_mae_std": float(np.std([m["mae"] for m in fold_metrics], ddof=1)),
                    "cv_smape_mean_pct": float(np.mean([m["smape_pct"] for m in fold_metrics])),
                    "elapsed_s": elapsed,
                }
            )

    summary = pd.DataFrame(summary_rows).sort_values(
        ["cv_rmse_mean", "cv_mae_mean", "cv_r2_mean"],
        ascending=[True, True, False],
    ).reset_index(drop=True)
    folds = pd.DataFrame(fold_rows)
    return summary, folds


def paired_feature_set_comparison(folds_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (target, model), g in folds_df.groupby(["target", "model"]):
        raw = g[g["feature_set"] == "RAW"].sort_values("fold")
        if len(raw) != N_FOLDS:
            continue
        raw_rmse = raw["rmse"].to_numpy()
        for comparator in ("PHYSICS", "PHYSICS_STRENGTH23"):
            other = g[g["feature_set"] == comparator].sort_values("fold")
            if len(other) != N_FOLDS:
                continue
            other_rmse = other["rmse"].to_numpy()
            diff = raw_rmse - other_rmse  # positivo = representação comparada melhor
            try:
                p = 1.0 if np.allclose(diff, 0) else float(
                    wilcoxon(raw_rmse, other_rmse, alternative="two-sided").pvalue
                )
            except Exception:
                p = np.nan
            rows.append(
                {
                    "target": target,
                    "model": model,
                    "baseline_feature_set": "RAW",
                    "comparison_feature_set": comparator,
                    "baseline_rmse_mean": float(raw_rmse.mean()),
                    "comparison_rmse_mean": float(other_rmse.mean()),
                    "rmse_reduction_pct": float(100.0 * (raw_rmse.mean() - other_rmse.mean()) / raw_rmse.mean())
                    if raw_rmse.mean() != 0 else np.nan,
                    "paired_wilcoxon_p": p,
                    "comparison_better_folds": int(np.sum(other_rmse < raw_rmse)),
                    "n_folds": N_FOLDS,
                }
            )
    return pd.DataFrame(rows)


# =============================================================================
# 7. TUNING DO VENCEDOR NO TREINO
# =============================================================================


def tuning_space(model_name: str) -> dict:
    spaces = {
        "Ridge": {
            "model__alpha": [1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1000.0],
        },
        "RandomForest": {
            "model__n_estimators": [300, 500, 700],
            "model__max_depth": [None, 8, 12, 16, 24],
            "model__min_samples_leaf": [1, 2, 4, 6],
            "model__max_features": [0.5, 0.7, 0.9, 1.0],
        },
        "ExtraTrees": {
            "model__n_estimators": [300, 500, 700],
            "model__max_depth": [None, 8, 12, 16, 24],
            "model__min_samples_leaf": [1, 2, 4, 6],
            "model__max_features": [0.5, 0.7, 0.9, 1.0],
        },
        "GradientBoosting": {
            "model__n_estimators": [200, 350, 500, 700],
            "model__learning_rate": [0.015, 0.03, 0.05, 0.08],
            "model__max_depth": [2, 3, 4],
            "model__min_samples_leaf": [2, 4, 8, 12],
            "model__subsample": [0.7, 0.85, 1.0],
        },
        "HistGradientBoosting": {
            "model__learning_rate": [0.02, 0.04, 0.06, 0.1],
            "model__max_iter": [200, 350, 500, 700],
            "model__max_leaf_nodes": [15, 31, 63],
            "model__min_samples_leaf": [10, 20, 30, 40],
            "model__l2_regularization": [0.0, 0.01, 0.1, 1.0],
        },
        "MLP": {
            "model__hidden_layer_sizes": [(64, 32), (128, 64), (128, 64, 32), (256, 128, 64)],
            "model__alpha": [1e-5, 1e-4, 1e-3, 1e-2],
            "model__learning_rate_init": [3e-4, 1e-3, 3e-3],
        },
        "XGBoost": {
            "model__n_estimators": [300, 500, 700, 900],
            "model__max_depth": [3, 4, 5, 6, 8],
            "model__learning_rate": [0.02, 0.04, 0.06, 0.1],
            "model__subsample": [0.7, 0.85, 1.0],
            "model__colsample_bytree": [0.7, 0.85, 1.0],
            "model__reg_lambda": [0.1, 1.0, 5.0, 10.0],
            "model__min_child_weight": [1, 3, 5, 10],
        },
    }
    return spaces.get(model_name, {})


def effective_tuning_iterations(space: dict) -> int:
    """Avoid sampling more times than a finite parameter grid contains."""
    if not space:
        return 0
    total = int(np.prod([len(values) for values in space.values()], dtype=np.int64))
    return min(TUNING_ITERATIONS, total)


def tune_all_candidates(
    train_df: pd.DataFrame,
    target: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[tuple[str, str], Pipeline]]:
    """Optimize and compare every pre-specified candidate on training data only.

    RandomizedSearchCV performs 5-fold tuning for each representation-model
    pair.  The returned estimators are refit on all training cases by sklearn.
    A second fixed-parameter 5-fold pass records R², RMSE and MAE for reporting
    only.  Neither operation reads calibration or test records.
    """
    X = train_df[RAW_INPUT_COLUMNS].copy()
    y = train_df[target].to_numpy(dtype=float)
    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    factories = model_factories()
    summary_rows: list[dict] = []
    fold_rows: list[dict] = []
    fitted: dict[tuple[str, str], Pipeline] = {}

    for feature_set in feature_sets_for_target(target):
        for model_name, factory in factories.items():
            space = tuning_space(model_name)
            n_iter = effective_tuning_iterations(space)
            print(
                f"\n[{target}] OTIMIZAÇÃO FORMAL | {feature_set} | {model_name} "
                f"| candidatos={n_iter} | folds={N_FOLDS}",
                flush=True,
            )
            base = make_pipeline(feature_set, factory())
            started = time.time()
            if n_iter:
                search = RandomizedSearchCV(
                    estimator=base,
                    param_distributions=space,
                    n_iter=n_iter,
                    scoring="neg_root_mean_squared_error",
                    cv=kf,
                    refit=True,
                    random_state=RANDOM_STATE,
                    n_jobs=1,
                    return_train_score=False,
                    verbose=1,
                )
                search.fit(X, y)
                best_estimator = search.best_estimator_
                best_params = search.best_params_
                tuning_cv_rmse = float(-search.best_score_)
                tuning_cv_std = float(search.cv_results_["std_test_score"][search.best_index_])
            else:
                best_estimator = base.fit(X, y)
                best_params = {}
                tuning_cv_rmse = np.nan
                tuning_cv_std = np.nan

            # Repeated with the frozen winning hyperparameters for complete
            # fold-level metrics.  This is still an internal training-CV report.
            candidate_folds = []
            for fold, (tr_idx, va_idx) in enumerate(kf.split(X), start=1):
                pipe = clone(best_estimator)
                t0 = time.time()
                pipe.fit(X.iloc[tr_idx], y[tr_idx])
                prediction = np.maximum(0.0, pipe.predict(X.iloc[va_idx]))
                m = regression_metrics(y[va_idx], prediction)
                m.update({
                    "target": target,
                    "feature_set": feature_set,
                    "model": model_name,
                    "fold": fold,
                    "n_train_fold": len(tr_idx),
                    "n_val_fold": len(va_idx),
                    "fit_predict_s": time.time() - t0,
                    "tuning_best_params": json.dumps(best_params, sort_keys=True, default=str),
                })
                candidate_folds.append(m)
                fold_rows.append(m)

            summary_rows.append({
                "target": target,
                "feature_set": feature_set,
                "model": model_name,
                "cv_r2_mean": float(np.mean([m["r2"] for m in candidate_folds])),
                "cv_r2_std": float(np.std([m["r2"] for m in candidate_folds], ddof=1)),
                "cv_rmse_mean": float(np.mean([m["rmse"] for m in candidate_folds])),
                "cv_rmse_std": float(np.std([m["rmse"] for m in candidate_folds], ddof=1)),
                "cv_mae_mean": float(np.mean([m["mae"] for m in candidate_folds])),
                "cv_mae_std": float(np.std([m["mae"] for m in candidate_folds], ddof=1)),
                "cv_smape_mean_pct": float(np.mean([m["smape_pct"] for m in candidate_folds])),
                "tuning_method": "RandomizedSearchCV",
                "tuning_scoring": "neg_root_mean_squared_error",
                "tuning_n_iter": n_iter,
                "tuning_best_cv_rmse": tuning_cv_rmse,
                "tuning_best_cv_rmse_std": tuning_cv_std,
                "best_params": json.dumps(best_params, sort_keys=True, default=str),
                "elapsed_s": time.time() - started,
            })
            fitted[(feature_set, model_name)] = best_estimator

    summary = pd.DataFrame(summary_rows).sort_values(
        ["tuning_best_cv_rmse", "cv_rmse_mean", "cv_mae_mean"],
        ascending=[True, True, True],
    ).reset_index(drop=True)
    return summary, pd.DataFrame(fold_rows), fitted


def tune_winner(
    train_df: pd.DataFrame,
    target: str,
    feature_set: str,
    model_name: str,
):
    factories = model_factories()
    base = make_pipeline(feature_set, factories[model_name]())

    if not TUNE_WINNER:
        base.fit(train_df[RAW_INPUT_COLUMNS], train_df[target].to_numpy(dtype=float))
        return base, {"tuning": False, "best_params": {}}

    space = tuning_space(model_name)
    if not space:
        base.fit(train_df[RAW_INPUT_COLUMNS], train_df[target].to_numpy(dtype=float))
        return base, {"tuning": False, "best_params": {}}

    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    search = RandomizedSearchCV(
        estimator=base,
        param_distributions=space,
        n_iter=TUNING_ITERATIONS,
        scoring="neg_root_mean_squared_error",
        cv=kf,
        refit=True,
        random_state=RANDOM_STATE,
        n_jobs=1,  # evita nested parallelism; modelos de árvores/XGB usam CPU internamente
        return_train_score=False,
        verbose=1,
    )
    search.fit(train_df[RAW_INPUT_COLUMNS], train_df[target].to_numpy(dtype=float))

    info = {
        "tuning": True,
        "best_params": search.best_params_,
        "best_cv_rmse": float(-search.best_score_),
        "n_iter": TUNING_ITERATIONS,
    }
    return search.best_estimator_, info


# =============================================================================
# 8. CALIBRAÇÃO CONFORMAL E TESTE FINAL
# =============================================================================


def calibrate_and_test(
    fitted_pipeline: Pipeline,
    calibration_df: pd.DataFrame,
    test_df: pd.DataFrame,
    target: str,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    X_cal = calibration_df[RAW_INPUT_COLUMNS]
    y_cal = calibration_df[target].to_numpy(dtype=float)
    pred_cal = np.maximum(0.0, fitted_pipeline.predict(X_cal))
    abs_resid = np.abs(y_cal - pred_cal)
    q = conformal_quantile(abs_resid, CONFORMAL_ALPHA)

    X_test = test_df[RAW_INPUT_COLUMNS]
    y_test = test_df[target].to_numpy(dtype=float)
    inference_start = time.perf_counter()
    pred_test = np.maximum(0.0, fitted_pipeline.predict(X_test))
    inference_seconds_total = time.perf_counter() - inference_start

    lower = np.maximum(0.0, pred_test - q)
    upper = pred_test + q
    covered = (y_test >= lower) & (y_test <= upper)

    metrics = regression_metrics(y_test, pred_test)
    metrics.update(
        {
            "conformal_alpha": CONFORMAL_ALPHA,
            "conformal_nominal_coverage": 1.0 - CONFORMAL_ALPHA,
            "conformal_q_abs": float(q),
            "test_empirical_coverage": float(np.mean(covered)),
            "test_mean_interval_width": float(np.mean(upper - lower)),
            "n_calibration": int(len(calibration_df)),
            "n_test": int(len(test_df)),
            "test_inference_seconds_total": float(inference_seconds_total),
            "test_inference_seconds_per_case": float(inference_seconds_total / len(test_df)),
        }
    )

    pred_cal_df = pd.DataFrame(
        {
            "id_simulacao": calibration_df["id_simulacao"].to_numpy(),
            "target": target,
            "y_true": y_cal,
            "y_pred": pred_cal,
            "abs_residual": abs_resid,
        }
    )

    pred_test_df = pd.DataFrame(
        {
            "id_simulacao": test_df["id_simulacao"].to_numpy(),
            "target": target,
            "y_true": y_test,
            "y_pred": pred_test,
            "pi_lower_95": lower,
            "pi_upper_95": upper,
            "covered_95": covered.astype(int),
            "residual": y_test - pred_test,
        }
    )

    return metrics, pred_cal_df, pred_test_df


# =============================================================================
# 9. FIGURAS DO MODELO FINAL
# =============================================================================


def save_parity_plot(pred_df: pd.DataFrame, target: str) -> None:
    y = pred_df["y_true"].to_numpy(dtype=float)
    p = pred_df["y_pred"].to_numpy(dtype=float)
    lo = min(float(y.min()), float(p.min()))
    hi = max(float(y.max()), float(p.max()))

    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    ax.scatter(y, p, alpha=0.65, s=24)
    ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.2)
    ax.set_xlabel("Valor simulado")
    ax.set_ylabel("Valor predito")
    ax.set_title(f"Parity plot — {target}")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / f"parity_{target}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_residual_plot(pred_df: pd.DataFrame, target: str) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.scatter(pred_df["y_pred"], pred_df["residual"], alpha=0.65, s=24)
    ax.axhline(0.0, linestyle="--", linewidth=1.2)
    ax.set_xlabel("Valor predito")
    ax.set_ylabel("Resíduo (simulado - predito)")
    ax.set_title(f"Resíduos — {target}")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / f"residuals_{target}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# 10. TREINO FINAL POR TARGET
# =============================================================================


def train_final_target(
    target: str,
    train_df: pd.DataFrame,
    cal_df: pd.DataFrame,
    test_df: pd.DataFrame,
    benchmark_summary: pd.DataFrame,
    tuned_candidates: dict[tuple[str, str], Pipeline] | None = None,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    ranking = benchmark_summary[benchmark_summary["target"] == target].copy()
    ranking_columns = ["cv_rmse_mean", "cv_mae_mean", "cv_r2_mean"]
    if "tuning_best_cv_rmse" in ranking.columns:
        ranking_columns = ["tuning_best_cv_rmse", "cv_rmse_mean", "cv_mae_mean", "cv_r2_mean"]
    ranking_ascending = [True, True, False]
    if "tuning_best_cv_rmse" in ranking.columns:
        ranking_ascending = [True, True, True, False]
    ranking = ranking.sort_values(ranking_columns, ascending=ranking_ascending)
    if ranking.empty:
        raise RuntimeError(f"Sem benchmark válido para {target}.")

    winner = ranking.iloc[0]
    feature_set = str(winner["feature_set"])
    model_name = str(winner["model"])

    print("\n" + "=" * 78)
    print(f"FINALIZANDO TARGET: {target}")
    print(f"Vencedor CV: {feature_set} + {model_name}")
    print(f"CV RMSE: {winner['cv_rmse_mean']:.6g} ± {winner['cv_rmse_std']:.6g}")
    print("=" * 78)

    candidate_key = (feature_set, model_name)
    if tuned_candidates is not None and candidate_key in tuned_candidates:
        fitted = tuned_candidates[candidate_key]
        row = winner.to_dict()
        tuning_info = {
            "tuning": True,
            "tuning_method": "RandomizedSearchCV",
            "tuning_scope": "all pre-specified candidates, training partition only",
            "best_params": json.loads(str(row.get("best_params", "{}"))),
            "best_cv_rmse": float(row["tuning_best_cv_rmse"]),
            "n_iter": int(row["tuning_n_iter"]),
        }
    else:
        fitted, tuning_info = tune_winner(train_df, target, feature_set, model_name)
    metrics, pred_cal, pred_test = calibrate_and_test(fitted, cal_df, test_df, target)

    bundle = {
        "script_version": SCRIPT_VERSION,
        "created_utc": utc_now(),
        "target": target,
        "feature_set": feature_set,
        "model_name": model_name,
        "pipeline": fitted,
        "conformal_alpha": CONFORMAL_ALPHA,
        "conformal_q_abs": metrics["conformal_q_abs"],
        "input_columns": RAW_INPUT_COLUMNS,
        "material_order_constraint": "each of the seven materials exactly once",
        "layer_order": "layer 1 = operator/inner face; layer 7 = source-impact/outer face",
        "training_ids": train_df["id_simulacao"].astype(int).tolist(),
        "calibration_ids": cal_df["id_simulacao"].astype(int).tolist(),
        "test_ids": test_df["id_simulacao"].astype(int).tolist(),
        "tuning": tuning_info,
        "test_metrics": metrics,
    }

    file_name = "surrogate_dose.joblib" if target == TARGET_DOSE else "surrogate_penetration.joblib"
    joblib.dump(bundle, MODEL_DIR / file_name, compress=3)

    metadata = {k: v for k, v in bundle.items() if k != "pipeline"}
    dump_json(metadata, MODEL_DIR / file_name.replace(".joblib", "_metadata.json"))

    pred_cal.to_csv(OUTPUT_DIR / f"calibration_predictions_{target}.csv", index=False)
    pred_test.to_csv(OUTPUT_DIR / f"test_predictions_{target}.csv", index=False)
    save_parity_plot(pred_test, target)
    save_residual_plot(pred_test, target)

    print(
        f"TESTE FINAL | R²={metrics['r2']:.4f} | RMSE={metrics['rmse']:.6g} | "
        f"MAE={metrics['mae']:.6g} | cobertura95={metrics['test_empirical_coverage']:.3f}"
    )

    return metadata, pred_cal, pred_test


# =============================================================================
# 11. MAIN
# =============================================================================


def main() -> None:
    np.random.seed(RANDOM_STATE)

    print("=" * 88)
    print("FASE 3 — SURROGATE MODELING: DOSE MCNP + PENETRAÇÃO LS-DYNA")
    print("=" * 88)
    print(f"Script: {SCRIPT_VERSION}")
    print(f"RUN_MODE: {RUN_MODE}")
    print("Targets: dose efetiva em 5 min + max_penetration_mm")
    print("Massa/custo: determinísticos; excluídos do ML")

    print(f"\nMCNP ML:    {DOSE_DATASET_PATH}")
    print(f"LS-DYNA ML: {PENETRATION_DATASET_PATH}")

    # ==========================================================================
    # DATASETS ESPECÍFICOS POR TARGET
    # ==========================================================================

    dose_dataset, dose_audit = load_consolidated_target_dataset(
        path=DOSE_DATASET_PATH,
        target=TARGET_DOSE,
        expected_rows=EXPECTED_CAMPAIGN_CASES,
    )

    penetration_dataset, lsdyna_audit = load_consolidated_target_dataset(
        path=PENETRATION_DATASET_PATH,
        target=TARGET_PEN,
        expected_rows=EXPECTED_LSDYNA_VALID_CASES,
    )

    enforce_run_mode(
        dose_dataset=dose_dataset,
        penetration_dataset=penetration_dataset,
        lsdyna_audit=lsdyna_audit,
    )

    print("\nDATASETS CIENTÍFICOS:")
    print(f"   • Dose MCNP:           {len(dose_dataset)} casos")
    print(f"   • Penetração LS-DYNA: {len(penetration_dataset)} casos")

    print(
        f"   • LS-DYNA válidos: "
        f"{lsdyna_audit['rows_valid']}/"
        f"{lsdyna_audit['rows_total']}"
    )
    print(
        f"   • Status balístico: "
        f"{lsdyna_audit['ballistic_status_counts']}"
    )
    print(
        "   • Perfuração: não incluída no recorte ML consolidado "
        "(não é necessária para o treino do alvo de penetração)"
    )

    if RUN_MODE.upper() == "SMOKE_TEST":
        dose_dataset = dose_dataset.sample(
            n=min(SMOKE_SAMPLE_SIZE, len(dose_dataset)),
            random_state=RANDOM_STATE,
        ).sort_values("id_simulacao").reset_index(drop=True)
        penetration_dataset = penetration_dataset.sample(
            n=min(SMOKE_SAMPLE_SIZE, len(penetration_dataset)),
            random_state=RANDOM_STATE,
        ).sort_values("id_simulacao").reset_index(drop=True)
        print(
            f"   • SMOKE_TEST: amostra determinística de {len(dose_dataset)} "
            "casos por target; artefatos separados dos resultados finais."
        )

    # ==========================================================================
    # PERSISTÊNCIA DOS DATASETS ESPECÍFICOS
    # ==========================================================================

    dose_dataset.to_csv(
        OUTPUT_DIR / "dataset_dose_fase3.csv",
        index=False,
    )

    penetration_dataset.to_csv(
        OUTPUT_DIR / "dataset_penetracao_fase3.csv",
        index=False,
    )

    # ==========================================================================
    # SPLITS ESPECÍFICOS POR TARGET
    # ==========================================================================

    dose_splits = make_splits(dose_dataset)

    penetration_splits = make_splits(
        penetration_dataset
    )

    save_split_ids(
        dose_splits,
        OUTPUT_DIR / "split_ids_dose.csv",
    )

    save_split_ids(
        penetration_splits,
        OUTPUT_DIR / "split_ids_penetration.csv",
    )

    # ==========================================================================
    # DATAFRAMES DE TREINO / CALIBRAÇÃO / TESTE
    # ==========================================================================

    dose_train_df = subset_by_ids(
        dose_dataset,
        dose_splits["train"],
    )
    dose_cal_df = subset_by_ids(
        dose_dataset,
        dose_splits["calibration"],
    )
    dose_test_df = subset_by_ids(
        dose_dataset,
        dose_splits["test"],
    )

    pen_train_df = subset_by_ids(
        penetration_dataset,
        penetration_splits["train"],
    )
    pen_cal_df = subset_by_ids(
        penetration_dataset,
        penetration_splits["calibration"],
    )
    pen_test_df = subset_by_ids(
        penetration_dataset,
        penetration_splits["test"],
    )

    print("\nSPLIT — DOSE:")
    print(
        f"   treino={len(dose_train_df)} "
        f"({len(dose_train_df)/len(dose_dataset):.1%}) | "
        f"calibração={len(dose_cal_df)} "
        f"({len(dose_cal_df)/len(dose_dataset):.1%}) | "
        f"teste={len(dose_test_df)} "
        f"({len(dose_test_df)/len(dose_dataset):.1%})"
    )

    print("\nSPLIT — PENETRAÇÃO:")
    print(
        f"   treino={len(pen_train_df)} "
        f"({len(pen_train_df)/len(penetration_dataset):.1%}) | "
        f"calibração={len(pen_cal_df)} "
        f"({len(pen_cal_df)/len(penetration_dataset):.1%}) | "
        f"teste={len(pen_test_df)} "
        f"({len(pen_test_df)/len(penetration_dataset):.1%})"
    )

    # ==========================================================================
    # CONTRATO TARGET -> DATASET/SPLITS
    # ==========================================================================

    target_data = {
        TARGET_DOSE: {
            "dataset": dose_dataset,
            "train": dose_train_df,
            "calibration": dose_cal_df,
            "test": dose_test_df,
            "splits": dose_splits,
            "split_file": "split_ids_dose.csv",
        },

        TARGET_PEN: {
            "dataset": penetration_dataset,
            "train": pen_train_df,
            "calibration": pen_cal_df,
            "test": pen_test_df,
            "splits": penetration_splits,
            "split_file": "split_ids_penetration.csv",
        },
    }

    # ==========================================================================
    # MANIFESTO
    # ==========================================================================

    manifest = {
        "script_version": SCRIPT_VERSION,
        "created_utc": utc_now(),
        "run_mode": RUN_MODE,
        "random_state": RANDOM_STATE,
        "n_folds": N_FOLDS,
        "candidate_optimization": {
            "enabled": TUNE_ALL_CANDIDATES,
            "method": "RandomizedSearchCV",
            "scoring": "neg_root_mean_squared_error",
            "n_iter_requested": TUNING_ITERATIONS,
            "scope": "all pre-specified model-representation candidates; training partition only",
        },

        "datasets": {
            TARGET_DOSE: {
                "n_total": len(dose_dataset),
                "train": len(dose_train_df),
                "calibration": len(dose_cal_df),
                "test": len(dose_test_df),
                "split_file": "split_ids_dose.csv",
            },

            TARGET_PEN: {
                "n_total": len(penetration_dataset),
                "train": len(pen_train_df),
                "calibration": len(pen_cal_df),
                "test": len(pen_test_df),
                "split_file": "split_ids_penetration.csv",
            },
        },

        "targets": TARGETS,
        "raw_input_columns": RAW_INPUT_COLUMNS,
        "physics_features": physics_feature_names(),

        "files": {
            "mcnp_ml_dataset": dose_audit,
            "lsdyna_ml_dataset": lsdyna_audit,
        },

        "lsdyna_audit": lsdyna_audit,

        "mcnp_dose_origin": "dose_efetiva_mSv_seg * 300 s; consolidado Fase 2",

        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": getattr(
                sklearn,
                "__version__",
                "unknown",
            ),
            "xgboost": (
                getattr(xgb, "__version__", None)
                if XGBOOST_AVAILABLE
                else None
            ),
        },
    }

    dump_json(
        manifest,
        OUTPUT_DIR / "manifesto_fase3.json",
    )

    # Otimização e comparação de todos os candidatos, exclusivamente no treino.
    all_summary = []
    all_folds = []
    all_tuned_candidates: dict[str, dict[tuple[str, str], Pipeline]] = {}

    for target in TARGETS:

        train_df_target = target_data[target]["train"]

        print("\n" + "=" * 88)
        print(f"SELEÇÃO DE CANDIDATOS: {target}")
        print(
            f"Casos disponíveis no treinamento: "
            f"{len(train_df_target)}"
        )
        print("=" * 88)

        if TUNE_ALL_CANDIDATES:
            summary, folds, tuned = tune_all_candidates(train_df_target, target)
            all_tuned_candidates[target] = tuned
        else:
            summary, folds = benchmark_target(train_df_target, target)
            all_tuned_candidates[target] = {}

        all_summary.append(summary)
        all_folds.append(folds)

    benchmark_summary = pd.concat(all_summary, ignore_index=True)
    fold_results = pd.concat(all_folds, ignore_index=True)

    benchmark_summary.to_csv(OUTPUT_DIR / "benchmark_cv_summary.csv", index=False)
    fold_results.to_csv(OUTPUT_DIR / "benchmark_cv_folds.csv", index=False)

    comparison = paired_feature_set_comparison(fold_results)
    comparison.to_csv(OUTPUT_DIR / "raw_vs_physics_paired_comparison.csv", index=False)

    # Treino, calibração conformal e teste final dos vencedores.
    final_metadata = {}

    for target in TARGETS:

        target_train_df = target_data[target]["train"]
        target_cal_df = target_data[target]["calibration"]
        target_test_df = target_data[target]["test"]

        print("\n" + "=" * 88)
        print(f"MODELO FINAL: {target}")
        print("=" * 88)

        metadata, _, _ = train_final_target(
            target,
            target_train_df,
            target_cal_df,
            target_test_df,
            benchmark_summary,
            all_tuned_candidates[target],
        )

        final_metadata[target] = metadata

    dump_json(final_metadata, OUTPUT_DIR / "final_model_summary.json")
    table_dir = generate_article_tables(OUTPUT_DIR, tuning_space, final_metadata)

    print("\n" + "=" * 88)
    print("FASE 3 CONCLUÍDA")
    print("=" * 88)
    print(f"Resultados: {OUTPUT_DIR}")
    print(f"Modelos:    {MODEL_DIR}")
    print("Arquivos principais:")
    print("  - benchmark_cv_summary.csv")
    print("  - benchmark_cv_folds.csv")
    print("  - raw_vs_physics_paired_comparison.csv")
    print("  - split_ids_dose.csv")
    print("  - split_ids_penetration.csv")
    print("  - dataset_dose_fase3.csv")
    print("  - dataset_penetracao_fase3.csv")
    print("  - surrogate_dose.joblib")
    print("  - surrogate_penetration.joblib")
    print(f"  - tabelas_artigo ({table_dir})")
    print("\nIMPORTANTE: Fase 4 deve receber SOMENTE as 14 entradas brutas; os bundles")
    print("salvos já incorporam feature engineering e preprocessing dentro do pipeline.")


def configure_cli() -> None:
    """Configura um teste operacional sem tocar nos artefatos científicos finais."""
    global RUN_MODE, N_FOLDS, TUNE_ALL_CANDIDATES, TUNE_WINNER, OUTPUT_DIR, MODEL_DIR, FIGURE_DIR
    parser = argparse.ArgumentParser(description="Fase 3 — metamodelagem da tese")
    parser.add_argument(
        "--smoke-test",
        "--debug",
        action="store_true",
        dest="smoke_test",
        help="teste curto em amostra de 160 casos; não gera bundles finais",
    )
    args = parser.parse_args()
    if args.smoke_test:
        RUN_MODE = "SMOKE_TEST"
        N_FOLDS = 2
        TUNE_ALL_CANDIDATES = False
        TUNE_WINNER = False
        OUTPUT_DIR = BASE_DIR / "analise_desempenho_modelos (Fase 3) - MCNP63 - smoke"
        MODEL_DIR = BASE_DIR / "modelos_treinados (Fase 3) - MCNP63 - smoke"
        FIGURE_DIR = OUTPUT_DIR / "figuras"
        for directory in (OUTPUT_DIR, MODEL_DIR, FIGURE_DIR):
            directory.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    configure_cli()
    main()
