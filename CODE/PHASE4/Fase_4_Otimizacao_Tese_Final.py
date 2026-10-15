# -*- coding: utf-8 -*-
"""
FASE 4 — OTIMIZAÇÃO MULTIOBJETIVO DEFINITIVA DA TESE
=====================================================

Problema
--------
Encontrar blindagens semicilíndricas de sete camadas que minimizem simultaneamente:

    f1(x) = dose efetiva prevista em 5 min [mSv]
    f2(x) = massa total determinística [kg]
    f3(x) = custo total determinístico [unidade monetária]

sujeitas a:

    1) cada um dos sete materiais aparece exatamente uma vez;
    2) 0.2 <= t_i <= 2.0 cm;
    3) sum(t_i) = 6.0 cm;
    4) restrição balística conservadora:
           max_penetration_pred_mm + q_conformal_95 <= 60 mm

A dose e a penetração são obtidas pelos bundles salvos pela Fase 3:
    modelos_treinados (Fase 3) - MCNP63_DEDF_AP/surrogate_dose.joblib
    modelos_treinados (Fase 3) - MCNP63_DEDF_AP/surrogate_penetration.joblib

Massa e custo NÃO são surrogates.

Metodologia de otimização
-------------------------
- NSGA-II de variáveis mistas, sem mutação adaptativa/restart.
- Permutação de materiais: Order Crossover (OX) + swap mutation.
- Espessuras: Simulated Binary Crossover (SBX) + polynomial mutation.
- Após operadores contínuos, projeção no simplex com caixa:
      sum(t_i)=6, 0.2<=t_i<=2.0.
- Resolução de projeto: 0.001 cm, coerente com o DoE executado.
- Constraint-domination pelas regras de viabilidade de Deb:
      viável > inviável;
      entre inviáveis, menor violação é preferida;
      entre viáveis, dominância de Pareto.
- Seleção ambiental elitista:
      pais + filhos -> non-dominated sorting -> crowding -> N sobreviventes.
- Execução final multiseed para avaliar estabilidade do otimizador.
- Hypervolume 3D com ponto de referência FIXO derivado do DoE.
- Seleção de compromisso TOPSIS-IQR (igual peso nos três objetivos).
- Exportação de candidatos representativos para validação MCNP/LS-DYNA.

IMPORTANTE
----------
Este script NÃO roda MCNP nem LS-DYNA. Ele apenas carrega os surrogates já
validados na Fase 3, executa a otimização e exporta candidatos para posterior
verificação de alta fidelidade.

Camadas:
    camada 1 = face interna / operador
    camada 7 = face externa / fonte Co-60 e impacto do FSP
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import itertools
import json
import math
import os
import pickle
import platform
import sys
import time
from typing import Any

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn

# Os bundles serializados da Fase 3 referenciam este módulo. A Fase 4 deve
# importar a mesma versão, inclusive PHYSICS_STRENGTH23, antes do unpickle.
PHASE3_FEATURE_DIR = Path(__file__).resolve().parent.parent / "Fase 3"
if str(PHASE3_FEATURE_DIR) not in sys.path:
    sys.path.insert(0, str(PHASE3_FEATURE_DIR))
from feature_engineering_final import (
    MATERIALS,
    MATERIAL_COLS,
    THICKNESS_COLS,
    validate_design_dataframe,
)


# =============================================================================
# 1. CONFIGURAÇÃO CIENTÍFICA CONGELADA
# =============================================================================

SCRIPT_VERSION = "FASE4-THESIS-CONSTRAINED-MIXED-NSGA2-v1.5-MCNP63-DEDF-AP-PARETO-METRICS"

RUN_MODE = "FINAL"  # "FINAL" ou "SMOKE_TEST"

# Execução final: 10 réplicas independentes, orçamento fixo por réplica.
FINAL_N_RUNS = 10
FINAL_POPULATION_SIZE = 300
FINAL_NUM_GENERATIONS = 200

# Teste rápido apenas para validar a pipeline depois que os .joblib existirem.
SMOKE_N_RUNS = 1
SMOKE_POPULATION_SIZE = 60
SMOKE_NUM_GENERATIONS = 15

RUN_SEEDS = [42, 137, 233, 349, 457, 563, 677, 809, 929, 1049]

# Operadores genéticos fixos — sem adaptação/restart.
CROSSOVER_PROB = 0.90
PERM_SWAP_MUTATION_PROB = 0.20
THICKNESS_GENE_MUTATION_PROB = 1.0 / 7.0
SBX_ETA = 15.0
POLYNOMIAL_MUTATION_ETA = 20.0

# Restrições geométricas.
N_LAYERS = 7
MIN_THICKNESS_CM = 0.2
MAX_THICKNESS_CM = 2.0
TOTAL_THICKNESS_CM = 6.0
THICKNESS_RESOLUTION_CM = 0.001

INNER_RADIUS_CM = 35.0
SHIELD_HEIGHT_CM = 70.0

# Restrição balística.
BALLISTIC_LIMIT_MM = 60.0
CONSTRAINT_TOL_MM = 1e-9

# Métricas / checkpoint.
METRIC_INTERVAL = 5
CHECKPOINT_INTERVAL = 50
RESUME = True
FORCE_NEW = False

# Pós-processamento.
REFERENCE_MARGIN = 0.10
TOPSIS_WEIGHTS = np.array([1.0, 1.0, 1.0], dtype=float)
N_HIGH_FIDELITY_CANDIDATES = 5
DOMAIN_DISTANCE_REPORT = True

BASE_DIR = Path(__file__).resolve().parent
PHASE3_DIR = BASE_DIR.parent / "Fase 3"
MODEL_DIR = PHASE3_DIR / "modelos_treinados (Fase 3) - MCNP63_DEDF_AP"
DOSE_MODEL_PATH = MODEL_DIR / "surrogate_dose.joblib"
PEN_MODEL_PATH = MODEL_DIR / "surrogate_penetration.joblib"

DOE_CANDIDATES = [
    BASE_DIR.parent / "Fase 2" / "DATASET_CONJUNTO_MCNP63_DEDF_AP_LSDYNA_FINAL_20261002" / "DATASET_MESTRE_CONJUNTO_MCNP63_DEDF_AP_LSDYNA_2000_FINAL.csv",
    BASE_DIR / "simulacao_configs (Fase 1)" / "configs_simulacao_reais_CONSOLIDADO.csv",
    BASE_DIR / "simulacao_configs (Fase 1)" / "configs_simulacao_reais(2).csv",
    BASE_DIR / "simulacao_configs (Fase 1)" / "configs_simulacao_reais.csv",
    BASE_DIR / "configs_simulacao_reais(2).csv",
    BASE_DIR / "configs_simulacao_reais.csv",
]

RESULTS_ROOT = BASE_DIR / "resultados_otimizacao (Fase 4) - MCNP63_DEDF_AP"

EXPECTED_DOSE_TARGET = "dose_efetiva_5min_mSv"
EXPECTED_PEN_TARGET = "max_penetration_mm"
EXPECTED_INPUT_COLUMNS = list(MATERIAL_COLS) + list(THICKNESS_COLS)
TRAINING_SKLEARN_VERSION = "1.4.2"

# As 5040 permutações são pequenas o bastante para serem pré-computadas.
ALL_PERMUTATIONS = np.asarray(
    list(itertools.permutations(range(N_LAYERS))),
    dtype=np.int16,
)
MATERIAL_TO_INDEX = {m: i for i, m in enumerate(MATERIALS)}


# =============================================================================
# 2. PROPRIEDADES DETERMINÍSTICAS — MESMA BASE DA FASE 1 FINAL
# =============================================================================

BASE_COST_PER_KG = {
    "SiC": 0.90,
    "UHMWPE": 1.33,
    "Al": 3.38,
    "B4C": 10.67,
    "Fe": 0.46,
    "W": 245.20,
    "Epoxy": 2.14,
    "Bi2O3": 21.83,
    "HDPE": 1.03,
}

PURE_DENSITY = {
    "SiC": 3.21,
    "UHMWPE": 0.96,
    "Al": 2.70,
    "Fe": 7.87,
    "B4C": 2.52,
    "W": 19.30,
    "Epoxy": 1.15,
    "Bi2O3": 8.90,
    "HDPE": 0.96,
}


def _rho_mass_fraction(components: dict[str, float]) -> float:
    return 1.0 / sum(w / PURE_DENSITY[k] for k, w in components.items())


def _rho_volume_fraction(components: dict[str, float]) -> float:
    return sum(v * PURE_DENSITY[k] for k, v in components.items())


def _cost_mass_fraction(components: dict[str, float]) -> float:
    return sum(w * BASE_COST_PER_KG[k] for k, w in components.items())


def _cost_volume_fraction(components: dict[str, float]) -> float:
    masses_g = {k: v * PURE_DENSITY[k] for k, v in components.items()}
    total_mass_g = sum(masses_g.values())
    return sum(
        (masses_g[k] / total_mass_g) * BASE_COST_PER_KG[k]
        for k in components
    )


MATERIAL_PROPERTIES = {
    "SiC": {
        "density": PURE_DENSITY["SiC"],
        "cost_kg": BASE_COST_PER_KG["SiC"],
    },
    "UHMWPE": {
        "density": PURE_DENSITY["UHMWPE"],
        "cost_kg": BASE_COST_PER_KG["UHMWPE"],
    },
    "Al": {
        "density": PURE_DENSITY["Al"],
        "cost_kg": BASE_COST_PER_KG["Al"],
    },
    "Fe": {
        "density": PURE_DENSITY["Fe"],
        "cost_kg": BASE_COST_PER_KG["Fe"],
    },
    "B4C-HDPE": {
        "density": _rho_volume_fraction({"B4C": 0.30, "HDPE": 0.70}),
        "cost_kg": _cost_volume_fraction({"B4C": 0.30, "HDPE": 0.70}),
    },
    "W-Epoxy": {
        "density": _rho_mass_fraction({"W": 0.80, "Epoxy": 0.20}),
        "cost_kg": _cost_mass_fraction({"W": 0.80, "Epoxy": 0.20}),
    },
    "Bi2O3-HDPE": {
        "density": _rho_mass_fraction({"Bi2O3": 0.70, "HDPE": 0.30}),
        "cost_kg": _cost_mass_fraction({"Bi2O3": 0.70, "HDPE": 0.30}),
    },
}

DENSITY_VECTOR = np.asarray(
    [MATERIAL_PROPERTIES[m]["density"] for m in MATERIALS],
    dtype=float,
)
COST_VECTOR = np.asarray(
    [MATERIAL_PROPERTIES[m]["cost_kg"] for m in MATERIALS],
    dtype=float,
)

# Congelamento da densidade adotada na tese: 4.643 g/cm³ (arredondado).
assert abs(MATERIAL_PROPERTIES["W-Epoxy"]["density"] - 4.643) < 0.001
assert len(MATERIALS) == N_LAYERS == 7


# =============================================================================
# 3. UTILITÁRIOS
# =============================================================================


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def dump_json(data: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    def _default(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, Path):
            return str(obj)
        raise TypeError(f"Objeto não serializável: {type(obj)}")

    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False, default=_default)


def resolve_first_existing(candidates: list[Path], label: str) -> Path:
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(
        f"{label} não encontrado. Caminhos procurados:\n"
        + "\n".join(f"  - {p}" for p in candidates)
    )


def atomic_pickle_dump(obj: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("wb") as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, path)


def stable_signature(payload: dict) -> str:
    text = json.dumps(payload, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def selected_run_parameters() -> tuple[int, int, int]:
    mode = RUN_MODE.upper()
    if mode == "FINAL":
        return FINAL_N_RUNS, FINAL_POPULATION_SIZE, FINAL_NUM_GENERATIONS
    if mode == "SMOKE_TEST":
        return SMOKE_N_RUNS, SMOKE_POPULATION_SIZE, SMOKE_NUM_GENERATIONS
    raise ValueError("RUN_MODE deve ser 'FINAL' ou 'SMOKE_TEST'.")


# =============================================================================
# 4. CARREGAMENTO E CONTRATO DOS SURROGATES
# =============================================================================


def load_surrogate_bundle(path: Path, expected_target: str) -> dict:
    if not path.exists():
        raise FileNotFoundError(
            f"Bundle ausente: {path}\n"
            "Execute primeiro a Fase 3 FINAL para gerar os dois surrogates."
        )

    bundle = joblib.load(path)
    if not isinstance(bundle, dict):
        raise TypeError(f"{path.name}: esperado dict, recebido {type(bundle)}.")

    required = {
        "target",
        "pipeline",
        "conformal_q_abs",
        "input_columns",
        "feature_set",
        "model_name",
        "test_metrics",
    }
    missing = sorted(required - set(bundle))
    if missing:
        raise ValueError(f"{path.name}: chaves ausentes: {missing}")

    if str(bundle["target"]) != expected_target:
        raise ValueError(
            f"{path.name}: target='{bundle['target']}', esperado='{expected_target}'."
        )

    if list(bundle["input_columns"]) != EXPECTED_INPUT_COLUMNS:
        raise ValueError(
            f"{path.name}: contrato de inputs incompatível com a Fase 4."
        )

    q = float(bundle["conformal_q_abs"])
    if not np.isfinite(q) or q < 0:
        raise ValueError(f"{path.name}: conformal_q_abs inválido: {q}")

    if not hasattr(bundle["pipeline"], "predict"):
        raise TypeError(f"{path.name}: pipeline não possui predict().")

    return bundle


# =============================================================================
# 5. REPRESENTAÇÃO DO PROJETO / VALIDAÇÃO
# =============================================================================


def project_box_simplex(
    v: np.ndarray,
    total: float = TOTAL_THICKNESS_CM,
    lo: float = MIN_THICKNESS_CM,
    hi: float = MAX_THICKNESS_CM,
) -> np.ndarray:
    """
    Projeção Euclidiana em:
        {x | sum(x)=total, lo<=x_i<=hi}
    via busca bissecional do multiplicador de Lagrange.
    """
    v = np.asarray(v, dtype=float).copy()

    n = len(v)
    if not (n * lo <= total <= n * hi):
        raise ValueError("Simplex com caixa inviável.")

    # x_i = clip(v_i - lambda, lo, hi)
    low = float(np.min(v - hi) - abs(total) - 1.0)
    high = float(np.max(v - lo) + abs(total) + 1.0)

    for _ in range(100):
        mid = 0.5 * (low + high)
        x = np.clip(v - mid, lo, hi)
        s = float(x.sum())
        if s > total:
            low = mid
        else:
            high = mid

    x = np.clip(v - 0.5 * (low + high), lo, hi)

    # Correção residual numérica sem sair dos limites.
    residual = total - float(x.sum())
    for _ in range(20):
        if abs(residual) <= 1e-12:
            break
        if residual > 0:
            slack = hi - x
        else:
            slack = x - lo
        idx = int(np.argmax(slack))
        step = math.copysign(min(abs(residual), float(slack[idx])), residual)
        x[idx] += step
        residual = total - float(x.sum())

    if not np.isclose(x.sum(), total, atol=1e-9, rtol=0.0):
        raise RuntimeError("Falha na projeção: soma != 6 cm.")
    if np.any(x < lo - 1e-10) or np.any(x > hi + 1e-10):
        raise RuntimeError("Falha na projeção: bounds violados.")

    return x


def quantize_thickness(x: np.ndarray) -> np.ndarray:
    """
    Quantiza para 0.001 cm e mantém soma exatamente 6.000 cm em unidades inteiras.
    """
    x = project_box_simplex(np.asarray(x, dtype=float))

    scale = int(round(1.0 / THICKNESS_RESOLUTION_CM))
    min_u = int(round(MIN_THICKNESS_CM * scale))
    max_u = int(round(MAX_THICKNESS_CM * scale))
    target_u = int(round(TOTAL_THICKNESS_CM * scale))

    raw = x * scale
    units = np.floor(raw + 1e-12).astype(int)
    units = np.clip(units, min_u, max_u)
    frac = raw - np.floor(raw)

    delta = target_u - int(units.sum())

    if delta > 0:
        order = np.argsort(frac)[::-1]
        while delta > 0:
            progressed = False
            for idx in order:
                if units[idx] < max_u:
                    units[idx] += 1
                    delta -= 1
                    progressed = True
                    if delta == 0:
                        break
            if not progressed:
                raise RuntimeError("Falha ao quantizar: sem slack superior.")

    elif delta < 0:
        order = np.argsort(frac)
        while delta < 0:
            progressed = False
            for idx in order:
                if units[idx] > min_u:
                    units[idx] -= 1
                    delta += 1
                    progressed = True
                    if delta == 0:
                        break
            if not progressed:
                raise RuntimeError("Falha ao quantizar: sem slack inferior.")

    out = units.astype(float) / scale

    if not np.isclose(out.sum(), TOTAL_THICKNESS_CM, atol=1e-12, rtol=0.0):
        raise RuntimeError("Quantização não preservou soma.")
    if np.any(out < MIN_THICKNESS_CM) or np.any(out > MAX_THICKNESS_CM):
        raise RuntimeError("Quantização violou bounds.")

    return out


def validate_population(perms: np.ndarray, thickness: np.ndarray) -> None:
    perms = np.asarray(perms)
    thickness = np.asarray(thickness, dtype=float)

    if perms.ndim != 2 or perms.shape[1] != N_LAYERS:
        raise ValueError("perms deve ter shape (N,7).")
    if thickness.shape != perms.shape:
        raise ValueError("thickness deve ter o mesmo shape de perms.")

    expected = np.arange(N_LAYERS)
    for i, row in enumerate(perms):
        if not np.array_equal(np.sort(row.astype(int)), expected):
            raise ValueError(f"Indivíduo {i}: materiais não formam permutação 0..6.")

    if np.any(thickness < MIN_THICKNESS_CM - 1e-12):
        raise ValueError("Espessura abaixo do mínimo.")
    if np.any(thickness > MAX_THICKNESS_CM + 1e-12):
        raise ValueError("Espessura acima do máximo.")
    if not np.allclose(
        thickness.sum(axis=1),
        TOTAL_THICKNESS_CM,
        atol=1e-10,
        rtol=0.0,
    ):
        raise ValueError("Soma de espessuras diferente de 6 cm.")


def design_dataframe(perms: np.ndarray, thickness: np.ndarray) -> pd.DataFrame:
    validate_population(perms, thickness)
    n = len(perms)
    data: dict[str, Any] = {}

    for i, col in enumerate(MATERIAL_COLS):
        data[col] = np.asarray(MATERIALS, dtype=object)[perms[:, i].astype(int)]

    for i, col in enumerate(THICKNESS_COLS):
        data[col] = thickness[:, i].astype(float)

    df = pd.DataFrame(data, index=np.arange(n))
    return df[EXPECTED_INPUT_COLUMNS]


# =============================================================================
# 6. MASSA E CUSTO EXATOS
# =============================================================================


def deterministic_mass_cost(
    perms: np.ndarray,
    thickness: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    perms = np.asarray(perms, dtype=int)
    thickness = np.asarray(thickness, dtype=float)
    validate_population(perms, thickness)

    # Raios de cada camada.
    cumulative_before = np.cumsum(thickness, axis=1) - thickness
    r_inner = INNER_RADIUS_CM + cumulative_before
    r_outer = r_inner + thickness

    volume_cm3 = (
        0.5
        * math.pi
        * (r_outer**2 - r_inner**2)
        * SHIELD_HEIGHT_CM
    )

    rho = DENSITY_VECTOR[perms]
    cost_kg = COST_VECTOR[perms]

    layer_mass_kg = rho * volume_cm3 / 1000.0
    total_mass_kg = layer_mass_kg.sum(axis=1)
    total_cost = (layer_mass_kg * cost_kg).sum(axis=1)

    return total_mass_kg, total_cost


# =============================================================================
# 7. AVALIAÇÃO VETORIZADA DOS SURROGATES
# =============================================================================


@dataclass
class Evaluator:
    dose_bundle: dict
    pen_bundle: dict

    def __post_init__(self):
        self.dose_pipeline = self.dose_bundle["pipeline"]
        self.pen_pipeline = self.pen_bundle["pipeline"]
        self.q_dose = float(self.dose_bundle["conformal_q_abs"])
        self.q_pen = float(self.pen_bundle["conformal_q_abs"])

    def evaluate(
        self,
        perms: np.ndarray,
        thickness: np.ndarray,
    ) -> dict[str, np.ndarray]:
        df = design_dataframe(perms, thickness)

        dose = np.asarray(self.dose_pipeline.predict(df), dtype=float).ravel()
        pen = np.asarray(self.pen_pipeline.predict(df), dtype=float).ravel()

        if len(dose) != len(df) or len(pen) != len(df):
            raise RuntimeError("Surrogate retornou dimensão incompatível.")

        dose = np.maximum(0.0, dose)
        pen = np.maximum(0.0, pen)

        mass, cost = deterministic_mass_cost(perms, thickness)

        dose_low = np.maximum(0.0, dose - self.q_dose)
        dose_up = dose + self.q_dose

        pen_low = np.maximum(0.0, pen - self.q_pen)
        pen_up = pen + self.q_pen

        cv = np.maximum(0.0, pen_up - BALLISTIC_LIMIT_MM)
        feasible = cv <= CONSTRAINT_TOL_MM

        F = np.column_stack([dose, mass, cost])

        arrays = {
            "F": F,
            "dose": dose,
            "dose_low": dose_low,
            "dose_up": dose_up,
            "mass": mass,
            "cost": cost,
            "pen": pen,
            "pen_low": pen_low,
            "pen_up": pen_up,
            "cv": cv,
            "feasible": feasible.astype(bool),
        }

        for k, arr in arrays.items():
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"Avaliação produziu NaN/Inf em {k}.")

        return arrays


def eval_take(evaluation: dict[str, np.ndarray], idx) -> dict[str, np.ndarray]:
    return {k: v[idx] for k, v in evaluation.items()}


def eval_concat(a: dict[str, np.ndarray], b: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {k: np.concatenate([a[k], b[k]], axis=0) for k in a}


# =============================================================================
# 8. INICIALIZAÇÃO
# =============================================================================


def sample_constrained_dirichlet(rng: np.random.Generator, n: int) -> np.ndarray:
    accepted = []
    available = TOTAL_THICKNESS_CM - N_LAYERS * MIN_THICKNESS_CM

    while sum(len(x) for x in accepted) < n:
        batch_n = max(64, 2 * (n - sum(len(x) for x in accepted)))
        frac = rng.dirichlet(np.ones(N_LAYERS), size=batch_n)
        x = MIN_THICKNESS_CM + available * frac
        mask = np.all(x <= MAX_THICKNESS_CM + 1e-12, axis=1)
        if np.any(mask):
            accepted.append(x[mask])

    out = np.vstack(accepted)[:n]
    out = np.vstack([quantize_thickness(row) for row in out])
    return out


def initialize_population(
    rng: np.random.Generator,
    n: int,
) -> tuple[np.ndarray, np.ndarray]:
    if n <= len(ALL_PERMUTATIONS):
        idx = rng.choice(len(ALL_PERMUTATIONS), size=n, replace=False)
    else:
        idx = rng.choice(len(ALL_PERMUTATIONS), size=n, replace=True)

    perms = ALL_PERMUTATIONS[idx].copy()
    thickness = sample_constrained_dirichlet(rng, n)
    validate_population(perms, thickness)
    return perms, thickness


# =============================================================================
# 9. OPERADORES GENÉTICOS DE VARIÁVEIS MISTAS
# =============================================================================


def order_crossover(
    p1: np.ndarray,
    p2: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Order Crossover (OX) preservando permutação."""
    n = len(p1)
    a, b = sorted(rng.choice(n, size=2, replace=False))
    b += 1  # intervalo [a,b)

    def make_child(base, donor):
        child = np.full(n, -1, dtype=np.int16)
        child[a:b] = base[a:b]
        used = set(int(x) for x in child[a:b])

        donor_order = [int(x) for x in donor if int(x) not in used]
        positions = list(range(b, n)) + list(range(0, a))
        for pos, gene in zip(positions, donor_order):
            child[pos] = gene
        return child

    c1 = make_child(p1, p2)
    c2 = make_child(p2, p1)

    if not np.array_equal(np.sort(c1), np.arange(n)):
        raise RuntimeError("OX produziu permutação inválida.")
    if not np.array_equal(np.sort(c2), np.arange(n)):
        raise RuntimeError("OX produziu permutação inválida.")

    return c1, c2


def bounded_sbx_pair(
    x1: np.ndarray,
    x2: np.ndarray,
    rng: np.random.Generator,
    eta: float = SBX_ETA,
) -> tuple[np.ndarray, np.ndarray]:
    """
    SBX limitado por gene, seguido por projeção/quantização da restrição de soma.
    """
    c1 = x1.astype(float).copy()
    c2 = x2.astype(float).copy()

    for i in range(N_LAYERS):
        if rng.random() > 0.5:
            continue

        y1 = float(min(x1[i], x2[i]))
        y2 = float(max(x1[i], x2[i]))

        if abs(y1 - y2) <= 1e-14:
            continue

        yl = MIN_THICKNESS_CM
        yu = MAX_THICKNESS_CM
        rand = float(rng.random())

        beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1))
        alpha = 2.0 - beta ** (-(eta + 1.0))
        if rand <= 1.0 / alpha:
            betaq = (rand * alpha) ** (1.0 / (eta + 1.0))
        else:
            betaq = (1.0 / (2.0 - rand * alpha)) ** (1.0 / (eta + 1.0))

        child1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1))

        beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1))
        alpha = 2.0 - beta ** (-(eta + 1.0))
        if rand <= 1.0 / alpha:
            betaq = (rand * alpha) ** (1.0 / (eta + 1.0))
        else:
            betaq = (1.0 / (2.0 - rand * alpha)) ** (1.0 / (eta + 1.0))

        child2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1))

        child1 = float(np.clip(child1, yl, yu))
        child2 = float(np.clip(child2, yl, yu))

        if rng.random() <= 0.5:
            c1[i], c2[i] = child2, child1
        else:
            c1[i], c2[i] = child1, child2

    c1 = quantize_thickness(c1)
    c2 = quantize_thickness(c2)
    return c1, c2


def polynomial_mutation(
    x: np.ndarray,
    rng: np.random.Generator,
    eta: float = POLYNOMIAL_MUTATION_ETA,
    gene_prob: float = THICKNESS_GENE_MUTATION_PROB,
) -> np.ndarray:
    y = x.astype(float).copy()
    yl = MIN_THICKNESS_CM
    yu = MAX_THICKNESS_CM

    for i in range(N_LAYERS):
        if rng.random() >= gene_prob:
            continue

        value = float(y[i])
        delta1 = (value - yl) / (yu - yl)
        delta2 = (yu - value) / (yu - yl)
        rand = float(rng.random())
        mut_pow = 1.0 / (eta + 1.0)

        if rand < 0.5:
            xy = 1.0 - delta1
            val = 2.0 * rand + (1.0 - 2.0 * rand) * (xy ** (eta + 1.0))
            deltaq = val**mut_pow - 1.0
        else:
            xy = 1.0 - delta2
            val = (
                2.0 * (1.0 - rand)
                + 2.0 * (rand - 0.5) * (xy ** (eta + 1.0))
            )
            deltaq = 1.0 - val**mut_pow

        value += deltaq * (yu - yl)
        y[i] = np.clip(value, yl, yu)

    return quantize_thickness(y)


def mutate_permutation(
    perm: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    out = perm.copy()
    if rng.random() < PERM_SWAP_MUTATION_PROB:
        i, j = rng.choice(N_LAYERS, size=2, replace=False)
        out[i], out[j] = out[j], out[i]
    return out


def make_offspring(
    perms: np.ndarray,
    thickness: np.ndarray,
    ranks: np.ndarray,
    crowding: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    n = len(perms)

    def tournament() -> int:
        a, b = rng.choice(n, size=2, replace=False)
        if ranks[a] < ranks[b]:
            return int(a)
        if ranks[b] < ranks[a]:
            return int(b)
        if crowding[a] > crowding[b]:
            return int(a)
        if crowding[b] > crowding[a]:
            return int(b)
        return int(a if rng.random() < 0.5 else b)

    child_perms = []
    child_thickness = []

    while len(child_perms) < n:
        ia = tournament()
        ib = tournament()

        p1p, p2p = perms[ia], perms[ib]
        p1t, p2t = thickness[ia], thickness[ib]

        if rng.random() < CROSSOVER_PROB:
            c1p, c2p = order_crossover(p1p, p2p, rng)
            c1t, c2t = bounded_sbx_pair(p1t, p2t, rng)
        else:
            c1p, c2p = p1p.copy(), p2p.copy()
            c1t, c2t = p1t.copy(), p2t.copy()

        c1p = mutate_permutation(c1p, rng)
        c2p = mutate_permutation(c2p, rng)
        c1t = polynomial_mutation(c1t, rng)
        c2t = polynomial_mutation(c2t, rng)

        child_perms.extend([c1p, c2p])
        child_thickness.extend([c1t, c2t])

    cp = np.asarray(child_perms[:n], dtype=np.int16)
    ct = np.asarray(child_thickness[:n], dtype=float)
    validate_population(cp, ct)
    return cp, ct


# =============================================================================
# 10. CONSTRAINT-DOMINATION / NSGA-II
# =============================================================================


def dominance_matrix(F: np.ndarray, cv: np.ndarray) -> np.ndarray:
    """
    Matriz D onde D[i,j]=True se i domina j segundo regras de viabilidade de Deb.
    """
    F = np.asarray(F, dtype=float)
    cv = np.asarray(cv, dtype=float).ravel()

    feasible = cv <= CONSTRAINT_TOL_MM
    fi = feasible[:, None]
    fj = feasible[None, :]

    A = F[:, None, :]
    B = F[None, :, :]

    pareto = np.all(A <= B + 1e-12, axis=2) & np.any(A < B - 1e-12, axis=2)

    better_infeasible = (
        (~fi)
        & (~fj)
        & (cv[:, None] < cv[None, :] - 1e-12)
    )

    dom = (
        (fi & (~fj))
        | (fi & fj & pareto)
        | better_infeasible
    )

    np.fill_diagonal(dom, False)
    return dom


def fast_non_dominated_sort(
    F: np.ndarray,
    cv: np.ndarray,
) -> tuple[list[np.ndarray], np.ndarray]:
    """
    Ordenação não-dominada vetorizada.
    """
    n = len(F)
    dom = dominance_matrix(F, cv)
    dominated_count = dom.sum(axis=0).astype(int)

    remaining = np.ones(n, dtype=bool)
    ranks = np.full(n, -1, dtype=int)
    fronts: list[np.ndarray] = []
    rank = 0

    while np.any(remaining):
        current = np.where(remaining & (dominated_count == 0))[0]

        if len(current) == 0:
            # Não deveria ocorrer em uma relação de dominância válida.
            current = np.where(remaining)[0]

        fronts.append(current)
        ranks[current] = rank
        remaining[current] = False

        decrement = dom[current].sum(axis=0).astype(int)
        dominated_count = dominated_count - decrement
        dominated_count[~remaining] = -1

        rank += 1

    return fronts, ranks


def crowding_distance_for_front(
    F: np.ndarray,
    front: np.ndarray,
) -> np.ndarray:
    """
    Crowding distance NSGA-II normalizada por amplitude de cada objetivo.
    Retorna distâncias na ordem dos índices de 'front'.
    """
    front = np.asarray(front, dtype=int)
    n = len(front)

    if n == 0:
        return np.asarray([], dtype=float)
    if n <= 2:
        return np.full(n, np.inf, dtype=float)

    values = F[front]
    distances = np.zeros(n, dtype=float)

    for j in range(F.shape[1]):
        order = np.argsort(values[:, j], kind="mergesort")
        v = values[order, j]
        span = float(v[-1] - v[0])

        if span <= 1e-15:
            continue

        distances[order[0]] = np.inf
        distances[order[-1]] = np.inf

        for k in range(1, n - 1):
            idx = order[k]
            if not np.isinf(distances[idx]):
                distances[idx] += (v[k + 1] - v[k - 1]) / span

    return distances


def ranks_and_crowding(
    F: np.ndarray,
    cv: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    fronts, ranks = fast_non_dominated_sort(F, cv)
    crowding = np.zeros(len(F), dtype=float)

    for front in fronts:
        d = crowding_distance_for_front(F, front)
        crowding[front] = d

    return ranks, crowding, fronts


def environmental_selection(
    perms: np.ndarray,
    thickness: np.ndarray,
    evaluation: dict[str, np.ndarray],
    n_survive: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    F = evaluation["F"]
    cv = evaluation["cv"]
    fronts, _ = fast_non_dominated_sort(F, cv)

    selected: list[int] = []

    for front in fronts:
        if len(selected) + len(front) <= n_survive:
            selected.extend(front.tolist())
            continue

        slots = n_survive - len(selected)
        if slots <= 0:
            break

        d = crowding_distance_for_front(F, front)
        order = np.argsort(d, kind="mergesort")[::-1]
        selected.extend(front[order[:slots]].tolist())
        break

    if len(selected) != n_survive:
        raise RuntimeError(
            f"Seleção ambiental retornou {len(selected)}/{n_survive} indivíduos."
        )

    idx = np.asarray(selected, dtype=int)
    return perms[idx], thickness[idx], eval_take(evaluation, idx)


# =============================================================================
# 11. HYPERVOLUME 3D EXATO
# =============================================================================


def pure_pareto_indices(F: np.ndarray) -> np.ndarray:
    """Índices não-dominados para um conjunto já viável."""
    F = np.asarray(F, dtype=float)
    if len(F) == 0:
        return np.asarray([], dtype=int)
    zeros = np.zeros(len(F), dtype=float)
    fronts, _ = fast_non_dominated_sort(F, zeros)
    return fronts[0]


def _union_area_2d(points_yz: np.ndarray, ref_yz: np.ndarray) -> float:
    """
    Área exata da união de retângulos [y_i, ref_y] x [z_i, ref_z].
    """
    if len(points_yz) == 0:
        return 0.0

    ry, rz = map(float, ref_yz)
    p = points_yz[
        (points_yz[:, 0] < ry)
        & (points_yz[:, 1] < rz)
    ]
    if len(p) == 0:
        return 0.0

    ys = np.sort(np.unique(np.concatenate([p[:, 0], [ry]])))
    area = 0.0
    z_best = rz

    for i in range(len(ys) - 1):
        y = ys[i]
        active = p[p[:, 0] <= y + 1e-15]
        if len(active):
            z_best = min(z_best, float(active[:, 1].min()))
        width = float(ys[i + 1] - y)
        height = max(0.0, rz - z_best)
        area += width * height

    return float(area)


def hypervolume_3d(points: np.ndarray, reference: np.ndarray) -> float:
    """
    Hypervolume exato para três objetivos de minimização.
    """
    points = np.asarray(points, dtype=float)
    reference = np.asarray(reference, dtype=float)

    if points.size == 0:
        return 0.0
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("hypervolume_3d requer array (N,3).")

    mask = np.all(points < reference[None, :], axis=1)
    p = points[mask]
    if len(p) == 0:
        return 0.0

    p = p[pure_pareto_indices(p)]

    xs = np.sort(np.unique(np.concatenate([p[:, 0], [reference[0]]])))
    hv = 0.0

    for i in range(len(xs) - 1):
        x = xs[i]
        active = p[p[:, 0] <= x + 1e-15]
        if len(active) == 0:
            continue

        dx = float(xs[i + 1] - x)
        area = _union_area_2d(active[:, 1:3], reference[1:3])
        hv += dx * area

    return float(hv)


def build_hv_anchor(F_doe: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ideal = np.min(F_doe, axis=0)
    nadir_anchor = np.max(F_doe, axis=0)
    span = np.maximum(nadir_anchor - ideal, np.maximum(np.abs(nadir_anchor), 1.0) * 1e-6)
    reference = nadir_anchor + REFERENCE_MARGIN * span
    return ideal, reference


def normalized_hypervolume(
    F: np.ndarray,
    ideal_anchor: np.ndarray,
    reference_anchor: np.ndarray,
) -> float:
    if len(F) == 0:
        return 0.0

    denom = reference_anchor - ideal_anchor
    scaled = (F - ideal_anchor[None, :]) / denom[None, :]
    scaled = np.clip(scaled, 0.0, 1.0)
    return hypervolume_3d(scaled, np.ones(3, dtype=float))


# =============================================================================
# 12. MÉTRICAS DE CONVERGÊNCIA
# =============================================================================


def spacing_metric(F: np.ndarray) -> float:
    """
    Schott-like spacing no espaço objetivo normalizado.
    Menor = espaçamento mais uniforme. Apenas diagnóstico.
    """
    F = np.asarray(F, dtype=float)
    if len(F) < 3:
        return float("nan")

    lo = np.min(F, axis=0)
    hi = np.max(F, axis=0)
    span = np.where((hi - lo) > 1e-15, hi - lo, 1.0)
    z = (F - lo) / span

    diff = z[:, None, :] - z[None, :, :]
    dist = np.linalg.norm(diff, axis=2)
    np.fill_diagonal(dist, np.inf)
    nearest = np.min(dist, axis=1)

    return float(np.std(nearest, ddof=1))


def normalize_objectives(
    F: np.ndarray,
    ideal_anchor: np.ndarray,
    reference_anchor: np.ndarray,
) -> np.ndarray:
    """Normaliza objetivos na âncora fixa derivada do DoE."""
    F = np.asarray(F, dtype=float)
    denom = reference_anchor - ideal_anchor
    return np.clip(
        (F - ideal_anchor[None, :]) / denom[None, :],
        0.0,
        1.0,
    )


def inverted_generational_distance(
    approximation: np.ndarray,
    reference_front: np.ndarray,
) -> float:
    """
    IGD no espaço objetivo já normalizado. Menor = melhor.

    A frente de referência é uma aproximação empírica (a união Pareto das
    réplicas), não uma frente analítica verdadeira desconhecida.
    """
    approximation = np.asarray(approximation, dtype=float)
    reference_front = np.asarray(reference_front, dtype=float)
    if len(approximation) == 0 or len(reference_front) == 0:
        return float("nan")

    distances = np.linalg.norm(
        reference_front[:, None, :] - approximation[None, :, :],
        axis=2,
    )
    return float(np.mean(np.min(distances, axis=1)))


def generalized_spread(
    approximation: np.ndarray,
    reference_front: np.ndarray,
) -> float:
    """
    Generalized spread de Deb et al. para M objetivos. Menor = melhor.

    Mede tanto a cobertura dos extremos da frente de referência quanto a
    uniformidade dos vizinhos mais próximos. É preferido ao spread 2D nesta
    otimização de três objetivos.
    """
    approximation = np.unique(np.asarray(approximation, dtype=float), axis=0)
    reference_front = np.asarray(reference_front, dtype=float)
    if len(approximation) < 2 or len(reference_front) == 0:
        return float("nan")

    pair_dist = np.linalg.norm(
        approximation[:, None, :] - approximation[None, :, :],
        axis=2,
    )
    np.fill_diagonal(pair_dist, np.inf)
    nearest = np.min(pair_dist, axis=1)
    mean_nearest = float(np.mean(nearest))

    # Um extremo da frente de referência por objetivo de minimização.
    extreme_reference = reference_front[
        np.argmin(reference_front, axis=0)
    ]
    extreme_distances = np.linalg.norm(
        extreme_reference[:, None, :] - approximation[None, :, :],
        axis=2,
    )
    edge_error = float(np.sum(np.min(extreme_distances, axis=1)))
    denominator = edge_error + len(approximation) * mean_nearest
    if denominator <= 1e-15:
        return 0.0
    return float(
        (edge_error + np.sum(np.abs(nearest - mean_nearest))) / denominator
    )


def generation_metrics(
    generation: int,
    evaluation: dict[str, np.ndarray],
    ideal_anchor: np.ndarray,
    reference_anchor: np.ndarray,
) -> dict:
    feasible_idx = np.where(evaluation["feasible"])[0]

    if len(feasible_idx):
        Ff = evaluation["F"][feasible_idx]
        pidx = pure_pareto_indices(Ff)
        front = Ff[pidx]

        hv_raw = hypervolume_3d(front, reference_anchor)
        hv_norm = normalized_hypervolume(front, ideal_anchor, reference_anchor)
        spacing = spacing_metric(front)
        pareto_size = len(front)
    else:
        hv_raw = 0.0
        hv_norm = 0.0
        spacing = float("nan")
        pareto_size = 0

    return {
        "generation": int(generation),
        "feasible_count": int(len(feasible_idx)),
        "feasible_fraction": float(len(feasible_idx) / len(evaluation["F"])),
        "pareto_size": int(pareto_size),
        "hv_raw": float(hv_raw),
        "hv_normalized": float(hv_norm),
        "spacing_metric": spacing,
        "min_constraint_violation_mm": float(np.min(evaluation["cv"])),
        "median_constraint_violation_mm": float(np.median(evaluation["cv"])),
    }


# =============================================================================
# 13. CHECKPOINT / EXECUÇÃO DE UMA RÉPLICA
# =============================================================================


def make_run_signature(
    seed: int,
    pop_size: int,
    num_generations: int,
    model_hashes: dict,
) -> tuple[dict, str]:
    payload = {
        "script_version": SCRIPT_VERSION,
        "seed": seed,
        "population_size": pop_size,
        "num_generations": num_generations,
        "checkpoint_interval": CHECKPOINT_INTERVAL,
        "crossover_prob": CROSSOVER_PROB,
        "perm_swap_mutation_prob": PERM_SWAP_MUTATION_PROB,
        "thickness_gene_mutation_prob": THICKNESS_GENE_MUTATION_PROB,
        "sbx_eta": SBX_ETA,
        "poly_mut_eta": POLYNOMIAL_MUTATION_ETA,
        "bounds_cm": [MIN_THICKNESS_CM, MAX_THICKNESS_CM],
        "total_thickness_cm": TOTAL_THICKNESS_CM,
        "resolution_cm": THICKNESS_RESOLUTION_CM,
        "ballistic_limit_mm": BALLISTIC_LIMIT_MM,
        "model_hashes": model_hashes,
    }
    return payload, stable_signature(payload)


def save_checkpoint(
    path: Path,
    signature: str,
    generation_completed: int,
    perms: np.ndarray,
    thickness: np.ndarray,
    evaluation: dict[str, np.ndarray],
    history: list[dict],
    rng: np.random.Generator,
) -> None:
    payload = {
        "signature": signature,
        "generation_completed": generation_completed,
        "perms": perms,
        "thickness": thickness,
        "evaluation": evaluation,
        "history": history,
        "rng_state": rng.bit_generator.state,
    }
    atomic_pickle_dump(payload, path)


def load_checkpoint(path: Path, signature: str):
    with path.open("rb") as f:
        data = pickle.load(f)

    if data.get("signature") != signature:
        raise ValueError(
            "Checkpoint pertence a outra configuração/modelo. "
            "Use FORCE_NEW=True ou remova o checkpoint."
        )
    return data


def dataframe_from_population(
    perms: np.ndarray,
    thickness: np.ndarray,
    evaluation: dict[str, np.ndarray],
    run_id: int | None = None,
    seed: int | None = None,
) -> pd.DataFrame:
    rows = []

    for k in range(len(perms)):
        row: dict[str, Any] = {}

        if run_id is not None:
            row["run_id"] = int(run_id)
        if seed is not None:
            row["seed"] = int(seed)

        for i in range(N_LAYERS):
            row[f"material_camada_{i+1}"] = MATERIALS[int(perms[k, i])]
            row[f"espessura_camada_{i+1}"] = round(float(thickness[k, i]), 3)

        row.update(
            {
                "dose_pred_5min_mSv": float(evaluation["dose"][k]),
                "dose_pi_lower95_mSv": float(evaluation["dose_low"][k]),
                "dose_pi_upper95_mSv": float(evaluation["dose_up"][k]),
                "massa_total_kg": float(evaluation["mass"][k]),
                "custo_total_unidade": float(evaluation["cost"][k]),
                "max_penetration_pred_mm": float(evaluation["pen"][k]),
                "penetration_pi_lower95_mm": float(evaluation["pen_low"][k]),
                "penetration_pi_upper95_mm": float(evaluation["pen_up"][k]),
                "ballistic_constraint_violation_mm": float(evaluation["cv"][k]),
                "ballistic_feasible_95PI": int(bool(evaluation["feasible"][k])),
            }
        )
        rows.append(row)

    return pd.DataFrame(rows)


def extract_feasible_pareto(
    perms: np.ndarray,
    thickness: np.ndarray,
    evaluation: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    feasible = np.where(evaluation["feasible"])[0]
    if len(feasible) == 0:
        return (
            np.empty((0, N_LAYERS), dtype=np.int16),
            np.empty((0, N_LAYERS), dtype=float),
            eval_take(evaluation, np.asarray([], dtype=int)),
        )

    Ff = evaluation["F"][feasible]
    local = pure_pareto_indices(Ff)
    idx = feasible[local]

    return perms[idx], thickness[idx], eval_take(evaluation, idx)


def run_one_seed(
    run_id: int,
    seed: int,
    pop_size: int,
    num_generations: int,
    evaluator: Evaluator,
    ideal_anchor: np.ndarray,
    reference_anchor: np.ndarray,
    model_hashes: dict,
    output_dir: Path,
) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    run_dir = output_dir / f"run_{run_id:02d}_seed_{seed}"
    run_dir.mkdir(parents=True, exist_ok=True)

    checkpoint_path = run_dir / "checkpoint.pkl"
    complete_path = run_dir / "run_complete.json"
    pareto_path = run_dir / "pareto_final.csv"
    history_path = run_dir / "convergence_history.csv"

    sig_payload, signature = make_run_signature(
        seed, pop_size, num_generations, model_hashes
    )

    # Execução já concluída e compatível.
    if RESUME and not FORCE_NEW and complete_path.exists() and pareto_path.exists():
        with complete_path.open("r", encoding="utf-8") as f:
            done = json.load(f)
        if done.get("signature") == signature:
            print(f"Run {run_id:02d}: já concluído -> reutilizando resultados.")
            pareto_df = pd.read_csv(pareto_path)
            history_df = pd.read_csv(history_path)
            return pareto_df, done["summary"], history_df

    rng = np.random.default_rng(seed)
    history: list[dict] = []

    start_time = time.time()

    if RESUME and not FORCE_NEW and checkpoint_path.exists():
        data = load_checkpoint(checkpoint_path, signature)
        perms = data["perms"]
        thickness = data["thickness"]
        evaluation = data["evaluation"]
        history = data["history"]
        rng.bit_generator.state = data["rng_state"]
        start_gen = int(data["generation_completed"]) + 1
        print(
            f"Run {run_id:02d} seed={seed}: retomando da geração "
            f"{start_gen}/{num_generations}."
        )
    else:
        perms, thickness = initialize_population(rng, pop_size)
        evaluation = evaluator.evaluate(perms, thickness)
        start_gen = 1

        m0 = generation_metrics(
            0, evaluation, ideal_anchor, reference_anchor
        )
        history.append(m0)

        save_checkpoint(
            checkpoint_path,
            signature,
            0,
            perms,
            thickness,
            evaluation,
            history,
            rng,
        )

        print(
            f"Run {run_id:02d} seed={seed}: geração 0 | "
            f"viáveis={m0['feasible_count']}/{pop_size} | "
            f"HVn={m0['hv_normalized']:.6f}"
        )

    for generation in range(start_gen, num_generations + 1):
        ranks, crowding, _ = ranks_and_crowding(
            evaluation["F"], evaluation["cv"]
        )

        child_perm, child_thick = make_offspring(
            perms, thickness, ranks, crowding, rng
        )
        child_eval = evaluator.evaluate(child_perm, child_thick)

        combined_perm = np.vstack([perms, child_perm])
        combined_thick = np.vstack([thickness, child_thick])
        combined_eval = eval_concat(evaluation, child_eval)

        perms, thickness, evaluation = environmental_selection(
            combined_perm,
            combined_thick,
            combined_eval,
            pop_size,
        )

        if generation % METRIC_INTERVAL == 0 or generation == num_generations:
            m = generation_metrics(
                generation,
                evaluation,
                ideal_anchor,
                reference_anchor,
            )
            history.append(m)

            print(
                f"Run {run_id:02d} | gen {generation:03d}/{num_generations} | "
                f"viáveis={m['feasible_count']:3d}/{pop_size} | "
                f"Pareto={m['pareto_size']:3d} | "
                f"HVn={m['hv_normalized']:.6f} | "
                f"CVmin={m['min_constraint_violation_mm']:.4f} mm"
            )

        if (
            generation % CHECKPOINT_INTERVAL == 0
            or generation == num_generations
        ):
            save_checkpoint(
                checkpoint_path,
                signature,
                generation,
                perms,
                thickness,
                evaluation,
                history,
                rng,
            )
            pd.DataFrame(history).to_csv(history_path, index=False)

    runtime_s = time.time() - start_time

    p_perm, p_thick, p_eval = extract_feasible_pareto(
        perms, thickness, evaluation
    )

    if len(p_perm):
        pareto_df = dataframe_from_population(
            p_perm,
            p_thick,
            p_eval,
            run_id=run_id,
            seed=seed,
        )
        pareto_df.to_csv(pareto_path, index=False)
    else:
        pareto_df = pd.DataFrame()
        # Salvar os menos inviáveis para diagnóstico.
        best_idx = np.argsort(evaluation["cv"])[: min(20, len(perms))]
        dataframe_from_population(
            perms[best_idx],
            thickness[best_idx],
            eval_take(evaluation, best_idx),
            run_id=run_id,
            seed=seed,
        ).to_csv(run_dir / "best_infeasible.csv", index=False)

    final_metrics = generation_metrics(
        num_generations,
        evaluation,
        ideal_anchor,
        reference_anchor,
    )

    summary = {
        "run_id": run_id,
        "seed": seed,
        "runtime_s": runtime_s,
        "population_size": pop_size,
        "num_generations": num_generations,
        "final_feasible_count": final_metrics["feasible_count"],
        "final_feasible_fraction": final_metrics["feasible_fraction"],
        "final_pareto_size": final_metrics["pareto_size"],
        "final_hv_raw": final_metrics["hv_raw"],
        "final_hv_normalized": final_metrics["hv_normalized"],
        "final_spacing_metric": final_metrics["spacing_metric"],
        "min_constraint_violation_mm": final_metrics["min_constraint_violation_mm"],
    }

    complete = {
        "signature": signature,
        "signature_payload": sig_payload,
        "completed_utc": utc_now(),
        "summary": summary,
    }
    dump_json(complete, complete_path)
    pd.DataFrame(history).to_csv(history_path, index=False)

    return pareto_df, summary, pd.DataFrame(history)


# =============================================================================
# 14. DOE COMO ÂNCORA / BASELINE INTERNO
# =============================================================================


def load_doe(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "sim_id" in df.columns and "id_simulacao" not in df.columns:
        df = df.rename(columns={"sim_id": "id_simulacao"})
    thickness_cm = {
        f"espessura_camada_{i}_cm": f"espessura_camada_{i}"
        for i in range(1, N_LAYERS + 1)
    }
    if all(source in df.columns for source in thickness_cm):
        df = df.rename(columns=thickness_cm)

    if "id_simulacao" not in df.columns:
        raise ValueError("DoE sem id_simulacao.")

    keep = ["id_simulacao"] + EXPECTED_INPUT_COLUMNS
    missing = [c for c in keep if c not in df.columns]
    if missing:
        raise ValueError(f"DoE: colunas ausentes: {missing}")

    out = df[keep].copy()
    out["id_simulacao"] = pd.to_numeric(
        out["id_simulacao"], errors="raise"
    ).astype(int)

    if out["id_simulacao"].duplicated().any():
        raise ValueError("DoE contém IDs duplicados.")

    for c in THICKNESS_COLS:
        out[c] = pd.to_numeric(out[c], errors="raise").astype(float)

    validate_design_dataframe(out[EXPECTED_INPUT_COLUMNS])
    return out


def arrays_from_doe(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    perms = np.zeros((len(df), N_LAYERS), dtype=np.int16)
    thickness = df[THICKNESS_COLS].to_numpy(dtype=float)

    for i, c in enumerate(MATERIAL_COLS):
        vals = df[c].astype(str).map(MATERIAL_TO_INDEX)
        if vals.isna().any():
            bad = df.loc[vals.isna(), c].unique().tolist()
            raise ValueError(f"Material desconhecido no DoE: {bad}")
        perms[:, i] = vals.to_numpy(dtype=np.int16)

    # O DoE já está em 0.001 cm; não alteramos as configurações.
    validate_population(perms, thickness)
    return perms, thickness


def evaluate_doe_baseline(
    doe: pd.DataFrame,
    evaluator: Evaluator,
    output_dir: Path,
) -> tuple[
    np.ndarray,
    np.ndarray,
    pd.DataFrame,
    np.ndarray,
    np.ndarray,
    float,
]:
    perms, thickness = arrays_from_doe(doe)
    ev = evaluator.evaluate(perms, thickness)

    full = dataframe_from_population(perms, thickness, ev)
    full.insert(0, "id_simulacao", doe["id_simulacao"].to_numpy())
    full.to_csv(output_dir / "doe_evaluado_pelos_surrogates.csv", index=False)

    ideal_anchor, reference_anchor = build_hv_anchor(ev["F"])

    p_perm, p_thick, p_eval = extract_feasible_pareto(
        perms, thickness, ev
    )
    if len(p_perm):
        pareto_doe = dataframe_from_population(p_perm, p_thick, p_eval)
        pareto_doe.to_csv(output_dir / "pareto_doe_surrogate_baseline.csv", index=False)
        hv_doe = normalized_hypervolume(
            p_eval["F"], ideal_anchor, reference_anchor
        )
    else:
        pareto_doe = pd.DataFrame()
        hv_doe = 0.0

    return (
        ideal_anchor,
        reference_anchor,
        full,
        perms,
        thickness,
        hv_doe,
    )


# =============================================================================
# 15. UNIÃO MULTIRUN / TOPSIS / DISTÂNCIA AO DOE
# =============================================================================


def deduplicate_designs(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    keys = list(MATERIAL_COLS) + list(THICKNESS_COLS)
    out = df.copy()

    for c in THICKNESS_COLS:
        out[c] = out[c].astype(float).round(3)

    out = out.drop_duplicates(subset=keys, keep="first").reset_index(drop=True)
    return out


def nondominated_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    objective_cols = [
        "dose_pred_5min_mSv",
        "massa_total_kg",
        "custo_total_unidade",
    ]
    F = df[objective_cols].to_numpy(dtype=float)
    idx = pure_pareto_indices(F)
    return df.iloc[idx].reset_index(drop=True)


def topsis_iqr(
    F: np.ndarray,
    weights: np.ndarray = TOPSIS_WEIGHTS,
) -> tuple[int, np.ndarray]:
    """
    TOPSIS com normalização robusta por IQR, todos os objetivos de minimização.
    """
    F = np.asarray(F, dtype=float)
    if len(F) == 0:
        raise ValueError("TOPSIS recebeu frente vazia.")
    if len(F) == 1:
        return 0, np.asarray([1.0])

    q1 = np.percentile(F, 25, axis=0)
    q3 = np.percentile(F, 75, axis=0)
    iqr = q3 - q1
    iqr = np.where(iqr > 1e-15, iqr, 1.0)

    Z = (F - q1[None, :]) / iqr[None, :]

    w = np.asarray(weights, dtype=float)
    if len(w) != F.shape[1] or np.any(w < 0) or np.sum(w) <= 0:
        raise ValueError("Pesos TOPSIS inválidos.")
    w = w / w.sum()
    Z = Z * w[None, :]

    ideal = np.min(Z, axis=0)
    anti = np.max(Z, axis=0)

    d_ideal = np.linalg.norm(Z - ideal[None, :], axis=1)
    d_anti = np.linalg.norm(Z - anti[None, :], axis=1)

    score = d_anti / np.maximum(d_ideal + d_anti, 1e-15)
    best = int(np.argmax(score))
    return best, score


def add_domain_diagnostics(
    pareto_df: pd.DataFrame,
    doe_perms: np.ndarray,
    doe_thickness: np.ndarray,
) -> pd.DataFrame:
    if pareto_df.empty:
        return pareto_df.copy()

    out = pareto_df.copy()

    doe_t = np.asarray(doe_thickness, dtype=float)
    t_scale = MAX_THICKNESS_CM - MIN_THICKNESS_CM

    hamming_min = []
    thick_rmse_at_nearest = []
    mixed_min = []
    permutation_seen = []
    same_perm_thick_rmse = []

    for _, row in out.iterrows():
        p = np.asarray(
            [MATERIAL_TO_INDEX[str(row[c])] for c in MATERIAL_COLS],
            dtype=int,
        )
        t = np.asarray([float(row[c]) for c in THICKNESS_COLS], dtype=float)

        dmat = np.mean(doe_perms != p[None, :], axis=1)
        dth = np.sqrt(
            np.mean(((doe_t - t[None, :]) / t_scale) ** 2, axis=1)
        )
        dmixed = np.sqrt((dmat**2 + dth**2) / 2.0)

        j = int(np.argmin(dmixed))
        same = np.all(doe_perms == p[None, :], axis=1)

        hamming_min.append(float(dmat[j]))
        thick_rmse_at_nearest.append(float(dth[j] * t_scale))
        mixed_min.append(float(dmixed[j]))
        permutation_seen.append(int(np.any(same)))

        if np.any(same):
            rmse = np.sqrt(
                np.mean((doe_t[same] - t[None, :]) ** 2, axis=1)
            )
            same_perm_thick_rmse.append(float(np.min(rmse)))
        else:
            same_perm_thick_rmse.append(float("nan"))

    out["doe_min_material_hamming_fraction"] = hamming_min
    out["doe_nearest_thickness_rmse_cm"] = thick_rmse_at_nearest
    out["doe_mixed_distance"] = mixed_min
    out["material_permutation_seen_in_doe"] = permutation_seen
    out["same_permutation_min_thickness_rmse_cm"] = same_perm_thick_rmse

    return out


def select_validation_candidates(
    pareto_df: pd.DataFrame,
    n_select: int = N_HIGH_FIDELITY_CANDIDATES,
) -> pd.DataFrame:
    if pareto_df.empty:
        return pareto_df.copy()

    df = pareto_df.copy().reset_index(drop=True)
    objective_cols = [
        "dose_pred_5min_mSv",
        "massa_total_kg",
        "custo_total_unidade",
    ]
    F = df[objective_cols].to_numpy(dtype=float)

    best_topsis, scores = topsis_iqr(F)
    df["topsis_iqr_score"] = scores

    selected_labels: dict[int, list[str]] = {}

    def add(idx: int, label: str):
        selected_labels.setdefault(int(idx), []).append(label)

    add(best_topsis, "TOPSIS_IQR")
    add(int(np.argmin(F[:, 0])), "MIN_DOSE")
    add(int(np.argmin(F[:, 1])), "MIN_MASS")
    add(int(np.argmin(F[:, 2])), "MIN_COST")

    lo = np.min(F, axis=0)
    hi = np.max(F, axis=0)
    span = np.where((hi - lo) > 1e-15, hi - lo, 1.0)
    Z = (F - lo) / span

    while len(selected_labels) < min(n_select, len(df)):
        selected_idx = np.asarray(sorted(selected_labels), dtype=int)
        candidates = np.asarray(
            [i for i in range(len(df)) if i not in selected_labels],
            dtype=int,
        )
        if len(candidates) == 0:
            break

        d = np.linalg.norm(
            Z[candidates, None, :] - Z[selected_idx][None, :, :],
            axis=2,
        )
        min_d = np.min(d, axis=1)
        chosen = int(candidates[int(np.argmax(min_d))])
        add(chosen, "MAXIMIN_PARETO_COVERAGE")

    rows = []
    for idx in sorted(selected_labels):
        row = df.iloc[idx].copy()
        row["validation_selection_reason"] = ";".join(selected_labels[idx])
        rows.append(row)

    out = pd.DataFrame(rows).reset_index(drop=True)
    out.insert(0, "validation_candidate_id", np.arange(1, len(out) + 1))
    return out


# =============================================================================
# 16. FIGURAS
# =============================================================================


def plot_pareto_3d(df: pd.DataFrame, path: Path) -> None:
    if df.empty:
        return

    fig = plt.figure(figsize=(8.0, 6.5))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        df["dose_pred_5min_mSv"],
        df["massa_total_kg"],
        df["custo_total_unidade"],
        s=22,
        alpha=0.75,
    )
    ax.set_xlabel("Dose efetiva prevista em 5 min (mSv)")
    ax.set_ylabel("Massa (kg)")
    ax.set_zlabel("Custo")
    ax.set_title("Frente de Pareto — união multirun")
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_projection(
    df: pd.DataFrame,
    x: str,
    y: str,
    xlabel: str,
    ylabel: str,
    path: Path,
) -> None:
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(7.0, 5.5))
    ax.scatter(df[x], df[y], s=22, alpha=0.75)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_hv_convergence(histories: list[pd.DataFrame], path: Path) -> None:
    if not histories:
        return

    fig, ax = plt.subplots(figsize=(8.0, 5.5))

    for h in histories:
        if not h.empty:
            ax.plot(
                h["generation"],
                h["hv_normalized"],
                alpha=0.65,
            )

    ax.set_xlabel("Geração")
    ax.set_ylabel("Hypervolume normalizado")
    ax.set_title("Convergência multirun do NSGA-II")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# 17. RELATÓRIO FINAL
# =============================================================================


def final_multirun_analysis(
    pareto_runs: list[pd.DataFrame],
    run_summaries: list[dict],
    histories: list[pd.DataFrame],
    ideal_anchor: np.ndarray,
    reference_anchor: np.ndarray,
    hv_doe: float,
    doe_perms: np.ndarray,
    doe_thickness: np.ndarray,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    nonempty = [df for df in pareto_runs if df is not None and not df.empty]

    run_summary_df = pd.DataFrame(run_summaries)
    run_summary_df["hv_gain_vs_doe_normalized"] = (
        run_summary_df["final_hv_normalized"] - hv_doe
    )

    if not nonempty:
        raise RuntimeError(
            "Nenhuma réplica encontrou solução balisticamente viável sob o limite "
            "conformal. Inspecione surrogate_penetration, q_conformal e o limite de 60 mm."
        )

    union = pd.concat(nonempty, ignore_index=True)
    union = deduplicate_designs(union)

    # Todas já são viáveis em suas réplicas; extrair Pareto global da união.
    pareto = nondominated_dataframe(union)

    if DOMAIN_DISTANCE_REPORT:
        pareto = add_domain_diagnostics(
            pareto, doe_perms, doe_thickness
        )

    objective_cols = [
        "dose_pred_5min_mSv",
        "massa_total_kg",
        "custo_total_unidade",
    ]
    F = pareto[objective_cols].to_numpy(dtype=float)
    F_reference_normalized = normalize_objectives(
        F, ideal_anchor, reference_anchor
    )

    # Indicadores complementares por réplica. A união global não dominada é a
    # referência empírica comum, pois a frente Pareto verdadeira é desconhecida.
    quality_rows: list[dict] = []
    for run_id, run_pareto in enumerate(pareto_runs, start=1):
        row = {
            "scope": "run",
            "run_id": run_id,
            "reference_front": "global_pareto_union_empirical",
            "n_solutions": int(len(run_pareto)) if run_pareto is not None else 0,
        }
        if run_pareto is None or run_pareto.empty:
            row.update({
                "igd_normalized": float("nan"),
                "generalized_spread_normalized": float("nan"),
            })
        else:
            F_run = run_pareto[objective_cols].to_numpy(dtype=float)
            F_run_normalized = normalize_objectives(
                F_run, ideal_anchor, reference_anchor
            )
            row.update({
                "igd_normalized": inverted_generational_distance(
                    F_run_normalized, F_reference_normalized
                ),
                "generalized_spread_normalized": generalized_spread(
                    F_run_normalized, F_reference_normalized
                ),
            })
        quality_rows.append(row)

    quality_rows.append({
        "scope": "global_union",
        "run_id": "",
        "reference_front": "global_pareto_union_empirical",
        "n_solutions": int(len(pareto)),
        "igd_normalized": inverted_generational_distance(
            F_reference_normalized, F_reference_normalized
        ),
        "generalized_spread_normalized": generalized_spread(
            F_reference_normalized, F_reference_normalized
        ),
    })
    quality_df = pd.DataFrame(quality_rows)
    quality_df.to_csv(output_dir / "pareto_quality_multirun.csv", index=False)

    run_quality = quality_df[quality_df["scope"] == "run"].copy()
    run_summary_df = run_summary_df.merge(
        run_quality[
            ["run_id", "igd_normalized", "generalized_spread_normalized"]
        ],
        on="run_id",
        how="left",
        validate="one_to_one",
    )
    run_summary_df.to_csv(output_dir / "resumo_runs.csv", index=False)

    best_idx, scores = topsis_iqr(F)
    pareto["topsis_iqr_score"] = scores
    pareto["is_topsis_iqr"] = 0
    pareto.loc[best_idx, "is_topsis_iqr"] = 1

    pareto = pareto.sort_values(
        ["dose_pred_5min_mSv", "massa_total_kg", "custo_total_unidade"]
    ).reset_index(drop=True)
    pareto.insert(0, "pareto_id", np.arange(1, len(pareto) + 1))

    pareto.to_csv(output_dir / "pareto_union_multirun.csv", index=False)

    candidates = select_validation_candidates(
        pareto,
        n_select=N_HIGH_FIDELITY_CANDIDATES,
    )
    candidates.to_csv(
        output_dir / "candidatos_validacao_alta_fidelidade.csv",
        index=False,
    )

    hv_raw_union = hypervolume_3d(F, reference_anchor)
    hv_norm_union = normalized_hypervolume(
        F, ideal_anchor, reference_anchor
    )

    hvs = run_summary_df["final_hv_normalized"].to_numpy(dtype=float)

    summary = {
        "script_version": SCRIPT_VERSION,
        "completed_utc": utc_now(),
        "n_runs_completed": int(len(run_summary_df)),
        "n_global_pareto": int(len(pareto)),
        "hv_doe_surrogate_baseline_normalized": float(hv_doe),
        "hv_union_raw": float(hv_raw_union),
        "hv_union_normalized": float(hv_norm_union),
        "run_hv_normalized_mean": float(np.mean(hvs)),
        "run_hv_normalized_std": float(np.std(hvs, ddof=1)) if len(hvs) > 1 else 0.0,
        "run_hv_normalized_median": float(np.median(hvs)),
        "run_hv_normalized_iqr": float(
            np.percentile(hvs, 75) - np.percentile(hvs, 25)
        ),
        "pareto_quality_reference": (
            "global_pareto_union_empirical_normalized_by_fixed_doe_anchor"
        ),
        "run_igd_normalized_mean": float(
            run_quality["igd_normalized"].mean()
        ),
        "run_igd_normalized_std": float(
            run_quality["igd_normalized"].std(ddof=1)
        ) if len(run_quality) > 1 else 0.0,
        "run_generalized_spread_normalized_mean": float(
            run_quality["generalized_spread_normalized"].mean()
        ),
        "run_generalized_spread_normalized_std": float(
            run_quality["generalized_spread_normalized"].std(ddof=1)
        ) if len(run_quality) > 1 else 0.0,
        "global_union_igd_normalized": float(
            quality_df.loc[
                quality_df["scope"] == "global_union", "igd_normalized"
            ].iloc[0]
        ),
        "global_union_generalized_spread_normalized": float(
            quality_df.loc[
                quality_df["scope"] == "global_union",
                "generalized_spread_normalized",
            ].iloc[0]
        ),
        "ideal_anchor": ideal_anchor,
        "reference_anchor": reference_anchor,
        "ballistic_constraint": (
            "max_penetration_pred_mm + conformal_q_abs <= 60 mm"
        ),
        "topsis_weights": TOPSIS_WEIGHTS / TOPSIS_WEIGHTS.sum(),
        "high_fidelity_candidates": int(len(candidates)),
    }
    dump_json(summary, output_dir / "resumo_final_fase4.json")

    figures = output_dir / "figuras"
    figures.mkdir(exist_ok=True)

    plot_pareto_3d(
        pareto,
        figures / "pareto_3d_union.png",
    )
    plot_projection(
        pareto,
        "dose_pred_5min_mSv",
        "massa_total_kg",
        "Dose efetiva prevista em 5 min (mSv)",
        "Massa (kg)",
        figures / "pareto_dose_massa.png",
    )
    plot_projection(
        pareto,
        "dose_pred_5min_mSv",
        "custo_total_unidade",
        "Dose efetiva prevista em 5 min (mSv)",
        "Custo",
        figures / "pareto_dose_custo.png",
    )
    plot_projection(
        pareto,
        "massa_total_kg",
        "custo_total_unidade",
        "Massa (kg)",
        "Custo",
        figures / "pareto_massa_custo.png",
    )
    plot_hv_convergence(
        histories,
        figures / "convergencia_hypervolume_multirun.png",
    )

    return pareto, candidates


# =============================================================================
# 18. SELF-TESTS SEM MODELOS
# =============================================================================


def run_self_tests() -> None:
    print("Executando self-tests da Fase 4...")

    # Projeção/quantização.
    rng = np.random.default_rng(123)
    for _ in range(100):
        v = rng.normal(0.85, 0.5, size=7)
        x = quantize_thickness(v)
        assert np.isclose(x.sum(), 6.0, atol=1e-12)
        assert np.all(x >= 0.2)
        assert np.all(x <= 2.0)
        assert np.allclose(x * 1000, np.round(x * 1000))

    # OX.
    p1 = np.arange(7, dtype=np.int16)
    p2 = p1[::-1].copy()
    for _ in range(50):
        c1, c2 = order_crossover(p1, p2, rng)
        assert np.array_equal(np.sort(c1), np.arange(7))
        assert np.array_equal(np.sort(c2), np.arange(7))

    # Massa/custo.
    perms = np.vstack([p1, p2])
    thick = np.vstack([
        np.asarray([0.8, 0.8, 0.8, 0.8, 0.8, 1.0, 1.0]),
        np.asarray([1.0, 1.0, 0.8, 0.8, 0.8, 0.8, 0.8]),
    ])
    validate_population(perms, thick)
    mass, cost = deterministic_mass_cost(perms, thick)
    assert np.all(mass > 0)
    assert np.all(cost > 0)

    # Hypervolume de um ponto: (1-a)(1-b)(1-c).
    point = np.asarray([[0.2, 0.3, 0.4]])
    ref = np.asarray([1.0, 1.0, 1.0])
    hv = hypervolume_3d(point, ref)
    assert np.isclose(hv, 0.8 * 0.7 * 0.6, atol=1e-12)

    # Dois pontos: inclusão-exclusão.
    pts = np.asarray([[0.2, 0.5, 0.5], [0.5, 0.2, 0.5]])
    expected = (
        (0.8 * 0.5 * 0.5)
        + (0.5 * 0.8 * 0.5)
        - (0.5 * 0.5 * 0.5)
    )
    assert np.isclose(hypervolume_3d(pts, ref), expected, atol=1e-12)

    # Indicadores de qualidade da frente: referência idêntica implica IGD zero;
    # uma aproximação que contém a própria referência não pode piorar o IGD.
    reference_front = np.asarray([
        [0.1, 0.8, 0.8],
        [0.4, 0.4, 0.6],
        [0.8, 0.2, 0.3],
    ])
    assert np.isclose(
        inverted_generational_distance(reference_front, reference_front),
        0.0,
        atol=1e-12,
    )
    assert np.isclose(
        inverted_generational_distance(
            np.vstack([reference_front, [0.2, 0.2, 0.2]]),
            reference_front,
        ),
        0.0,
        atol=1e-12,
    )
    delta = generalized_spread(reference_front, reference_front)
    assert np.isfinite(delta) and delta >= 0.0

    # Constraint-domination.
    F = np.asarray([
        [1.0, 1.0, 1.0],   # viável
        [0.1, 0.1, 0.1],   # inviável, mesmo "melhor" em objetivos
        [2.0, 2.0, 2.0],   # viável e dominado pelo 0
    ])
    cv = np.asarray([0.0, 1.0, 0.0])
    fronts, ranks = fast_non_dominated_sort(F, cv)
    assert ranks[0] < ranks[2] < ranks[1]

    print("SELF-TESTS: PASS")


# =============================================================================
# 19. MAIN
# =============================================================================


def main() -> None:
    if sklearn.__version__ != TRAINING_SKLEARN_VERSION:
        raise RuntimeError(
            "Os bundles da Fase 3 foram treinados com scikit-learn "
            f"{TRAINING_SKLEARN_VERSION}, mas este interpretador usa "
            f"{sklearn.__version__}. Execute a Fase 4 com "
            r"C:\Users\Admin\anaconda3\envs\spyder_dl8\python.exe."
        )
    n_runs, pop_size, num_generations = selected_run_parameters()

    if n_runs > len(RUN_SEEDS):
        raise ValueError("RUN_SEEDS não contém seeds suficientes.")

    print("=" * 92)
    print("FASE 4 — NSGA-II RESTRITO / DOSE + MASSA + CUSTO / CONSTRAINT BALÍSTICA")
    print("=" * 92)
    print(f"Script: {SCRIPT_VERSION}")
    print(f"RUN_MODE: {RUN_MODE}")
    print(
        f"Orçamento: {n_runs} runs × {pop_size} indivíduos × "
        f"{num_generations} gerações"
    )
    print("Objetivos: minimizar dose 5 min, massa e custo")
    print(
        "Restrição balística: penetração prevista + q conformal 95% <= "
        f"{BALLISTIC_LIMIT_MM:.1f} mm"
    )
    print(
        "Decisão: 7 materiais únicos + 7 espessuras; "
        "0.2–2.0 cm; soma=6.0 cm"
    )

    # Contrato dos modelos.
    dose_bundle = load_surrogate_bundle(
        DOSE_MODEL_PATH, EXPECTED_DOSE_TARGET
    )
    pen_bundle = load_surrogate_bundle(
        PEN_MODEL_PATH, EXPECTED_PEN_TARGET
    )
    evaluator = Evaluator(dose_bundle, pen_bundle)

    print("\nSurrogates:")
    print(
        f"  Dose: {dose_bundle['feature_set']} + {dose_bundle['model_name']} | "
        f"q95={evaluator.q_dose:.6g} mSv"
    )
    print(
        f"  Pen.: {pen_bundle['feature_set']} + {pen_bundle['model_name']} | "
        f"q95={evaluator.q_pen:.6g} mm"
    )

    model_hashes = {
        "surrogate_dose": sha256_file(DOSE_MODEL_PATH),
        "surrogate_penetration": sha256_file(PEN_MODEL_PATH),
    }

    feature_path = PHASE3_FEATURE_DIR / "feature_engineering_final.py"
    if feature_path.exists():
        model_hashes["feature_engineering"] = sha256_file(feature_path)

    doe_path = resolve_first_existing(DOE_CANDIDATES, "CSV do DoE")
    doe = load_doe(doe_path)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = RESULTS_ROOT / f"fase4_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # DoE como baseline/âncora fixa.
    (
        ideal_anchor,
        reference_anchor,
        doe_eval_df,
        doe_perms,
        doe_thickness,
        hv_doe,
    ) = evaluate_doe_baseline(doe, evaluator, output_dir)

    print("\nÂncora de hypervolume derivada do DoE:")
    print(f"  Ideal:     {ideal_anchor}")
    print(f"  Referência:{reference_anchor}")
    print(f"  HVn DoE baseline: {hv_doe:.6f}")

    manifest = {
        "script_version": SCRIPT_VERSION,
        "created_utc": utc_now(),
        "run_mode": RUN_MODE,
        "n_runs": n_runs,
        "population_size": pop_size,
        "num_generations": num_generations,
        "run_seeds": RUN_SEEDS[:n_runs],
        "operators": {
            "permutation_crossover": "Order Crossover (OX)",
            "continuous_crossover": "bounded SBX + box-simplex projection",
            "permutation_mutation": "swap mutation",
            "continuous_mutation": "polynomial mutation + box-simplex projection",
            "crossover_probability": CROSSOVER_PROB,
            "perm_swap_mutation_probability": PERM_SWAP_MUTATION_PROB,
            "thickness_gene_mutation_probability": THICKNESS_GENE_MUTATION_PROB,
            "sbx_eta": SBX_ETA,
            "polynomial_mutation_eta": POLYNOMIAL_MUTATION_ETA,
            "adaptive_mutation": False,
            "adaptive_restart": False,
        },
        "design_constraints": {
            "materials": list(MATERIALS),
            "each_material_exactly_once": True,
            "min_thickness_cm": MIN_THICKNESS_CM,
            "max_thickness_cm": MAX_THICKNESS_CM,
            "total_thickness_cm": TOTAL_THICKNESS_CM,
            "resolution_cm": THICKNESS_RESOLUTION_CM,
        },
        "ballistic_constraint": {
            "limit_mm": BALLISTIC_LIMIT_MM,
            "definition": (
                "prediction + split-conformal absolute quantile <= limit"
            ),
            "q_conformal_mm": evaluator.q_pen,
        },
        "objectives": [
            "dose_pred_5min_mSv",
            "massa_total_kg (deterministic)",
            "custo_total_unidade (deterministic)",
        ],
        "models": {
            "dose": {
                "path": str(DOSE_MODEL_PATH),
                "sha256": model_hashes["surrogate_dose"],
                "feature_set": dose_bundle["feature_set"],
                "model_name": dose_bundle["model_name"],
                "test_metrics": dose_bundle["test_metrics"],
                "conformal_q_abs": evaluator.q_dose,
            },
            "penetration": {
                "path": str(PEN_MODEL_PATH),
                "sha256": model_hashes["surrogate_penetration"],
                "feature_set": pen_bundle["feature_set"],
                "model_name": pen_bundle["model_name"],
                "test_metrics": pen_bundle["test_metrics"],
                "conformal_q_abs": evaluator.q_pen,
            },
        },
        "doe": {
            "path": str(doe_path),
            "sha256": sha256_file(doe_path),
            "n": len(doe),
            "hv_normalized_baseline": hv_doe,
        },
        "hv_anchor": {
            "ideal": ideal_anchor,
            "reference": reference_anchor,
            "reference_margin": REFERENCE_MARGIN,
        },
        "material_properties": MATERIAL_PROPERTIES,
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "joblib": getattr(joblib, "__version__", "unknown"),
        },
    }
    dump_json(manifest, output_dir / "manifesto_fase4.json")

    pareto_runs: list[pd.DataFrame] = []
    run_summaries: list[dict] = []
    histories: list[pd.DataFrame] = []

    for run_id, seed in enumerate(RUN_SEEDS[:n_runs], start=1):
        print("\n" + "-" * 92)
        print(f"INICIANDO RUN {run_id}/{n_runs} | seed={seed}")
        print("-" * 92)

        pareto_df, summary, history_df = run_one_seed(
            run_id=run_id,
            seed=seed,
            pop_size=pop_size,
            num_generations=num_generations,
            evaluator=evaluator,
            ideal_anchor=ideal_anchor,
            reference_anchor=reference_anchor,
            model_hashes=model_hashes,
            output_dir=output_dir,
        )

        pareto_runs.append(pareto_df)
        run_summaries.append(summary)
        history_df = history_df.copy()
        history_df.insert(0, "run_id", run_id)
        history_df.insert(1, "seed", seed)
        histories.append(history_df)

    all_hist = pd.concat(histories, ignore_index=True)
    all_hist.to_csv(output_dir / "convergencia_todos_runs.csv", index=False)

    pareto, candidates = final_multirun_analysis(
        pareto_runs=pareto_runs,
        run_summaries=run_summaries,
        histories=histories,
        ideal_anchor=ideal_anchor,
        reference_anchor=reference_anchor,
        hv_doe=hv_doe,
        doe_perms=doe_perms,
        doe_thickness=doe_thickness,
        output_dir=output_dir,
    )

    print("\n" + "=" * 92)
    print("FASE 4 CONCLUÍDA")
    print("=" * 92)
    print(f"Resultados: {output_dir}")
    print(f"Pareto global: {len(pareto)} soluções")
    print(f"Candidatos para validação de alta fidelidade: {len(candidates)}")
    print("Arquivos principais:")
    print("  - pareto_union_multirun.csv")
    print("  - candidatos_validacao_alta_fidelidade.csv")
    print("  - resumo_runs.csv")
    print("  - pareto_quality_multirun.csv")
    print("  - convergencia_todos_runs.csv")
    print("  - resumo_final_fase4.json")
    print("  - manifesto_fase4.json")
    print("\nNenhuma simulação MCNP/LS-DYNA foi executada pela Fase 4.")


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        run_self_tests()
    else:
        if "--smoke-test" in sys.argv:
            RUN_MODE = "SMOKE_TEST"
            RESULTS_ROOT = BASE_DIR / "resultados_otimizacao (Fase 4) - MCNP63_DEDF_AP - smoke"
            MODEL_DIR = PHASE3_DIR / "modelos_treinados (Fase 3) - MCNP63_DEDF_AP"
            DOSE_MODEL_PATH = MODEL_DIR / "surrogate_dose.joblib"
            PEN_MODEL_PATH = MODEL_DIR / "surrogate_penetration.joblib"
        main()
