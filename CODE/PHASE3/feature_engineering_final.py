# -*- coding: utf-8 -*-
"""
Feature engineering compartilhado — Fase 3 (surrogates) e Fase 4 (otimização)
===============================================================================
Projeto: blindagem semicilíndrica multicamada para RDD Co-60 + fragmento balístico.

Princípios metodológicos
------------------------
1. Nenhuma feature usa dose, penetração, perfuração, energia residual ou qualquer
   outra saída de MCNP/LS-DYNA.
2. Nenhuma coluna ausente é preenchida artificialmente. Dados inválidos abortam.
3. Cada configuração deve conter os sete materiais exatamente uma vez.
4. As sete espessuras devem estar entre 0,2 e 2,0 cm e somar 6,0 cm.
5. Camada 1 = face interna / operador; camada 7 = face externa / fonte e impacto.
6. Massa e custo não são targets de ML; permanecem cálculos determinísticos.

Conjuntos de entrada suportados
-------------------------------
RAW:
    7 materiais categóricos + 7 espessuras.
PHYSICS:
    RAW + features físicas determinísticas e auditáveis.
PHYSICS_STRENGTH23:
    PHYSICS + 23 descritores mecânicos pré-especificados, extraídos somente
    das propriedades constitutivas do modelo LS-DYNA congelado.

O módulo foi deliberadamente reduzido em relação à versão anterior: foram
removidas features constantes, propriedades balísticas arbitrárias, imputações
silenciosas, hashes/label-encodings artificiais e features derivadas de targets.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin


# =============================================================================
# 1. DOMÍNIO CONGELADO
# =============================================================================

N_LAYERS = 7
TOTAL_THICKNESS_CM = 6.0
MIN_LAYER_CM = 0.2
MAX_LAYER_CM = 2.0
INNER_RADIUS_CM = 35.0
OUTER_RADIUS_CM = 41.0

MATERIALS = (
    "SiC",
    "UHMWPE",
    "Al",
    "Fe",
    "B4C-HDPE",
    "W-Epoxy",
    "Bi2O3-HDPE",
)
MATERIAL_SET = set(MATERIALS)

MATERIAL_COLS = [f"material_camada_{i}" for i in range(1, N_LAYERS + 1)]
THICKNESS_COLS = [f"espessura_camada_{i}" for i in range(1, N_LAYERS + 1)]

# Densidades canônicas da tese (g/cm³).
DENSITY_G_CM3 = {
    "SiC": 3.210,
    "UHMWPE": 0.960,
    "Al": 2.700,
    "Fe": 7.870,
    "B4C-HDPE": 1.428,
    "W-Epoxy": 4.643,
    "Bi2O3-HDPE": 2.557,
}

# Coeficientes de atenuação linear aproximados para ~1,25 MeV (cm^-1),
# adotados no Artigo 1 e mantidos aqui como descritores physics-informed.
# São FEATURES auxiliares; o target continua sendo a resposta MCNP F4+FM4.
MU_CO60_CM_INV = {
    "SiC": 0.1826,
    "UHMWPE": 0.0624,
    "Al": 0.1484,
    "Fe": 0.4210,
    "B4C-HDPE": 0.0843,
    "W-Epoxy": 0.2635,
    "Bi2O3-HDPE": 0.1561,
}

# Z_eff adotado no Artigo 1 como descritor físico auxiliar.
Z_EFF = {
    "SiC": 10.8,
    "UHMWPE": 5.9,
    "Al": 13.0,
    "Fe": 26.0,
    "B4C-HDPE": 5.2,
    "W-Epoxy": 54.1,
    "Bi2O3-HDPE": 46.7,
}

# Classificação simples usada somente para descritores agregados.
HIGH_Z_MATERIALS = {m for m, z in Z_EFF.items() if z >= 26.0}
LOW_Z_MATERIALS = MATERIAL_SET - HIGH_Z_MATERIALS

# Propriedades transcritas do modelo mecânico LS-DYNA congelado.  G é o módulo
# de cisalhamento disponível no cartão constitutivo. strength_ref representa o
# parâmetro de início/limiar disponível: HEL para SiC, Johnson--Cook A para Al
# e Fe, e SIGY para os compósitos MAT_003.  Esses valores são descritores de
# entrada, não resultados de simulação e não tornam mecanismos de falha de
# famílias constitutivas distintas diretamente equivalentes.
FROZEN_SHEAR_MODULUS_GPA = {
    "SiC": 193.0,
    "UHMWPE": 21.33 / (2.0 * 1.27),
    "Al": 70.0 / (2.0 * 1.33),
    "Fe": 205.0 / (2.0 * 1.29),
    "B4C-HDPE": 1.1493 / (2.0 * 1.35),
    "W-Epoxy": 10.0 / (2.0 * 1.30),
    "Bi2O3-HDPE": 1.50 / (2.0 * 1.38),
}
FROZEN_STRENGTH_REFERENCE_GPA = {
    "SiC": 11.700,
    "UHMWPE": 0.320,
    "Al": 0.324,
    "Fe": 0.792,
    "B4C-HDPE": 0.0239,
    "W-Epoxy": 0.0942,
    "Bi2O3-HDPE": 0.0200,
}


# =============================================================================
# 2. VALIDAÇÃO
# =============================================================================


def _as_dataframe(X) -> pd.DataFrame:
    if isinstance(X, pd.DataFrame):
        return X.copy()
    raise TypeError(
        "ShieldFeatureEngineer exige pandas.DataFrame para preservar nomes de colunas."
    )


def validate_design_dataframe(
    df: pd.DataFrame,
    material_cols: Sequence[str] = MATERIAL_COLS,
    thickness_cols: Sequence[str] = THICKNESS_COLS,
    atol_sum: float = 1e-6,
) -> None:
    """Valida rigorosamente o domínio físico/combinatório do DoE."""
    required = list(material_cols) + list(thickness_cols)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Colunas obrigatórias ausentes: {missing}")

    if df[required].isna().any().any():
        bad = df[required].isna().sum()
        bad = bad[bad > 0].to_dict()
        raise ValueError(f"Valores ausentes nas features de entrada: {bad}")

    # Espessuras numéricas e dentro do domínio.
    t = df[list(thickness_cols)].apply(pd.to_numeric, errors="coerce")
    if t.isna().any().any():
        raise ValueError("Há espessuras não numéricas.")

    if ((t < MIN_LAYER_CM - 1e-12) | (t > MAX_LAYER_CM + 1e-12)).any().any():
        bad_rows = np.where(
            ((t < MIN_LAYER_CM - 1e-12) | (t > MAX_LAYER_CM + 1e-12)).any(axis=1)
        )[0][:10]
        raise ValueError(
            f"Espessura fora de [{MIN_LAYER_CM}, {MAX_LAYER_CM}] cm. "
            f"Primeiras linhas problemáticas: {bad_rows.tolist()}"
        )

    sums = t.sum(axis=1).to_numpy(dtype=float)
    if not np.allclose(sums, TOTAL_THICKNESS_CM, atol=atol_sum, rtol=0.0):
        bad = np.where(~np.isclose(sums, TOTAL_THICKNESS_CM, atol=atol_sum, rtol=0.0))[0][:10]
        raise ValueError(
            f"Soma das espessuras diferente de {TOTAL_THICKNESS_CM} cm. "
            f"Primeiras linhas problemáticas: {bad.tolist()}"
        )

    # Cada projeto deve ser uma permutação dos sete materiais.
    for pos, row in enumerate(df[list(material_cols)].itertuples(index=False, name=None)):
        row_set = set(row)
        if len(row) != N_LAYERS or len(row_set) != N_LAYERS or row_set != MATERIAL_SET:
            raise ValueError(
                f"Linha {pos}: materiais não formam uma permutação dos sete materiais. "
                f"Recebido: {list(row)}"
            )


# =============================================================================
# 3. NOMES DAS FEATURES FÍSICAS
# =============================================================================


def physics_feature_names() -> list[str]:
    names: list[str] = []

    # Propriedades posicionais por camada.
    for i in range(1, N_LAYERS + 1):
        names.extend(
            [
                f"rho_layer_{i}",
                f"mu_layer_{i}",
                f"z_eff_layer_{i}",
                f"areal_density_layer_{i}",
                f"mu_x_layer_{i}",
            ]
        )

    # Espessura associada a cada material (ordem-invariante quanto à quantidade,
    # mas preservamos também as features posicionais acima para aprender a ordem).
    names.extend([f"thick_{m}" for m in MATERIALS])

    # Física radiológica agregada.
    names.extend(
        [
            "mu_x_total",
            "transmission_beer_lambert",
            "z_eff_weighted_thickness",
            "z_rho_t_total",
            "high_z_weighted_thickness",
            "low_z_weighted_thickness",
            "high_low_z_ratio",
            "z_rho_t_x_mu_x",
            "mu_x_radial_moment",
            "mu_x_radial_centroid",
            "operator_face_mu",
            "operator_face_z_eff",
            "source_face_mu",
            "source_face_z_eff",
        ]
    )

    # Descritores balísticos/geométricos seguros.
    names.extend(
        [
            "areal_density_total",
            "areal_density_radial_moment",
            "areal_density_radial_centroid",
            "operator_face_density",
            "operator_face_thickness_cm",
            "impact_face_density",
            "impact_face_thickness_cm",
            "thickness_std",
            "thickness_min",
            "thickness_max",
            "thickness_range",
            "thickness_median",
            "thickness_max_min_ratio",
            "thickness_entropy",
        ]
    )

    for i in range(1, N_LAYERS):
        names.extend(
            [
                f"delta_thickness_{i}_{i+1}",
                f"delta_density_{i}_{i+1}",
                f"abs_delta_density_{i}_{i+1}",
            ]
        )

    return names


def ballistic_strength23_feature_names() -> list[str]:
    """Names of the pre-specified compact ballistic descriptor family.

    The 23 variables were selected in a train-only, equal-budget comparison.
    They summarize layerwise strength-thickness capacity, material-strength
    discontinuities at interfaces, positional strength references, and serial
    shear compliance.  No response variable is used.
    """
    names = [f"mech_areal_strength_{i}" for i in range(1, N_LAYERS + 1)]
    names.extend(["mech_areal_strength_total", "mech_areal_strength_centroid"])
    names.extend([f"mech_abs_strength_jump_{i}_{i + 1}" for i in range(1, N_LAYERS)])
    names.extend([f"mech_strength_ref_{i}" for i in range(1, N_LAYERS + 1)])
    names.append("mech_shear_compliance_series")
    if len(names) != 23:
        raise RuntimeError("Strength23 feature contract must contain exactly 23 variables.")
    return names


# =============================================================================
# 4. CÁLCULO DAS FEATURES
# =============================================================================


def _normalized_entropy(values: np.ndarray) -> np.ndarray:
    """Entropia de Shannon normalizada por log(7), linha a linha."""
    p = values / values.sum(axis=1, keepdims=True)
    eps = np.finfo(float).tiny
    h = -np.sum(p * np.log(np.clip(p, eps, None)), axis=1)
    return h / np.log(values.shape[1])


def build_physics_features(
    df: pd.DataFrame,
    material_cols: Sequence[str] = MATERIAL_COLS,
    thickness_cols: Sequence[str] = THICKNESS_COLS,
) -> pd.DataFrame:
    """Retorna somente as features físicas numéricas, com índice preservado."""
    validate_design_dataframe(df, material_cols, thickness_cols)

    mats = df[list(material_cols)].astype(str)
    th = df[list(thickness_cols)].astype(float).to_numpy()
    n = len(df)

    rho = np.zeros((n, N_LAYERS), dtype=float)
    mu = np.zeros((n, N_LAYERS), dtype=float)
    z = np.zeros((n, N_LAYERS), dtype=float)

    for i, c in enumerate(material_cols):
        rho[:, i] = mats[c].map(DENSITY_G_CM3).to_numpy(dtype=float)
        mu[:, i] = mats[c].map(MU_CO60_CM_INV).to_numpy(dtype=float)
        z[:, i] = mats[c].map(Z_EFF).to_numpy(dtype=float)

    if not (np.isfinite(rho).all() and np.isfinite(mu).all() and np.isfinite(z).all()):
        raise ValueError("Falha ao mapear propriedades físicas de algum material.")

    areal = rho * th                    # g/cm²
    mu_x = mu * th                      # adimensional

    out = pd.DataFrame(index=df.index)

    for i in range(N_LAYERS):
        k = i + 1
        out[f"rho_layer_{k}"] = rho[:, i]
        out[f"mu_layer_{k}"] = mu[:, i]
        out[f"z_eff_layer_{k}"] = z[:, i]
        out[f"areal_density_layer_{k}"] = areal[:, i]
        out[f"mu_x_layer_{k}"] = mu_x[:, i]

    # Espessura específica de cada material.
    for material in MATERIALS:
        vals = np.zeros(n, dtype=float)
        for i, c in enumerate(material_cols):
            mask = mats[c].to_numpy() == material
            vals[mask] = th[mask, i]
        out[f"thick_{material}"] = vals

    # Coordenada radial normalizada do centro de cada camada: 0 = operador,
    # 1 = face externa/fonte/impacto.
    cumulative_before = np.cumsum(th, axis=1) - th
    center_from_inner = cumulative_before + 0.5 * th
    radial_pos = center_from_inner / TOTAL_THICKNESS_CM

    mu_x_total = mu_x.sum(axis=1)
    z_rho_t_total = np.sum(z * areal, axis=1)
    z_weighted = np.sum(z * th, axis=1) / TOTAL_THICKNESS_CM

    high_mask = np.isin(mats.to_numpy(), list(HIGH_Z_MATERIALS))
    low_mask = ~high_mask
    high_z = np.sum(z * th * high_mask, axis=1)
    low_z = np.sum(z * th * low_mask, axis=1)

    out["mu_x_total"] = mu_x_total
    out["transmission_beer_lambert"] = np.exp(-mu_x_total)
    out["z_eff_weighted_thickness"] = z_weighted
    out["z_rho_t_total"] = z_rho_t_total
    out["high_z_weighted_thickness"] = high_z
    out["low_z_weighted_thickness"] = low_z
    out["high_low_z_ratio"] = high_z / np.maximum(low_z, 1e-12)
    out["z_rho_t_x_mu_x"] = z_rho_t_total * mu_x_total

    mu_moment = np.sum(mu_x * radial_pos, axis=1)
    out["mu_x_radial_moment"] = mu_moment
    out["mu_x_radial_centroid"] = mu_moment / np.maximum(mu_x_total, 1e-12)

    # Camada 1 = operador; camada 7 = fonte/impacto.
    out["operator_face_mu"] = mu[:, 0]
    out["operator_face_z_eff"] = z[:, 0]
    out["source_face_mu"] = mu[:, -1]
    out["source_face_z_eff"] = z[:, -1]

    areal_total = areal.sum(axis=1)
    areal_moment = np.sum(areal * radial_pos, axis=1)
    out["areal_density_total"] = areal_total
    out["areal_density_radial_moment"] = areal_moment
    out["areal_density_radial_centroid"] = areal_moment / np.maximum(areal_total, 1e-12)

    out["operator_face_density"] = rho[:, 0]
    out["operator_face_thickness_cm"] = th[:, 0]
    out["impact_face_density"] = rho[:, -1]
    out["impact_face_thickness_cm"] = th[:, -1]

    out["thickness_std"] = np.std(th, axis=1, ddof=0)
    out["thickness_min"] = np.min(th, axis=1)
    out["thickness_max"] = np.max(th, axis=1)
    out["thickness_range"] = out["thickness_max"] - out["thickness_min"]
    out["thickness_median"] = np.median(th, axis=1)
    out["thickness_max_min_ratio"] = out["thickness_max"] / np.maximum(out["thickness_min"], 1e-12)
    out["thickness_entropy"] = _normalized_entropy(th)

    for i in range(N_LAYERS - 1):
        k = i + 1
        j = i + 2
        out[f"delta_thickness_{k}_{j}"] = th[:, i + 1] - th[:, i]
        out[f"delta_density_{k}_{j}"] = rho[:, i + 1] - rho[:, i]
        out[f"abs_delta_density_{k}_{j}"] = np.abs(rho[:, i + 1] - rho[:, i])

    expected = physics_feature_names()
    missing = [c for c in expected if c not in out.columns]
    extra = [c for c in out.columns if c not in expected]
    if missing or extra:
        raise RuntimeError(f"Contrato de features inconsistente. missing={missing}, extra={extra}")

    if out.isna().any().any() or not np.isfinite(out.to_numpy(dtype=float)).all():
        raise ValueError("Feature engineering produziu NaN/Inf.")

    return out[expected]


def build_ballistic_strength23_features(
    df: pd.DataFrame,
    material_cols: Sequence[str] = MATERIAL_COLS,
    thickness_cols: Sequence[str] = THICKNESS_COLS,
) -> pd.DataFrame:
    """Build the compact frozen-model ballistic descriptor family.

    Layer order is preserved.  The impact-facing layer is seven and the
    operator-facing layer is one.  The calculation is deterministic and is
    intentionally independent of LS-DYNA output files and all targets.
    """
    validate_design_dataframe(df, material_cols, thickness_cols)
    mats = df[list(material_cols)].astype(str)
    th = df[list(thickness_cols)].to_numpy(dtype=float)
    n = len(df)
    shear = np.zeros((n, N_LAYERS), dtype=float)
    strength = np.zeros((n, N_LAYERS), dtype=float)
    for i, col in enumerate(material_cols):
        shear[:, i] = mats[col].map(FROZEN_SHEAR_MODULUS_GPA).to_numpy(dtype=float)
        strength[:, i] = mats[col].map(FROZEN_STRENGTH_REFERENCE_GPA).to_numpy(dtype=float)
    if not (np.isfinite(shear).all() and np.isfinite(strength).all()):
        raise ValueError("Falha ao mapear propriedades mecânicas congeladas.")

    areal_strength = strength * th
    radial = (np.cumsum(th, axis=1) - 0.5 * th) / TOTAL_THICKNESS_CM
    out = pd.DataFrame(index=df.index)
    for i in range(N_LAYERS):
        out[f"mech_areal_strength_{i + 1}"] = areal_strength[:, i]
    total = areal_strength.sum(axis=1)
    out["mech_areal_strength_total"] = total
    out["mech_areal_strength_centroid"] = (areal_strength * radial).sum(axis=1) / np.maximum(total, 1e-12)
    for i in range(N_LAYERS - 1):
        out[f"mech_abs_strength_jump_{i + 1}_{i + 2}"] = np.abs(strength[:, i + 1] - strength[:, i])
    for i in range(N_LAYERS):
        out[f"mech_strength_ref_{i + 1}"] = strength[:, i]
    out["mech_shear_compliance_series"] = (th / shear).sum(axis=1)

    expected = ballistic_strength23_feature_names()
    if list(out.columns) != expected:
        raise RuntimeError("Strength23 feature names are inconsistent with its contract.")
    if out.isna().any().any() or not np.isfinite(out.to_numpy(dtype=float)).all():
        raise ValueError("Strength23 feature engineering produziu NaN/Inf.")
    return out


# =============================================================================
# 5. TRANSFORMER SCIKIT-LEARN / JOBLIB
# =============================================================================


class ShieldFeatureEngineer(BaseEstimator, TransformerMixin):
    """
    Transformer determinístico e sem y.

    feature_set='RAW'                -> 14 entradas originais.
    feature_set='PHYSICS'            -> RAW + descritores físicos.
    feature_set='PHYSICS_STRENGTH23' -> PHYSICS + 23 descritores mecânicos
                                        do modelo LS-DYNA congelado.
    """

    def __init__(
        self,
        feature_set: str = "PHYSICS",
        material_cols: Sequence[str] | None = None,
        thickness_cols: Sequence[str] | None = None,
    ):
        # Não copie nem normalize estes parâmetros aqui. O clone() do scikit-learn
        # verifica que o construtor preserva a identidade dos objetos recebidos.
        self.feature_set = feature_set
        self.material_cols = material_cols if material_cols is not None else MATERIAL_COLS
        self.thickness_cols = thickness_cols if thickness_cols is not None else THICKNESS_COLS

    def fit(self, X, y=None):
        df = _as_dataframe(X)
        validate_design_dataframe(df, self.material_cols, self.thickness_cols)
        feature_set = str(self.feature_set).upper()
        if feature_set not in {"RAW", "PHYSICS", "PHYSICS_STRENGTH23"}:
            raise ValueError("feature_set deve ser 'RAW', 'PHYSICS' ou 'PHYSICS_STRENGTH23'.")
        self.feature_set_ = feature_set
        self.n_features_in_ = len(self.material_cols) + len(self.thickness_cols)
        self.feature_names_in_ = np.array(self.material_cols + self.thickness_cols, dtype=object)
        return self

    def transform(self, X):
        if not hasattr(self, "feature_set_"):
            raise RuntimeError("Transformer não ajustado. Execute fit() primeiro.")

        df = _as_dataframe(X)
        validate_design_dataframe(df, self.material_cols, self.thickness_cols)

        base = df[self.material_cols + self.thickness_cols].copy()
        for c in self.thickness_cols:
            base[c] = pd.to_numeric(base[c], errors="raise").astype(float)

        if self.feature_set_ == "RAW":
            return base

        engineered = build_physics_features(df, self.material_cols, self.thickness_cols)
        parts = [base, engineered]
        if self.feature_set_ == "PHYSICS_STRENGTH23":
            parts.append(build_ballistic_strength23_features(df, self.material_cols, self.thickness_cols))
        result = pd.concat(parts, axis=1)
        if result.columns.duplicated().any():
            raise RuntimeError("Features duplicadas após engenharia física.")
        return result

    def get_feature_names_out(self, input_features=None):
        base = list(self.material_cols) + list(self.thickness_cols)
        feature_set = getattr(self, "feature_set_", str(self.feature_set).upper())
        if feature_set == "RAW":
            return np.array(base, dtype=object)
        names = base + physics_feature_names()
        if feature_set == "PHYSICS_STRENGTH23":
            names += ballistic_strength23_feature_names()
        return np.array(names, dtype=object)


# Alias explícito para compatibilidade de importação entre Fase 3 e Fase 4.
PhysicsFeatureEngineer = ShieldFeatureEngineer


if __name__ == "__main__":
    print("feature_engineering.py carregado.")
    print(f"Materiais: {list(MATERIALS)}")
    print(f"Features RAW: {len(MATERIAL_COLS) + len(THICKNESS_COLS)}")
    print(f"Features físicas adicionais: {len(physics_feature_names())}")
