"""Publication-table generator for the final Phase-3 outputs.

The tables distinguish training-only model selection from the untouched final
test.  They are CSV sources for the EAAI manuscript, not typeset tables.
"""
from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd


ALGORITHM_CLASS = {
    "Ridge": "Linear regularized regression",
    "RandomForest": "Bagged decision-tree ensemble",
    "ExtraTrees": "Extremely randomized tree ensemble",
    "GradientBoosting": "Sequential gradient-boosted trees",
    "HistGradientBoosting": "Histogram-based gradient-boosted trees",
    "MLP": "Feed-forward neural network",
    "XGBoost": "Regularized gradient-boosted trees",
}


def _format_space(space: dict) -> str:
    if not space:
        return "No tuned hyperparameters"
    return "; ".join(f"{key.replace('model__', '')}: {list(values)}" for key, values in space.items())


def _feature_rows() -> list[dict]:
    return [
        {"feature_group": "Raw design variables", "feature": "Layer material", "formula": "m_i, i=1,...,7", "physical_rationale": "Discrete material sequence of the seven-layer shield."},
        {"feature_group": "Raw design variables", "feature": "Layer thickness", "formula": "t_i, i=1,...,7; sum(t_i)=6 cm", "physical_rationale": "Allocated radial thickness at each ordered layer."},
        {"feature_group": "Radiological physics", "feature": "Layer attenuation descriptor", "formula": "(mu/rho)_i rho_i t_i", "physical_rationale": "Approximate layerwise photon attenuation at the Co-60 energy scale."},
        {"feature_group": "Radiological physics", "feature": "Total attenuation descriptor", "formula": "sum_i mu_i t_i", "physical_rationale": "Aggregate attenuation along the radial transport path."},
        {"feature_group": "Radiological physics", "feature": "Transmission proxy", "formula": "exp[-sum_i(mu_i t_i)]", "physical_rationale": "Beer-Lambert reference descriptor, not a replacement for MCNP transport."},
        {"feature_group": "Radiological physics", "feature": "Effective atomic-number descriptors", "formula": "sum_i Z_eff,i t_i / sum_i t_i", "physical_rationale": "Captures material-composition variation relevant to photon interaction."},
        {"feature_group": "Geometry and mass", "feature": "Areal-density descriptors", "formula": "rho_i t_i; sum_i(rho_i t_i)", "physical_rationale": "Layerwise and total mass per unit area, retaining radial order."},
        {"feature_group": "Geometry and mass", "feature": "Face and radial-moment descriptors", "formula": "f(m_1,t_1,m_7,t_7); sum_i[w_i rho_i t_i]", "physical_rationale": "Distinguishes operator-facing and impact/source-facing material placement."},
        {"feature_group": "Geometry and mass", "feature": "Thickness-distribution descriptors", "formula": "std(t), min(t), max(t), H(t)", "physical_rationale": "Quantifies concentration and uniformity of the thickness allocation."},
        {"feature_group": "Frozen ballistic mechanics", "feature": "Layer areal strength", "formula": "S_i t_i", "physical_rationale": "Strength-thickness capacity using the frozen LS-DYNA strength reference."},
        {"feature_group": "Frozen ballistic mechanics", "feature": "Total and radial-centroid areal strength", "formula": "sum_i(S_i t_i); sum_i[w_i S_i t_i]/sum_i(S_i t_i)", "physical_rationale": "Summarizes resistance magnitude and its location along the impact path."},
        {"feature_group": "Frozen ballistic mechanics", "feature": "Strength discontinuity at interfaces", "formula": "abs(S_(i+1)-S_i)", "physical_rationale": "Represents ordered material mismatch at each of the six interfaces."},
        {"feature_group": "Frozen ballistic mechanics", "feature": "Serial shear compliance", "formula": "sum_i(t_i/G_i)", "physical_rationale": "Proxy for through-thickness shear deformability from the frozen constitutive cards."},
    ]


def generate_article_tables(output_dir: Path, tuning_space: callable, final_metadata: dict) -> Path:
    """Write four article-ready CSV tables and one full supplementary matrix."""
    table_dir = output_dir / "tabelas_artigo"
    table_dir.mkdir(parents=True, exist_ok=True)
    benchmark = pd.read_csv(output_dir / "benchmark_cv_summary.csv")

    algorithms = sorted(benchmark["model"].dropna().unique())
    pd.DataFrame([
        {
            "algorithm_class": ALGORITHM_CLASS.get(model, "Regression model"),
            "algorithm": model,
            "optimized_hyperparameters_search_space": _format_space(tuning_space(model)),
            "optimization_protocol": "RandomizedSearchCV, 5-fold CV, training partition only",
        }
        for model in algorithms
    ]).to_csv(table_dir / "table_1_algorithms_and_optimization.csv", index=False)

    pd.DataFrame(_feature_rows()).to_csv(
        table_dir / "table_2_feature_groups_and_physical_rationale.csv", index=False
    )

    rows = []
    for (target, model), group in benchmark.groupby(["target", "model"], sort=True):
        values = group.set_index("feature_set")
        raw = values.loc["RAW"] if "RAW" in values.index else None
        physics = values.loc["PHYSICS"] if "PHYSICS" in values.index else None
        strength = values.loc["PHYSICS_STRENGTH23"] if "PHYSICS_STRENGTH23" in values.index else None
        candidates = [x for x in (raw, physics, strength) if x is not None]
        best = max(candidates, key=lambda row: float(row["cv_r2_mean"]))
        rows.append({
            "target": target,
            "model": model,
            "raw_features_cv_r2": float(raw["cv_r2_mean"]) if raw is not None else np.nan,
            "physics_features_cv_r2": float(physics["cv_r2_mean"]) if physics is not None else np.nan,
            "physics_strength23_cv_r2": float(strength["cv_r2_mean"]) if strength is not None else np.nan,
            "best_representation": best.name,
            "best_cv_r2": float(best["cv_r2_mean"]),
            "delta_r2_best_minus_raw": float(best["cv_r2_mean"] - raw["cv_r2_mean"]) if raw is not None else np.nan,
            "selection_scope": "5-fold CV on training partition only",
        })
    pd.DataFrame(rows).to_csv(table_dir / "table_3_feature_representation_comparison.csv", index=False)

    final_rows = []
    for target, metadata in final_metadata.items():
        metrics = metadata["test_metrics"]
        final_rows.append({
            "target": target,
            "selected_model": metadata["model_name"],
            "selected_feature_set": metadata["feature_set"],
            "test_r2": metrics["r2"],
            "test_mae": metrics["mae"],
            "test_rmse": metrics["rmse"],
            "test_smape_pct": metrics["smape_pct"],
            "test_inference_seconds_per_case": metrics["test_inference_seconds_per_case"],
            "conformal_q_abs": metrics["conformal_q_abs"],
            "test_coverage_95": metrics["test_empirical_coverage"],
            "evaluation_scope": "Untouched held-out test partition",
        })
    pd.DataFrame(final_rows).to_csv(table_dir / "table_4_selected_models_untouched_test.csv", index=False)

    benchmark.to_csv(table_dir / "supplementary_all_optimized_candidates.csv", index=False)
    manifest = {
        "table_1": "Algorithm classes and RandomizedSearchCV search spaces.",
        "table_2": "Feature groups, formulas, and physical rationale.",
        "table_3": "Training-only feature-representation comparison.",
        "table_4": "Selected-surrogate results on the untouched test partition.",
        "supplementary": "Full matrix of optimized candidate combinations.",
    }
    (table_dir / "README_tables_article.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return table_dir
