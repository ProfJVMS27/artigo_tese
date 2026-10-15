"""Controlled supplementary optimizer study. Never launches MCNP/LS-DYNA.

Frozen surrogate comparison, not physical validation or optimizer tuning.
Run --help. Local checkpoint pickle files must never come from untrusted sources.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import pickle
import platform
import sys
import time

import numpy as np
import pandas as pd
import scipy
import xgboost
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon
import pymoo
from pymoo.algorithms.moo.nsga3 import ReferenceDirectionSurvival
from pymoo.core.population import Population
from pymoo.core.problem import Problem
from pymoo.util.ref_dirs import get_reference_directions
import Fase_4_Otimizacao_Tese_Final as f4

BASE = Path(__file__).resolve().parent
ROOT = BASE / 'comparacoes_otimizacao_20261012'
VERSION = '1.0'
OBJ = ['dose_pred_5min_mSv', 'massa_total_kg', 'custo_total_unidade']


def write_json(path, obj):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str), encoding='utf-8')
    os.replace(tmp, path)


def jobs():
    conditions = [('NSGA2', limit, margin) for limit in (60, 15, 20, 30)
                  for margin in (True, False)]
    conditions += [('RANDOM_ARCHIVE', 60, True), ('NSGA3', 60, True)]
    return [dict(method=m, limit_mm=l, margin=q, seed=s,
                 job=f'{m}_L{l}_' + ('Q' if q else 'POINT') + f'_seed{s}')
            for m, l, q in conditions for s in f4.RUN_SEEDS]


def configure(smoke):
    if pymoo.__version__ != '0.6.1.5':
        raise RuntimeError('Validated pymoo version is 0.6.1.5. Do not silently change it.')
    root = ROOT / ('smoke' if smoke else 'final')
    root.mkdir(parents=True, exist_ok=True)
    dose = f4.load_surrogate_bundle(f4.DOSE_MODEL_PATH, f4.EXPECTED_DOSE_TARGET)
    pen = f4.load_surrogate_bundle(f4.PEN_MODEL_PATH, f4.EXPECTED_PEN_TARGET)
    evaluator = f4.Evaluator(dose, pen)
    doe_path = f4.resolve_first_existing(f4.DOE_CANDIDATES, 'DoE')
    doe = f4.load_doe(doe_path)
    p, t = f4.arrays_from_doe(doe)
    ev = evaluator.evaluate(p, t)
    ideal, ref = f4.build_hv_anchor(ev['F'])
    paths = [Path(__file__), Path(f4.__file__), f4.DOSE_MODEL_PATH,
             f4.PEN_MODEL_PATH, doe_path,
             f4.PHASE3_FEATURE_DIR / 'feature_engineering_final.py']
    n, generations = (24, 3) if smoke else (300, 200)
    selected = jobs()
    if smoke:
        selected = [j for j in selected if j['seed'] == f4.RUN_SEEDS[0]]
    plan = dict(version=VERSION, smoke=smoke, population=n, generations=generations,
                evaluations_per_run=n*(generations+1), jobs=selected,
                hashes={str(p): f4.sha256_file(p) for p in paths},
                python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                sklearn=f4.sklearn.__version__, pymoo=pymoo.__version__, xgboost=xgboost.__version__,
                platform=platform.platform(), ideal=ideal.tolist(), reference=ref.tolist(),
                q_pen=evaluator.q_pen, checkpoint_interval=50, metrics_interval=5,
                reference_directions='Das-Dennis 3 objectives H=23: 300 directions (smoke H=5:21)',
                random_baseline='independent Dirichlet/permutation draws; bounded rank/crowding archive N=300',
                primary='final population feasible front HV, un-clipped fixed DoE affine normalization',
                secondary='IGD/generalized spread against scenario-specific empirical pooled front',
                no_HF_runs=True, posthoc_sensitivity=True)
    plan_path = root / 'plan_frozen.json'
    if plan_path.exists() and json.loads(plan_path.read_text(encoding='utf-8')) != plan:
        raise RuntimeError('Frozen plan/input/environment changed. Preserve this study and use a new version.')
    write_json(plan_path, plan)
    baseline = []
    for limit in (60,15,20,30):
        for margin in (True,False):
            e = apply_constraint(ev, dict(limit_mm=limit,margin=margin), evaluator.q_pen)
            ids = front(e)
            baseline.append(dict(limit_mm=limit,margin=margin,DOE_n=len(p),
                                 feasible_n=int(e['feasible'].sum()),
                                 hv=hv(e['F'][ids],ideal,ref),
                                 point_pass_margin_fail=int(np.sum((ev['pen']<=limit)&(ev['pen']+evaluator.q_pen>limit)))))
    pd.DataFrame(baseline).to_csv(root / 'DOE_baseline_by_scenario.csv',index=False)
    sig = hashlib.sha256(json.dumps(plan, sort_keys=True).encode()).hexdigest()
    write_json(root / 'preflight.json', dict(status='PASS', signature=sig,
                                          timestamp=f4.utc_now(), total_jobs=len(selected)))
    return root, plan, sig, evaluator, ideal, ref


def apply_constraint(ev, job, q):
    ev = {k: v.copy() for k, v in ev.items()}
    # Intervals remain the original calibrated intervals even in the point-only arm.
    ev['cv'] = np.maximum(0., ev['pen'] + (q if job['margin'] else 0.) - job['limit_mm'])
    ev['feasible'] = ev['cv'] <= f4.CONSTRAINT_TOL_MM
    return ev


def select_nsga3(p, t, ev, survival):
    problem = Problem(n_var=14, n_obj=3, n_ieq_constr=1)
    pop = Population.new(F=ev['F'], G=ev['cv'][:, None], source_idx=np.arange(len(p)))
    selected = survival.do(problem, pop, n_survive=len(p)//2)
    idx = selected.get('source_idx').astype(int)
    return p[idx], t[idx], f4.eval_take(ev, idx)


def front(ev):
    feasible = np.flatnonzero(ev['feasible'])
    if not len(feasible):
        return feasible
    return feasible[f4.pure_pareto_indices(ev['F'][feasible])]


def hv(F, ideal, ref):
    # Unlike the historical implementation, do not clip values below the DoE ideal.
    return f4.hypervolume_3d((F-ideal)/(ref-ideal), np.ones(3)) if len(F) else 0.


def run_job(root, plan, signature, evaluator, ideal, ref, job, stop_at=None):
    folder = root / job['job']
    folder.mkdir(exist_ok=True)
    done_path = folder / 'complete.json'
    if done_path.exists():
        done = json.loads(done_path.read_text())
        if done['signature'] != signature:
            raise RuntimeError('Completed job signature mismatch')
        print('SKIP complete ' + job['job'], flush=True)
        return
    checkpoint = folder / 'checkpoint.pkl'
    n, generations = plan['population'], plan['generations']
    rng = np.random.default_rng(job['seed'])
    np.random.seed(job['seed'])  # pymoo 0.6.1.5 niching uses the legacy RNG
    survival = ReferenceDirectionSurvival(get_reference_directions('das-dennis', 3,
                                          n_partitions=5 if plan['smoke'] else 23))
    history, elapsed, count, feasible_seen = [], 0., 0, 0
    if checkpoint.exists():
        with checkpoint.open('rb') as fh:
            state = pickle.load(fh)
        if state['signature'] != signature:
            raise RuntimeError('Checkpoint signature mismatch')
        p, t, ev = state['p'], state['t'], state['ev']
        history, elapsed, count = state['history'], state['elapsed'], state['count']
        feasible_seen = state['feasible_seen']
        rng.bit_generator.state = state['rng']
        np.random.set_state(state['legacy_rng'])
        survival = state['survival']
        start = state['generation'] + 1
    else:
        p, t = f4.initialize_population(rng, n)
        tic = time.perf_counter()
        ev = apply_constraint(evaluator.evaluate(p, t), job, evaluator.q_pen)
        elapsed += time.perf_counter()-tic
        count, feasible_seen, start = n, int(ev['feasible'].sum()), 0
        f4.design_dataframe(p, t).to_csv(folder / 'initial_population.csv', index=False)
    print('START/RESUME ' + job['job'] + f' generation={start}', flush=True)
    tic = time.perf_counter()
    for gen in range(start, generations+1):
        if gen:
            if job['method'] == 'RANDOM_ARCHIVE':
                cp, ct = f4.initialize_population(rng, n)
            else:
                if job['method'] == 'NSGA3':
                    # Standard NSGA-III feasibility tournament: CV then random.
                    ranks, crowd = ev['cv'], np.zeros(n)
                else:
                    ranks, crowd, _ = f4.ranks_and_crowding(ev['F'], ev['cv'])
                cp, ct = f4.make_offspring(p, t, ranks, crowd, rng)
            ce = apply_constraint(evaluator.evaluate(cp, ct), job, evaluator.q_pen)
            count += n
            feasible_seen += int(ce['feasible'].sum())
            pp, tt, ee = np.vstack([p, cp]), np.vstack([t, ct]), f4.eval_concat(ev, ce)
            if job['method'] == 'NSGA3':
                p, t, ev = select_nsga3(pp, tt, ee, survival)
            else:
                p, t, ev = f4.environmental_selection(pp, tt, ee, n)
            f4.validate_population(p, t)
        if gen % 5 == 0 or gen == generations:
            ids = front(ev)
            row = dict(generation=gen, n_evaluations=count, hv=hv(ev['F'][ids], ideal, ref),
                       feasible_fraction=float(ev['feasible'].mean()), front_size=len(ids),
                       sampled_feasible_fraction=feasible_seen/count,
                       near_limit_1mm=int(np.sum(np.abs(ev['pen']+(evaluator.q_pen if job['margin'] else 0)-job['limit_mm']) <= 1)),
                       elapsed_seconds=elapsed+time.perf_counter()-tic)
            history.append(row)
            pd.DataFrame(history).to_csv(folder / 'convergence.csv', index=False)
            print(f"{job['job']} gen={gen}/{generations} eval={count} HV={row['hv']:.6f} feasible={row['feasible_fraction']:.3f}", flush=True)
        if gen % 50 == 0 or gen == generations or gen == stop_at:
            state = dict(signature=signature, generation=gen, p=p, t=t, ev=ev,
                         history=history, elapsed=elapsed+time.perf_counter()-tic,
                         count=count, feasible_seen=feasible_seen, rng=rng.bit_generator.state,
                         legacy_rng=np.random.get_state(), survival=survival)
            f4.atomic_pickle_dump(state, checkpoint)
        if gen == stop_at and gen < generations:
            return
    assert count == n*(generations+1)
    df = f4.dataframe_from_population(p, t, ev)
    df = df.rename(columns={'ballistic_feasible_95PI': 'feasible_under_study_rule'})
    df.to_csv(folder / 'final_population.csv', index=False)
    df.iloc[front(ev)].to_csv(folder / 'pareto_final.csv', index=False)
    write_json(done_path, dict(signature=signature, **job, **history[-1],
                              completed_utc=f4.utc_now()))


def quality(A, R):
    if not len(A) or not len(R):
        return np.nan, np.nan
    tree = cKDTree(A)
    igd = float(tree.query(R)[0].mean())
    if len(A) < 2:
        return igd, np.nan
    nearest = tree.query(A, k=2)[0][:, 1]
    edge = tree.query(R[np.argmin(R, axis=0)])[0].sum()
    denom = edge + len(A)*nearest.mean()
    spread = (edge+np.abs(nearest-nearest.mean()).sum())/denom if denom else 0.
    return igd, float(spread)


def pooled_front(F):
    """Exact 3-D minimization skyline, O(n log n) memory-safe pooled reference."""
    if not len(F):
        return F
    F=np.unique(F,axis=0)  # lexicographically sorted x,y,z
    ys=np.unique(F[:,1]); tree=np.full(len(ys)+1,np.inf); keep=[]
    for row in F:
        pos=int(np.searchsorted(ys,row[1]))+1
        k=pos; best=np.inf
        while k:
            best=min(best,tree[k]);k-=k&-k
        if best>row[2]:keep.append(row)
        k=pos
        while k<len(tree):
            tree[k]=min(tree[k],row[2]);k+=k&-k
    return np.asarray(keep).reshape(-1,3)


def summarize(root, plan, sig, ideal, ref):
    rows, frames = [], {}
    missing = []
    for j in plan['jobs']:
        path = root / j['job'] / 'complete.json'
        if not path.exists():
            missing.append(j['job'])
            continue
        row = json.loads(path.read_text())
        if row['signature'] != sig:
            raise RuntimeError('Summary signature mismatch')
        rows.append(row)
        frames[j['job']] = pd.read_csv(path.parent / 'pareto_final.csv')
    if missing:
        write_json(root / 'incomplete_jobs.json', missing)
        raise RuntimeError(f'{len(missing)} jobs incomplete. Resume the same run command.')
    summary = pd.DataFrame(rows)
    for (limit, margin), group in summary.groupby(['limit_mm', 'margin']):
        pool = np.vstack([frames[j][OBJ].to_numpy() for j in group.job])
        pool = np.unique(pool, axis=0)
        R = pooled_front(pool)
        pd.DataFrame(R, columns=OBJ).to_csv(root / f'empirical_reference_L{limit}_margin{margin}.csv', index=False)
        for i, row in group.iterrows():
            A = np.unique(frames[row.job][OBJ].to_numpy(), axis=0)
            igd, spread = quality((A-ideal)/(ref-ideal), (R-ideal)/(ref-ideal))
            summary.loc[i, 'igd_empirical'] = igd
            summary.loc[i, 'generalized_spread'] = spread
    summary.to_csv(root / 'comparison_per_run.csv', index=False)
    summary.groupby(['method', 'limit_mm', 'margin'])[['hv', 'igd_empirical', 'generalized_spread', 'elapsed_seconds']].agg(['mean','std','median','min','max']).to_csv(root / 'comparison_summary.csv')
    pairs = []
    comparisons = [('NSGA2',60,True,m,60,True) for m in ('RANDOM_ARCHIVE','NSGA3')]
    comparisons += [('NSGA2',l,True,'NSGA2',l,False) for l in (60,15,20,30)]
    for ma, la, qa, mb, lb, qb in comparisons:
        a = summary.query('method == @ma and limit_mm == @la and margin == @qa').set_index('seed')
        b = summary.query('method == @mb and limit_mm == @lb and margin == @qb').set_index('seed')
        d = (a.hv-b.hv).to_numpy()
        pval = 1. if np.all(np.abs(d)<1e-14) else float(wilcoxon(d, method='auto').pvalue)
        rng = np.random.default_rng(20261012)
        boots = rng.choice(d, (10000, len(d)), replace=True).mean(axis=1)
        pairs.append(dict(A=f'{ma}_L{la}_Q{qa}', B=f'{mb}_L{lb}_Q{qb}', n_pairs=len(d),
                          mean_HV_A_minus_B=float(d.mean()), CI95_low=float(np.quantile(boots,.025)),
                          CI95_high=float(np.quantile(boots,.975)), p_wilcoxon=pval))
    order = np.argsort([p['p_wilcoxon'] for p in pairs])
    adjusted = 0.
    for rank, idx in enumerate(order):
        adjusted = max(adjusted, min(1., pairs[idx]['p_wilcoxon']*(len(pairs)-rank)))
        pairs[idx]['p_Holm_family6'] = adjusted
    pd.DataFrame(pairs).to_csv(root / 'paired_HV_comparisons.csv', index=False)
    ablation=[]
    boundaries=[]
    design_cols=f4.EXPECTED_INPUT_COLUMNS
    for limit in (60,15,20,30):
        for seed in sorted(summary.seed.unique()):
            jpoint=f'NSGA2_L{limit}_POINT_seed{seed}'
            jq=f'NSGA2_L{limit}_Q_seed{seed}'
            point, margin = frames[jpoint],frames[jq]
            a=set(map(tuple,point[design_cols].to_numpy()))
            b=set(map(tuple,margin[design_cols].to_numpy()))
            pen=point['max_penetration_pred_mm'].to_numpy()
            ablation.append(dict(limit_mm=limit,seed=seed,point_front_n=len(point),margin_front_n=len(margin),
                                 design_Jaccard=len(a&b)/len(a|b) if a|b else np.nan,
                                 point_front_rejected_by_margin=int(np.sum(pen+plan['q_pen']>limit+f4.CONSTRAINT_TOL_MM))))
            if len(margin):
                chosen=margin.assign(distance_to_limit_mm=limit-margin['penetration_pi_upper95_mm'])
                chosen=chosen.sort_values('distance_to_limit_mm').head(3).copy()
                chosen['sensitivity_limit_mm']=limit;chosen['seed']=seed
                boundaries.append(chosen)
    pd.DataFrame(ablation).to_csv(root / 'margin_ablation_per_seed.csv',index=False)
    if boundaries:
        pd.concat(boundaries).drop_duplicates(subset=design_cols+['sensitivity_limit_mm']).to_csv(
            root / 'near_boundary_candidates_NOT_HF_VALIDATED.csv',index=False)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8,5))
    for method in ('NSGA2','RANDOM_ARCHIVE','NSGA3'):
        group = summary.query('method == @method and limit_mm == 60 and margin == True')
        histories = [pd.read_csv(root / j / 'convergence.csv') for j in group.job]
        vals = np.array([h.hv for h in histories])
        x = histories[0].n_evaluations
        ax.plot(x, vals.mean(axis=0), label=method)
        ax.fill_between(x, vals.min(axis=0), vals.max(axis=0), alpha=.15)
    ax.set(xlabel='Surrogate evaluations', ylabel='Fixed-reference normalized hypervolume')
    ax.legend(); fig.tight_layout(); fig.savefig(root / 'equal_budget_convergence.pdf'); plt.close(fig)
    write_json(root / 'study_complete.json', dict(status='COMPLETE', jobs=len(rows),
               signature=sig, smoke=plan['smoke'], timestamp=f4.utc_now(),
               interpretation='Surrogate-space comparisons only. No experimental or HF evidence added.'))
    print('COMPLETE: ' + str(root), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['preflight','run','summarize'])
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()
    root, plan, sig, evaluator, ideal, ref = configure(args.smoke)
    print(f"PLAN: {len(plan['jobs'])} runs x {plan['evaluations_per_run']} evaluations. Output: {root}", flush=True)
    if args.action == 'preflight':
        return
    lock = root / 'RUNNING.lock'
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise RuntimeError(f'Existing lock {lock}. Confirm previous Python is stopped before removing only this lock.')
    os.write(descriptor, f'PID={os.getpid()} UTC={f4.utc_now()}'.encode()); os.close(descriptor)
    try:
        if args.action == 'run':
            for job in plan['jobs']:
                run_job(root, plan, sig, evaluator, ideal, ref, job)
        summarize(root, plan, sig, ideal, ref)
    finally:
        lock.unlink(missing_ok=True)


if __name__ == '__main__':
    main()
