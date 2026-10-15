"""Read existing solver histories only. No solver launches or input changes."""
import csv
import hashlib
import importlib.util
import json
import re
from pathlib import Path

ROOT = Path(r'C:\Users\Admin\Desktop\Scripts Tese 1')
STUDY = ROOT / 'Fase 5/VALIDACAO_ALTA_FIDELIDADE/fase4_20261003_043316_MCNP63_DEDF_AP_30/study_30'
OUT = Path(__file__).parent
parser = ROOT / 'Fase 2/LS-DYNA/summarize_lsdyna_study.py'
spec = importlib.util.spec_from_file_location('existing_parser', parser)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def tally(p):
    if not p.exists() and p.with_suffix('.txt').exists():
        p = p.with_suffix('.txt')
    text = p.read_text(encoding='latin-1')
    block = text.split('1tally fluctuation charts')[-1]
    assert re.search(r'tally\s+4\b', block), p
    rows = []
    for line in block.splitlines():
        vals = line.split()
        if len(vals) == 6 and vals[0].isdigit():
            try:
                rows.append(dict(zip(['nps','mean','r','vov','slope','fom'], map(float, vals))))
            except ValueError:
                pass
    assert rows and rows[-1]['nps'] == 100000000, p
    assert all(a['nps'] < b['nps'] for a,b in zip(rows,rows[1:])), p
    end = rows[-1]
    tail = [x for x in rows if x['nps'] >= 50000000]
    return {'path':str(p),'sha256':sha(p),'tfc':rows,
            'ten_checks_pass':bool(re.search(r'4\s+passed the 10 statistical checks',text)),
            'last_half_max_deviation_from_final_pct':max(abs(x['mean']/end['mean']-1)*100 for x in tail),
            'dose_5min_mSv_tfc_rounded':300*end['mean']}

def ballistic(folder):
    p = folder/'nodout'
    records = mod.nodout_records(mod.read_text(p))
    assert records, p
    b = mod.ballistic(records)
    peak = min(records,key=lambda x:x['front_x_mm'])
    at03 = [x for x in records if x['time_ms'] <= 0.300001]
    b.update(path=str(p),sha256=sha(p),n_records=len(records),
             peak_time_ms=peak['time_ms'],
             max_pen_at_030_mm=max(0,410-min(x['front_x_mm'] for x in at03)),
             penetration_gain_last_005_ms=max(0,min(x['front_x_mm'] for x in records if x['time_ms'] <= records[-1]['time_ms']-.05)-peak['front_x_mm']))
    k = folder/(folder.name+'.k')
    b['keyword_sha256'] = sha(k)
    return b

report = {'scope':'30 Phase-5 cases only, not the entire training campaign',
          'parser_sha256':sha(parser),'mcnp':{},'lsdyna':{},'endtim060':{}}
for n in range(9001,9031):
    report['mcnp'][str(n)] = tally(STUDY/f'mcnp/sim_{n}/output_{n}.o')
    report['lsdyna'][str(n)] = ballistic(STUDY/f'lsdyna/sim_{n}')
for n in (9016,9026):
    report['endtim060'][str(n)] = ballistic(STUDY/f'lsdyna_endtim060_sensitivity/sim_{n}_endtim060')
ordered = sorted(report['lsdyna'],key=lambda n:(report['lsdyna'][n]['max_penetration_mm'],int(n)))
report['mesh_selection'] = {'rule':'minimum, lower median (rank 15/30), maximum baseline penetration; ties by sim_id',
                          'ids':[ordered[0],ordered[14],ordered[-1]]}
m = report['mcnp']
ids = sorted(m,key=lambda n:(m[n]['dose_5min_mSv_tfc_rounded'],int(n)))
report['mcnp_selection'] = {'rule':'minimum dose, maximum dose, largest final R, deduplicated',
                          'ids':list(dict.fromkeys([ids[0],ids[-1],max(m,key=lambda n:m[n]['tfc'][-1]['r'])]))}
(OUT/'auditoria_historicos.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print('MCNP 10/10:',sum(x['ten_checks_pass'] for x in m.values()),'/30')
print('MCNP final R range:',min(x['tfc'][-1]['r'] for x in m.values()),max(x['tfc'][-1]['r'] for x in m.values()))
print('MCNP last-half maximum drift (%):',max(x['last_half_max_deviation_from_final_pct'] for x in m.values()))
print('Mesh selection:',report['mesh_selection'])
print('MCNP selection:',report['mcnp_selection'])
for n,b in report['endtim060'].items():
    print('Duration',n,b)
print(OUT/'auditoria_historicos.json')
