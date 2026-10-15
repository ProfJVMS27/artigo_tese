"""Sequential, streamed execution. Default is dry-run, --run explicitly starts solver.

Completed cases are skipped only with matching hashes and verified termination.
Interrupted/failed output is preserved and blocks automatic overwriting.
"""
import argparse
import hashlib
import importlib.util
import json
import re
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parent
EXE=Path(r'C:\Program Files\LS-DYNA Suite R16.1 Student\lsdyna\ls-dyna_smp_d_R16.1_180-gd50332dbe5_winx64_ifort190_sse2_studentversion.exe')

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def terminated(folder):
    return any('normaltermination' in re.sub(r'\s+','',p.read_text(encoding='latin-1').lower()) for p in [folder/'d3hsp',folder/'messag'] if p.exists())

def summarize(row):
    parser=ROOT.parent.parent/'Fase 2/LS-DYNA/summarize_lsdyna_study.py'
    spec=importlib.util.spec_from_file_location('parser',parser)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    folder=Path(row['keyword']).parent
    text=mod.read_text(folder/'nodout')
    result={}
    for label,ids in [('all_fsp_nodes',row['fsp_nodes']),('original_eight_corners',row['fsp_original_corner_nodes'])]:
        mod.NODE_IDS=set(ids)
        records=mod.nodout_records(text)
        if not records: raise RuntimeError(f'Missing complete nodout records: {folder}')
        b=mod.ballistic(records)
        if b['final_time_ms']<0.299: raise RuntimeError('Incomplete time history')
        result[label]=dict(b,n_records=len(records))
    if result['all_fsp_nodes']['n_records']!=result['original_eight_corners']['n_records']:
        raise RuntimeError('Incomplete FSP node history')
    er=mod.ratios(mod.read_text(folder/'glstat'))
    if not er: raise RuntimeError('Missing energy history')
    result['energy_ratio_final']=er[-1]
    result['max_energy_deviation_abs']=max(abs(x-1) for x in er)
    result['note']='Nodal mean velocity is not a mass-weighted COM velocity. Report both node definitions. Numerical convergence not certified by normal termination.'
    result['parser_sha256']=sha(parser)
    result['nodout_sha256']=sha(folder/'nodout')
    (folder/'postprocess.json').write_text(json.dumps(result,indent=2))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--run',action='store_true'); ap.add_argument('--ncpu',type=int,default=3)
    args=ap.parse_args()
    if args.ncpu<1: raise ValueError('ncpu must be positive')
    rows=json.loads((ROOT/'mesh_preparation_manifest.json').read_text())
    if not EXE.is_file(): raise FileNotFoundError(EXE)
    for row in rows:
        if sha(Path(row['keyword']))!=row['keyword_sha256']: raise RuntimeError('Study input hash changed')
        if sha(Path(row['source']))!=row['source_sha256']: raise RuntimeError('Original input hash changed')
    if not args.run:
        for r in rows: print('DRY RUN',r['case'],r['level'],r['nodes'],'nodes',r['elements'],'elements')
        print('No solver started. Use --run for sequential execution.'); return
    lock=ROOT/'MESH_RUNNING.lock'
    with lock.open('x') as f: f.write(datetime.now(timezone.utc).isoformat())
    try:
        for row in rows:
            p=Path(row['keyword']); folder=p.parent; done=folder/'completed.json'
            if done.exists():
                saved=json.loads(done.read_text())
                if saved['keyword_sha256']!=sha(p) or not terminated(folder): raise RuntimeError('Invalid completion marker')
                print('SKIP verified completed',folder.name,flush=True); continue
            if any((folder/n).exists() for n in ('d3hsp','messag','nodout','console.log')):
                raise RuntimeError(f'Existing incomplete/unverified output preserved: {folder}. Do not delete the lock while a solver is active.')
            if shutil.disk_usage(folder).free < 5*1024**3:
                raise RuntimeError('Less than 5 GiB free. Outputs preserved. Free disk space before continuing.')
            cmd=[str(EXE),f'i={p.name}',f'ncpu={args.ncpu}','memory=100m']
            print('START',folder.name,flush=True); start=time.monotonic()
            with (folder/'console.log').open('w',encoding='utf-8') as log:
                proc=subprocess.Popen(cmd,cwd=folder,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace',bufsize=1)
                for line in proc.stdout:
                    print(line,end='',flush=True); log.write(line); log.flush()
                code=proc.wait()
            if code!=0 or not terminated(folder): raise RuntimeError(f'Solver failure: {folder}, exit={code}')
            summarize(row)
            status=dict(keyword_sha256=sha(p),elapsed_s=time.monotonic()-start,return_code=code,utc=datetime.now(timezone.utc).isoformat())
            done.write_text(json.dumps(status,indent=2))
            with (ROOT/'execution_progress.jsonl').open('a') as f:
                f.write(json.dumps(dict(case=row['case'],level=row['level'],**status))+'\n')
            print('END',folder.name,'completed, not convergence-certified',flush=True)
    except BaseException:
        print('Stopped. Lock preserved. Confirm no solver is active before any manual recovery.',flush=True)
        raise
    else:
        lock.unlink()

if __name__=='__main__': main()
