"""Read-only solver forensics. Diagnostic estimates do not replace dataset labels."""
import hashlib
import json
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parent
NUMBER=r'[+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def analyze(keyword):
    folder=keyword.parent
    text=keyword.read_text()
    block=text.split('*ELEMENT_SOLID\n')[1].split('*')[0]
    elements=[list(map(int,l.split())) for l in block.splitlines() if l.strip() and not l.startswith('$')]
    fsp=[e for e in elements if e[1]==99]
    ids={i for e in fsp for i in e[2:]}
    adjacency={i:[e[0] for e in fsp if i in e[2:]] for i in ids}
    msg=(folder/'messag').read_text(encoding='latin1')
    deleted={int(i):float(t) for i,t in re.findall(r'node\s+number\s+(\d+) deleted at time\s+('+NUMBER+')',msg)}
    failed={int(i):float(t) for i,t in re.findall(r'solid element\s+(\d+)\s+failed at time\s+('+NUMBER+')',msg)}
    negatives={int(i):float(t) for i,t in re.findall(r'negative volume in solid element #\s*(\d+) cycle\s+\d+ time\s+('+NUMBER+')',msg)}
    for i,t in negatives.items(): failed[i]=min(failed.get(i,float('inf')),t)
    peaks={name:None for name in ('all_nodes','active_by_log')}
    t=None
    for line in (folder/'nodout').read_text(encoding='latin1').splitlines():
        match=re.search(r'\(\s*at time\s+('+NUMBER+')\s*\)',line)
        if match: t=float(match[1]); continue
        match=re.match(r'\s*(\d+)\s+',line)
        if t is None or not match or int(match[1]) not in ids: continue
        i=int(match[1]); v=[float(x) for x in re.findall(r'[+\-]?\d+\.\d+E[+\-]\d+',line)]
        if len(v)<12: continue
        active=deleted.get(i,float('inf'))>t and any(failed.get(e,float('inf'))>t for e in adjacency[i])
        for name,include in [('all_nodes',True),('active_by_log',active)]:
            if include and (peaks[name] is None or v[9]<peaks[name]['x_mm']):
                peaks[name]=dict(node=i,time_ms=t,x_mm=v[9],penetration_mm=410-v[9],node_deletion_time=deleted.get(i),adjacent_elements=adjacency[i],adjacent_failure_times={e:failed.get(e) for e in adjacency[i]})
    gl=(folder/'glstat').read_text(encoding='latin1')
    def series(label): return [float(x) for x in re.findall(r'^\s*'+re.escape(label)+r'\.+\s*('+NUMBER+')',gl,re.M)]
    ratio=series('total energy / initial energy'); sliding=series('sliding interface energy')
    final_mat=(folder/'matsum').read_text(encoding='latin1').split('mat.#=   99')[-1][:450]
    return dict(case=folder.name,peaks=peaks,fsp_elements=len(fsp),fsp_failed={e[0]:failed[e[0]] for e in fsp if e[0] in failed},fsp_negative_volume={e[0]:negatives[e[0]] for e in fsp if e[0] in negatives},fsp_deleted_nodes={i:t for i,t in deleted.items() if i in ids},energy_ratio_final=ratio[-1],max_energy_deviation=max(abs(x-1) for x in ratio),sliding_energy_final_J=sliding[-1],sliding_energy_min_J=min(sliding),final_fsp_matsum_excerpt=final_mat,files={n:dict(path=str(folder/n),sha256=sha(folder/n)) for n in ('messag','nodout','glstat','matsum')})

if __name__=='__main__':
    manifest=json.loads((ROOT/'mesh_preparation_manifest.json').read_text())
    paths=list(dict.fromkeys([row['source'] for row in manifest]))+[row['keyword'] for row in manifest]
    results=[analyze(Path(p)) for p in paths]
    out=dict(note='Active-by-log estimate assumes all deletion events are logged. No d3plot binary confirmation. No label replacement. Not convergence certification.',results=results)
    (ROOT/'diagnostico_extremos_20261015.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    for r in results:
        print(r['case'],'raw=',round(r['peaks']['all_nodes']['penetration_mm'],4),'active_by_log=',round(r['peaks']['active_by_log']['penetration_mm'],4),'FSP_negative_volume=',len(r['fsp_negative_volume']),'/',r['fsp_elements'],'E/E0=',r['energy_ratio_final'],'contact_J=',r['sliding_energy_final_J'],flush=True)
