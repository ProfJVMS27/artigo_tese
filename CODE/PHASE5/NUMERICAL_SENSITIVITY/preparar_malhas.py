"""Prepare isolated structured mesh studies from executed decks, never run LS-DYNA."""
import hashlib
import itertools
import json
import math
import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent/'VALIDACAO_ALTA_FIDELIDADE/fase4_20261003_043316_MCNP63_DEDF_AP_30/study_30/lsdyna'
ALLOWED = {'*NODE','*ELEMENT_SOLID','*INITIAL_VELOCITY_NODE','*BOUNDARY_SPC_NODE','*DATABASE_HISTORY_NODE'}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def blocks(text):
    return re.findall(r'^\*[^\n]*(?:\n(?!\*)[^\n]*)*',text,re.M)

def key(b):
    return b.splitlines()[0].strip()

def data(b):
    return [s for s in b.splitlines()[1:] if s.strip() and not s.startswith('$')]

def mesh(bs):
    nodes = {int(s[:8]):tuple(float(s[8+i*16:24+i*16]) for i in range(3)) for s in data(next(b for b in bs if key(b)=='*NODE'))}
    elements = [tuple(int(s[i:i+8]) for i in range(0,80,8)) for s in data(next(b for b in bs if key(b)=='*ELEMENT_SOLID'))]
    return nodes,elements

def quality(nodes, elems):
    xyz=np.array([nodes[n] for e in elems for n in e[2:]]).reshape(-1,8,3)
    signs=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
    volumes=np.zeros(len(elems)); minimum=float('inf')
    for pt in itertools.product((-1/math.sqrt(3),1/math.sqrt(3)),repeat=3):
        d=np.zeros((8,3))
        for j in range(3):
            other=[k for k in range(3) if k!=j]
            d[:,j]=signs[:,j]*np.prod(1+signs[:,other]*np.array(pt)[other],axis=1)/8
        det=np.linalg.det(np.einsum('eni,nj->eij',xyz,d))
        minimum=min(minimum,float(det.min())); volumes+=det
    assert minimum>0, minimum
    perpart={str(pid):float(sum(v for e,v in zip(elems,volumes) if e[1]==pid)) for pid in sorted({e[1] for e in elems})}
    return {'min_gauss_jacobian':minimum,'volume_mm3_by_part':perpart}

def refine_axis(a, lo, hi, split):
    out=[]
    for x,y in zip(a,a[1:]):
        n=split if x>=lo-1e-6 and y<=hi+1e-6 else 1
        out.extend(x+(y-x)*j/n for j in range(n))
    return out+[a[-1]]

def prepare(case, level):
    src=SOURCE/f'sim_{case}/sim_{case}.k'; original=src.read_text(); digest=sha(src)
    bs=blocks(original); nodes,elems=mesh(bs)
    fsp=[e for e in elems if e[1]==99]; assert len(fsp)==1
    original_corners=fsp[0][2:]; assert len(nodes)==19853 and len(elems)==17473
    # Recover exact structured grids from executed node order, not a mutable generator.
    nt=next(i for i in range(2,19846) if abs(nodes[i][2]-nodes[1][2])>1e-5)-1
    stride=next(i for i in range(2,19846) if abs(math.hypot(*nodes[i][:2])-350)>1e-4)-1
    nz=stride//nt; assert stride%nt==0 and 19845==15*stride
    theta=[math.atan2(nodes[i+1][1],nodes[i+1][0]) for i in range(nt)]
    z=[nodes[1+i*nt][2] for i in range(nz)]
    r=[math.hypot(*nodes[1+i*stride][:2]) for i in range(15)]
    pids=[elems[i*2*(nt-1)*(nz-1)][1] for i in range(7)]
    assert len(set(pids))==7
    # Validate source connectivity and layer mapping before rebuilding.
    expected=[]
    for ir in range(14):
        for iz in range(nz-1):
            for it in range(nt-1):
                a=1+ir*stride+iz*nt+it; b=a+stride
                expected.append((len(expected)+1,pids[ir//2],a,a+1,a+nt+1,a+nt,b,b+1,b+nt+1,b+nt))
    assert expected==elems[:-1]
    th=refine_axis(theta,-math.radians(5),math.radians(5),level+1)
    zz=refine_axis(z,1225,1275,level+1)
    nr=2+level
    rr=[r[2*j]+(r[2*j+2]-r[2*j])*i/nr for j in range(7) for i in range(nr)]+[r[-1]]
    nn={}; ee=[]
    for rad in rr:
        for height in zz:
            for angle in th:
                nn[len(nn)+1]=(rad*math.cos(angle),rad*math.sin(angle),height)
    nt2=len(th); stride2=len(th)*len(zz)
    for ir in range(len(rr)-1):
        for iz in range(len(zz)-1):
            for it in range(nt2-1):
                a=1+ir*stride2+iz*nt2+it; b=a+stride2
                ee.append((len(ee)+1,pids[ir//nr],a,a+1,a+nt2+1,a+nt2,b,b+1,b+nt2+1,b+nt2))
    fixed=[i for i,p in nn.items() if abs(p[2]-900)<1e-7]
    oldfixed=[int(s[:10]) for s in data(next(b for b in bs if key(b)=='*BOUNDARY_SPC_NODE'))]
    assert set(oldfixed)=={i for i,p in nodes.items() if abs(p[2]-900)<1e-7}
    corners=np.array([nodes[i] for i in original_corners]); low=corners.min(axis=0); high=corners.max(axis=0)
    f=level+1; start=len(nn)+1
    def fid(i,j,k): return start+(k*(f+1)+j)*(f+1)+i
    for k in range(f+1):
        for j in range(f+1):
            for i in range(f+1): nn[fid(i,j,k)]=tuple(low+(high-low)*np.array([i,j,k])/f)
    for k in range(f):
        for j in range(f):
            for i in range(f):
                ee.append((len(ee)+1,99,fid(i,j,k),fid(i+1,j,k),fid(i+1,j+1,k),fid(i,j+1,k),fid(i,j,k+1),fid(i+1,j,k+1),fid(i+1,j+1,k+1),fid(i,j+1,k+1)))
    fids=list(range(start,len(nn)+1))
    cornerids=[fid(i*f,j*f,k*f) for i,j,k in ((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))]
    assert len(nn)<128000 and len(ee)<128000
    velocity=data(next(b for b in bs if key(b)=='*INITIAL_VELOCITY_NODE'))
    assert len({s[10:] for s in velocity})==1
    assert {int(s[:10]) for s in velocity}==set(original_corners)
    spc=data(next(b for b in bs if key(b)=='*BOUNDARY_SPC_NODE'))
    assert len({s[10:] for s in spc})==1
    replacement={
        '*NODE':'*NODE\n'+''.join(f'{i:8d}'+''.join(f'{v:16.8E}' for v in p)+f'{0:8d}{0:8d}\n' for i,p in nn.items()),
        '*ELEMENT_SOLID':'*ELEMENT_SOLID\n'+''.join(''.join(f'{v:8d}' for v in e)+'\n' for e in ee),
        '*INITIAL_VELOCITY_NODE':'*INITIAL_VELOCITY_NODE\n'+''.join(f'{i:10d}'+velocity[0][10:]+'\n' for i in fids),
        '*BOUNDARY_SPC_NODE':'*BOUNDARY_SPC_NODE\n'+''.join(f'{i:10d}'+spc[0][10:]+'\n' for i in fixed),
        '*DATABASE_HISTORY_NODE':'*DATABASE_HISTORY_NODE\n'+''.join(''.join(f'{i:10d}' for i in fids[j:j+8])+'\n' for j in range(0,len(fids),8))}
    result=''.join(replacement.get(key(b),b.rstrip('\n')+'\n') for b in bs)
    rb=blocks(result)
    assert [(key(b),b.rstrip()) for b in bs if key(b) not in ALLOWED]==[(key(b),b.rstrip()) for b in rb if key(b) not in ALLOWED]
    parsedn,parsede=mesh(rb); assert parsede==ee
    qa=quality(parsedn,parsede); baseqa=quality(nodes,elems)
    volratio={pid:qa['volume_mm3_by_part'][pid]/v-1 for pid,v in baseqa['volume_mm3_by_part'].items()}
    assert max(abs(v) for v in volratio.values())<0.001
    assert abs(volratio['99'])<1e-6
    folder=ROOT/'malhas'/f'sim_{case}_M{level}'; folder.mkdir(parents=True,exist_ok=True)
    dest=folder/f'sim_{case}_M{level}.k'
    if dest.exists(): assert dest.read_text()==result, 'Refuse to overwrite different deck'
    else: dest.write_text(result,encoding='ascii')
    assert sha(src)==digest
    info=dict(case=case,level=level,source=str(src),source_sha256=digest,keyword=str(dest),keyword_sha256=sha(dest),nodes=len(nn),elements=len(ee),fsp_nodes=fids,fsp_original_corner_nodes=cornerids,radial_per_layer=nr,central_theta_deg=0.5/(level+1),central_z_mm=5/(level+1),fsp_subdivision=f,volume_relative_change_by_part=volratio,quality=qa,physics_blocks_unchanged=True)
    (folder/'mesh_manifest.json').write_text(json.dumps(info,indent=2))
    print(case,level,len(nn),len(ee),'QA PASS',flush=True)
    return info

if __name__=='__main__':
    rows=[prepare(case,level) for case in (9012,9015,9021) for level in (1,2)]
    (ROOT/'mesh_preparation_manifest.json').write_text(json.dumps(rows,indent=2))
