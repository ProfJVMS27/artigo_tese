"""Build two factor-isolation decks from verified original and M1 meshes."""
import json
from pathlib import Path
import preparar_malhas as p

out=p.ROOT.parent/'CONVERGENCIA_9012_ISOLADO'
base=next(r for r in json.loads((p.ROOT/'mesh_preparation_manifest.json').read_text()) if r['case']==9012 and r['level']==1)
orig=Path(base['source']); fine=Path(base['keyword'])
assert p.sha(orig)==base['source_sha256'] and p.sha(fine)==base['keyword_sha256']
ob=p.blocks(orig.read_text()); fb=p.blocks(fine.read_text())
rows=[]
for label,sb,pb in [('BLINDAGEM_M1_FSP_ORIGINAL',fb,ob),('BLINDAGEM_ORIGINAL_FSP_M1',ob,fb)]:
    sn,se=p.mesh(sb); fn,fe=p.mesh(pb)
    se=[e for e in se if e[1]!=99]; fe=[e for e in fe if e[1]==99]
    sids={i for e in se for i in e[2:]}; fids=sorted({i for e in fe for i in e[2:]})
    nn={i:sn[i] for i in sorted(sids)}
    mapping={i:max(sids)+1+j for j,i in enumerate(fids)}
    nn.update({mapping[i]:fn[i] for i in fids})
    ee=list(se)+[(len(se)+1+j,99,*(mapping[i] for i in e[2:])) for j,e in enumerate(fe)]
    corners=[]
    on,oe=p.mesh(ob)
    oldcorners=next(e[2:] for e in oe if e[1]==99)
    for i in oldcorners:
        match=[mapping[j] for j in fids if max(abs(a-b) for a,b in zip(fn[j],on[i]))<1e-6]
        assert len(match)==1; corners.extend(match)
    velocity=p.data(next(b for b in ob if p.key(b)=='*INITIAL_VELOCITY_NODE'))[0][10:]
    replacement={
      '*NODE':'*NODE\n'+''.join(f'{i:8d}'+''.join(f'{v:16.8E}' for v in xyz)+f'{0:8d}{0:8d}\n' for i,xyz in nn.items()),
      '*ELEMENT_SOLID':'*ELEMENT_SOLID\n'+''.join(''.join(f'{v:8d}' for v in e)+'\n' for e in ee),
      '*INITIAL_VELOCITY_NODE':'*INITIAL_VELOCITY_NODE\n'+''.join(f'{i:10d}'+velocity+'\n' for i in mapping.values()),
      '*BOUNDARY_SPC_NODE':next(b for b in sb if p.key(b)=='*BOUNDARY_SPC_NODE').rstrip()+'\n',
      '*DATABASE_HISTORY_NODE':'*DATABASE_HISTORY_NODE\n'+''.join(''.join(f'{i:10d}' for i in list(mapping.values())[j:j+8])+'\n' for j in range(0,len(mapping),8))}
    text=''.join(replacement.get(p.key(b),b.rstrip('\n')+'\n') for b in ob)
    bs=p.blocks(text); n,e=p.mesh(bs)
    assert e==ee and len(n)==len(nn) and len(set(n))==len(nn)
    assert {i for elem in e for i in elem[2:]}==set(n)
    assert [b.rstrip() for b in bs if p.key(b) not in p.ALLOWED]==[b.rstrip() for b in ob if p.key(b) not in p.ALLOWED]
    # Each constituent's coordinates/connectivity are preserved, only FSP IDs change.
    assert {i:n[i] for i in sids}=={i:sn[i] for i in sids}
    assert all(n[mapping[i]]==fn[i] for i in fids)
    assert len(n)<128000 and len(e)<128000
    qa=p.quality(n,e)
    folder=out/label;folder.mkdir(parents=True,exist_ok=True);dest=folder/'sim_9012.k'
    if dest.exists(): assert dest.read_text()==text,'Refusing overwrite'
    else: dest.write_text(text,encoding='ascii')
    row=dict(case=9012,level=label,source=str(orig),source_sha256=p.sha(orig),keyword=str(dest),keyword_sha256=p.sha(dest),fsp_nodes=list(mapping.values()),fsp_original_corner_nodes=corners,nodes=len(n),elements=len(e),quality=qa)
    rows.append(row);print(label,len(n),len(e),'static PASS',flush=True)
(out/'mesh_preparation_manifest.json').write_text(json.dumps(rows,indent=2))
assert p.sha(orig)==base['source_sha256'] and p.sha(fine)==base['keyword_sha256']
