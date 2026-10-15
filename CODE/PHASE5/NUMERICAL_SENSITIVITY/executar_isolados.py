"""Run only two isolated 9012 studies. Without --run, dry-run only."""
import json
import sys
from pathlib import Path

root=Path(__file__).resolve().parent
sys.path.insert(0,str(root.parent/'CONVERGENCIA_20261014'))
import executar_malhas as runner
from diagnosticar_extremos import analyze

runner.ROOT=root
old_summary=runner.summarize
def summarize(row):
    old_summary(row)
    folder=Path(row['keyword']).parent
    report=analyze(Path(row['keyword']))
    (folder/'diagnostico_atividade.json').write_text(json.dumps(report,indent=2))
    p=folder/'postprocess.json'; data=json.loads(p.read_text())
    data['warning']='Raw nodal extrema may include deleted nodes. Consult diagnostico_atividade.json. Active-by-log values are diagnostic, not validated labels.'
    p.write_text(json.dumps(data,indent=2))
    print('DIAGNOSTIC active-by-log penetration:',report['peaks']['active_by_log']['penetration_mm'],flush=True)
runner.summarize=summarize
if __name__=='__main__': runner.main()
