"""Kondo MPFR oracle accuracy, bounded field grids and export contracts."""
import csv
from decimal import Decimal,localcontext
import json,os
from pathlib import Path
import subprocess,sys,tempfile
program,enabled=sys.argv[1:]
enabled=enabled.upper() in ('ON','TRUE','1')
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',UNI20_COLOR='never')
def run(args,code=0):
    p=subprocess.run([program,*map(str,args)],capture_output=True,text=True,env=env,timeout=240)
    assert p.returncode==code,(args,p.returncode,p.stdout,p.stderr)
    return p.stdout
def table(args,bits,code=0):
    d=json.loads(run([*args,'--format','json'],code))
    t=next(iter(d['tables'].values())) if 'tables' in d else d
    assert int(t['metadata']['Working precision (bits)'])==bits
    for row in t['rows']:
        for cell in row:
            if isinstance(cell,dict) and 'precision' in cell:assert cell['precision']==bits,cell
    return t
def rows(t):
    return [dict(zip([c['id'] for c in t['columns']],[c['value'] if isinstance(c,dict) else c for c in row])) for row in t['rows']]
base=['--field','1','--scale','1']
if enabled:
    refs={'1': {'magnetization': '0.19211164621093772781603846463983704837166745414312923905831121436958852701905427577840832527934417128241503293479403595707395', 'energy': '-0.1067922423030353039230318678651420002441806867445409557659118577044758869637202097108773746557594663999694233177803816046091'}, '2': {'magnetization': '0.2751636379106730468713672485815389261523262796848405221986864233858923994227302345525603596294153041830138422145228546957692', 'energy': '-0.34629170027457867780906807893015178590470848543255021208916271877973380318861316866770553677958007446373773403295083008090387'}}
    with localcontext() as ctx:
        ctx.prec=140
        for field in ('1','2'):
            previous=Decimal(1)
            for bits in ((80,160,256,333) if field=='1' else (80,160)):
                r=rows(table(['--field',field,'--scale','1','--precision',f'bits-{bits}'],bits))[0]
                assert r['status']=='converged'
                em=abs(Decimal(r['magnetization'])-Decimal(refs[field]['magnetization']))
                ee=abs(Decimal(r['energy_change'])-Decimal(refs[field]['energy']))/Decimal(field)
                assert em<=Decimal(r['magnetization_error']) and ee<=Decimal(r['scaled_energy_error']),(bits,r,em,ee)
                error=max(em,ee)
                assert error<previous*Decimal('1e-14'),(field,bits,error,previous)
                previous=error
            assert previous<Decimal('1e-92' if field=='1' else '1e-42')
        # Non-dyadic fractions, crossing zero and the low/high-field crossover.
        grid=['--field','-2','--field-end','2','--points','10','--scale','1','--precision','bits-128','--tolerance','1e-10']
        serial=table(grid,128)
        assert serial['rows']==table([*grid,'--threads','2','--no-retain'],128)['rows']
        for i,r in enumerate(rows(serial)):
            assert abs(Decimal(r['field'])-(-2+Decimal(4)*i/9))<Decimal('1e-36')
        # More rows than the admission window; all cheap failed solves are drained.
        failed=['--field','-2','--field-end','2','--points','130','--scale','1','--precision','bits-128','--max-series-terms','0']
        a=table(failed,128,2)
        assert a['rows']==table([*failed,'--threads','2','--no-retain'],128,2)['rows']
        assert all(r['magnetization'] is None and r['energy_change'] is None and r['status']=='series_limit' for r in rows(a))
        for option,status in [('--max-lobes','tail_limit'),('--max-evaluations','quadrature_limit')]:
            r=rows(table(['--field','2','--scale','1','--precision','bits-160',option,'0'],160,2))[0]
            assert r['status']==status and r['magnetization'] is None
        literal='1.0000000000000000000000000000000000000000000000000000001'
        t=table(['--field',literal,'--scale','2','--precision','digits-100','--tolerance','1e-20'],333)
        assert abs(Decimal(rows(t)[0]['field'])-Decimal(literal))<Decimal('1e-95')
    with tempfile.TemporaryDirectory() as folder:
        args=[*base,'--precision','bits-160']
        expected=rows(table(args,160))
        for ext,sep in [('csv',','),('tsv','\t')]:
            path=Path(folder)/('response.'+ext)
            run([*args,f'--{ext}-table',f'response={path}','--quiet','--no-retain'])
            actual=list(csv.DictReader([line for line in path.read_text().splitlines() if not line.startswith('#')],delimiter=sep))
            for a,e in zip(actual,expected):
                for key,value in e.items():
                    assert a[key]==(str(value).lower() if isinstance(value,bool) else str(value)),(key,a[key],value)
            assert len(actual)==len(expected)
with tempfile.TemporaryDirectory() as folder:
    path=Path(folder)/'preserve.json'
    path.write_text('sentinel')
    invalid=[['--precision',x] for x in ('bits-0','digits--1','bits-128junk')]
    if enabled:
        invalid += [['--precision','bits-128','--tolerance',x] for x in ('0','-1','nan','inf')]
        invalid += [['--precision','bits-128','--max-series-terms','4097']]
    else:invalid += [['--precision',x] for x in ('bits-160','digits-100')]
    for args in invalid:
        run([*base,*args,'--json',path,'--force'],1)
        assert path.read_text()=='sentinel'
print('Kondo runtime precision checks passed')
