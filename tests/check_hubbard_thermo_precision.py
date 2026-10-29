"""Explicit Hubbard precision, CPU/provider solves, failures and full exports."""
import csv
from decimal import Decimal, localcontext
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

lines, edges, arbitrary = sys.argv[1:]
enabled = lambda x: x.upper() in ('ON','TRUE','1')
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')

def run(program, args, code=0):
    p = subprocess.run([program,*args],text=True,capture_output=True,env=env,timeout=300)
    assert p.returncode == code, (args,p.returncode,p.stdout,p.stderr)
    return p

def table(program, args, code=0):
    doc=json.loads(run(program,[*args,'--format','json'],code).stdout)
    return next(iter(doc['tables'].values())) if 'tables' in doc else doc

def rows(t):
    return [dict(zip([c['id'] for c in t['columns']], (v['value'] if isinstance(v,dict) else v for v in r))) for r in t['rows']]

base=['--u','4','--precision','bits-128']
if not enabled(arbitrary):
    for program in (lines,edges):
        assert 'UNI20_ENABLE_MPREAL=ON' in run(program,[*base,'--momentum','0'],1).stderr
    raise SystemExit(0)

momentum='0.147435530616030589362559903539606986354219786134238672485504679036370550127469923860951005'
energy=Decimal('0.180195419367165283255675569344877169275747982857939066888904265913419968700554935899809673')
with localcontext() as context:
    context.prec=110
    errors=[]
    for bits in (128,256):
        t=table(lines,['--u','4','--branch','spinon','--momentum',momentum,'--precision',f'bits-{bits}'])
        assert int(t['metadata']['Working precision (bits)']) == bits
        r=rows(t)[0];assert r['status']=='converged'
        error=abs(Decimal(r['energy'])-energy);errors.append(error)
        assert error < Decimal('1e-32' if bits==128 else '1e-70'),r
    assert errors[1] < errors[0]*Decimal('1e-30')

    # Non-dyadic grid fractions remain exact rationals until multiplied by pi.
    # Check their eventual working precision and value, as well as complete exports.
    cases=[(lines,[*base,'--points','4','--tolerance','1e-10']),
           (edges,[*base,'--points','4','--tolerance','1e-10'])]
    for program,args in cases:
        serial=table(program,args)
        for r in serial['rows']:
            for cell in r:
                if isinstance(cell,dict) and 'precision' in cell:
                    assert cell['precision']==128,cell
        for i,r in enumerate(rows(serial)):
            fraction=Decimal(i % 4)/3
            expected=(fraction if r.get('branch')=='spinon' else 2*fraction-1)
            assert abs(Decimal(r['p_over_pi' if program==lines else 'momentum_over_pi'])-expected)<Decimal('1e-36'),r
        parallel=table(program,[*args,'--threads','2'])
        assert serial['columns']==parallel['columns'] and serial['rows']==parallel['rows']
        streamed=table(program,[*args,'--threads','2','--no-retain'])
        assert serial['rows']==streamed['rows']
        for fmt,separator in [('csv',','),('tsv','\t')]:
            out=run(program,[*args,'--format',fmt]).stdout
            parsed=list(csv.DictReader(io.StringIO('\n'.join(s for s in out.splitlines() if not s.startswith('#'))),delimiter=separator))
            assert len(parsed)==len(serial['rows'])
            for actual,expected in zip(parsed,rows(serial)):
                assert set(actual)==set(expected)
                for key,value in expected.items():
                    if value is None: assert actual[key]=='',key
                    elif isinstance(value,bool): assert actual[key]==str(value).lower(),(key,actual[key],value)
                    else:
                        try: assert Decimal(actual[key])==Decimal(str(value)),(key,actual[key],value)
                        except ArithmeticError: assert actual[key]==str(value),(key,actual[key],value)
    bad=table(lines,[*base,'--branch','holon','--momentum','1','--max-evaluations','0'],2)
    assert rows(bad)[0]['energy'] is None and rows(bad)[0]['symmetric_energy'] is None
    bad=table(edges,[*base,'--momentum','1','--max-evaluations','0'],2)
    assert rows(bad)[0]['lower_energy'] is None and rows(bad)[0]['upper_energy'] is None
    tiny=table(lines,['--u','4','--branch','holon','--momentum','1','--precision','bits-333','--tolerance','1e-80'])
    assert Decimal(tiny['metadata']['Relative tolerance'])==Decimal('1e-80') or abs(Decimal(tiny['metadata']['Relative tolerance'])/Decimal('1e-80')-1)<Decimal('1e-95')
    doped=[*base,'--density','0.324321270322984098050102851457799211925925445763299821413665672162111298405245789656083022','--branch','spinon','--momentum','0','--u','40']
    # Replace the unique interaction token rather than pass a repeated scalar option.
    doped=doped[:-2];doped[1]='40'
    t=table(lines,doped)
    assert rows(t)[0]['status']=='converged'
    assert abs(Decimal(t['metadata']['Fermi rapidity Q'])-1)<Decimal('1e-32')

with tempfile.TemporaryDirectory() as directory:
    path=Path(directory)/'preserve.json'
    for program in (lines,edges):
        for extra in [['--precision','bits-0'],['--tolerance','nan'],['--u','-1']]:
            path.write_text('sentinel')
            args=['--u','4','--momentum','0','--precision','bits-128']
            for flag,value in zip(extra[::2],extra[1::2]):
                if flag in args: args[args.index(flag)+1]=value
                else: args += [flag,value]
            run(program,[*args,'--json',str(path),'--force'],1)
            assert path.read_text()=='sentinel'
print('Hubbard MPFR precision and complete export contracts passed')
