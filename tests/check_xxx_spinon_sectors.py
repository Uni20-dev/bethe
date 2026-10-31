"""Fresh selected-sector spectra, physical sums, failures and export contracts."""
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, mpfr = sys.argv[1:]
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')

def invoke(args, expected=0):
    p = subprocess.run([program, *map(str, args)], capture_output=True, text=True, env=env, timeout=180)
    assert p.returncode == expected, (args, p.returncode, p.stdout[-4000:], p.stderr)
    return p.stdout

def changed(args, overrides):
    names=set(overrides[::2])
    result=[]
    i=0
    while i<len(args):
        if args[i] in names:
            i+=2
        else:
            result.append(args[i]); i+=1
    return [*result,*overrides]

def document(args, expected=0):
    return json.loads(invoke([*args, '--format', 'json'], expected))

def rows(doc, table='spectrum'):
    t = doc['tables'][table]
    names = [c['id'] for c in t['columns']]
    return [dict(zip(names, [cell["value"] if isinstance(cell,dict) else cell for cell in row])) for row in t["rows"]]

def close(a, b, relative='2e-14'):
    with localcontext() as ctx:
        ctx.prec = 130
        assert abs(Decimal(a)-Decimal(b)) <= Decimal(relative)*max(abs(Decimal(a)), abs(Decimal(b)), Decimal('1e-100')), (a,b)

base = ['--sectors', '2,4', '--momentum', '1.5', '--omega', '1.8', '--precision', 'long-double',
        '--four-order', '2', '--four-cutoff', '2','--diagnostics']
doc = document(base)
r = {row['sector']: row for row in rows(doc)}
assert set(r) == {'2','4','sum:2+4'}
assert r['4']['status'] == 'estimated'
assert Decimal(r['4']['density']) > 0
close(r['sum:2+4']['density'], Decimal(r['2']['density'])+Decimal(r['4']['density']))
assert set(doc['tables']) == {'spectrum','continuum','diagnostics'}
diag = {row['sector']: row for row in rows(doc,'diagnostics')}
assert diag['4']['quadrature_difference'] is None
assert diag['4']['sampling_error'] is None
assert Decimal(diag['4']['contour_diagnostic']) >= 0
assert Decimal(diag['4']['arithmetic_diagnostic']) >= 0
parallel = document([*base, '--threads', '2', '--no-retain'])
for name in doc['tables']:
    assert doc['tables'][name]['rows'] == parallel['tables'][name]['rows']
raised = document([*base,'--channel','raising'])
for a,b in zip(rows(doc), rows(raised)):
    close(Decimal(a['density'])*2,b['density'])
scaled = document(changed(base,['--exchange','2','--omega','3.6']))
for a,b in zip(rows(doc),rows(scaled)):
    close(a['density'],Decimal(b['density'])*2)
failed = document([*base,'--four-max-nodes','1'],2)
f = {row['sector']:row for row in rows(failed)}
assert f['2']['density'] is not None
assert f['4']['density'] is None and f['4']['status']=='numerical_failure'
assert f['sum:2+4']['density'] is None and f['sum:2+4']['status']=='incomplete'
refined = document([*base,'--four-refine-order','3','--four-extended-cutoff','3'])
d = {row['sector']:row for row in rows(refined,'diagnostics')}['4']
assert d['quadrature_difference'] is not None and d['cutoff_change'] is not None

with tempfile.TemporaryDirectory() as tmp:
    path=Path(tmp)/'output.json'
    invoke([*base,'--json',path,'--quiet'])
    original=path.read_bytes()
    for bad in [['--sectors','4,4'],['--sectors','8'],['--sectors','2,'],['--four-order','1'],
                ['--four-cutoff','nan'],['--four-refine-order','1','--four-extended-cutoff','3'],
                ['--four-refine-order','4'],['--table','replicates'],['--table','failed'],
                ['--density-seconds','0'],['--density-terms','1'],
                ['--exchange','1e-4900','--omega','1e4900'],
                ['--exchange','1e4900','--omega','1e-4900']]:
        invoke([*changed(base,bad),'--json',path,'--force'],1)
        assert path.read_bytes()==original
    diagnostics=Path(tmp)/'diagnostics.csv'
    invoke([*base,'--quiet','--csv-table',f'diagnostics={diagnostics}'])
    assert 'contour_diagnostic' in diagnostics.read_text()

# Physical zero in all sectors needs no expensive density evaluations.
zero=document(['--sectors','2,4,6','--momentum','0','--omega','2','--precision','long-double','--replicate-output'])
assert all(Decimal(row['density'])==0 for row in rows(zero))
assert rows(zero,'replicates') == []

if mpfr.upper() in ('ON','TRUE','1'):
    six=['--sectors','2,4,6','--momentum','1.57079632679489661923132169163975','--omega','2',
         '--precision','bits-384','--four-order','2','--four-cutoff','2','--replicates','3','--samples','8','--seed','42',
         '--density-relative-tolerance','1e-5','--threads','2','--diagnostics','--replicate-output']
    doc=document(six)
    r={row['sector']:row for row in rows(doc)}
    assert all(r[str(n)]['density'] is not None for n in (2,4,6))
    assert Decimal(r['6']['density'])>0
    close(r['sum:2+4+6']['density'],sum(Decimal(r[str(n)]['density']) for n in (2,4,6)))
    reps=rows(doc,'replicates'); assert len(reps)==3
    assert all(row['failed_sample'] is None for row in reps)
    assert [int(row['replicate']) for row in reps]==[0,1,2]
    d={row['sector']:row for row in rows(doc,'diagnostics')}
    assert d['6']['sampling_error'] is not None
    assert d['6']['quadrature_difference'] is None
print('Selected spinon spectra, partial sums and export contracts passed')
