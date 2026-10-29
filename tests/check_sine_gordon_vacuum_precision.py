"""Sine-Gordon vacuum explicit precision, complete exports and honest failure."""
import csv
from decimal import Decimal,localcontext
import json,os
from pathlib import Path
import subprocess,sys,tempfile
program,available=sys.argv[1:]
enabled=available.upper() in ('ON','TRUE','1')
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',UNI20_COLOR='never')

def run(args,code=0):
 p=subprocess.run([program,*args],capture_output=True,text=True,env=env,timeout=120)
 assert p.returncode==code,(args,p.returncode,p.stdout,p.stderr)
 return p.stdout

def result(args,bits,code=0):
 doc=json.loads(run([*args,'--format','json'],code))
 assert doc['status']=='complete'
 t=doc['tables']['vacuum']
 assert int(t['metadata']['Working precision (bits)'])==bits
 assert t['summary']['Outcome']==('success' if code==0 else 'partial')
 for row in t['rows']:
  for col,cell in zip(t['columns'],row):
   if col['type']=='real' and cell is not None:
    assert isinstance(cell,dict),(col,cell)
    assert cell['precision']==bits,(col,cell)
 return t

def value(c):return c['value'] if isinstance(c,dict) else c
help=run(['--help']); assert 'bits-N' in help and 'digits-N' in help
base=['--length','1','--p','1']
if enabled:
 reference=Decimal('-0.34560421614102583595905220720551356251474577009985950635758577036317931727427046595187736252422296779127574676')
 previous=Decimal(1)
 for bits in (80,160,256,333):
  args=[*base,'--precision',f'bits-{bits}','--max-intervals','8192']
  t=result(args,bits)
  row=dict(zip([c['id'] for c in t['columns']],t['rows'][0]))
  with localcontext() as c:
   c.prec=120
   error=abs(Decimal(value(row['scaling_function']))-reference)
   assert error<previous*Decimal('1e-15'),(bits,error,previous)
   previous=error
  assert row['converged'] is True and row['status']=='converged'
  failed=result([*base,'--precision',f'bits-{bits}','--max-cutoffs','0'],bits,2)
  assert failed['rows'][0][:3]==[None,None,None]
 for bits in (80,160):
  args=[*base,'--precision',f'bits-{bits}']
  direct=result(args,bits)
  with tempfile.TemporaryDirectory() as folder:
   root=Path(folder)
   for sep,flag,suffix in ((',', '--csv-table','csv'),('\t','--tsv-table','tsv')):
    path=root/('vacuum.'+suffix)
    t=result([*args,flag,f'vacuum={path}'],bits)
    lines=path.read_text().splitlines(); assert any(line.startswith('#') for line in lines)
    rows=list(csv.reader((l for l in lines if not l.startswith('#')),delimiter=sep))
    assert rows[0]==[c['id'] for c in t['columns']]
    expected=[[str(value(c)).lower() if isinstance(value(c),bool) else '' if c is None else str(value(c)) for c in row] for row in t['rows']]
    assert rows[1:]==expected
   file=root/'vacuum.json'
   run([*args,'--json',str(file),'--no-retain','--quiet'])
   assert json.loads(file.read_text())['tables']['vacuum']['rows']==direct['rows']
 parsed=result(['--mass','1.23456789012345678901234567890123456789','--length','1','--p','2',
   '--precision','digits-100','--tolerance','1e-80','--max-kernel-evaluations','0'],333,2)
 with localcontext() as c:
  c.prec=110
  assert abs(Decimal(parsed['metadata']['Soliton mass'])-Decimal('1.23456789012345678901234567890123456789'))<Decimal('1e-95')
  assert abs(Decimal(parsed['metadata']['Absolute Y tolerance per state'])/Decimal('1e-80')-1)<Decimal('1e-95')
else:
 assert not run([*base,'--precision','bits-160'],1)
with tempfile.TemporaryDirectory() as folder:
 path=Path(folder)/'keep.json'
 for args in (['--length','-1','--p','1'],['--length','1','--p','0'],[*base,'--tolerance','nan'],[*base,'--max-intervals','1'],[*base,'--contour-shift','2']):
  for precision in (['fp64','bits-160'] if enabled else ['fp64']):
   path.write_text('keep\n')
   run([*args,'--precision',precision,'--json',str(path),'--force'],1)
   assert path.read_text()=='keep\n'
print('Sine-Gordon vacuum MPFR precision, accuracy, exports and failures passed')
