"""Independent real-axis Dirac vacuum integral; no shifted NLIE or Bethe roots."""
import json,sys
import mpmath as mp
mp.mp.dps = int(sys.argv[1]) if len(sys.argv)>1 else 110
rows=[]
for text in ['0.25','1','5']:
 u=mp.mpf(text)
 cutoff=mp.acosh((mp.mp.dps+30)*mp.log(10)/u)
 f=lambda x:mp.cosh(x)*mp.log1p(mp.exp(-u*mp.cosh(x)))
 y=-2*u/mp.pi*mp.quad(f,[0,cutoff/4,cutoff/2,3*cutoff/4,cutoff])
 # log(1+exp(-u cosh x)) < exp(-u cosh x); cosh x <= exp x.
 # For x>=B, cosh x>=exp x/2, giving this analytic omitted tail bound.
 tail=4/mp.pi*mp.exp(-u*mp.exp(cutoff)/2)
 rows.append(dict(u=text,y=mp.nstr(y,mp.mp.dps),tail_bound=mp.nstr(tail,mp.mp.dps)))
print(json.dumps(dict(digits=mp.mp.dps,fixtures=rows),indent=2))
