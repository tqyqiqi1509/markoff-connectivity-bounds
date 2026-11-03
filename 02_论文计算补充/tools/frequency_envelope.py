"""Generate exactly the 15104 full-distribution envelope cells on [0,32]."""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing,json,hashlib,math,time
from flint import arb,acb,ctx
ctx.prec=192
ETA=F(1,18);EPS=F(63,1000);DELTA=F(829,9180);PRIME_CUTOFF=4099
def aa(x):
    x=F(x)
    return arb(x.numerator)/x.denominator

def pack(x):
    assert x.is_finite()
    return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display_only':str(x)}

def primes_below(stop):
    sieve=bytearray(b'\x01')*stop; sieve[:2]=b'\x00\x00'
    for q in range(2,math.isqrt(stop-1)+1):
        if sieve[q]:
            for m in range(q*q,stop,q): sieve[m]=0
    return [q for q in range(2,stop) if sieve[q]]

def local_data():
    data=[]; rows=[]; count=0
    for q in primes_below(PRIME_CUTOFF):
        logq=arb(q).log(); r=(-aa(ETA)*logq).exp(); cq=arb(q-1)/q
        epsw=(-aa(EPS)*logq).exp(); infinite=1+cq*r/(1-r)
        apow=arb(1); rpow=arb(1); sumjr=arb(0); weights=[]; a=0
        while True:
            a+=1; apow*=epsw; rpow*=r; sumjr+=a*rpow
            if apow*infinite<1: break
            weights.append((apow,(apow*cq*logq*sumjr).upper()))
            count+=1; assert a<1000
        rows.append({'q':q,'first_omitted_exponent':a,'tail_guard':pack(apow*infinite)})
        data.append((q,logq,r,cq,weights))
    lq=arb(PRIME_CUTOFF).log()
    guard=(-aa(EPS)*lq).exp()*(1+(-aa(ETA)*lq).exp())
    assert guard<1
    assert len(data)==564 and count==1355
    return data,rows,guard,count

def frequency_chunk(job):
    offset,denominator,start,length=job
    midpoint0=aa(F(2*start+1,2*denominator)); radius=aa(F(1,2*denominator))
    step=aa(F(1,denominator)); products=[arb(1) for _ in range(length)]
    # Exact-grid recurrence removes repeated trigonometric calls only.
    # Arb propagation encloses every exact r*exp(i*midpoint*log(q)).
    for q,logq,r,cq,weights in _LOCAL:
        sine,cosine=(logq*midpoint0).sin_cos(); unit=acb(cosine,sine)*r
        sine,cosine=(logq*step).sin_cos(); rotation=acb(cosine,sine)
        for k in range(length):
            power=acb(1); partial=acb(1); best=arb(1)
            for weight,derivative in weights:
                power*=unit; partial+=cq*power
                candidate=(weight*partial.abs_upper()+derivative*radius).upper()
                best=best.max(candidate).upper()
            products[k]*=best
            if k+1<length: unit*=rotation
    result=[]
    for product in products:
        value=(product*_DENOMINATOR_INVERSE).upper().min(arb(1)).upper()
        assert 0<value<=1
        result.append(str(value.fmpq()))
    return offset,result
def worker_init():
 global _LOCAL,_DENOMINATOR_INVERSE
 ctx.prec=192;_LOCAL=local_data()[0]
 _DENOMINATOR_INVERSE=(-arb(14000).log()-aa(DELTA-EPS)*aa(F(10**153,2)).log()).exp()
def main():
 start=time.monotonic();local,tails,guard,count=local_data();jobs=[];offset=0
 for den,first,last in [(16384,0,4096),(1024,256,8192),(128,1024,4096)]:
  for i in range(first,last,256):jobs.append((offset+i-first,den,i,min(256,last-i)))
  offset+=last-first
 assert offset==15104;values=[None]*offset;done=0
 with ProcessPoolExecutor(max_workers=4,initializer=worker_init,mp_context=multiprocessing.get_context('fork')) as pool:
  for future in as_completed([pool.submit(frequency_chunk,j) for j in jobs]):
   i,v=future.result();values[i:i+len(v)]=v;done+=len(v)
   if done%2048==0:print('FREQUENCY_CELLS',done,flush=True)
 assert all(x is not None for x in values)
 out={'status':'PASS_FULL_FREQUENCY_ENVELOPE','frequency_upper_values':values,'frequency_cells':offset,'prime_count':len(local),'nonzero_exponents':count,'local_tails':tails,'prime_tail_guard':pack(guard),'grid_segments':[[0,'1/4',16384,4096],['1/4',8,1024,7936],[8,32,128,3072]],'resource_threshold':'14000*(10^153/2)^(12533/459000)','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start}
 (Path(__file__).resolve().parent.parent/'certificates/FREQUENCY_ENVELOPE.json').write_text(json.dumps(out,indent=2));print('PASS_FREQUENCY_ENVELOPE',offset,flush=True)
if __name__=='__main__':main()
