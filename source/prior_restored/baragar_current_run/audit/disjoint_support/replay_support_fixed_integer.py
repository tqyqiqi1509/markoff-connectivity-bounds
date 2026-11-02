#!/usr/bin/env python3
"""Separate all-support certificate replay using only standard-library integers.
No imports from the producer. Reads all node bytes, rebuilds all integral
corrections and runs all allowed-prime sequences independently from scratch.
Arb is not required here. Analytic node enclosures remain explicitly supplied
by the producer and its accepted predecessor, whose bytes/hashes are checked.
"""
import argparse,gzip,hashlib,json,struct,sys
from pathlib import Path
sys.set_int_max_str_digits(0)
BITS=160;S=1<<BITS

def cd(n,d):return -((-n)//d)
def mul(a,b):return cd(a*b,S)
def dyad(d):
    m,e=int(d['mantissa']),int(d['exponent'])
    return m*S*(1<<e) if e>=0 else cd(m*S,1<<(-e))
def rational(d):return cd(int(d['numerator'])*S,int(d['denominator']))
def read_nodes(path,codec,panels,blocks,key):
    count=len(codec.unpack(bytes(codec.size)));sums=[0]*count;bs=[0]*count;bi=0;h=hashlib.sha256()
    with gzip.open(path,'rb') as f:
        for i in range(panels):
            b=f.read(codec.size);assert len(b)==codec.size;h.update(b)
            row=codec.unpack(b)
            for j,v in enumerate(row):sums[j]+=v;bs[j]+=v
            if i==blocks[bi]['last_node']:
                assert bs==blocks[bi][key];bs=[0]*count;bi+=1
        assert f.read(1)==b''
    assert bi==len(blocks)
    return sums,h.hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',default='baragar_run/audit/disjoint_support/support_169_32768.json');ap.add_argument('--output',default='baragar_run/audit/disjoint_support/replay_169_fixed_integer.json');a=ap.parse_args()
    path=Path(a.input);d=json.loads(path.read_text());info=d['integration'];src=Path(info['source_file']);old=json.loads(src.read_text())
    assert hashlib.sha256(src.read_bytes()).hexdigest()==info['source_sha256']
    params=old['parameters']['active_primes'];n=len(params);k=len(info['small_primes']);start=k+1;N=info['panels'];scale=1<<info['node_scale_bits']
    oldinfo=old['integration'];oldpath=src.parent/oldinfo['node_file'];newpath=path.parent/info['node_file']
    assert hashlib.sha256(oldpath.read_bytes()).hexdigest()==oldinfo['node_gzip_sha256']
    assert hashlib.sha256(newpath.read_bytes()).hexdigest()==info['gzip_sha256']
    oldsum,h=read_nodes(oldpath,struct.Struct('<'+'Q'*(n+1)),N,oldinfo['blocks'],'tail_sums')
    assert h==oldinfo['node_raw_sha256']==info['source_node_sha256']
    assert oldsum==[int(t['node_sum_integer']) for t in oldinfo['tails']]
    subsum,h=read_nodes(newpath,struct.Struct('<'+'Q'*(1<<k)),N,info['blocks'],'subset_sums')
    assert h==info['raw_sha256'] and subsum==[int(t['node_sum_integer']) for t in info['subset_integrals']]
    q=[p['q'] for p in params];w=[dyad(p['w_upper']) for p in params];r=[dyad(p['r_upper']) for p in params]
    K=[dyad(p['K_upper']) for p in params];L=[dyad(p['L_upper']) for p in params]
    Delta=cd(rational(info['pi_upper']),2*N)
    KK=[S]*(n+1);FF=[S]*(n+1);WW=[S]*(n+1)
    for j in range(n-1,-1,-1):
        KK[j]=mul(K[j],KK[j+1]);FF[j]=mul(S+mul(L[j],Delta),FF[j+1]);WW[j]=min(KK[j],cd(oldsum[j]*FF[j],N*scale))
    subW=[]
    for mask in range(1<<k):
        cap=KK[start];cor=FF[start]
        for j in range(1,start):
            if mask&(1<<(j-1)):
                cap=mul(cap,K[j]);cor=mul(cor,S+mul(L[j],Delta))
        subW.append(min(cap,cd(subsum[mask]*cor,N*scale)))
    allC=[];oneC=[];details=[]
    for mask in range(1<<k):
        seq=[0]+[j for j in range(1,start) if mask&(1<<(j-1))]+list(range(start,n))
        tails=[S]*(len(seq)+1);waves=[]
        for t,j in enumerate(seq):
            if j>=start:waves.append(WW[j])
            elif j>0:
                suffix=sum(1<<(h-1) for h in seq if j<=h<start);waves.append(subW[suffix])
            else:waves.append(0) # W starting with 2 is never queried.
        for t in range(len(seq)-1,-1,-1):tails[t]=mul(K[seq[t]],tails[t+1])
        C=S;one=None;maxrows=[]
        for t in range(len(seq)-1,-1,-1):
            j=seq[t];best=C;best_a=0;wa=ra=S;aexp=0
            while True:
                aexp+=1;wa=mul(wa,w[j]);ra=mul(ra,r[j])
                if mul(wa,tails[t+1]+C)<=C:break
                J=tails[t+1];E=S
                for u in range(t+1,len(seq)):
                    h=1
                    while q[seq[u]]**h<q[j]**aexp:h+=1
                    J=min(J,mul(E,h*waves[u]));E=mul(E,K[seq[u]])
                val=mul(wa,J+mul(ra,C))
                if j==0 and aexp==1:one=val
                if val>best:best,best_a=val,aexp
            assert one is not None or j!=0
            C=best;maxrows.append([q[j],best_a,aexp])
        allC.append(C);oneC.append(one);details.append({'mask':mask,'C_all_scaled':str(C),'C_one_scaled':str(one),'maxima_and_stops':maxrows})
    pairs=[allC[m]+oneC[((1<<k)-1)^m] for m in range(1<<k)]
    maximum=max(pairs);C=cd(maximum,S);worst=pairs.index(maximum)
    assert C==d['recursion']['joint_integer_upper']
    P=10**d['endpoint']['power'];T=int(d['endpoint']['T']);KB=181265
    lhs=(216*(2*KB+3))**900*27**420*C**6660*(P+1)**673;rhs=(2*P)**900
    checks={'same_integer_constant':True,'main_exact':lhs<rhs,'root_exact':(T-1)**44<(P+1)**5<=T**44,
        'safe8':36864*2541933**4*2541935**2*T**6<P,'cage':256*2541933**4*T**4<P}
    assert all(checks.values())
    result={'status':'PASSED','input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'fixed_bits':BITS,
            'all_old_values_read':N*(n+1),'all_new_values_read':N*(1<<k),'all_support_sequences':1<<k,
            'joint_integer_upper':C,'worst_high_v2_mask':worst,'joint_upper_scaled':str(maximum),
            'checks':checks,'details':details,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='details'},indent=2))
if __name__=='__main__':main()
