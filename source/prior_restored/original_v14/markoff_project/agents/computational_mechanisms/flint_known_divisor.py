#!/usr/bin/env python3
"""Call the stable FLINT C API for a *proved* determinant with known divisor.

The public Python binding does not expose this useful preconditioner.
Only documented exported C functions are used; no Python-FLINT object
internals or matrix-struct offsets are accessed.  Matrix storage is an
oversized, long-aligned opaque buffer, initialized/cleared by FLINT itself.
Integers and matrices are owned entirely by this one system library.

API and exactness contract:
  https://flintlib.org/download/flint-3.1.0.pdf
  fmpz_mat_det_modular_given_divisor(det,A,positive_divisor,proved=1).
The caller must prove that the supplied divisor divides det(A).  A minimal
denominator of an exact rational solve supplies such a divisor.
"""
import ctypes,ctypes.util,hashlib,json,math,os,pathlib,sys,time
sys.set_int_max_str_digits(0)
_name=os.environ.get('MARKOFF_FLINT_LIBRARY') or ctypes.util.find_library('flint')
if not _name:raise ImportError('A system FLINT shared library is required for this optional optimization')
_lib=ctypes.CDLL(_name)
P=ctypes.c_void_p;L=ctypes.c_long;I=ctypes.c_int
for name,args,result in [
 ('fmpz_init',[P],None),('fmpz_clear',[P],None),('fmpz_set_str',[P,ctypes.c_char_p,I],I),('fmpz_get_str',[P,I,P],P),
 ('fmpz_mat_init',[P,L,L],None),('fmpz_mat_clear',[P],None),('fmpz_mat_entry',[P,L,L],P),
 ('fmpz_mat_nrows',[P],L),('fmpz_mat_ncols',[P],L),('fmpz_mat_det_modular_given_divisor',[P,P,P,I],None)]:
 fn=getattr(_lib,name);fn.argtypes=args;fn.restype=result


def determinant_given_divisor(rows,divisor):
    n=len(rows);assert divisor>0 and all(len(r)==n for r in rows)
    # Current public struct contains a few pointers/lengths.  This buffer
    # deliberately exceeds its size and avoids assuming its field layout.
    matrix=(L*64)();den=L();det=L();_lib.fmpz_init(ctypes.byref(den));_lib.fmpz_init(ctypes.byref(det));_lib.fmpz_mat_init(ctypes.byref(matrix),n,n)
    try:
        assert _lib.fmpz_mat_nrows(ctypes.byref(matrix))==n and _lib.fmpz_mat_ncols(ctypes.byref(matrix))==n
        assert _lib.fmpz_set_str(ctypes.byref(den),str(divisor).encode(),10)==0
        maxbits=0
        for i,row in enumerate(rows):
            for j,x in enumerate(row):
                x=int(x);maxbits=max(maxbits,abs(x).bit_length());entry=_lib.fmpz_mat_entry(ctypes.byref(matrix),i,j)
                assert _lib.fmpz_set_str(entry,str(x).encode(),10)==0
        _lib.fmpz_mat_det_modular_given_divisor(ctypes.byref(det),ctypes.byref(matrix),ctypes.byref(den),1)
        # |det| <= n! maxentry^n; log10(2)<1/3 gives a rigorous, deliberately
        # loose decimal-buffer length without a floating-point size estimate.
        capacity=(n*maxbits+2)//3+n*len(str(max(1,n)))+8
        out=ctypes.create_string_buffer(max(32,capacity));_lib.fmpz_get_str(out,10,ctypes.byref(det));return int(out.value)
    finally:
        _lib.fmpz_mat_clear(ctypes.byref(matrix));_lib.fmpz_clear(ctypes.byref(det));_lib.fmpz_clear(ctypes.byref(den))


if __name__=='__main__':
    assert determinant_given_divisor([[2,3],[5,7]],1)==-1
    assert determinant_given_divisor([[2,0],[0,3]],6)==6
    root=pathlib.Path(__file__).resolve().parent
    for arg in sys.argv[1:]:
        d=int(arg);c=json.loads((root/f'three_block_d{d}_kernel.json').read_text());s=json.loads((root/c['source']).read_text());selected=[r[:3*d-2] for r,m in zip(s['rows'],s['row_metadata']) if m['n']<=3*d]
        restored=list(map(int,c['restored_content_factors']));A=[[x*f for x in row[:-1]] for row,f in zip(selected,restored)];w=list(map(int,c['primitive_kernel']));den=int(c['last_coordinate']);assert w[-1]==den and math.gcd(*w)==1
        assert all(sum(x*y for x,y in zip(row,w))==0 for row in selected)
        st=time.time();value=determinant_given_divisor(A,den);assert str(value)==c['first_square_minor_det']
        print(json.dumps({'d':d,'library':_name,'proved':1,'determinant_matches_existing':True,'seconds':time.time()-st,'determinant_digits':len(str(abs(value))),'given_divisor_digits':len(str(den))}),flush=True)
