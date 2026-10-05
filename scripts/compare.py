#!/usr/bin/env python3
"""Same fixed Jacobi sweeps: Python loops, NumPy slices, and student C++."""
import argparse
import math
from pathlib import Path
import statistics
import subprocess
import time
import numpy as np

p=argparse.ArgumentParser()
p.add_argument("--n",type=int,default=128)
p.add_argument("--sweeps",type=int,default=100)
p.add_argument("--repeats",type=int,default=3)
p.add_argument("--binary",type=Path,help="Optional path to a built C++ bench executable.")
args=p.parse_args()
n=args.n; sweeps=args.sweeps
if min(n,sweeps,args.repeats)<1: raise SystemExit("Arguments must be positive.")
h=1/(n+1); d=4/h**2; c=-1/h**2
xx,yy=np.meshgrid(np.arange(1,n+1)*h,np.arange(1,n+1)*h)
b=2*math.pi**2*np.sin(math.pi*xx)*np.sin(math.pi*yy)

def python_loops():
    rhs=b.ravel().tolist(); x=[0.]*(n*n); y=x.copy()
    start=time.perf_counter()
    for _ in range(sweeps):
        for j in range(n):
            for i in range(n):
                q=j*n+i; t=rhs[q]
                if i: t-=c*x[q-1]
                if i+1<n: t-=c*x[q+1]
                if j: t-=c*x[q-n]
                if j+1<n: t-=c*x[q+n]
                y[q]=t/d
        x,y=y,x
    return time.perf_counter()-start,sum(x)

def numpy_slices():
    x=np.zeros((n,n)); y=np.empty_like(x)
    start=time.perf_counter()
    for _ in range(sweeps):
        y[:]=b
        y[:,1:]-=c*x[:,:-1]; y[:,:-1]-=c*x[:,1:]
        y[1:,:]-=c*x[:-1,:]; y[:-1,:]-=c*x[1:,:]
        y/=d; x,y=y,x
    return time.perf_counter()-start,float(x.sum())

def cpp():
    root=Path(__file__).resolve().parents[1]
    output=subprocess.check_output([str(args.binary or root/"build/bench"),str(n),str(sweeps)],text=True)
    return tuple(map(float,output.split()))

answers={}
for name,fun in (("Python loops",python_loops),("NumPy slices",numpy_slices),("C++",cpp)):
    samples=[fun() for _ in range(args.repeats)]
    elapsed=statistics.median(t for t,_ in samples); checksum=samples[-1][1]
    answers[name]=(elapsed,checksum)
    print(f"{name:14s} median_seconds={elapsed:.8g} sum={checksum:.12g}")
for name,(_,checksum) in answers.items():
    if not math.isfinite(checksum) or not math.isclose(checksum,answers["NumPy slices"][1],rel_tol=1e-10,abs_tol=1e-12):
        raise SystemExit(f"{name}: checksum mismatch; verify the Jacobi implementation first.")
print(f"Python loops / C++: {answers['Python loops'][0]/answers['C++'][0]:.3f}x")
print(f"NumPy slices / C++: {answers['NumPy slices'][0]/answers['C++'][0]:.3f}x")
print("Fixed sweeps, zero initial guess, double precision; initial array setup and process startup excluded; temporaries inside sweeps included.")
