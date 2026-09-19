import numpy as np
rng=np.random.default_rng(7); N=200000; R=5.0
p=rng.normal(size=(N,3)); p/=np.linalg.norm(p,axis=1)[:,None]
p*=R*rng.random(N)[:,None]**(1/3); xs,ys,zs=p.T

H={'X':lambda x,y,z:x,'Y':lambda x,y,z:y,'Z':lambda x,y,z:z,
 'Z2':lambda x,y,z:2*z**2-x**2-y**2,'X2-Y2':lambda x,y,z:x**2-y**2,'XY':lambda x,y,z:x*y,
 'ZX':lambda x,y,z:z*x,'ZY':lambda x,y,z:z*y,
 'Z3':lambda x,y,z:z*(2*z**2-3*x**2-3*y**2),'X3m1':lambda x,y,z:x*(4*z**2-x**2-y**2),
 'Y3m1':lambda x,y,z:y*(4*z**2-x**2-y**2),'Z(X2-Y2)':lambda x,y,z:z*(x**2-y**2),
 'ZXY':lambda x,y,z:z*x*y,'X3':lambda x,y,z:x*(x**2-3*y**2),'Y3':lambda x,y,z:y*(3*x**2-y**2)}
names=list(H); norm={k:H[k](xs,ys,zs).std() for k in names}   # нормировка на сфере R=5

def sample(b,a=1.0,M=120000):
    t=2*np.pi*rng.random(M); rr=a*np.sqrt(rng.random(M))
    return b*(2*rng.random(M)-1), rr*np.cos(t), rr*np.sin(t)

for b in (5.0,3.0):
    xc,yc,zc=sample(b)
    A=np.column_stack([H[k](xc,yc,zc)/norm[k] for k in names])
    A-=A.mean(0)                       # константа линию не уширяет
    A/=np.sqrt(len(xc))
    U,S,Vt=np.linalg.svd(A,full_matrices=False)
    print(f"\n===== столб образца 2b = {2*b:.0f} мм =====")
    print(" вклад в ширину (σ по образцу, ед. σ на сфере R=5):")
    for k in names: print(f"   {k:9s} {np.std(H[k](xc,yc,zc)/norm[k]):6.2f}")
    print(" сингулярные числа (какие комбинации вообще различимы):")
    print("   "+"  ".join(f"{s/S[0]:.3f}" for s in S))
    n=(S/S[0]>0.05).sum(); print(f"   значимых направлений (>5% от старшего): {n} из {len(names)}")
    for i in range(min(6,len(S))):
        w=Vt[i]; idx=np.argsort(-abs(w))[:3]
        comp=" + ".join(f"{w[j]:+.2f}·{names[j]}" for j in idx)
        print(f"   мода {i+1}  s={S[i]/S[0]:.3f}   {comp}")
