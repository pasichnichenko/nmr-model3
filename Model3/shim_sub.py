import numpy as np
rng=np.random.default_rng(3); N=200000; R=5.0
p=rng.normal(size=(N,3)); p/=np.linalg.norm(p,axis=1)[:,None]; p*=R*rng.random(N)[:,None]**(1/3)
xs,ys,zs=p.T
H={'X':lambda x,y,z:x,'Y':lambda x,y,z:y,'Z':lambda x,y,z:z,
 'Z2':lambda x,y,z:2*z**2-x**2-y**2,'X2-Y2':lambda x,y,z:x**2-y**2,'XY':lambda x,y,z:x*y,
 'ZX':lambda x,y,z:z*x,'X3':lambda x,y,z:x*(x**2-3*y**2),'X3m1':lambda x,y,z:x*(4*z**2-x**2-y**2)}
norm={k:H[k](xs,ys,zs).std() for k in H}
M=200000; b,a=5.0,1.0
t=2*np.pi*rng.random(M); rr=a*np.sqrt(rng.random(M))
xc=b*(2*rng.random(M)-1); yc=rr*np.cos(t); zc=rr*np.sin(t)
def v(k):
    u=H[k](xc,yc,zc)/norm[k]; return u-u.mean()
def resid(target,tools):
    y=v(target); A=np.column_stack([v(t_) for t_ in tools])
    c,*_=np.linalg.lstsq(A,y,rcond=None)
    return (y-A@c).std()/y.std()
print("остаток после гашения одной гармоники другой (доля от исходного вклада):")
for tgt,tl in [('X2-Y2',['Z2']),('Z2',['X2-Y2']),('X3',['X']),('X3',['X','Z2']),
               ('X3m1',['X3']),('XY',['ZX']),('X3',['X','X2-Y2']),
               ('Z2',['X']),('X2-Y2',['X','Z2'])]:
    print(f"   {tgt:6s} гасим {str(tl):22s} -> остаток {resid(tgt,tl)*100:5.1f} %")
