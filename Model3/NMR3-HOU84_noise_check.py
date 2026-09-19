import numpy as np

def analyze(f0_tank, Q, wL, Z0, Rp_list, f_scan=None):
    # coil: L with reactance wL at w0(tank); r = wL/Q
    w0 = 2*np.pi*f0_tank
    L  = wL/w0
    r  = wL/Q
    C  = 1.0/(w0**2*L)
    Cp = np.sqrt(C/(Q*w0*Z0))          # Eq.[6]
    return L, r, C, Cp, w0

f0 = 60e6; Q=100; wL=50.0; Z0=50.0
L,r,C,Cp,w0 = analyze(f0,Q,wL,Z0,None)
print(f"L={L*1e9:.1f} nH  r={r:.3f} Ohm  C={C*1e12:.2f} pF  C'={Cp*1e12:.3f} pF  1/(w0C')={1/(w0*Cp):.1f} Ohm")

def node_Y(w, Rp):
    """Return: Zsrc seen by R' at node Y (with R' removed), H_coil = VY/e_s, H_R = VY/e_R'."""
    ZL = r + 1j*w*L
    ZC = 1.0/(1j*w*C)
    ZCp= 1.0/(1j*w*Cp)
    # Thevenin looking from Y back (R' removed): tank impedance in series with C'
    Ztank = ZL*ZC/(ZL+ZC)
    Zsrc  = Ztank + ZCp
    # signal transfer: e_s in series with coil branch
    # Node X: current divider. Easier: solve nodal.
    # Unknowns Vx, Vy.  e_s drives through (r + jwL) from ground.
    # (Vx - e_s)/ZL + Vx/ZC + (Vx-Vy)/ZCp = 0
    # (Vy-Vx)/ZCp + Vy/Rp = 0
    A = np.array([[1/ZL+1/ZC+1/ZCp, -1/ZCp],
                  [-1/ZCp, 1/ZCp+1/Rp]])
    b = np.array([1/ZL, 0.0])           # e_s = 1
    V = np.linalg.solve(A,b)
    H_coil = V[1]
    # noise of R': EMF in series with R'
    b2 = np.array([0.0, 1/Rp])
    V2 = np.linalg.solve(A,b2)
    H_R = V2[1]
    return Zsrc, H_coil, H_R

# find operating frequency where input impedance at Y is real, for R'=50
from scipy.optimize import brentq
def imagZin(f, Rp):
    w=2*np.pi*f
    ZL=r+1j*w*L; ZC=1/(1j*w*C); ZCp=1/(1j*w*Cp)
    Ztank=ZL*ZC/(ZL+ZC)
    return (Ztank+ZCp).imag
fm = brentq(lambda f: imagZin(f,50), 0.80*f0, 0.999*f0)
wm = 2*np.pi*fm
print(f"matching freq f_m = {fm/1e6:.3f} MHz  (f0={f0/1e6} MHz)  Zsrc = {(node_Y(wm,50)[0]):.2f}")

print("\nR'      Zsrc(f_m)          SNRpenalty_dB   (=10log10(1+ |H_R|^2 R' / (|H_c|^2 r)))")
for Rp in [50,100,200,500,1000,2000,5000]:
    Zs,Hc,HR = node_Y(wm,Rp)
    deg = 1 + (abs(HR)**2*Rp)/(abs(Hc)**2*r)
    print(f"{Rp:6.0f}  {Zs.real:8.2f}{Zs.imag:+8.2f}j   {10*np.log10(deg):8.3f}   simple 1+Rs/R': {10*np.log10(1+50/Rp):.3f}")
