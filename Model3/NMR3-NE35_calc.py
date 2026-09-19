import numpy as np
f0=21.7e6; w0=2*np.pi*f0
T=318.0; kB=1.380649e-23; q=1.602177e-19; fourkT=4*kB*T
print(f"f0={f0/1e6} MHz  T={T} K  4kT={fourkT:.4e}\n")

def NF(Rs,Ea,Ia): return 10*np.log10(1+(Ea**2+(Ia*Rs)**2)/(fourkT*Rs))

print("=== 0. How much amplifier noise can the Hoult architecture actually tolerate? ===")
print("   (source = Q*w0*L, its own thermal noise sets the bar)")
for Rs in [2900,5000,6600,26880]:
    eth=np.sqrt(fourkT*Rs)
    print(f"  Rs={Rs:6.0f} Ohm  thermal noise = {eth*1e9:5.2f} nV/sqrt(Hz)")
    for NFt in [0.5,1.0,2.0]:
        F=10**(NFt/10); need=np.sqrt((F-1)*fourkT*Rs)
        print(f"      for NF <= {NFt:.1f} dB the amplifier may contribute up to {need*1e9:5.2f} nV/sqrt(Hz) total")

print("\n=== 1. GAT1 back-fitted from Hoult Fig.5(a) (5 MHz, eyeballed points) ===")
# points read off: NF 3.2 dB at 1k, 1.2 dB at 5k, 1.35 dB at 10k, 1.8 dB at 20k
pts=[(1000,3.2),(3000,1.5),(5000,1.2),(10000,1.35),(20000,1.8)]
from scipy.optimize import least_squares
def resid(p):
    Ea,Ia=np.exp(p)
    return [NF(R,Ea,Ia)-v for R,v in pts]
sol=least_squares(resid,[np.log(4e-9),np.log(1e-12)])
Ea_g,Ia_g=np.exp(sol.x)
print(f"  fit: E_a = {Ea_g*1e9:.2f} nV/sqrt(Hz), I_a = {Ia_g*1e12:.2f} pA/sqrt(Hz), R0 = {Ea_g/Ia_g/1000:.1f} kOhm")
print(f"       F_min = {10*np.log10(1+2*Ea_g*Ia_g/fourkT):.2f} dB")
print(f"       equivalent gate leakage if I_a were pure shot noise: {(Ia_g**2/(2*q))*1e6:.2f} uA")

print("\n=== 2. Device comparison at 21.7 MHz, T=318 K ===")
gmJ=0.017; EaJ=np.sqrt(fourkT*(2/3)/gmJ)
IaJ=w0*4e-12*np.sqrt(fourkT*(4/3)/(5*gmJ))
gmN=0.080
print(f"  J309   : gm=17 mS  -> E_a(white) = {EaJ*1e9:.2f} nV/sqrt(Hz); I_a(induced gate) = {IaJ*1e12:.3f} pA/sqrt(Hz)")
for gamma,lab in [(1.0,'gamma=1'),(2.0,'gamma=2 (hot electron)')]:
    print(f"  NE3509 : gm=80 mS  -> E_a(white,{lab}) = {np.sqrt(fourkT*gamma/gmN)*1e9:.2f} nV/sqrt(Hz)")
EaN=np.sqrt(fourkT*1.5/gmN)
print(f"  NE3509 : take E_a(white) = {EaN*1e9:.2f} nV/sqrt(Hz) (gamma=1.5)")
print("  NE3509 1/f voltage noise at 21.7 MHz for assumed corner f_c:")
for fc in [50e6,200e6,1e9,5e9]:
    print(f"      f_c={fc/1e6:6.0f} MHz -> E_a(21.7 MHz) = {EaN*np.sqrt(1+fc/f0)*1e9:5.2f} nV/sqrt(Hz)")

print("\n=== 3. Gate leakage shot noise ===")
for Ig,lab in [(1e-11,'J309 typ  10 pA'),(1e-9,'J309 max   1 nA'),
               (5e-9,'NE3509 hypothetical 5 nA @ Vgs=-0.26 V'),
               (5e-8,'NE3509 hypothetical 50 nA'),
               (5e-7,'NE3509 IGSO TYP  0.5 uA @ Vgs=-3 V'),
               (1e-5,'NE3509 IGSO MAX   10 uA @ Vgs=-3 V')]:
    print(f"  I_G = {lab:42s} -> i_shot = {np.sqrt(2*q*Ig)*1e12:8.3f} pA/sqrt(Hz)")

print("\n=== 4. NF(Rs) at 21.7 MHz ===")
cases=[("J309  (Ea=0.83nV, Ia=0.286pA est.)",EaJ,IaJ),
       ("NE3509, I_G=5 nA   (Ea=0.51nV)",EaN,np.sqrt(2*q*5e-9)),
       ("NE3509, I_G=50 nA",EaN,np.sqrt(2*q*5e-8)),
       ("NE3509, I_G=0.5 uA (IGSO typ)",EaN,np.sqrt(2*q*5e-7)),
       ("NE3509, I_G=10 uA  (IGSO max)",EaN,np.sqrt(2*q*1e-5)),
       ("GAT1 back-fit (Hoult, 5 MHz)",Ea_g,Ia_g)]
hdr=f"{'case':38s}" + "".join(f"{r/1000:>9.1f}k" for r in [2900,6600,15000,26880])
print(hdr); print("-"*len(hdr))
for lab,Ea,Ia in cases:
    print(f"{lab:38s}" + "".join(f"{NF(r,Ea,Ia):10.2f}" for r in [2900,6600,15000,26880]))
print(f"{'MAR-6SM+ (measured, 50 ohm system)':38s}" + f"{'  ---':>10}"*2 + f"{'2.27 dB':>10}"*0 + "   [2.27 dB at 50 ohm]")

print("\n=== 5. Headroom between protection clamp and destruction ===")
for lab,Vbr,Vop in [("GAT1 (Hoult)",12.0,None),("J309",25.0,-1.2),("NE3509M04",3.0,-0.26)]:
    print(f"  {lab:14s} BV = {Vbr:5.1f} V ; clamp at 0.7 V -> margin {Vbr/0.7:5.1f}x = {20*np.log10(Vbr/0.7):5.1f} dB"
          + (f" ; operating V_GS = {Vop:+.2f} V -> {abs(Vop):.2f} V to forward conduction" if Vop else ""))

print("\n=== 6. NE3509M04 operating point ===")
for Idss,Voff in [(30,-0.25),(45,-0.5),(60,-0.75)]:
    Vgs=Voff*(1-np.sqrt(10/Idss))
    print(f"  IDSS={Idss} mA, VGS(off)={Voff} V -> V_GS at I_D=10 mA = {Vgs:+.3f} V")
print("  J309 for comparison: V_GS at I_D = 10 mA spans -0.087 .. -1.69 V")

print("\n=== 7. Magnetic worst case, M04 package vs SOT-23 ===")
Bs=1.3; mu0=4*np.pi*1e-7; B0=0.51
for mass_mg,lab in [(2.5,'M04 flat-lead, ~2.5 mg lead metal'),(25,'SOT-23, ~25 mg')]:
    V=mass_mg*1e-6/8100.0; m=Bs*V/mu0
    row=f"  {lab:36s}"
    for r in [0.02,0.03,0.05,0.10]:
        B=mu0*m/(2*np.pi*r**3)
        row+=f"  r={r*1000:3.0f}mm: grad {3*0.006/r*B/B0*1e6:7.3f} ppm"
    print(row)

print("\n=== 8. Dissipation ===")
print(f"  NE3509M04 at VDS=2 V, ID=10 mA -> {2*0.010*1000:.0f} mW  ({100*2*0.010/0.150:.0f} % of Ptot=150 mW)")
print(f"  J309      at VDS=5 V, ID=10 mA -> {5*0.010*1000:.0f} mW")
print(f"  MAR-6SM+                        -> 56 mW")
