import numpy as np
f0=21.7e6; w0=2*np.pi*f0
T=318.0; k=1.380649e-23; q=1.602177e-19; fourkT=4*k*T
print(f"f0={f0/1e6} MHz  T={T} K (+45 C)  4kT={fourkT:.4e} V^2/(Hz*Ohm)\n")

print("=== A. E_a (thermal channel floor, gamma=2/3) ===")
for gm in [0.010,0.017,0.024,0.030]:
    Ea=np.sqrt(fourkT*(2/3)/gm)
    print(f"  gm={gm*1000:5.1f} mS  ->  E_a = {Ea*1e9:5.3f} nV/sqrt(Hz)")

print("\n=== B. I_a contributions at 21.7 MHz ===")
gm=0.017; Cgs=4e-12; delta=4/3
i_ind=w0*Cgs*np.sqrt(fourkT*delta/(5*gm))
print(f"  induced gate noise (van der Ziel, Cgs=4pF, gm=17mS): {i_ind*1e12:.3f} pA/sqrt(Hz)")
for Ig,lab in [(1e-11,'10 pA (datasheet TYP text)'),(1e-9,'1 nA (datasheet MAX)')]:
    print(f"  shot noise, I_G={lab:28s}: {np.sqrt(2*q*Ig)*1e12:.4f} pA/sqrt(Hz)")
Ia=np.sqrt(i_ind**2+(np.sqrt(2*q*1e-11))**2)
Ea=np.sqrt(fourkT*(2/3)/gm)
R0=Ea/Ia
print(f"  => I_a(total) = {Ia*1e12:.3f} pA/sqrt(Hz);  E_a = {Ea*1e9:.3f} nV/sqrt(Hz)")
print(f"  => R0 = E_a/I_a = {R0/1000:.2f} kOhm")
Fmin=1+2*Ea*Ia/fourkT
print(f"  => F_min = {Fmin:.4f} = {10*np.log10(Fmin):.3f} dB")

def NF(Rs,Ea,Ia):
    F=1+(Ea**2+ (Ia*Rs)**2)/(fourkT*Rs); return 10*np.log10(F)

print("\n=== C. NF vs source impedance (= Q*w0*L when gate sits directly on the probe) ===")
print(f"{'Q*w0*L':>10}{'S=Rs/R0':>9}{'NF nominal':>12}{'NF if I_a x3':>14}{'NF if I_a x10':>15}")
for Rs in [2000,2893,6600,10000,15000,26880,40000]:
    print(f"{Rs:10.0f}{Rs/R0:9.2f}{NF(Rs,Ea,Ia):12.3f}{NF(Rs,Ea,3*Ia):14.3f}{NF(Rs,Ea,10*Ia):15.3f}")
print("  (MAR-6SM+ at +45 C, 20 MHz, 16 mA:  NF = 2.27 dB  [NMR3-MAR6-006])")

print("\n=== D. Required open-loop gain  |a-1| ===")
print("  direct gate connection:  |a-1| = D*C/(Q*C')")
for wL,N in [(55,'N~10'),(224,'N~20')]:
    L=wL/w0; C=1/(w0**2*L); Q=120; D=10
    for Cp in [1e-12,2e-12,3e-12]:
        print(f"   {N}  wL={wL:3.0f}  C={C*1e12:5.1f} pF  C'={Cp*1e12:.0f} pF -> |a-1| = {D*C/(Q*Cp):6.2f}   (Q*w0*L={Q*wL/1000:.1f} k)")
print("\n  after lossless noise-match to R0:  |a-1| = D/(w0*R0*C')   -- independent of coil!")
for Cp in [0.5e-12,1e-12,2e-12,3e-12,5e-12]:
    print(f"   C'={Cp*1e12:4.1f} pF -> |a-1| = {10/(w0*R0*Cp):6.2f}")

print("\n=== E. Cgs as a fraction of probe tuning capacitance ===")
for wL,N in [(55,'N~10'),(224,'N~20')]:
    L=wL/w0; C=1/(w0**2*L)
    for Cin in [4e-12,5e-12]:
        print(f"   {N}: C={C*1e12:5.1f} pF, C_gs={Cin*1e12:.0f} pF -> {100*Cin/C:5.1f} % of C; retune shift {0.5*Cin/C*1e6:6.0f} ppm = {0.5*Cin/C*f0/1e3:6.0f} kHz")

print("\n=== F. Dissipation ===")
for Vds,Id in [(5,0.010),(10,0.010),(10,0.020),(10,0.030)]:
    print(f"   VDS={Vds:2.0f} V, ID={Id*1000:2.0f} mA -> {Vds*Id*1000:5.0f} mW   ({100*Vds*Id/0.360:.0f} % of PD_max=360 mW)")
print("   MAR-6SM+: 16 mA x 3.5 V = 56 mW")

print("\n=== G. Worst-case magnetic perturbation of a ferromagnetic lead frame ===")
Bs=1.3; mu0=4*np.pi*1e-7; B0=0.51
for mass_mg,lab in [(25,'~25 mg lead frame (SOT-23, alloy 42)'),(60,'~60 mg (TO-92, alloy 42)')]:
    V=mass_mg*1e-6/8100.0   # m^3, rho=8.1 g/cm3
    m=Bs*V/mu0
    print(f"  {lab}: V={V*1e9:.2f} mm^3, m={m*1e3:.2f} mA*m^2")
    for r in [0.03,0.05,0.10,0.15,0.20]:
        B=mu0*m/(2*np.pi*r**3)
        print(f"     r={r*1000:4.0f} mm -> dB={B*1e6:8.2f} uT = {B/B0*1e6:8.2f} ppm of 0.51 T ; gradient over 6 mm ~ {3*0.006/r*B/B0*1e6:7.3f} ppm")
