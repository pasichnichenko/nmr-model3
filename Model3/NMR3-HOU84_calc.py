import numpy as np
f0=21.7e6; w0=2*np.pi*f0; Z0=50.0
print(f"f0={f0/1e6} MHz  w0={w0:.4e} rad/s\n")
hdr=f"{'wL,Om':>6}{'L,uH':>7}{'Q':>5}{'Ropt,Om':>9}{'R\",Om':>8}{'Qeff':>7}{'tau_unl':>9}{'tau_50':>8}{'tau_H':>8}{'6t_unl':>8}{'6t_50':>8}{'6t_H':>7}{'gain':>7}{'dSNR':>7}{'conv':>7}"
print(hdr); print('-'*len(hdr))
for wL in [50,75,100,150,224]:
    L=wL/w0
    for Q in [80,120,150]:
        Ropt=np.sqrt(Q*Z0*wL); Rpp=2*Ropt; Qeff=2*np.sqrt(Q*Z0/wL)
        t_unl=2*Q/w0; t_50=Q/w0; t_H=2*Qeff/w0
        pen=10*np.log10(1+Qeff/(2*Q)); conv=10*np.log10(Q/Qeff)
        print(f"{wL:6.0f}{L*1e6:7.3f}{Q:5.0f}{Ropt:9.0f}{Rpp:8.0f}{Qeff:7.1f}"
              f"{t_unl*1e6:9.2f}{t_50*1e6:8.2f}{t_H*1e6:8.3f}"
              f"{6*t_unl*1e6:8.1f}{6*t_50*1e6:8.1f}{6*t_H*1e6:7.2f}"
              f"{t_50/t_H:7.2f}{pen:7.2f}{conv:7.2f}")
print("\ntau in us; 6t = 6 tau (project criterion); gain = tau_50/tau_H;")
print("dSNR = 10log10(1+Qeff/2Q) dB penalty of Hoult damping; conv = 10log10(Q/Qeff) dB penalty of classic series damping")

print("\n--- 4-pole optimum (Appendix), Qeff=3.1, x=1, y=0.63 ---")
print(f"{'wL':>5}{'Q':>5}{'Ropt':>7}{'L,uH':>8}{'r,Om':>7}{'tau4,ns':>9}{'6t4,us':>8}{'SRF@2pF':>9}{'X_L/X_Cstray':>13}")
for wL in [50,100,224]:
    for Q in [120]:
        Ropt=np.sqrt(Q*Z0*wL); Lp=Ropt/w0; rp=0.63*Ropt
        t4=2*3.125/w0
        srf=1/(2*np.pi*np.sqrt(Lp*2e-12))
        ratio=(1/(w0*2e-12))/Ropt
        print(f"{wL:5.0f}{Q:5.0f}{Ropt:7.0f}{Lp*1e6:8.2f}{rp:7.0f}{t4*1e9:9.1f}{6*t4*1e6:8.2f}{srf/1e6:9.1f}{ratio:13.1f}")

print("\n--- cable lengths at 21.7 MHz ---")
c=3e8
for name,vf in [("RG-58/RG-316 PE (0.66)",0.66),("PTFE RG-142/400 (0.695)",0.695),("foam PE (0.82)",0.82)]:
    lam=c/f0*vf
    print(f"{name:26s} lambda={lam:6.2f} m   lambda/2={lam/2:5.2f} m   lambda/4={lam/4:5.2f} m")

print("\n--- first-order phase error vs dead time, 10 ppm span at 21.7 MHz ---")
span=10*21.7   # Hz
for td in [5e-6,10e-6,20e-6,50e-6,100e-6]:
    print(f"t_dead={td*1e6:6.1f} us -> phi1 = {360*span*td:6.2f} deg across 10 ppm ; signal loss at T2*=29.3 ms: {100*(1-np.exp(-td/29.3e-3)):.3f} %")
