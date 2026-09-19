import numpy as np
f0=21.7e6; w0=2*np.pi*f0
print(f"f0 = {f0/1e6} MHz, w0 = {w0:.4e} rad/s\n")

print("=== 1. Quench (phase-inverted) pulse: null time ===")
print("V(t) = -k*V0 + (1+k)V0*exp(-t/tau)  ->  t_null = tau*ln((1+k)/k)")
for k in [0.5,0.75,1.0,1.5,2.0,3.0]:
    print(f"  k(quench/main amplitude) = {k:4.2f}   t_null = {np.log((1+k)/k):.4f} tau")

print("\n=== 2. tau during TX (probe loaded by 50 ohm TX) and quench duration ===")
print(f"{'Q_free':>7}{'Q_eff(TX)':>11}{'tau,us':>9}{'t_null,ns':>11}{'ticks@10ns':>12}")
for Q in [80,120,150]:
    for load,name in [(2,'50R match')]:
        Qe=Q/load; tau=2*Qe/w0; tn=tau*np.log(2)
        print(f"{Q:7.0f}{Qe:11.0f}{tau*1e6:9.3f}{tn*1e9:11.1f}{tn/10e-9:12.1f}   ({name})")
print("  mismatched TX (Q_eff=12, sigma=5):")
Qe=12; tau=2*Qe/w0; tn=tau*np.log(2)
print(f"{'-':>7}{Qe:11.0f}{tau*1e6:9.3f}{tn*1e9:11.1f}{tn/10e-9:12.1f}")

print("\n=== 3. Quench residual budget (residual/V0) ===")
tau=2*60/w0   # Q_eff=60 during TX, Q_free=120
print(f"  reference tau = {tau*1e9:.0f} ns")
rows=[("timing grain 10 ns (F411 @100 MHz)", 10e-9/tau),
      ("timing grain 7.69 ns (AD9959 SYSCLK 130 MHz)", 7.69e-9/tau),
      ("ASF amplitude trim, 10 bit  (dt = tau/2048)", 1/2048),
      ("phase word 14 bit (0.022 deg)", np.deg2rad(360/2**14)),
      ("TX amplitude reproducibility 1 %", 0.01),
      ("TX amplitude reproducibility 0.1 %", 0.001),
      ("TX amplitude reproducibility 0.01 %", 0.0001)]
for n,r in rows:
    print(f"  {n:46s} residual = {r:.2e}  = {-20*np.log10(r):5.1f} dB  = {np.log(1/r):4.1f} tau saved")

print("\n=== 4. Feedback damping: required open-loop gain  a = 1 + j*D*C/(Q*C') ===")
print(f"{'wL,ohm':>7}{'L,uH':>7}{'C,pF':>7}{'Q':>5}{'Q*wL,kohm':>11}{'D':>4}{'C-prime,pF':>11}{'|a-1|':>8}{'Zin,ohm':>9}")
for wL in [55,224]:
    L=wL/w0; C=1/(w0**2*L)
    for Q in [120]:
        for D in [10]:
            for Cp in [0.5,1.0,2.0]:
                a1=D*C/(Q*Cp*1e-12); Zin=1/(w0*Cp*1e-12*a1)
                print(f"{wL:7.0f}{L*1e6:7.3f}{C*1e12:7.1f}{Q:5.0f}{Q*wL/1000:11.2f}{D:4.0f}{Cp:11.1f}{a1:8.2f}{Zin:9.0f}")

print("\n=== 5. Detuning from feedback phase error theta ===")
print("  C_eq = C'*|a-1|*sin(theta) = (C/Q)*D*sin(theta);  df/f = -C_eq/(2C)")
for wL in [55,224]:
    L=wL/w0; C=1/(w0**2*L); Q=120; D=10
    for th in [1,3,10]:
        Ceq=(C/Q)*D*np.sin(np.deg2rad(th))
        print(f"  wL={wL:3.0f} ohm, C={C*1e12:5.1f} pF, theta={th:2.0f} deg -> C_eq={Ceq*1e12:5.2f} pF, df/f={-Ceq/(2*C)*1e6:7.0f} ppm = {-Ceq/(2*C)*f0/1e3:7.1f} kHz")
        
print("\n=== 6. Bandwidth vs NF trade (Hoult eq 20/21) ===")
print("  S=(F-1)/(Fopt-1);  Q*w0*L = S*R0;  dw = sqrt(S^2-1)*w0/Q")
Fopt=10**0.1
for NFdB in [1.0,1.5,2.0,3.0]:
    F=10**(NFdB/10); S=(F-1)/(Fopt-1)
    if S<1: print(f"  NF target {NFdB} dB : below Fopt"); continue
    dw=np.sqrt(S**2-1)
    print(f"  NF target {NFdB:3.1f} dB -> S={S:5.2f}, Q*w0*L={S:5.2f}*R0, dw={dw:5.2f}*w0/Q = {dw*f0/120/1e3:7.1f} kHz (Q=120)")
print(f"  our required band: 10 ppm span = {10*21.7:.0f} Hz; probe band w0/Q at Q=120 = {f0/120/1e3:.0f} kHz")

print("\n=== 7. Transmitter power for faster rise/fall (eq 26): P x sigma/2 ===")
print(f"{'sigma':>6}{'Q_eff':>7}{'tau,us':>9}{'P_mismatch,W':>14}{'P_damping,W':>13}")
P0=0.154
for s in [1,2,5,10,20]:
    Qe=120/(2*s); tau=2*Qe/w0
    print(f"{s:6.0f}{Qe:7.1f}{tau*1e6:9.3f}{P0*s/2:14.3f}{P0*2*s:13.3f}")
