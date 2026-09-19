#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NMR3-MISP — численная проверка выводов Mispelter/Lupu/Briguet (2015)
применительно к параметрам проекта NMR3 (21.7 МГц, 1H).

Формулы и номера уравнений — по 2-му изданию книги.
"""
import math

f_MHz = 21.7
f = f_MHz * 1e6
w = 2 * math.pi * f
Z0 = 50.0
GAMMA_H = 42.57748e6          # Гц/Тл (Приложение A, табл. A.1)

# ---------------------------------------------------------------- 1. Соленоид
def L_terman(r_mm, l_mm, n):
    """Eq. (2.85), приближение Нагаоки. L в нГн, размеры в мм, r — РАДИУС."""
    return 39.4 * r_mm**2 * n**2 / (9 * r_mm + 10 * l_mm)

def psi(D, l):
    """Eq. (2.91), Callendar/Medhurst, для оптимального шага и n>3."""
    return 1.0 / (1.03 + 0.4 * D / l)

def Q_medhurst(D_mm, l_mm, f_MHz):
    """Eq. (2.90): Q = 7.5 * D * psi * sqrt(f_MHz)."""
    return 7.5 * D_mm * psi(D_mm, l_mm) * math.sqrt(f_MHz)

def lam_srf(alpha, l_wire_mm):
    """Eq. (2.92): lambda_SRF = alpha_SRF * l_wire."""
    return alpha * l_wire_mm

print("=" * 74)
print("1. СОЛЕНОИД 21.7 МГц — геометрия, L, Q, число витков")
print("=" * 74)
print(f"{'D,мм':>5} {'l,мм':>5} {'l/D':>5} {'psi':>6} {'Q_Medh':>7} "
      f"{'n(X=100Ом)':>10} {'L,нГн':>7} {'r,Ом':>6} {'tau,мкс':>8}")
for D, l in [(8, 8), (10, 10), (10, 14), (12, 12), (12, 17), (14, 14), (16, 16)]:
    Q = Q_medhurst(D, l, f_MHz)
    L_target = 100.0 / w * 1e9           # нГн для X_L = 100 Ом
    # подбор n
    n = math.sqrt(L_target / L_terman(D / 2, l, 1))
    L = L_terman(D / 2, l, n)
    r = (L * 1e-9 * w) / Q
    tau = 2 * Q / w * 1e6
    print(f"{D:5.0f} {l:5.0f} {l/D:5.2f} {psi(D,l):6.3f} {Q:7.1f} "
          f"{n:10.1f} {L:7.0f} {r:6.3f} {tau:8.2f}")

# ---------------------------------------------------------------- 2. Кабель
def alpha_dB_per_m(a, b, f_MHz):
    """Eq. (2.59) / табл. 2.4: alpha = a*sqrt(f) + b*f, дБ/м."""
    return a * math.sqrt(f_MHz) + b * f_MHz

CABLES = {                    # (a, b) из табл. 2.4
    "RG316": (0.027, 10e-5),
    "RG223": (0.013, 8e-5),
    "RG402": (0.012, 3e-5),
    "RG401": (0.007, 2.5e-5),
    "RG214": (0.006, 7e-5),
}

def efficiency(ZL, length_m, a, b, f_MHz):
    """Eq. (2.71): PL/Pin = (1-|G|^2) / (A - |G|^2/A)."""
    G = abs((ZL - Z0) / (ZL + Z0))
    ML = alpha_dB_per_m(a, b, f_MHz) * length_m
    A = 10 ** (ML / 10.0)
    E = (1 - G**2) / (A - G**2 / A)
    return E, G, ML

print()
print("=" * 74)
print("2. ПОТЕРИ В ОТРЕЗКЕ RG402 МЕЖДУ КАТУШКОЙ И БЛОКОМ СОГЛАСОВАНИЯ")
print("   (Eq. 2.71; катушка X_L=100 Ом, Q=250 -> r=0.40 Ом, 21.7 МГц)")
print("=" * 74)
Qc = 250.0
XL = 100.0
r_coil = XL / Qc
Z_raw = complex(r_coil, XL)              # катушка БЕЗ преднастройки
Z_pre = complex(r_coil, 0.0)             # катушка предварительно настроена в серии
print(f"|Gamma| без преднастройки = {abs((Z_raw-Z0)/(Z_raw+Z0)):.5f}")
print(f"|Gamma| с преднастройкой  = {abs((Z_pre-Z0)/(Z_pre+Z0)):.5f}")
print()
print(f"{'L,см':>6} | {'RG402 сырая':>14} {'RG402 предн.':>14} "
      f"{'RG401 сырая':>14} {'RG316 сырая':>14}")
for Lcm in (5, 10, 15, 20, 30, 50, 100, 150):
    row = f"{Lcm:6.0f} |"
    for cab, Z in (("RG402", Z_raw), ("RG402", Z_pre),
                   ("RG401", Z_raw), ("RG316", Z_raw)):
        a, b = CABLES[cab]
        E, G, ML = efficiency(Z, Lcm / 100, a, b, f_MHz)
        row += f" {E:6.3f} ({10*math.log10(E):+5.2f}дБ)"
    print(row)

print()
print("   Для сравнения — та же схема на 200 МГц (пример из книги, рис. 4.18):")
for Lcm in (10,):
    a, b = CABLES["RG316"]
    E, G, ML = efficiency(complex(1.0, 100.0), Lcm/100, a, b, 200.0)
    print(f"   RG316, 10 см, Q=100, 200 МГц: E={E:.3f} ({10*math.log10(E):+.2f} дБ)"
          f"  [книга: измерено 2–3 дБ]")

# ---------------------------------------------------------------- 3. B1 и 90°
print()
print("=" * 74)
print("3. B1 И ДЛИТЕЛЬНОСТЬ 90° (рис. 10.21: I=sqrt(2*P*Q*w*Ceq),")
print("   B1+ = 2*pi*n*I/sqrt(d^2+l^2) гаусс; Ceq = C для соленоида)")
print("=" * 74)
print(f"{'P,мВт':>6} {'D,мм':>5} {'l,мм':>5} {'n':>4} {'Q':>5} {'I,А':>6} "
      f"{'B1+,Гс':>7} {'nu1,кГц':>8} {'t90,мкс':>8} {'Vpk,В':>7}")
for P in (0.05, 0.1, 0.5, 1.0):
    for D, l in [(10, 10), (12, 17)]:
        Q = Q_medhurst(D, l, f_MHz)
        L_nH = 100.0 / w * 1e9
        n = math.sqrt(L_nH / L_terman(D / 2, l, 1))
        L = L_nH * 1e-9
        C = 1.0 / (L * w**2)
        I = math.sqrt(2 * P * Q * w * C)
        B1_G = 2 * math.pi * n * I / math.sqrt(D**2 + l**2)
        B1_T = B1_G * 1e-4
        nu1 = GAMMA_H * B1_T
        t90 = 1 / (4 * nu1) * 1e6
        Vpk = I / (C * w)
        print(f"{P*1000:6.0f} {D:5.0f} {l:5.0f} {n:4.1f} {Q:5.0f} {I:6.3f} "
              f"{B1_G:7.3f} {nu1/1000:8.2f} {t90:8.1f} {Vpk:7.1f}")

# ---------------------------------------------------------------- 4. Звон
print()
print("=" * 74)
print("4. ПОСТОЯННАЯ ЗВОНА tau=2Q/w (Eq. 4.85) И МЁРТВОЕ ВРЕМЯ")
print("=" * 74)
print(f"{'Q':>5} {'tau,мкс':>9} {'t(1кВ->1мкВ)=20.7tau':>22} "
      f"{'t(70В->0.5мкВ)':>16}")
for Q in (50, 100, 150, 200, 250, 300):
    tau = 2 * Q / w
    t_book = math.log(1e9) * tau           # 1 кВ -> 1 мкВ
    t_nmr3 = math.log(70 / 0.5e-6) * tau   # 70 В -> 0.5 мкВ
    print(f"{Q:5.0f} {tau*1e6:9.2f} {t_book*1e6:22.1f} {t_nmr3*1e6:16.1f}")

# ------------------------------------------------- 5. Скин-слой и сопротивление
print()
print("=" * 74)
print("5. СКИН-СЛОЙ И r_AC ПРОВОДА (Eq. 2.25, 2.27)")
print("=" * 74)
delta = 66 / math.sqrt(f_MHz)
print(f"delta(Cu) @ {f_MHz} МГц = {delta:.2f} мкм")
print(f"{'d,мм':>6} {'r_AC,мОм/м':>12} {'d/delta':>9}")
for d in (0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5):
    rac = 83.2 * math.sqrt(f_MHz) / d
    print(f"{d:6.2f} {rac:12.1f} {d/(delta*1e-3):9.0f}")

# ------------------------------------------------- 6. Q_matched, ошибка вдвое
print()
print("=" * 74)
print("6. Q ИЗ S11: Q_matched = Q/2 (Eq. 10.23) — типичная ошибка вдвое")
print("=" * 74)
for Q in (150, 250):
    Qm = Q / 2
    print(f"истинный Q={Q:3.0f} -> по -3дБ ширине S11 измерится {Qm:5.1f}; "
          f"полоса df = f/Qm = {f/Qm/1e3:6.1f} кГц")
