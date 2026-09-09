#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Серия JCAMP-DX: этанол/DMSO-d6 @ 21.7 МГц при разных ширинах линии.

Отличие от первого набора: шум задан ОДИН на всю серию во ВРЕМЕННОЙ области.
Это физически правильно — для одного образца и одного NS шум приёмника не
зависит от однородности поля. Площадь линии сохраняется, поэтому высота пика
падает ~1/LW, и широкие линии проигрывают не только в разрешении, но и в SNR.
Нормировка шума по высоте пика (как в первом наборе) этот эффект скрывает.
"""
import math
import numpy as np
from scipy.signal import find_peaks

import make_ethanol_jcamp as G

BF1 = G.BF1_MHZ
SW  = G.SW_HZ
TD  = G.TD_COMPLEX
CARRIER = 0.0
SNR_REF_AT_LW = 3.0      # опорная ширина линии, Гц
SNR_REF       = 250.0    # peak SNR по центральной линии CH3 при опорной ширине

SERIES = [1.20, 2.00, 3.00, 4.34, 6.51, 10.85, 21.70]


def fid_noiseless(lw):
    dw = 1.0 / SW
    t = np.arange(TD) * dw
    f = np.zeros(TD, dtype=complex)
    for nu, amp in G.build_lines(CARRIER):
        f += amp * np.exp(2j * np.pi * nu * t)
    return t, f * np.exp(-t / (1.0 / (np.pi * lw)))


def ft(fid, lb=0.3, si=32768):
    n = len(fid); t = np.arange(n) / SW
    w = np.exp(-np.pi * lb * t); w[0] *= 0.5
    sp = np.fft.fftshift(np.fft.fft(fid * w, si))
    fr = np.fft.fftshift(np.fft.fftfreq(si, 1.0 / SW))
    return fr / BF1, sp          # ppm (несущая на 0 ppm)


# --- калибровка общего уровня шума по опорной ширине -----------------------
_, ref = fid_noiseless(SNR_REF_AT_LW)
ppm, sp = ft(ref)
peak_ref = np.abs(sp.real).max()
# sigma спектра = sigma_t * sqrt(TD) для каждой квадратуры
# эмпирическая калибровка: измеряем реальный sigma спектра от единичного sigma_t
_rng0 = np.random.default_rng(1)
_n = _rng0.normal(0,1,TD) + 1j*_rng0.normal(0,1,TD)
_pp, _sn = ft(_n)
_gain = _sn.real.std()                      # sigma_spec при sigma_t = 1
SIGMA_T = (peak_ref / SNR_REF) / _gain

rng = np.random.default_rng(20260905)


def analyse(ppm, y, centre, half, J_ppm, noise=None):
    """число линий и глубина провала в окне вокруг мультиплета"""
    m = (ppm > centre - half) & (ppm < centre + half)
    x, yy = ppm[m], y[m]
    prom = max(yy.max() * 0.02, 5.0 * noise) if noise else yy.max() * 0.02
    pk, _ = find_peaks(yy, prominence=prom, height=5.0 * noise if noise else None)
    n = len(pk)
    dip = None
    if n >= 3:
        pos = np.sort(x[pk])
        lo, hi = pos[len(pos) // 2 - 1], pos[len(pos) // 2]
        seg = (x > lo) & (x < hi)
        if seg.any():
            v = yy[seg].min()
            o = min(yy[np.argmin(abs(x - lo))], yy[np.argmin(abs(x - hi))])
            dip = 1 - v / o
    return n, dip


rows = []
files = []
for lw in SERIES:
    t, f = fid_noiseless(lw)
    f = f + rng.normal(0, SIGMA_T, TD) + 1j * rng.normal(0, SIGMA_T, TD)

    lw_ppm = lw / BF1
    fname = "etoh_21p7MHz_LW%05.2fHz_%05.3fppm.dx" % (lw, lw_ppm)
    fname = fname.replace(".", "p", 2).replace("p dx", ".dx")
    fname = "etoh_21p7MHz_LW%sHz_%sppm.dx" % (
        ("%.2f" % lw).replace(".", "p"), ("%.3f" % lw_ppm).replace(".", "p"))

    meta = dict(title="Ethanol in DMSO-d6, 21.7 MHz, LW %.2f Hz = %.3f ppm" % (lw, lw_ppm),
                bf1=BF1, carrier_ppm=CARRIER, sw_hz=SW, ns=G.NS,
                sample="synthetic, common time-domain noise, LW %.2f Hz" % lw)
    G.write_jcamp(fname, t, f, meta)
    files.append(fname)

    ppm, sp = ft(f)
    y = sp.real
    noise = y[(ppm > 8) & (ppm < 11)].std()
    snr = y.max() / noise

    n3, d3 = analyse(ppm, y, 1.06, 0.75, 7 / BF1, noise)      # CH3 триплет
    n2, d2 = analyse(ppm, y, 3.44, 0.95, 7 / BF1, noise)      # CH2 q x d, 8 линий
    n1, d1 = analyse(ppm, y, 4.63, 0.55, 5 / BF1, noise)      # OH триплет

    # разделение групп: минимум между CH2 и OH
    seg = (ppm > 3.9) & (ppm < 4.3)
    sep = 1 - y[seg].min() / min(y[(ppm > 3.2) & (ppm < 3.7)].max(),
                                 y[(ppm > 4.4) & (ppm < 4.9)].max())
    rows.append((lw, lw_ppm, snr, n3, d3, n2, n1, d1, sep))

print("  LW,Гц  LW,ppm   SNR |  CH3: линий провал | CH2 линий | OH: линий провал | CH2/OH разд.")
for lw, lp, snr, n3, d3, n2, n1, d1, sep in rows:
    f3 = "%3.0f%%" % (100 * d3) if d3 else "  — "
    f1 = "%3.0f%%" % (100 * d1) if d1 else "  — "
    print("%6.2f  %6.3f %5.0f |     %d    %s     |     %d     |    %d    %s     |   %3.0f%%"
          % (lw, lp, snr, n3, f3, n2, n1, f1, 100 * sep))

print("\nФайлы:")
for f in files:
    print(" ", f)
