#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NMR3-FMT / прототип экспортёра
==============================
Синтез ССИ (FID) 1H этанола в DMSO-d6 для поля 0.51 Тл (21.7 МГц)
и запись в JCAMP-DX 5.01, DATA CLASS = NTUPLES, DATA TYPE = NMR FID.

Цель файла — проверка ТРАКТА "массив + параметры -> файл, читаемый Mnova",
а не физическая точность спектра.

Химические сдвиги: Gottlieb 1997, столбец DMSO-d6 (документ NMR3-GOTT):
    CH3  1.06 ppm   t, J = 7.0 Гц      3H
    CH2  3.44 ppm   q(J=7.0) x d(J=5.0) 2H
    OH   4.63 ppm   t, J = 5.0 Гц      1H
"""

import math
import numpy as np
from datetime import date

# ---------------------------------------------------------------- параметры

BF1_MHZ      = 21.700000     # базовая частота, МГц  (= SF, частота 0 ppm)
CARRIER_PPM  = 2.85          # положение несущей O1 на шкале ppm
SW_HZ        = 500.0         # ширина спектра, Гц  (23.04 ppm)
TD_COMPLEX   = 2048          # комплексных отсчётов  -> TD = 4096
NS           = 2000          # число накоплений (только в метаданные)
SNR_PEAK     = 400.0         # отношение сигнал/шум по самой высокой линии
CONJUGATE    = False         # True -> зеркалирует ось частот (см. README)
PH0_DEG      = 0.0           # искусственная ошибка фазы 0-го порядка
SEED         = 20260904

# спиновая система: (ppm, число протонов, [(J, число соседей), ...])
SPIN_SYSTEM = [
    ("CH3", 1.06, 3.0, [(7.0, 2)]),
    ("CH2", 3.44, 2.0, [(7.0, 3), (5.0, 1)]),
    ("OH",  4.63, 1.0, [(5.0, 2)]),
]

# --------------------------------------------------------------- мультиплеты

def multiplet(j_list):
    """Свёртка биномиальных расщеплений -> [(смещение, вес), ...], сумма весов = 1."""
    offsets = np.array([0.0])
    weights = np.array([1.0])
    for j, n in j_list:
        w = np.array([math.comb(n, k) for k in range(n + 1)], dtype=float)
        w /= w.sum()
        o = (np.arange(n + 1) - n / 2.0) * j
        offsets = (offsets[:, None] + o[None, :]).ravel()
        weights = (weights[:, None] * w[None, :]).ravel()
    return offsets, weights


def build_lines(carrier_ppm=None):
    """Список (частота_Гц_от_несущей, амплитуда) для всех линий."""
    if carrier_ppm is None:
        carrier_ppm = CARRIER_PPM
    lines = []
    for name, ppm, nh, jl in SPIN_SYSTEM:
        nu0 = (ppm - carrier_ppm) * BF1_MHZ          # Гц от несущей
        off, w = multiplet(jl)
        for o, wi in zip(off, w):
            lines.append((nu0 + o, nh * wi))
        # площадь группы == nh, т.к. сумма весов = 1
    return lines


def synth_fid(lw_hz, rng, carrier_ppm=None):
    """Комплексный ССИ + шум."""
    if carrier_ppm is None:
        carrier_ppm = CARRIER_PPM
    dw = 1.0 / SW_HZ
    t = np.arange(TD_COMPLEX) * dw
    t2star = 1.0 / (np.pi * lw_hz)

    fid = np.zeros(TD_COMPLEX, dtype=complex)
    for nu, amp in build_lines(carrier_ppm):
        fid += amp * np.exp(2j * np.pi * nu * t)
    fid *= np.exp(-t / t2star)

    if CONJUGATE:
        fid = np.conj(fid)
    if PH0_DEG:
        fid *= np.exp(1j * np.deg2rad(PH0_DEG))

    # шум: нормируем по высоте самой сильной линии в спектре
    spec = np.fft.fftshift(np.fft.fft(fid, 16384))
    peak = np.abs(spec).max()
    sigma_spec = peak / SNR_PEAK
    # пересчёт дисперсии спектр -> временная область
    sigma_t = sigma_spec / np.sqrt(16384 / 2.0)
    fid += rng.normal(0, sigma_t, TD_COMPLEX) + 1j * rng.normal(0, sigma_t, TD_COMPLEX)
    return t, fid, t2star


# --------------------------------------------------------------- JCAMP-DX

def _table(x, y_int, xfactor, per_line=8):
    """Секция (X++(Y..Y)) в форме AFFN."""
    out = []
    n = len(y_int)
    for i in range(0, n, per_line):
        chunk = y_int[i:i + per_line]
        out.append("%.6f " % (x[i] / xfactor) + " ".join("%d" % v for v in chunk))
    return "\n".join(out)


def write_jcamp(path, t, fid, meta):
    n = len(fid)
    re, im = fid.real, fid.imag

    # целочисленное кодирование с общим FACTOR (как делает TopSpin)
    amax = max(np.abs(re).max(), np.abs(im).max())
    yfactor = amax / (2 ** 28)
    re_i = np.rint(re / yfactor).astype(np.int64)
    im_i = np.rint(im / yfactor).astype(np.int64)

    xfactor = 1.0
    dw = 1.0 / meta["sw_hz"]
    last_t = (n - 1) * dw

    sfo1 = meta["bf1"] * (1.0 + meta["carrier_ppm"] * 1e-6)
    o1_hz = meta["carrier_ppm"] * meta["bf1"]
    sw_ppm = meta["sw_hz"] / meta["bf1"]

    H = []
    a = H.append
    a("##TITLE= %s" % meta["title"])
    a("##JCAMP-DX= 5.01")
    a("##DATA TYPE= NMR FID")
    a("##DATA CLASS= NTUPLES")
    a("##ORIGIN= DIY-NMR Model 3 (simulated)")
    a("##OWNER= NMR3 project")
    a("##LONGDATE= %s" % date.today().isoformat())
    a("##SPECTROMETER/DATA SYSTEM= DIY-NMR Model 3 / synthetic")
    a("##.OBSERVE FREQUENCY= %.8f" % sfo1)
    a("##.OBSERVE NUCLEUS= ^1H")
    a("##.ACQUISITION MODE= SIMULTANEOUS")
    a("##.AVERAGES= %d" % meta["ns"])
    a("##.DIGITISER RES= 20")
    a("##.ZERO FILL= 0")
    a("##.SOLVENT NAME= DMSO-D6")
    a("##SAMPLE DESCRIPTION= %s" % meta["sample"])
    # приватные метки Bruker — их Mnova разбирает надёжнее стандартных
    a("##$BF1= %.8f" % meta["bf1"])
    a("##$SFO1= %.8f" % sfo1)
    a("##$SF= %.8f" % meta["bf1"])
    a("##$O1= %.4f" % o1_hz)
    a("##$SW_h= %.4f" % meta["sw_hz"])
    a("##$SW= %.6f" % sw_ppm)
    a("##$TD= %d" % (2 * n))
    a("##$NUC1= <1H>")
    a("##$AQ_mod= 3")
    a("##$DIGMOD= 0")
    a("##$DSPFVS= 0")
    a("##$DECIM= 1")
    a("##$GRPDLY= 0")
    a("##$NS= %d" % meta["ns"])
    a("##$TE= 318.0")
    a("##$PULPROG= <zg>")
    a("##$BYTORDA= 0")

    a("##NTUPLES= NMR FID")
    a("##VAR_NAME= TIME, SPECTRUM/REAL, SPECTRUM/IMAG, PAGE NUMBER")
    a("##SYMBOL= X, R, I, N")
    a("##VAR_TYPE= INDEPENDENT, DEPENDENT, DEPENDENT, PAGE")
    a("##VAR_FORM= AFFN, AFFN, AFFN, AFFN")
    a("##VAR_DIM= %d, %d, %d, 2" % (n, n, n))
    a("##UNITS= SECONDS, ARBITRARY UNITS, ARBITRARY UNITS, ")
    a("##FIRST= 0, %.8g, %.8g, 1" % (re_i[0]*yfactor, im_i[0]*yfactor))
    a("##LAST= %.8f, %.8g, %.8g, 2" % (last_t, re_i[-1]*yfactor, im_i[-1]*yfactor))
    a("##MIN= 0, %.8g, %.8g, 1" % (re_i.min()*yfactor, im_i.min()*yfactor))
    a("##MAX= %.8f, %.8g, %.8g, 2" % (last_t, re_i.max()*yfactor, im_i.max()*yfactor))
    a("##FACTOR= %.10g, %.10g, %.10g, 1" % (xfactor, yfactor, yfactor))

    a("##PAGE= N=1")
    a("##DATA TABLE= (X++(R..R)), XYDATA")
    a(_table(t, re_i, xfactor))

    a("##PAGE= N=2")
    a("##DATA TABLE= (X++(I..I)), XYDATA")
    a(_table(t, im_i, xfactor))

    a("##END NTUPLES= NMR FID")
    a("##END=")

    with open(path, "w", encoding="ascii", newline="\r\n") as f:
        f.write("\n".join(H) + "\n")
    return path


# --------------------------------------------------------- обратный разбор

def read_jcamp(path):
    """Минимальный ридер — только для самопроверки записанного файла."""
    txt = open(path, "r", encoding="ascii").read().replace("\r\n", "\n")
    meta, pages, cur = {}, {}, None
    for line in txt.split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith("##"):
            key, _, val = line[2:].partition("=")
            key = key.strip().upper()
            if key == "PAGE":
                cur = val.strip()
                pages[cur] = []
            elif key == "DATA TABLE":
                pass
            elif key.startswith("END"):
                cur = None
            else:
                meta[key] = val.strip()
                cur = None if key != "PAGE" else cur
            continue
        if cur is not None:
            pages[cur].append(line)

    factors = [float(v) for v in meta["FACTOR"].split(",")]
    xf, yf = factors[0], factors[1]

    def page_to_array(key):
        vals = []
        for line in pages[key]:
            parts = line.split()
            vals.extend(int(v) for v in parts[1:])
        return np.array(vals, dtype=float) * yf

    r = page_to_array("N=1")
    i = page_to_array("N=2")
    n = int(meta["VAR_DIM"].split(",")[0])
    dw = float(meta["$SW_H"]) ** -1
    return meta, np.arange(n) * dw, r[:n] + 1j * i[:n]


# --------------------------------------------------------------------- main

def process(fid, sfo1, bf1, sw_hz, si=32768):
    """FT в стиле Mnova: экспоненциальное окно, zero-fill, magnitude."""
    n = len(fid)
    t = np.arange(n) / sw_hz
    w = np.exp(-np.pi * 0.3 * t)                # LB = 0.3 Гц
    spec = np.fft.fftshift(np.fft.fft(fid * w, si))
    f = np.fft.fftshift(np.fft.fftfreq(si, 1.0 / sw_hz))   # Гц от несущей
    ppm = (sfo1 + f * 1e-6) / bf1 * 1e6 - 1e6
    return ppm, spec


if __name__ == "__main__":
    import sys, json
    rng = np.random.default_rng(SEED)

    cases = [
        # (файл, ширина линии Гц, несущая ppm, заголовок)
        ("ethanol_21p7MHz_lw1p2Hz.dx", 1.2, 0.0,
         "Ethanol in DMSO-d6, 21.7 MHz, LW 1.2 Hz (Branch B, multiplets resolved)"),
        ("ethanol_21p7MHz_lw8Hz.dx", 8.0, 0.0,
         "Ethanol in DMSO-d6, 21.7 MHz, LW 8 Hz (Branch A, singlets only)"),
        ("ethanol_21p7MHz_o1test.dx", 1.2, 2.85,
         "Ethanol in DMSO-d6, 21.7 MHz, LW 1.2 Hz, carrier at 2.85 ppm (referencing test)"),
    ]

    report = {}
    for fname, lw, carrier, title in cases:
        t, fid, t2s = synth_fid(lw, rng, carrier_ppm=carrier)
        meta = dict(title=title, bf1=BF1_MHZ, carrier_ppm=carrier,
                    sw_hz=SW_HZ, ns=NS,
                    sample="Ethanol / DMSO-d6, synthetic, LW %.1f Hz, O1 at %.2f ppm"
                           % (lw, carrier))
        write_jcamp(fname, t, fid, meta)

        m2, t2, fid2 = read_jcamp(fname)
        err = np.abs(fid2 - fid).max() / np.abs(fid).max()
        import os
        report[fname] = dict(linewidth_hz=lw, carrier_ppm=carrier,
                             T2star_s=round(t2s, 4),
                             roundtrip_rel_err=float(err),
                             td=int(m2["$TD"]), sfo1=m2["$SFO1"],
                             size_bytes=os.path.getsize(fname))
        print("%-30s LW=%4.1f Hz  O1=%.2f ppm  T2*=%.3f s  round-trip=%.1e  %6d B"
              % (fname, lw, carrier, t2s, err, os.path.getsize(fname)))

    json.dump(report, open("generation_report.json", "w"), indent=2)
