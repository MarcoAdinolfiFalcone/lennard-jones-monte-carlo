#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt

def read_two_col(path):
    # file: N  logW
    x, y = [], []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            x.append(int(float(parts[0])))
            y.append(float(parts[1]))
    return np.array(x, dtype=int), np.array(y, dtype=float)

def read_bad_windows(path):
    # accettiamo diversi formati: righe con "n" oppure "n n+1"
    bad = set()
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            try:
                n = int(float(parts[0]))
                bad.add(n)
            except:
                pass
    return bad

def main():
    # ESEMPIO:
    #   python3 plot_logW_muRef.py ../data/SUS/muVT_T1p40_V343_muRef_m3p50
    import sys
    if len(sys.argv) < 2:
        print("Usage: python3 plot_logW_muRef.py <run_dir>")
        raise SystemExit(1)

    run_dir = sys.argv[1]
    logw_path = os.path.join(run_dir, "logW_from_windows.dat")
    if not os.path.isfile(logw_path):
        raise FileNotFoundError(f"Missing: {logw_path}")

    N, logW = read_two_col(logw_path)

    # prova a leggere bad windows (se presenti)
    bad = set()
    for cand in ["bad_windows.dat", "bad_windows_from_windows.dat"]:
        p = os.path.join(run_dir, cand)
        if os.path.isfile(p):
            bad = read_bad_windows(p)
            bad_file_used = cand
            break
    else:
        bad_file_used = None

    # ordina
    idx = np.argsort(N)
    N = N[idx]
    logW = logW[idx]

    # plot
    fig = plt.figure()
    plt.plot(N, logW, marker="o", linestyle="None", markersize=3)

    # evidenzia N "bad" se disponibili
    if bad:
        mask = np.array([n in bad for n in N], dtype=bool)
        if np.any(mask):
            plt.plot(N[mask], logW[mask], marker="x", linestyle="None",
                     markersize=6, mew=1.5, label="bad window (edge/poor stats)")

    plt.xlabel(r"$N$")
    plt.ylabel(r"$\log W(N;\mu_{\mathrm{ref}}^*)$ (gauge-fixed)")
    plt.title(r"SUS reconstruction at $T^*=1.4$ (run at $\mu_{\mathrm{ref}}^*$)")
    plt.grid(True, alpha=0.3)

    if bad_file_used is not None:
        plt.legend(loc="best")

    # output in figure folder di Overleaf (se esiste) altrimenti in run_dir
    out1 = os.path.join(run_dir, "logW_muRef.png")
    plt.tight_layout()
    plt.savefig(out1, dpi=300)
    print("Saved:", out1)

if __name__ == "__main__":
    main()
