#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt

def read_eos_mu_rho_P(path):
    # columns: mu rho P Nmean lnXi
    data = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if (not s) or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            mu  = float(parts[0])
            rho = float(parts[1])
            P   = float(parts[2])
            data.append((mu, rho, P))
    if not data:
        raise RuntimeError(f"No data read from {path}")
    arr = np.array(data, dtype=float)
    mu, rho, P = arr[:,0], arr[:,1], arr[:,2]
    # sort by mu just in case
    idx = np.argsort(mu)
    return mu[idx], rho[idx], P[idx]

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", help="e.g. data/SUS/muVT_T1p40_V343_muRef_m3p50")
    ap.add_argument("--infile", default="EOS_mu_rho_P.dat")
    ap.add_argument("--out", default=None, help="output png (default: <run_dir>/P_vs_mu.png)")
    ap.add_argument("--title", default=r"$\mu VT$ + SUS at $T^*=1.4,\ V^*=343$")
    args = ap.parse_args()

    in_path = os.path.join(args.run_dir, args.infile)
    out_path = args.out if args.out else os.path.join(args.run_dir, "P_vs_mu.png")

    mu, rho, P = read_eos_mu_rho_P(in_path)

    plt.figure()
    plt.plot(mu, P, marker="o", linestyle="None")
    plt.grid(True, alpha=0.35)
    plt.xlabel(r"$\mu^*$")
    plt.ylabel(r"$P^*(\mu^*)$")
    plt.title(args.title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    print(f"Saved: {out_path}")

if __name__ == "__main__":
    main()
