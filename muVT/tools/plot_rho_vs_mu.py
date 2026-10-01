#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt

def read_eos_mu_file(path):
    # columns: mu  rho  P  Nmean  lnXi
    data = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if (not s) or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 5:
                continue
            mu, rho, P, Nmean, lnXi = map(float, parts[:5])
            data.append((mu, rho, P, Nmean, lnXi))
    if not data:
        raise RuntimeError(f"No data read from {path}")
    return np.array(data)

def main(run_dir):
    infile = os.path.join(run_dir, "EOS_mu_rho_P.dat")
    outpng = os.path.join(run_dir, "rho_vs_mu.png")

    arr = read_eos_mu_file(infile)
    mu   = arr[:, 0]
    rho  = arr[:, 1]

    # sort by mu (robust)
    idx = np.argsort(mu)
    mu, rho = mu[idx], rho[idx]

    plt.figure()
    plt.plot(mu, rho, marker="o", linewidth=1.0, markersize=3.5)
    plt.xlabel(r"$\mu^*$")
    plt.ylabel(r"$\rho^*(\mu^*)$")
    plt.title(r"$\mu VT$ + SUS at $T^*=1.4$, $V^*=343$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    print(f"Saved: {outpng}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python3 plot_rho_vs_mu.py <run_dir>")
        sys.exit(1)
    main(sys.argv[1])
