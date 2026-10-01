#!/usr/bin/env python3
# ==========================================================
#  plot_eos_sus.py
#  Visualizza i file prodotti da EOS_from_SUS_muVT.cpp
#
#  Input (nella stessa cartella o via --run_dir):
#    - EOS_mu_rho_P.dat
#    - (opzionale) W_reweighted_mu_*.dat
#
#  Output:
#    - EOS_P_vs_rho.png
#    - EOS_rho_vs_mu.png
#    - EOS_P_vs_mu.png
#    - (opzionale) W_reweighted_examples.png
#
#  Uso:
#    python3 plot_eos_sus.py --run_dir data/SUS/muVT_T1p40_V343_muRef_m3p50
#    python3 plot_eos_sus.py --run_dir ... --plot_W --max_W_files 6
# ==========================================================

import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt


def load_eos(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File non trovato: {path}")
    # ignora commenti '#'
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 5:
        raise ValueError(
            f"Formato inatteso in {path}. Atteso >=5 colonne: mu rho P Nmean lnXi. "
            f"Trovate {data.shape[1]} colonne."
        )

    mu = data[:, 0]
    rho = data[:, 1]
    P = data[:, 2]
    Nmean = data[:, 3]
    lnXi = data[:, 4]
    return mu, rho, P, Nmean, lnXi


def maybe_sort_by(x, *arrays):
    idx = np.argsort(x)
    return (x[idx],) + tuple(a[idx] for a in arrays)


def savefig(path, dpi=200):
    plt.tight_layout()
    plt.savefig(path, dpi=dpi)
    print(f"[saved] {path}")


def plot_eos(mu, rho, P, Nmean, lnXi, outdir):
    # Ordina per mu per i grafici vs mu
    mu_s, rho_s, P_s, Nmean_s, lnXi_s = maybe_sort_by(mu, rho, P, Nmean, lnXi)

    # -----------------------
    # 1) P vs rho (EOS)
    # -----------------------
    plt.figure()
    plt.plot(rho, P, "o")
    plt.xlabel(r"$\rho = \langle N\rangle/V$")
    plt.ylabel(r"$P$")
    plt.title("EOS from SUS (muVT): P vs rho")
    savefig(os.path.join(outdir, "EOS_P_vs_rho.png"))

    # -----------------------
    # 2) rho vs mu
    # -----------------------
    plt.figure()
    plt.plot(mu_s, rho_s, "o-")
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"$\rho$")
    plt.title("SUS (muVT): rho(mu)")
    savefig(os.path.join(outdir, "EOS_rho_vs_mu.png"))

    # -----------------------
    # 3) P vs mu
    # -----------------------
    plt.figure()
    plt.plot(mu_s, P_s, "o-")
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"$P$")
    plt.title("SUS (muVT): P(mu)")
    savefig(os.path.join(outdir, "EOS_P_vs_mu.png"))

    # -----------------------
    # 4) (opzionale) lnXi vs mu
    # -----------------------
    plt.figure()
    plt.plot(mu_s, lnXi_s, "o-")
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"$\ln \Xi$")
    plt.title("SUS (muVT): lnXi(mu)")
    savefig(os.path.join(outdir, "EOS_lnXi_vs_mu.png"))


def load_W_file(path):
    # Formato:
    #  # ...
    #  # columns: N  W(N)
    #  0  ...
    #  1  ...
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError(f"Formato inatteso in {path}: attese 2 colonne (N, W).")
    N = data[:, 0].astype(int)
    W = data[:, 1]
    return N, W


def plot_W_examples(run_dir, outdir, max_files=6):
    pattern = os.path.join(run_dir, "W_reweighted_mu_*.dat")
    files = sorted(glob.glob(pattern))
    if not files:
        print("[info] Nessun file W_reweighted_mu_*.dat trovato. Salto plot W(N).")
        return

    files = files[:max_files]
    plt.figure()

    for f in files:
        N, W = load_W_file(f)
        label = os.path.basename(f).replace("W_reweighted_mu_", "").replace(".dat", "")
        plt.plot(N, W, "-", label=label)

    plt.xlabel("N")
    plt.ylabel("W(N)")
    plt.title("Reweighted distributions W(N; mu) (examples)")
    plt.legend(fontsize=8)
    savefig(os.path.join(outdir, "W_reweighted_examples.png"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--run_dir",
        type=str,
        default=".",
        help="Cartella che contiene EOS_mu_rho_P.dat (e opzionalmente W_reweighted_mu_*.dat)",
    )
    ap.add_argument(
        "--eos_file",
        type=str,
        default="EOS_mu_rho_P.dat",
        help="Nome del file EOS (default: EOS_mu_rho_P.dat)",
    )
    ap.add_argument(
        "--plot_W",
        action="store_true",
        help="Se presente, cerca e plotta anche W_reweighted_mu_*.dat",
    )
    ap.add_argument(
        "--max_W_files",
        type=int,
        default=6,
        help="Numero massimo di file W_reweighted da plottare (default: 6)",
    )
    args = ap.parse_args()

    run_dir = args.run_dir
    eos_path = os.path.join(run_dir, args.eos_file)

    mu, rho, P, Nmean, lnXi = load_eos(eos_path)

    # Crea output nella stessa cartella run_dir
    outdir = run_dir

    plot_eos(mu, rho, P, Nmean, lnXi, outdir)

    if args.plot_W:
        plot_W_examples(run_dir, outdir, max_files=args.max_W_files)

    print("[done]")


if __name__ == "__main__":
    main()
