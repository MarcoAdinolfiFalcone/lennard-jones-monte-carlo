import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

TSTAR = 1.4
EOS_FILE = Path("data/EOS_NPT_T1p4.dat")
OUT_FIG  = Path("data/EOS_NPT_T1p4_with_ideal.png")

def load_eos(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Non trovo {path}. Hai creato prima il file EOS?")
    # colonne: P rho U/N V L (con header commentato #)
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise RuntimeError(f"{path} non ha abbastanza colonne. Attese almeno P e rho.")
    P   = data[:, 0]
    rho = data[:, 1]
    return rho, P

def main():
    rho, P = load_eos(EOS_FILE)

    # Ordina per rho (solo per la linea ideale e per limiti assi puliti)
    idx = np.argsort(rho)
    rho = rho[idx]
    P = P[idx]

    # Linea gas ideale: P* = T* rho*
    rho_line = np.linspace(rho.min(), rho.max(), 300)
    P_ideal  = TSTAR * rho_line

    plt.figure()
    # tuoi punti: SOLO scatter (nessuna interpolazione/linea tra punti)
    plt.scatter(rho, P, label="NPT MC (punti)", s=40)

    # riferimento gas ideale
    plt.plot(rho_line, P_ideal, linestyle="--", label=rf"Gas ideale: $P^*=T^*\rho^*$ (T*={TSTAR})")

    plt.xlabel(r"$\rho^* = N/V$")
    plt.ylabel(r"$P^*$")
    plt.title(rf"EOS NPT (T*={TSTAR}) vs gas ideale")
    plt.grid(True)
    plt.legend()

    plt.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    print(f"[OK] salvato: {OUT_FIG}")

    # Se vuoi anche visualizzarlo a schermo, scommenta:
    # plt.show()

if __name__ == "__main__":
    main()
