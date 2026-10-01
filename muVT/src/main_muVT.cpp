// ======================================================
//  main_muVT.cpp  —  Scan su mu per scegliere mu_ref
//  Output ordinato in:
//    data/<run_name>/mu_<tag>/
//      - muvt_observables.dat
//      - histN_raw.dat
//      - histN_norm.dat
//      - summary.txt
//    + data/<run_name>/scan_summary.dat
//
//  NOTE:
//   - usa MuVTSystem::run_simulation(out_prefix, histN)
//   - out_prefix viene passato come "data/.../mu_xxx/out"
//     così i file scritti dal system finiscono già nella cartella giusta.
// ======================================================

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>
#include <cmath>
#include <filesystem>

#include "parameters_muVT.hpp"
#include "system_muVT.hpp"

namespace fs = std::filesystem;

// ------------------------------------------------------
// Utility: trasforma mu in tag file-safe
//   -3.50 -> m3p50
// ------------------------------------------------------
static std::string mu_to_tag(double mu)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(2) << mu;
    std::string s = os.str();
    for (char& c : s) {
        if (c == '.') c = 'p';
        if (c == '-') c = 'm';
    }
    return s;
}

// ------------------------------------------------------
// Crea cartella (e genitori) se non esiste
// ------------------------------------------------------
static void ensure_dir_or_throw(const fs::path& p)
{
    std::error_code ec;
    if (fs::exists(p, ec)) return;
    if (!fs::create_directories(p, ec)) {
        throw std::runtime_error("Impossibile creare directory: " + p.string());
    }
}

// ------------------------------------------------------
// Scrive istogramma grezzo H(N)
// ------------------------------------------------------
static void write_histogram_raw(const fs::path& fname,
                                const std::vector<long>& hist_N)
{
    std::ofstream f(fname, std::ios::out | std::ios::trunc);
    if (!f) {
        std::cerr << "Errore: impossibile aprire " << fname << " per scrittura.\n";
        return;
    }
    f << "# N  H(N)\n";
    for (std::size_t n = 0; n < hist_N.size(); ++n)
        f << n << " " << hist_N[n] << "\n";
}

// ------------------------------------------------------
// Scrive istogramma normalizzato P(N)
// ------------------------------------------------------
static void write_histogram_norm(const fs::path& fname,
                                 const std::vector<long>& hist_N,
                                 double V)
{
    std::ofstream f(fname, std::ios::out | std::ios::trunc);
    if (!f) {
        std::cerr << "Errore: impossibile aprire " << fname << " per scrittura.\n";
        return;
    }

    long total = std::accumulate(hist_N.begin(), hist_N.end(), 0L);

    f << "# N   counts   P(N)   rho=N/V\n";
    for (std::size_t N = 0; N < hist_N.size(); ++N) {
        long c = hist_N[N];
        double P   = (total > 0) ? (double(c) / double(total)) : 0.0;
        double rho = (V > 0.0) ? (double(N) / V) : 0.0;
        f << N << "  " << c << "  " << std::setprecision(12) << P << "  " << rho << "\n";
    }
}

// ------------------------------------------------------
// Momenti da H(N): <N>, <N^2>, Var(N), <rho>
// ------------------------------------------------------
static void moments_from_hist(const std::vector<long>& hist_N,
                              double V,
                              double& Nmean,
                              double& N2mean,
                              double& Nvar,
                              double& rhomean)
{
    long total = std::accumulate(hist_N.begin(), hist_N.end(), 0L);
    if (total <= 0) {
        Nmean = N2mean = Nvar = rhomean = 0.0;
        return;
    }

    double sumN  = 0.0;
    double sumN2 = 0.0;
    for (std::size_t n = 0; n < hist_N.size(); ++n) {
        const double c = double(hist_N[n]);
        sumN  += c * double(n);
        sumN2 += c * double(n) * double(n);
    }

    Nmean  = sumN  / double(total);
    N2mean = sumN2 / double(total);
    Nvar   = N2mean - Nmean*Nmean;
    rhomean = (V > 0.0) ? (Nmean / V) : 0.0;
}

int main(int argc, char* argv[])
{
    // Nome run: data/<run_name>/
    // default: scan_muVT
    std::string run_name = "scan_muVT";
    if (argc >= 2) run_name = argv[1];

    // Base dir
    fs::path base_dir = fs::path("data") / run_name;

    try {
        ensure_dir_or_throw(base_dir);
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << "\n";
        return 1;
    }

    // -----------------------------
    // Parametri base
    // -----------------------------
    MuVTParams pars;

    // Box: fisso volume (cubo). Esempio: V=343 -> L=7
    pars.V_box   = 343.0;
    pars.rho_init = 0.3;  // ignorato se V_box>0
    pars.N_init  = 50;
    pars.set_box_from_init();

    // Thermo
    pars.T  = 1.4;
    pars.rc = 2.5;

    // MC
    pars.totsteps    = 400'000;
    pars.eqstps      = 100'000;
    pars.savemeasure = 50;
    pars.outstps     = 10'000;
    pars.deltra      = 0.2;

    pars.p_tra = 0.50;
    pars.p_ins = 0.25;
    pars.p_del = 0.25;
    pars.normalize_move_probs();

    // RNG
    pars.seed = 12345;

    // Lista mu da testare
    std::vector<double> mu_list = { -6.0, -5.0, -4.0, -3.5, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0 };

    // File summary globale
    {
        std::ofstream gsum(base_dir / "scan_summary.dat", std::ios::out | std::ios::trunc);
        gsum << "# mu   <N>   <rho>   Var(N)   V\n";
    }

    std::cout << "=== Scan muVT (LJ) per scegliere mu_ref ===\n";
    std::cout << "Output base: " << base_dir << "\n";
    std::cout << "T=" << pars.T << " rc=" << pars.rc << " V_box=" << pars.V_box
              << " => L=" << std::cbrt(pars.V_box) << "\n";
    std::cout << "totsteps=" << pars.totsteps << " eqstps=" << pars.eqstps
              << " savemeasure=" << pars.savemeasure << " outstps=" << pars.outstps << "\n";
    std::cout << "p_tra=" << pars.p_tra << " p_ins=" << pars.p_ins << " p_del=" << pars.p_del << "\n\n";

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Target: rho in [0.05,0.30] -> N in [0.05 V, 0.30 V] con sigma=1\n\n";

    for (double mu : mu_list)
    {
        MuVTParams p = pars;
        p.mu = mu;
        p.update_activity();

        std::string tag = mu_to_tag(mu);
        fs::path mu_dir = base_dir / ("mu_" + tag);

        try {
            ensure_dir_or_throw(mu_dir);
        } catch (const std::exception& e) {
            std::cerr << "ERRORE: " << e.what() << "\n";
            return 1;
        }

        // out_prefix: i file prodotti dal system finiscono direttamente qui
        std::string out_prefix = (mu_dir / "out").string();

        std::cout << "----------------------------------------\n";
        std::cout << "Run mu=" << mu << " -> " << mu_dir << "\n";

        MuVTSystem sys(p);
        sys.print_info();

        std::vector<long> hist_N;
        sys.run_simulation(out_prefix, hist_N);

        // Il system scrive già:
        //   <out_prefix>_muvt_observables.dat
        //   <out_prefix>_histN.dat  (grezzo)
        // Noi creiamo copie/nomi più chiari e aggiungiamo norm+summary.

        // Copia/riscrivi H(N) grezzo
        write_histogram_raw(mu_dir / "histN_raw.dat", hist_N);

        // P(N)
        write_histogram_norm(mu_dir / "histN_norm.dat", hist_N, sys.get_V());

        // Momenti
        double Nmean=0, N2mean=0, Nvar=0, rhomean=0;
        moments_from_hist(hist_N, sys.get_V(), Nmean, N2mean, Nvar, rhomean);

        // Summary per questo mu
        {
            std::ofstream f(mu_dir / "summary.txt", std::ios::out | std::ios::trunc);
            f << std::setprecision(12);
            f << "mu = " << mu << "\n";
            f << "T  = " << p.T << "\n";
            f << "V  = " << sys.get_V() << "\n";
            f << "L  = " << sys.get_L() << "\n";
            f << "<N>    = " << Nmean << "\n";
            f << "Var(N) = " << Nvar << "\n";
            f << "<rho>  = " << rhomean << "\n";
            f << "files:\n";
            f << "  out_muvt_observables.dat (grezzo dal system)\n";
            f << "  out_histN.dat            (grezzo dal system)\n";
            f << "  muvt_observables.dat      (alias/copia se vuoi farla a mano)\n";
            f << "  histN_raw.dat\n";
            f << "  histN_norm.dat\n";
        }

        // Summary globale
        {
            std::ofstream gsum(base_dir / "scan_summary.dat", std::ios::out | std::ios::app);
            gsum << std::fixed << std::setprecision(6)
                 << mu << "  "
                 << Nmean << "  "
                 << rhomean << "  "
                 << Nvar << "  "
                 << sys.get_V() << "\n";
        }

        std::cout << ">>> <N>=" << std::setprecision(3) << Nmean
                  << "  <rho>=" << rhomean
                  << "  Var(N)=" << Nvar << "\n";
    }

    std::cout << "\nScan completata.\n";
    std::cout << "Output in: " << base_dir << "\n";
    std::cout << "Riepilogo: " << (base_dir / "scan_summary.dat") << "\n";
    return 0;
}
