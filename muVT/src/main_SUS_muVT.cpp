// ======================================================
//  main_SUS_muVT.cpp — Successive Umbrella Sampling (µVT)
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <filesystem>
#include <algorithm>

#include "parameters_muVT.hpp"
#include "system_muVT.hpp"

namespace fs = std::filesystem;

// ------------------------------------------------------
// Struct risultato finestra SUS [n,n+1]
// ------------------------------------------------------
struct SUSWindowResult {
    int    n = 0;
    long   Hn = 0;
    long   Hn1 = 0;
    double ratio = 0.0;   // Hn1/Hn
};

// ------------------------------------------------------
// Utility: double -> tag file-safe (m3p50, p1p40, etc.)
// ------------------------------------------------------
static std::string to_tag(double x)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(2) << x;
    std::string s = os.str();
    for (char& c : s) {
        if (c == '.') c = 'p';
        if (c == '-') c = 'm';
    }
    return s;
}

// ------------------------------------------------------
// Utility: crea directory se non esiste
// ------------------------------------------------------
static void ensure_dir(const fs::path& p)
{
    std::error_code ec;
    if (fs::exists(p, ec)) return;
    if (!fs::create_directories(p, ec)) {
        throw std::runtime_error("Impossibile creare directory: " + p.string());
    }
}

// ------------------------------------------------------
// Esegue UNA finestra SUS con vincolo N ∈ {n, n+1}
// Usa le tue probabilità pars.p_tra / p_ins / p_del e deltra
//
// Nota: per evitare di dover toccare rng globale o firme strane,
// qui facciamo un seed diverso per finestra: pars.seed + 100000*n
// ------------------------------------------------------
static SUSWindowResult run_SUS_window_singlecpp(const MuVTParams& base_pars,
                                               int n,
                                               long n_steps,
                                               long n_equil,
                                               const fs::path& window_dir,
                                               bool verbose)
{
    MuVTParams pars = base_pars;

    // finestra [n,n+1]
    pars.N_init = n;

    // seed diverso per finestra (riproducibile)
    if (pars.seed >= 0) pars.seed = pars.seed + 100000 * n;

    MuVTSystem sys(pars);

    // log finestra minimale
    std::ofstream flog(window_dir / "window_log.txt", std::ios::out | std::ios::trunc);
    flog << "# SUS window n=" << n << " (N in {" << n << "," << (n+1) << "})\n";
    flog << "# steps=" << n_steps << " equil=" << n_equil << "\n";

    long Hn  = 0;
    long Hn1 = 0;

    for (long step = 0; step < n_steps; ++step) {

        // Scegli tipo mossa (come mc_step, ma con vincolo finestra)
        double r = rng.ranf();

        if (r < pars.p_tra) {
            // traslazione sempre permessa
            sys.mc_move_particle(pars.deltra);
        } else if (r < pars.p_tra + pars.p_ins) {
            // inserzione: permessa solo se Ncur == n
            int Ncur = sys.get_N();
            if (Ncur < n + 1) {
                sys.mc_insert_particle();
            } else {
                // rifiuto a priori
            }
        } else {
            // cancellazione: permessa solo se Ncur == n+1
            int Ncur = sys.get_N();
            if (Ncur > n) {
                sys.mc_delete_particle();
            } else {
                // rifiuto a priori
            }
        }

        // Misure dopo equil
        if (step >= n_equil) {
            int Ncur = sys.get_N();
            if (Ncur == n)      ++Hn;
            else if (Ncur == n+1) ++Hn1;
            else {
                // non dovrebbe accadere
                flog << "WARNING step=" << step << " N=" << Ncur << " fuori finestra\n";
            }
        }

        if (verbose && step>0 && (step % 50000 == 0)) {
            std::cout << "  [n=" << n << "] step " << step
                      << "  N=" << sys.get_N()
                      << "  Hn=" << Hn << " Hn1=" << Hn1 << "\n";
        }
    }

    double ratio = 0.0;
    if (Hn > 0 && Hn1 > 0) ratio = double(Hn1) / double(Hn);

    // salva un mini summary
    std::ofstream fsu(window_dir / "window_summary.dat", std::ios::out | std::ios::trunc);
    fsu << "# n  Hn  Hn1  ratio\n";
    fsu << n << " " << Hn << " " << Hn1 << " " << std::setprecision(15) << ratio << "\n";

    SUSWindowResult res;
    res.n = n;
    res.Hn = Hn;
    res.Hn1 = Hn1;
    res.ratio = ratio;
    return res;
}

// ------------------------------------------------------
// Ricostruisce logW da ratios e salva logW e W_norm
// ------------------------------------------------------
static void build_logW_and_save(const fs::path& run_dir,
                                const std::vector<SUSWindowResult>& res,
                                int nmax)
{
    std::vector<double> logW(nmax + 1, 0.0); // logW(0)=0 gauge

    std::ofstream bad(run_dir / "bad_windows.dat", std::ios::out | std::ios::trunc);
    bad << "# n  Hn  Hn1  ratio  NOTE\n";

    for (int n = 0; n < nmax; ++n) {
        const auto& r = res[n];
        if (r.Hn <= 0 || r.Hn1 <= 0 || !(r.ratio > 0.0)) {
            bad << n << " " << r.Hn << " " << r.Hn1 << " "
                << std::setprecision(15) << r.ratio << "  BAD_STATS\n";
            logW[n+1] = logW[n]; // placeholder
        } else {
            logW[n+1] = logW[n] + std::log(r.ratio);
        }
    }

    // salva logW
    {
        std::ofstream f(run_dir / "logW_muRef.dat", std::ios::out | std::ios::trunc);
        f << "# N  logW(N)   (gauge: logW(0)=0)\n";
        for (int N = 0; N <= nmax; ++N) {
            f << N << " " << std::setprecision(15) << logW[N] << "\n";
        }
    }

    // normalizza W
    double maxlog = *std::max_element(logW.begin(), logW.end());
    std::vector<double> W(nmax + 1, 0.0);

    double Z = 0.0;
    for (int N = 0; N <= nmax; ++N) {
        W[N] = std::exp(logW[N] - maxlog);
        Z += W[N];
    }
    if (Z > 0.0) {
        for (int N = 0; N <= nmax; ++N) W[N] /= Z;
    }

    {
        std::ofstream f(run_dir / "W_muRef_norm.dat", std::ios::out | std::ios::trunc);
        f << "# N  W_muRef_norm(N)   (sum=1)\n";
        for (int N = 0; N <= nmax; ++N) {
            f << N << " " << std::setprecision(15) << W[N] << "\n";
        }
    }
}

// ------------------------------------------------------
// Scrive meta
// ------------------------------------------------------
static void write_meta(const fs::path& run_dir,
                       const MuVTParams& pars,
                       int nmax,
                       long n_steps,
                       long n_equil)
{
    std::ofstream f(run_dir / "meta.txt", std::ios::out | std::ios::trunc);
    if (!f) throw std::runtime_error("Impossibile scrivere meta.txt");

    f << "RUN META (SUS µVT)\n";
    f << "------------------\n";
    f << "T       = " << pars.T << "\n";
    f << "mu_ref  = " << pars.mu << "\n";
    f << "beta    = " << pars.beta() << "\n";
    f << "z_ref   = " << pars.z << "\n";
    f << "sigma   = " << pars.sigma << "\n";
    f << "epsilon = " << pars.epsilon << "\n";
    f << "rc      = " << pars.rc << "\n";
    f << "V_box   = " << pars.V_box << "\n";
    f << "L       = " << pars.L(0) << "\n";
    f << "V=L^3   = " << (pars.L(0)*pars.L(0)*pars.L(0)) << "\n";
    f << "\nMove probs:\n";
    f << "p_tra   = " << pars.p_tra << "\n";
    f << "p_ins   = " << pars.p_ins << "\n";
    f << "p_del   = " << pars.p_del << "\n";
    f << "deltra  = " << pars.deltra << "\n";
    f << "\nSUS:\n";
    f << "nmax           = " << nmax << "\n";
    f << "steps/window   = " << n_steps << "\n";
    f << "equil/window   = " << n_equil << "\n";
    f << "seed(base)     = " << pars.seed << " (se >=0, shift per finestra)\n";
}

// ------------------------------------------------------
// MAIN
// ------------------------------------------------------
int main(int argc, char* argv[])
{
    try {
        fs::path out_root = "data/SUS";
        if (argc >= 2) out_root = fs::path(argv[1]);

        MuVTParams pars;

        // ---- thermo/LJ
        pars.T  = 1.4;
        pars.mu = -3.5;      // mu_ref scelto dallo scan
        pars.epsilon = 1.0;
        pars.sigma   = 1.0;
        pars.rc      = 2.5;

        // ---- box
        pars.V_box    = 343.0;  // L=7
        pars.rho_init = 0.3;
        pars.N_init   = 0;      // sovrascritto a n
        pars.set_box_from_init();

        // ---- moves
        pars.p_tra = 0.50;
        pars.p_ins = 0.25;
        pars.p_del = 0.25;
        pars.normalize_move_probs();
        pars.deltra = 0.2;

        // ---- RNG
        pars.seed = 12345; // >=0 riproducibile; <0 random
        pars.update_activity();

        // ---- SUS setup
        const double V = pars.V_box > 0.0 ? pars.V_box : pars.volume();
        const double rho_max = 130.0 / V;   // ~0.379
        const int nmax = int(std::floor(rho_max * V + 1e-12)); // ~102 con V=343

        const long n_steps = 300'000;
        const long n_equil = 50'000;

        // ---- run_id
        const std::string run_id =
            "muVT_T" + to_tag(pars.T) +
            "_V" + std::to_string(int(std::llround(V))) +
            "_muRef_" + to_tag(pars.mu) +
            "_Nmax" + std::to_string(nmax);

        fs::path run_dir = out_root / run_id;
        fs::path win_dir = run_dir / "windows";

        ensure_dir(win_dir);

        write_meta(run_dir, pars, nmax, n_steps, n_equil);

        // ratios globale
        std::ofstream fr(run_dir / "ratios.dat", std::ios::out | std::ios::trunc);
        fr << "# n  Hn  Hn1  ratio=Hn1/Hn\n";

        std::cout << "=== SUS µVT (single .cpp) ===\n";
        std::cout << "Run dir: " << run_dir.string() << "\n";
        std::cout << "T=" << pars.T << " mu_ref=" << pars.mu << " V=" << V
                  << " rho_max=" << rho_max << " => nmax=" << nmax << "\n";
        std::cout << "steps/window=" << n_steps << " equil=" << n_equil << "\n\n";

        std::vector<SUSWindowResult> results;
        results.reserve(nmax);

        for (int n = 0; n < nmax; ++n) {
            // cartella finestra
            std::ostringstream wname;
            wname << "window_" << std::setw(4) << std::setfill('0') << n;
            fs::path wdir = win_dir / wname.str();
            ensure_dir(wdir);

            bool verbose = false;
            SUSWindowResult r = run_SUS_window_singlecpp(pars, n, n_steps, n_equil, wdir, verbose);

            results.push_back(r);
            fr << r.n << " " << r.Hn << " " << r.Hn1 << " "
               << std::setprecision(15) << r.ratio << "\n";

            if (n % 10 == 0) {
                std::cout << "[done] n=" << n << " Hn=" << r.Hn << " Hn1=" << r.Hn1
                          << " ratio=" << std::setprecision(6) << r.ratio << "\n";
            }
        }
        fr.close();

        build_logW_and_save(run_dir, results, nmax);

        std::cout << "\nSUS completato.\n";
        std::cout << "File principali:\n";
        std::cout << "  " << (run_dir / "meta.txt").string() << "\n";
        std::cout << "  " << (run_dir / "ratios.dat").string() << "\n";
        std::cout << "  " << (run_dir / "logW_muRef.dat").string() << "\n";
        std::cout << "  " << (run_dir / "W_muRef_norm.dat").string() << "\n";
        std::cout << "  " << (run_dir / "bad_windows.dat").string() << "\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
