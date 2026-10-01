// ======================================================
// main_test.cpp
// Una singola simulazione MC NPT (NPTSystem) a P fissata
// Output:
//   - frames_NPT/frame_XXXX.xyz
//   - npt_observables.dat  (# sweep U V rho)
//   - npt_summary.dat      (medie su produzione)
// ======================================================

#include "parameters__npt.hpp"
#include "system__npt.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <filesystem>   // C++17

// ------------------------------------------------------
// nome file frame_0000.xyz, frame_0001.xyz, ...
// ------------------------------------------------------
static std::string make_frame_filename(int frame_id)
{
    std::ostringstream oss;
    oss << "frames_NPT/frame_"
        << std::setw(4) << std::setfill('0') << frame_id
        << ".xyz";
    return oss.str();
}

// ------------------------------------------------------
// Scrive snapshot .xyz (OVITO)
// ------------------------------------------------------
static void write_xyz(const NPTSystem& sys, const std::string& filename)
{
    const auto& part = sys.get_particles();
    const int N = static_cast<int>(part.size());
    const Box& box = sys.get_box();

    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Errore: impossibile aprire " << filename << " per scrittura\n";
        return;
    }

    out << N << "\n";
    out << "L = " << box.L << "\n";

    for (int i = 0; i < N; ++i) {
        const auto& r = part[i].r;
        out << "Ar " << r(0) << " " << r(1) << " " << r(2) << "\n";
    }
}

// ------------------------------------------------------
int main()
{
    // ==============================
    // PARAMETRI (1 run a P fissata)
    // ==============================
    NPTParams pars;
    pars.N   = 343;     // 7^3 reticolo SC pieno
    pars.rho = 0.20;    // SOLO per inizializzazione box (scelta libera)
    pars.T   = 1.4;     // progetto
    pars.P   = 0.10;    // <-- QUI scegli la pressione della singola run

    // Monte Carlo (usa i campi nuovi nel tuo parameters__npt.hpp)
    pars.delta_move   = 0.40;
    pars.delta_lnV    = 0.02;
    pars.n_sweeps_eq  = 10000;
    pars.n_sweeps_prod= 30000;
    pars.print_every  = 2000;

    // frames
    const int n_frames_desired = 10;

    // ==============================
    // PREPARA CARTELLA FRAMES
    // ==============================
    std::filesystem::create_directories("frames_NPT");

    // ==============================
    // COSTRUISCI SISTEMA
    // ==============================
    NPTSystem sys(pars);
    sys.print_info();

    const int N = pars.N;
    const long n_total = pars.n_sweeps_eq + pars.n_sweeps_prod;

    // ==============================
    // FRAME 0 (config iniziale)
    // ==============================
    int frame_id = 0;
    write_xyz(sys, make_frame_filename(frame_id));
    frame_id++;

    // stride per salvare ~n_frames_desired totali (incluso frame 0)
    int frame_stride = 1;
    if (n_frames_desired > 1) {
        frame_stride = static_cast<int>(n_total / (n_frames_desired - 1));
        if (frame_stride < 1) frame_stride = 1;
    }
    long next_frame_sweep = frame_stride; // salva a sweep=frame_stride, 2*stride, ...

    // ==============================
    // FILE OSSERVABILI
    // ==============================
    std::ofstream obs("data/npt_observables.dat");
    if (!obs) {
        std::cerr << "Errore: impossibile aprire npt_observables.dat\n";
        return 1;
    }
    obs << "# sweep   U          V          rho\n";
    obs << std::scientific << std::setprecision(10);

    // ==============================
    // ACCETTANZE (manuali, come main vecchio)
    // ==============================
    long long acc_moves = 0, tot_moves = 0;
    long long acc_vol   = 0, tot_vol   = 0;

    // ==============================
    // MEDIE SU PRODUZIONE (semplici)
    // ==============================
    long long n_prod_samples = 0;
    double sum_rho = 0.0;
    double sum_Un  = 0.0;
    double sum_V   = 0.0;

    std::cout << "\n--- Inizio MC NPT (single P) ---\n";
    std::cout << "N=" << pars.N << "  T=" << pars.T << "  P=" << pars.P << "  rho_init=" << pars.rho << "\n";
    std::cout << "delta_move=" << pars.delta_move << "  delta_lnV=" << pars.delta_lnV << "\n";
    std::cout << "eq=" << pars.n_sweeps_eq << "  prod=" << pars.n_sweeps_prod << "  total=" << n_total << "\n\n";

    // ==============================
    // LOOP MC
    // ==============================
    for (long sweep = 1; sweep <= n_total; ++sweep) {

        // 1 sweep = N mosse particella + 1 mossa volume
        for (int m = 0; m < N; ++m) {
            if (sys.mc_move_particle(pars.delta_move)) acc_moves++;
            tot_moves++;
        }

        if (sys.mc_move_volume(pars.delta_lnV)) acc_vol++;
        tot_vol++;

        // osservabili istantanei
        const double U   = sys.get_energy();
        const double V   = sys.get_box().volume();
        const double rho = static_cast<double>(N) / V;

        obs << sweep << "  " << U << "  " << V << "  " << rho << "\n";

        // medie solo in produzione
        if (sweep > pars.n_sweeps_eq) {
            sum_rho += rho;
            sum_Un  += U / static_cast<double>(N);
            sum_V   += V;
            n_prod_samples++;
        }

        // frame xyz
        if (sweep == next_frame_sweep && frame_id < n_frames_desired) {
            write_xyz(sys, make_frame_filename(frame_id));
            frame_id++;
            next_frame_sweep += frame_stride;
        }

        // stampa
        if (pars.print_every > 0 && (sweep % pars.print_every == 0)) {
            const double acc_m = (tot_moves > 0) ? double(acc_moves) / double(tot_moves) : 0.0;
            const double acc_v = (tot_vol   > 0) ? double(acc_vol)   / double(tot_vol)   : 0.0;

            std::cout << "sweep " << sweep << "/" << n_total
                      << "  rho=" << std::fixed << std::setprecision(6) << rho
                      << "  U/N=" << (U / N)
                      << "  V=" << V
                      << "  acc(move)=" << std::setprecision(4) << acc_m
                      << "  acc(vol)="  << std::setprecision(4) << acc_v
                      << std::scientific << std::setprecision(10)
                      << "\n";
        }

        if (sweep == pars.n_sweeps_eq) {
            std::cout << "--- Fine equilibrazione, inizio produzione ---\n";
        }
    }

    // ==============================
    // SUMMARY (medie semplici su produzione)
    // ==============================
    const double acc_m_tot = (tot_moves > 0) ? double(acc_moves) / double(tot_moves) : 0.0;
    const double acc_v_tot = (tot_vol   > 0) ? double(acc_vol)   / double(tot_vol)   : 0.0;

    double rho_mean = 0.0, Un_mean = 0.0, V_mean = 0.0;
    if (n_prod_samples > 0) {
        rho_mean = sum_rho / double(n_prod_samples);
        Un_mean  = sum_Un  / double(n_prod_samples);
        V_mean   = sum_V   / double(n_prod_samples);
    }

    std::ofstream sum("data/npt_summary.dat");
    if (!sum) {
        std::cerr << "Errore: impossibile aprire npt_summary.dat\n";
        return 1;
    }

    sum << "# Single NPT run summary (simple production averages)\n";
    sum << "# N=" << pars.N << "  T=" << pars.T << "  P=" << pars.P << "  rho_init=" << pars.rho << "\n";
    sum << "# eq_sweeps=" << pars.n_sweeps_eq << "  prod_sweeps=" << pars.n_sweeps_prod
        << "  prod_samples=" << n_prod_samples << "\n";
    sum << "# delta_move=" << pars.delta_move << "  delta_lnV=" << pars.delta_lnV << "\n";
    sum << "# acc_move=" << acc_m_tot << "  acc_vol=" << acc_v_tot << "\n";
    sum << "# <rho>  <V>  <U/N>\n";
    sum << std::scientific << std::setprecision(10)
        << rho_mean << "  " << V_mean << "  " << Un_mean << "\n";

    std::cout << "\n=== Fine MC NPT (single P) ===\n";
    std::cout << "P=" << pars.P
              << "  <rho>=" << rho_mean
              << "  <U/N>=" << Un_mean
              << "  acc(move)=" << acc_m_tot
              << "  acc(vol)="  << acc_v_tot
              << "  prod_samples=" << n_prod_samples
              << "\n";
    std::cout << "Wrote: npt_observables.dat, npt_summary.dat, frames_NPT/\n";

    return 0;
}
