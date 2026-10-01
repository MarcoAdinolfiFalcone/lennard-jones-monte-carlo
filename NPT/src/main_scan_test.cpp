// ======================================================
// main_EOS_NPT.cpp
//
// Esegue una serie di simulazioni MC NPT (T fissata) per una lista di pressioni,
// e salva tutto in modo ordinato dentro ./data/
//
// Richieste utente:
//  - obs_every = 1 (scrivi osservabili ad ogni sweep)
//  - box_every = 1 (tenta box move ad ogni sweep)
//  - salva anche L in observables
//  - per ogni pressione crea cartella dedicata con .dat + frames
//
// Output:
//  - data/EOS_NPT_T1p4.dat          (tabella P  <rho> <L> <V> <U/N> acc_move acc_vol ...)
//  - data/P0p07/ ...               (files per singola pressione)
//  - data/P0p08/ ...
//  - ...
// ======================================================

#include "parameters__npt.hpp"
#include "system__npt.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <filesystem>

// ------------------------------
// helper: formatta P in stile "P0p07" (0.07 -> P0p07)
// ------------------------------
static std::string pressure_tag(double P)
{
    // P tra 0.00 e 9.99 circa: convertiamo in centesimi (2 decimali)
    // 0.07 -> 7, 0.30 -> 30
    int ip = static_cast<int>(std::round(P * 100.0));
    std::ostringstream oss;
    oss << "P0p" << std::setw(2) << std::setfill('0') << ip;
    return oss.str();
}

// ------------------------------
// helper: frame filename
// ------------------------------
static std::string make_frame_filename(const std::string& frames_dir, int frame_id)
{
    std::ostringstream oss;
    oss << frames_dir << "/frame_"
        << std::setw(4) << std::setfill('0') << frame_id
        << ".xyz";
    return oss.str();
}

// ------------------------------
// write xyz snapshot (OVITO)
// ------------------------------
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
    out << "L = " << std::setprecision(12) << box.L << "\n";

    out << std::scientific << std::setprecision(12);
    for (int i = 0; i < N; ++i) {
        const auto& r = part[i].r;
        out << "Ar " << r(0) << " " << r(1) << " " << r(2) << "\n";
    }
}

int main()
{
    // ==============================
    // Lista pressioni richiesta
    // ==============================
    const std::vector<double> pressures = {
        0.07, 0.08, 0.09, 0.10, 0.12, 0.14, 0.16, 0.18,
        0.20, 0.22, 0.24, 0.26, 0.28, 0.30
    };

    // ==============================
    // Parametri comuni (tutti i run)
    // ==============================
    NPTParams pars;
    pars.N   = 343;
    pars.rho = 0.20;   // solo inizializzazione
    pars.T   = 1.4;

    pars.delta_move     = 0.40;
    pars.delta_lnV      = 0.02;
    pars.n_sweeps_eq    = 10000;
    pars.n_sweeps_prod  = 30000;
    pars.print_every    = 2000;

    // Richieste: obs_every=1, box_every=1
    const int obs_every = 1;
    const int box_every = 1;

    // frames per run
    const int n_frames_desired = 10;

    // ==============================
    // Cartella base output
    // ==============================
    const std::string base_dir = "data";
    std::filesystem::create_directories(base_dir);

    // file EOS (appendiamo una riga per P)
    const std::string eos_file = base_dir + "/EOS_NPT_T1p4.dat";
    std::ofstream eos(eos_file, std::ios::out | std::ios::trunc);
    if (!eos) {
        std::cerr << "Errore: impossibile aprire " << eos_file << "\n";
        return 1;
    }

    eos << "# EOS NPT (MC) Lennard-Jones reduced units\n";
    eos << "# T = " << pars.T << "   N = " << pars.N << "   rho_init = " << pars.rho << "\n";
    eos << "# eq_sweeps = " << pars.n_sweeps_eq << "   prod_sweeps = " << pars.n_sweeps_prod << "\n";
    eos << "# delta_move = " << pars.delta_move << "   delta_lnV = " << pars.delta_lnV << "\n";
    eos << "# columns: P  <rho>  <L>  <V>  <U/N>  acc_move  acc_vol  prod_samples\n";
    eos << std::scientific << std::setprecision(12);

    // ==============================
    // Loop sulle pressioni
    // ==============================
    for (double P : pressures) {

        pars.P = P;
        const std::string tag = pressure_tag(P);

        // cartella run: data/P0p07/
        const std::string run_dir    = base_dir + "/" + tag;
        const std::string frames_dir = run_dir + "/frames";
        std::filesystem::create_directories(frames_dir);

        // filenames per questo run
        const std::string obs_file = run_dir + "/npt_observables_" + tag + ".dat";
        const std::string sum_file = run_dir + "/npt_summary_"     + tag + ".dat";

        std::cout << "\n====================================================\n";
        std::cout << "RUN NPT: " << tag << "   (P=" << P << ", T=" << pars.T << ", N=" << pars.N << ")\n";
        std::cout << "Output in: " << run_dir << "\n";
        std::cout << "====================================================\n";

        // costruisci sistema (nuovo per ogni P)
        NPTSystem sys(pars);
        sys.print_info();

        const long n_total = pars.n_sweeps_eq + pars.n_sweeps_prod;

        // frame 0
        int frame_id = 0;
        write_xyz(sys, make_frame_filename(frames_dir, frame_id));
        frame_id++;

        long frame_stride = 1;
        if (n_frames_desired > 1) {
            frame_stride = n_total / (n_frames_desired - 1);
            if (frame_stride < 1) frame_stride = 1;
        }
        long next_frame_sweep = frame_stride;

        // apri file osservabili
        std::ofstream obs(obs_file, std::ios::out | std::ios::trunc);
        if (!obs) {
            std::cerr << "Errore: impossibile aprire " << obs_file << "\n";
            return 1;
        }

        obs << "# sweep  U_tot  V  L  rho  acc_move  acc_vol\n";
        obs << std::scientific << std::setprecision(12);

        // contatori acceptance
        long long acc_moves = 0, tot_moves = 0;
        long long acc_vol   = 0, tot_vol   = 0;

        // medie su produzione
        long long n_prod_samples = 0;
        double sum_rho = 0.0;
        double sum_V   = 0.0;
        double sum_L   = 0.0;
        double sum_Un  = 0.0;

        // loop MC
        for (long sweep = 1; sweep <= n_total; ++sweep) {

            // N mosse particella
            for (int m = 0; m < pars.N; ++m) {
                if (sys.mc_move_particle(pars.delta_move)) acc_moves++;
                tot_moves++;
            }

            // box move (box_every=1 -> sempre)
            bool did_box = false;
            if (box_every > 0 && (sweep % box_every == 0)) {
                did_box = true;
                if (sys.mc_move_volume(pars.delta_lnV)) acc_vol++;
                tot_vol++;
            }

            // osservabili istantanei
            const double U   = sys.get_energy();
            const double V   = sys.get_box().volume();
            const double L   = sys.get_box().L;
            const double rho = static_cast<double>(pars.N) / V;

            const double acc_m = (tot_moves > 0) ? double(acc_moves) / double(tot_moves) : 0.0;
            const double acc_v = (tot_vol   > 0) ? double(acc_vol)   / double(tot_vol)   : 0.0;

            // obs_every=1 -> sempre
            if (obs_every > 0 && (sweep % obs_every == 0)) {
                obs << sweep << "  " << U << "  " << V << "  " << L << "  " << rho
                    << "  " << acc_m << "  " << acc_v << "\n";
            }

            // medie (solo produzione)
            if (sweep > pars.n_sweeps_eq) {
                sum_rho += rho;
                sum_V   += V;
                sum_L   += L;
                sum_Un  += U / double(pars.N);
                n_prod_samples++;
            }

            // frames
            if (sweep == next_frame_sweep && frame_id < n_frames_desired) {
                write_xyz(sys, make_frame_filename(frames_dir, frame_id));
                frame_id++;
                next_frame_sweep += frame_stride;
            }

            // stampa su schermo
            if (pars.print_every > 0 && (sweep % pars.print_every == 0)) {
                std::cout << std::fixed << std::setprecision(6);
                std::cout << "sweep " << sweep << "/" << n_total
                          << "  rho=" << rho
                          << "  L="   << L
                          << "  V="   << V
                          << "  U/N=" << (U / double(pars.N))
                          << "  acc(move)=" << std::setprecision(4) << acc_m
                          << "  acc(vol)="  << std::setprecision(4) << acc_v;
                if (!did_box) std::cout << "  [no box]";
                std::cout << "\n";
                std::cout << std::scientific << std::setprecision(12);
            }

            if (sweep == pars.n_sweeps_eq) {
                std::cout << "--- Fine equilibrazione, inizio produzione ---\n";
            }
        }

        // medie finali
        const double acc_m_tot = (tot_moves > 0) ? double(acc_moves) / double(tot_moves) : 0.0;
        const double acc_v_tot = (tot_vol   > 0) ? double(acc_vol)   / double(tot_vol)   : 0.0;

        double rho_mean = 0.0, V_mean = 0.0, L_mean = 0.0, Un_mean = 0.0;
        if (n_prod_samples > 0) {
            rho_mean = sum_rho / double(n_prod_samples);
            V_mean   = sum_V   / double(n_prod_samples);
            L_mean   = sum_L   / double(n_prod_samples);
            Un_mean  = sum_Un  / double(n_prod_samples);
        }

        // summary file per questa P
        std::ofstream sum(sum_file, std::ios::out | std::ios::trunc);
        if (!sum) {
            std::cerr << "Errore: impossibile aprire " << sum_file << "\n";
            return 1;
        }

        sum << "# Single NPT run summary (production averages)\n";
        sum << "# " << tag << "  P=" << pars.P << "  T=" << pars.T << "  N=" << pars.N
            << "  rho_init=" << pars.rho << "\n";
        sum << "# eq_sweeps=" << pars.n_sweeps_eq
            << "  prod_sweeps=" << pars.n_sweeps_prod
            << "  prod_samples=" << n_prod_samples << "\n";
        sum << "# delta_move=" << pars.delta_move
            << "  delta_lnV=" << pars.delta_lnV
            << "  obs_every=" << obs_every
            << "  box_every=" << box_every << "\n";
        sum << "# acc_move=" << acc_m_tot << "  acc_vol=" << acc_v_tot << "\n";
        sum << "# <rho>  <V>  <L>  <U/N>\n";
        sum << std::scientific << std::setprecision(12)
            << rho_mean << "  " << V_mean << "  " << L_mean << "  " << Un_mean << "\n";

        // append a EOS
        eos << pars.P << "  " << rho_mean << "  " << L_mean << "  " << V_mean << "  "
            << Un_mean << "  " << acc_m_tot << "  " << acc_v_tot << "  "
            << n_prod_samples << "\n";
        eos.flush();

        std::cout << "\n=== FINE RUN " << tag << " ===\n";
        std::cout << "P=" << std::scientific << std::setprecision(12) << pars.P
                  << "  <rho>=" << rho_mean
                  << "  <L>="   << L_mean
                  << "  <U/N>=" << Un_mean
                  << "  acc(move)=" << acc_m_tot
                  << "  acc(vol)="  << acc_v_tot
                  << "  prod_samples=" << n_prod_samples
                  << "\n";
        std::cout << "Wrote: " << obs_file << ", " << sum_file << ", " << frames_dir << "/\n";
    }

    std::cout << "\n====================================================\n";
    std::cout << "TUTTO FINITO. EOS salvata in: " << eos_file << "\n";
    std::cout << "====================================================\n";

    return 0;
}
