// ======================================================
// NPT_system.hpp  — Header per il sistema Monte Carlo NPT
// (versione "stile prof" + compatibile col tuo codice attuale)
// ======================================================

#ifndef NPT_SYSTEM_HPP
#define NPT_SYSTEM_HPP

#include "pvector.hpp"
#include "parameters__npt.hpp"
#include "box.hpp"
#include "particle_npt.hpp"

#include <vector>
#include <string>
#include <iosfwd>

// Classe che rappresenta lo stato del sistema NPT:
// - parametri (N, T, P, rho)
// - box (L, volume, PBC)
// - coordinate delle particelle
// - energia potenziale U
class NPTSystem {
public:
    // --------------------------------------------------
    // COSTRUTTORE
    // --------------------------------------------------
    explicit NPTSystem(const NPTParams& params_in);

    // --------------------------------------------------
    // ACCESSO AI PARAMETRI E STRUTTURE
    // --------------------------------------------------
    const NPTParams& get_params() const { return pars; }
    const Box& get_box() const { return sim_box; }

    const std::vector<Particle>& get_particles() const { return part; }
    std::vector<Particle>& particles() { return part; }  // accesso non-const

    // --------------------------------------------------
    // OUTPUT / DEBUG
    // --------------------------------------------------
    void print_info() const;

    // --------------------------------------------------
    // ENERGIA
    // --------------------------------------------------
    double get_energy() const { return U_current; }

    // Calcola da zero l'energia totale U e aggiorna U_current
    double compute_total_energy();

    // --------------------------------------------------
    // MONTE CARLO (mosse base, come nel tuo codice)
    // --------------------------------------------------
    // 1) Mossa di particella (NVT)
    bool mc_move_particle(double delta);

    // 2) Mossa di volume (NPT isotropo)
    //    delta_lnV = ampiezza max per variazione di log-volume
    bool mc_move_volume(double delta_lnV);

    // --------------------------------------------------
    // MONTE CARLO
    void run(long n_sweeps_eq,
             long n_sweeps_prod,
             double delta_move,
             double delta_lnV,
             int print_every = 0);

    // Uno sweep di particelle (tipicamente N tentativi)
    void move_NVT(double delta_move);

    // Un tentativo di mossa di volume
    void move_box(double delta_lnV);

    // Acceptance ratio (utile da stampare)
    double acc_particle() const;
    double acc_volume() const;

    // Contatori (se ti servono)
    long long attempts_particle() const { return att_p; }
    long long accepts_particle()  const { return acc_p; }
    long long attempts_volume()   const { return att_V; }
    long long accepts_volume()    const { return acc_V; }

private:
    // --------------------------------------------------
    // STATO DEL SISTEMA
    // --------------------------------------------------
    NPTParams             pars;        // parametri della simulazione
    Box                   sim_box;     // scatola cubica con PBC
    std::vector<Particle> part;        // particelle
    double                U_current;   // energia potenziale attuale

    // --------------------------------------------------
    // CONTATORI MC (stile prof)
    // --------------------------------------------------
    long long att_p = 0, acc_p = 0; // particelle
    long long att_V = 0, acc_V = 0; // volume

    // --------------------------------------------------
    // STORE/RESTORE
    // --------------------------------------------------
    // particle move
    int    stored_i = -1;
    pvec3d stored_pos_i;

    // volume move
    double stored_L = 0.0;
    double stored_U = 0.0;
    std::vector<pvec3d> stored_all_pos; // buffer

    void store_particle(int i);
    void restore_particle();

    void store_box();
    void restore_box();

    // --------------------------------------------------
    // METROPOLIS HELPERS 
    // --------------------------------------------------
    // alpha = log(pi_new / pi_old)
    bool   acc(double alpha) const;

    // NVT: alpha = -beta * dU
    double alpha_particle(double dU) const;

    // NPT: alpha = -beta*(dU + P*dV) + N*ln(Vnew/Vold)
    double alpha_box(double dU, double V_old, double V_new) const;

    // --------------------------------------------------
    // METODI INTERNI
    // --------------------------------------------------
    // Inizializzazione delle posizioni (reticolo SC)
    void init_lattice();

    // Energia di interazione tra particelle i e j
    double pair_energy(int i, int j) const;

    // Potenziale Lennard-Jones con shift al cutoff rc (rc fissato nel .cpp)
    double lj_shifted(double r2) const;
};

#endif // NPT_SYSTEM_HPP
