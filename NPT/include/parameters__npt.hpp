// ======================================================
// parameters__npt.hpp — Parametri per simulazione MC NPT
// ======================================================

#ifndef NPT_PARAMS_HPP
#define NPT_PARAMS_HPP

struct NPTParams
{
    // ==============================
    //  SISTEMA / TERMODINAMICA
    // ==============================
    int    N;     // numero di particelle
    double T;     // temperatura
    double P;     // pressione esterna
    double rho;   // densità target per configurazione iniziale
    double rc;    // cutoff LJ (se vuoi usarlo anche nei potenziali)

    // ==============================
    //  MONTE CARLO (NUOVI NOMI)  <-- USA QUESTI
    // ==============================
    long   n_sweeps_eq;     // sweep di equilibrazione
    long   n_sweeps_prod;   // sweep di produzione
    long   print_every;     // stampa ogni quanti sweep (0 = no stampa)

    double delta_move;      // ampiezza max spostamento per coordinata (particle move)
    double delta_lnV;       // ampiezza max variazione di ln(V) (box move)

    // ==============================
    //  RNG
    // ==============================
    int    seed;            // -1 = random

    // ==============================
    //  COMPATIBILITÀ (VECCHI NOMI)
    //  Non usarli nei nuovi main: servono solo se hai vecchio codice.
    // ==============================
    long   n_steps;         // (legacy) passi totali (non usato nel nuovo run a sweep)
    long   n_equil;         // (legacy) passi equil (non usato nel nuovo run a sweep)
    long   measure_every;   // (legacy) misura ogni tot passi

    double max_disp;        // (legacy) uguale a delta_move
    double max_dV;          // (legacy) ATTENZIONE: nel tuo codice la mossa è su ln(V), quindi max_dV = delta_lnV

    // ==============================
    //  COSTRUTTORE DEFAULT
    // ==============================
    NPTParams()
        : N(343),
          T(1.0),
          P(0.5),
          rho(0.8),
          rc(2.5),
          n_sweeps_eq(2000),
          n_sweeps_prod(8000),
          print_every(500),
          delta_move(0.2),
          delta_lnV(0.01),
          seed(-1),

          // legacy: valori coerenti con i nuovi
          n_steps(1'000'000),
          n_equil(100'000),
          measure_every(1'000),
          max_disp(delta_move),
          max_dV(delta_lnV)
    {}
};

#endif // NPT_PARAMS_HPP
