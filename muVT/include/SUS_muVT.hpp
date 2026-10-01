#ifndef SUS_MUVT_HPP
#define SUS_MUVT_HPP

#include <string>
#include "parameters_muVT.hpp"   // basta questo nel .hpp

// Risultato di una singola finestra SUS [n, n+1]
struct SUSWindowResult {
    int    n   = 0;   // indice della finestra: [n, n+1]
    long   Hn  = 0;   // tempo passato in N = n
    long   Hn1 = 0;   // tempo passato in N = n+1
    double ratio = 0.0; // Hn1 / Hn (stima W(n+1)/W(n))
};

/**
 * Esegue una singola finestra SUS in ensemble µVT,
 * vincolando il numero di particelle a N ∈ {n, n+1}.
 *
 * - base_pars: parametri µVT (T, mu_ref, V_box/rho_init, p_tra/p_ins/p_del, deltra, ecc.)
 *              N_init verrà sovrascritto con n.
 * - n: indice finestra (N ∈ {n, n+1})
 * - n_steps: passi MC totali
 * - n_equil: passi MC di equilibrazione (non contano in Hn, Hn1)
 * - out_dir: se non vuoto salva file ordinati (cartella window_XXX)
 * - verbose: stampa info su cout
 */
SUSWindowResult run_SUS_window(const MuVTParams& base_pars,
                               int n,
                               long n_steps,
                               long n_equil,
                               const std::string& out_dir = "",
                               bool verbose = false);

#endif // SUS_MUVT_HPP
