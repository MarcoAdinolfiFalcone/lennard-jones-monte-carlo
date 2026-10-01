// ======================================================
//  system__npt.cpp — implementazione di NPTSystem
//  (MC NPT isotropo, stile prof, fisica invariata)
// ======================================================

#include "system__npt.hpp"
#include <cmath>
#include <iostream>

// ---------------------------
// Parametri Lennard-Jones (unità ridotte)
// ---------------------------
namespace {
    constexpr double epsilon = 1.0;
    constexpr double sigma   = 1.0;
    constexpr double rc      = 2.5 * sigma;
    constexpr double rc2     = rc * rc;
}

// ======================================================
//  COSTRUTTORE
// ======================================================
NPTSystem::NPTSystem(const NPTParams& params_in)
    : pars(params_in), sim_box(), part(), U_current(0.0)
{
    // Volume iniziale da densità target: V = N / rho
    double V = static_cast<double>(pars.N) / pars.rho;
    double L = std::cbrt(V);
    sim_box = Box(L);

    // Alloco N particelle
    part.resize(pars.N);

    // Buffer store/restore box (allocato UNA volta, riusato ad ogni box move)
    stored_all_pos.resize(part.size());

    // Posizioni iniziali su reticolo cubico semplice (SC)
    init_lattice();

    // Energia iniziale del sistema
    U_current = compute_total_energy();
}

// ======================================================
//  Stampa info principali
// ======================================================
void NPTSystem::print_info() const
{
    std::cout << "=== NPT system ===\n";
    std::cout << "  N   = " << pars.N << "\n";
    std::cout << "  T   = " << pars.T << "\n";
    std::cout << "  P   = " << pars.P << "\n";
    std::cout << "  rho(target)   = " << pars.rho << "\n";
    std::cout << "  L   = " << sim_box.L << "\n";
    std::cout << "  V   = " << sim_box.volume() << "\n";
    std::cout << "  rho(computed) = " << sim_box.density(pars.N) << "\n";
    std::cout << "  U   = " << U_current << "\n";
    std::cout << "  U/N = " << U_current / pars.N << "\n";
}

// ======================================================
//  Inizializza posizioni su reticolo SC
// ======================================================
void NPTSystem::init_lattice()
{
    int N = pars.N;

    // numero di celle per lato: n3^3 >= N
    int n3 = static_cast<int>(std::ceil(std::cbrt(static_cast<double>(N))));
    double spacing = sim_box.L / n3;

    int idx = 0;
    for (int ix = 0; ix < n3 && idx < N; ++ix) {
        for (int iy = 0; iy < n3 && idx < N; ++iy) {
            for (int iz = 0; iz < n3 && idx < N; ++iz) {

                pvec3d r;
                r[0] = (ix + 0.5) * spacing;
                r[1] = (iy + 0.5) * spacing;
                r[2] = (iz + 0.5) * spacing;

                sim_box.apply_pbc(r);
                part[idx].r = r;

                ++idx;
                if (idx >= N) break;
            }
        }
    }
}

// ======================================================
//  Potenziale Lennard-Jones shiftato (funzione di r^2)
// ======================================================
double NPTSystem::lj_shifted(double r2) const
{
    if (r2 >= rc2) return 0.0;  // oltre cutoff → nulla

    double sig2 = sigma * sigma;

    double inv_r2  = sig2 / r2;                // (σ/r)^2
    double inv_r6  = inv_r2 * inv_r2 * inv_r2; // (σ/r)^6
    double inv_r12 = inv_r6 * inv_r6;          // (σ/r)^12

    double u = 4.0 * epsilon * (inv_r12 - inv_r6);

    // shift: sottraggo u(rc)
    double inv_rc2  = sig2 / rc2;
    double inv_rc6  = inv_rc2 * inv_rc2 * inv_rc2;
    double inv_rc12 = inv_rc6 * inv_rc6;
    double u_rc     = 4.0 * epsilon * (inv_rc12 - inv_rc6);

    return u - u_rc;
}

// ======================================================
//  Energia di interazione tra due particelle i e j
// ======================================================
double NPTSystem::pair_energy(int i, int j) const
{
    const pvec3d& ri = part[i].r;
    const pvec3d& rj = part[j].r;

    // distanza a immagine minima
    pvec3d dr = sim_box.delta_pbc(ri, rj);
    double r2 = dr * dr;

    return lj_shifted(r2);
}

// ======================================================
//  Energia totale (somma su tutte le coppie i<j)
// ======================================================
double NPTSystem::compute_total_energy()
{
    double U = 0.0;
    int N = pars.N;

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            U += pair_energy(i, j);
        }
    }

    U_current = U;
    return U;
}

// ======================================================
//  Metropolis helpers (stile prof)
// ======================================================
bool NPTSystem::acc(double alpha) const
{
    if (alpha >= 0.0) return true;
    // accetta se log(r) < alpha
    return (std::log(rng.ranf()) < alpha);
}

double NPTSystem::alpha_particle(double dU) const
{
    double beta = 1.0 / pars.T;
    return -beta * dU;
}

double NPTSystem::alpha_box(double dU, double V_old, double V_new) const
{
    double beta = 1.0 / pars.T;
    double dV   = V_new - V_old;
    return -beta * (dU + pars.P * dV)
           + pars.N * std::log(V_new / V_old);
}

// ======================================================
//  Store/restore (stile prof)
// ======================================================
void NPTSystem::store_particle(int i)
{
    stored_i     = i;
    stored_pos_i = part[i].r;
}

void NPTSystem::restore_particle()
{
    if (stored_i >= 0) part[stored_i].r = stored_pos_i;
}

void NPTSystem::store_box()
{
    stored_L = sim_box.L;
    stored_U = U_current;

    for (std::size_t i = 0; i < part.size(); ++i) {
        stored_all_pos[i] = part[i].r;
    }
}

void NPTSystem::restore_box()
{
    sim_box.L = stored_L;
    U_current = stored_U;

    for (std::size_t i = 0; i < part.size(); ++i) {
        part[i].r = stored_all_pos[i];
    }
}

// ======================================================
//  Mossa MC di traslazione (NVT)
// ======================================================
bool NPTSystem::mc_move_particle(double delta)
{
    int N = pars.N;
    if (N == 0) return false;

    // 1) scegli particella a caso
    int i = static_cast<int>(rng.ranf() * N);
    if (i >= N) i = N - 1;

    // 2) energia locale vecchia (contributi che coinvolgono i)
    double U_old_i = 0.0;
    for (int j = 0; j < N; ++j) {
        if (j == i) continue;
        U_old_i += pair_energy(i, j);
    }

    // 3) store posizione
    store_particle(i);

    // 4) proponi spostamento dr ∈ [-delta,delta]^3
    pvec3d dr;
    for (int k = 0; k < 3; ++k) {
        dr[k] = (2.0 * rng.ranf() - 1.0) * delta;
    }

    part[i].r += dr;
    sim_box.apply_pbc(part[i].r);

    // 5) energia locale nuova
    double U_new_i = 0.0;
    for (int j = 0; j < N; ++j) {
        if (j == i) continue;
        U_new_i += pair_energy(i, j);
    }

    double dU = U_new_i - U_old_i;

    // 6) Metropolis
    bool accepted = acc(alpha_particle(dU));

    if (accepted) {
        U_current += dU;
        return true;
    } else {
        restore_particle();
        return false;
    }
}

// ======================================================
//  Mossa MC di volume isotropa (NPT)
// ======================================================
bool NPTSystem::mc_move_volume(double delta_lnV)
{
    // Stato attuale
    double V_old = sim_box.volume();
    double L_old = sim_box.L;
    double U_old = U_current;

    // store completo (pos + L + U)
    store_box();

    // 1) Proponi lnV_new = lnV_old + xi * delta_lnV, xi∈[-1,1]
    double xi = 2.0 * rng.ranf() - 1.0;
    double lnV_old = std::log(V_old);
    double lnV_new = lnV_old + xi * delta_lnV;
    double V_new   = std::exp(lnV_new);
    double L_new   = std::cbrt(V_new);

    double scale = L_new / L_old;

    // 2) applica nuovo box e scala tutte le coordinate
    sim_box.L = L_new;
    for (auto& p : part) {
        p.r *= scale;
        sim_box.apply_pbc(p.r);
    }

    // 3) ricalcolo energia totale (come nel tuo codice)
    double U_new = compute_total_energy();
    double dU    = U_new - U_old;

    // 4) Metropolis NPT
    bool accepted = acc(alpha_box(dU, V_old, V_new));

    if (accepted) {
        U_current = U_new;
        return true;
    } else {
        restore_box();
        return false;
    }
}

// ======================================================
//  Loop MC (stile prof) + contatori acceptance
// ======================================================
void NPTSystem::move_NVT(double delta_move)
{
    for (int k = 0; k < pars.N; ++k) {
        ++att_p;
        if (mc_move_particle(delta_move)) ++acc_p;
    }
}

void NPTSystem::move_box(double delta_lnV)
{
    ++att_V;
    if (mc_move_volume(delta_lnV)) ++acc_V;
}

double NPTSystem::acc_particle() const
{
    return (att_p > 0) ? static_cast<double>(acc_p) / static_cast<double>(att_p) : 0.0;
}

double NPTSystem::acc_volume() const
{
    return (att_V > 0) ? static_cast<double>(acc_V) / static_cast<double>(att_V) : 0.0;
}

void NPTSystem::run(long n_sweeps_eq,
                    long n_sweeps_prod,
                    double delta_move,
                    double delta_lnV,
                    int print_every)
{
    long n_tot = n_sweeps_eq + n_sweeps_prod;

    for (long sweep = 1; sweep <= n_tot; ++sweep) {

        move_NVT(delta_move);
        move_box(delta_lnV);

        if (print_every > 0 && (sweep % print_every == 0)) {
            double V = sim_box.volume();
            double rho = static_cast<double>(pars.N) / V;

            std::cout << "Sweep " << sweep
                      << "  U/N=" << (U_current / static_cast<double>(pars.N))
                      << "  V="   << V
                      << "  rho=" << rho
                      << "  acc_p=" << acc_particle()
                      << "  acc_V=" << acc_volume()
                      << "\n";
        }
    }
}
