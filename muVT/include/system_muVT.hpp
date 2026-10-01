#ifndef SYSTEM_MUVT_HPP
#define SYSTEM_MUVT_HPP

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "parameters_muVT.hpp"
#include "particle_muVT.hpp"
#include "box.hpp"
#include "pvector.hpp"
#include "randnumgen.hpp"

// ======================================================
//   MuVTSystem — Lennard-Jones in ensemble µVT
//   (versione coerente col tuo Box [0,L) e delta_pbc())
// ======================================================
class MuVTSystem {
private:
    MuVTParams pars;

    Box box;                    // box cubico [0,L)
    std::vector<Particle> part; // posizioni

    double beta;                // 1/T
    double V;                   // L^3
    double U_current;           // energia totale

    // LJ cache
    double epsilon;
    double sigma;
    double rc;
    double rc2;
    double vshift;              // V_LJ(rc)

private:
    inline void sync_box_from_pars()
    {
        box.L = pars.L(0);
        V = box.volume();
    }

    inline pvec3d random_position_uniform()
    {
        // uniforme in [0,L)
        pvec3d r = { rng.ranf()*box.L, rng.ranf()*box.L, rng.ranf()*box.L };
        box.apply_pbc(r);
        return r;
    }

public:
    explicit MuVTSystem(const MuVTParams& p)
        : pars(p),
          box(),
          part(),
          beta(0.0),
          V(0.0),
          U_current(0.0),
          epsilon(p.epsilon),
          sigma(p.sigma),
          rc(p.rc),
          rc2(p.rc*p.rc),
          vshift(0.0)
    {
        // RNG
        if (pars.seed < 0) rng.rseed();
        else rng.seed((unsigned)pars.seed);

        // init box + probs + activity
        pars.set_box_from_init();
        pars.normalize_move_probs();
        pars.update_activity();

        sync_box_from_pars();

        beta = pars.beta();

        // init N e vettore particelle
        pars.Np = pars.N_init;
        part.assign(pars.Np, Particle{});

        // shift LJ
        const double invr2  = (sigma*sigma) / rc2;
        const double invr6  = invr2*invr2*invr2;
        const double invr12 = invr6*invr6;
        vshift = 4.0*epsilon*(invr12 - invr6);

        // lattice + energia
        init_lattice();
        U_current = compute_total_energy();
    }

    void print_info() const
    {
        std::cout << "=== MuVT system ===\n";
        std::cout << "Np=" << pars.Np << "  (N_init=" << pars.N_init << ")\n";
        std::cout << "T=" << pars.T << "  beta=" << pars.beta() << "\n";
        std::cout << "mu=" << pars.mu << "  z=" << pars.z << "\n";
        std::cout << "L=" << box.L << "  V=" << box.volume() << "  rho=" << box.density(pars.Np) << "\n";
        std::cout << "sigma=" << sigma << "  epsilon=" << epsilon << "  rc=" << rc << "\n";
        std::cout << "U=" << U_current << "  U/N=" << (pars.Np>0 ? U_current/double(pars.Np) : 0.0) << "\n";
        std::cout << "p_tra=" << pars.p_tra << "  p_ins=" << pars.p_ins << "  p_del=" << pars.p_del << "\n";
        std::cout << "deltra=" << pars.deltra << "\n";
    }

    void init_lattice()
    {
        if (pars.Np <= 0) return;

        const int Nloc = pars.Np;
        const int n3   = (int)std::ceil(std::cbrt((double)Nloc));
        const double a = box.L / n3;

        int cc = 0;
        for (int ix=0; ix<n3 && cc<Nloc; ++ix)
            for (int iy=0; iy<n3 && cc<Nloc; ++iy)
                for (int iz=0; iz<n3 && cc<Nloc; ++iz)
                {
                    part[cc].r = { (ix+0.5)*a, (iy+0.5)*a, (iz+0.5)*a };
                    box.apply_pbc(part[cc].r);
                    ++cc;
                }
    }

    double lj_shifted(double r2) const
    {
        if (r2 >= rc2) return 0.0;
        const double invr2  = (sigma*sigma) / r2;
        const double invr6  = invr2*invr2*invr2;
        const double invr12 = invr6*invr6;
        return 4.0*epsilon*(invr12 - invr6) - vshift;
    }

    double pair_energy(int i, int j) const
    {
        const pvec3d dr = box.delta_pbc(part[i].r, part[j].r);
        const double r2 = dr*dr;
        return lj_shifted(r2);
    }

    double compute_total_energy()
    {
        double U = 0.0;
        for (int i=0; i<pars.Np; ++i)
            for (int j=i+1; j<pars.Np; ++j)
                U += pair_energy(i,j);

        U_current = U;
        return U;
    }

    double delta_energy_insert(const pvec3d& r_new) const
    {
        double dU = 0.0;
        for (int j=0; j<pars.Np; ++j)
        {
            const pvec3d dr = box.delta_pbc(r_new, part[j].r);
            dU += lj_shifted(dr*dr);
        }
        return dU;
    }

    double delta_energy_delete(int index) const
    {
        double e = 0.0;
        for (int j=0; j<pars.Np; ++j)
        {
            if (j==index) continue;
            const pvec3d dr = box.delta_pbc(part[index].r, part[j].r);
            e += lj_shifted(dr*dr);
        }
        return e; // energia rimossa; ΔU=-e
    }

    double delta_energy_move(int i, const pvec3d& r_trial) const
    {
        double e_old = 0.0, e_new = 0.0;
        for (int j=0; j<pars.Np; ++j)
        {
            if (j==i) continue;

            const pvec3d dro = box.delta_pbc(part[i].r, part[j].r);
            e_old += lj_shifted(dro*dro);

            const pvec3d drn = box.delta_pbc(r_trial, part[j].r);
            e_new += lj_shifted(drn*drn);
        }
        return (e_new - e_old);
    }

    bool mc_move_particle(double delta)
    {
        if (pars.Np <= 0) return false;

        const int i = (int)(rng.ranf() * pars.Np);

        pvec3d r_trial = part[i].r + pvec3d{
            delta * 2.0 * (rng.ranf()-0.5),
            delta * 2.0 * (rng.ranf()-0.5),
            delta * 2.0 * (rng.ranf()-0.5)
        };
        box.apply_pbc(r_trial);

        const double dU = delta_energy_move(i, r_trial);
        const double xi = rng.ranf();

        if (dU <= 0.0 || xi < std::exp(-beta*dU))
        {
            part[i].r = r_trial;
            U_current += dU;
            return true;
        }
        return false;
    }

    bool mc_insert_particle()
    {
        const pvec3d r_new = random_position_uniform();
        const double dU = delta_energy_insert(r_new);

        // A_ins = min[1, (p_del/p_ins) * (z V/(N+1)) * exp(-beta dU)]
        const double pref = (pars.p_del / pars.p_ins) * (pars.z * V / (double(pars.Np) + 1.0));
        const double A = pref * std::exp(-beta*dU);

        const double xi = rng.ranf();
        if (xi < (A < 1.0 ? A : 1.0))
        {
            Particle p; p.r = r_new;
            part.push_back(p);
            ++pars.Np;
            U_current += dU;
            return true;
        }
        return false;
    }

    bool mc_delete_particle()
    {
        if (pars.Np <= 0) return false;

        const int i = (int)(rng.ranf() * pars.Np);
        const double e_removed = delta_energy_delete(i);

        // A_del = min[1, (p_ins/p_del) * (N/(zV)) * exp(+beta e_removed)]
        const double pref = (pars.p_ins / pars.p_del) * (double(pars.Np) / (pars.z * V));
        const double A = pref * std::exp(+beta*e_removed);

        const double xi = rng.ranf();
        if (xi < (A < 1.0 ? A : 1.0))
        {
            part[i] = part[pars.Np - 1];
            part.pop_back();
            --pars.Np;
            U_current -= e_removed;
            return true;
        }
        return false;
    }

    void mc_step()
    {
        const double r = rng.ranf();
        if (r < pars.p_tra)      mc_move_particle(pars.deltra);
        else if (r < pars.p_tra + pars.p_ins) mc_insert_particle();
        else                    mc_delete_particle();
    }

    void run_simulation(const std::string& out_prefix,
                        std::vector<long>& hist_N)
    {
        {
            std::ofstream f(out_prefix + "_muvt_observables.dat", std::ios::out | std::ios::trunc);
            f << "# step  N  V  rho  U  U_per_particle\n";
        }

        const int Nmax_hist = std::max(4*pars.N_init, 500);
        hist_N.assign(Nmax_hist + 1, 0);

        for (long t=0; t<pars.totsteps; ++t)
        {
            mc_step();

            const bool prod = (t >= pars.eqstps);
            if (prod)
            {
                if (pars.savemeasure > 0 && (t % pars.savemeasure == 0))
                {
                    std::ofstream f(out_prefix + "_muvt_observables.dat", std::ios::out | std::ios::app);
                    const double rho = box.density(pars.Np);
                    f << t << " " << pars.Np << " " << V << " " << rho << " "
                      << U_current << " " << (pars.Np>0 ? U_current/double(pars.Np) : 0.0) << "\n";
                }

                if (pars.Np >= 0 && pars.Np <= Nmax_hist) hist_N[pars.Np]++;

                if (t > 0 && pars.outstps > 0 && (t % pars.outstps == 0))
                {
                    std::cout << "[step " << t << "] N=" << pars.Np
                              << " rho=" << box.density(pars.Np)
                              << " U/N=" << (pars.Np>0 ? U_current/double(pars.Np) : 0.0)
                              << "\n";
                }
            }
        }

        {
            std::ofstream f(out_prefix + "_histN.dat", std::ios::out | std::ios::trunc);
            f << "# N  H(N)\n";
            for (int n=0; n<=Nmax_hist; ++n) f << n << " " << hist_N[n] << "\n";
        }
    }

    void write_xyz(const std::string& filename) const
    {
        std::ofstream out(filename);
        if (!out) return;

        out << pars.Np << "\n";
        out << "L=" << box.L << "\n";
        for (int i=0; i<pars.Np; ++i)
        {
            out << "Ar "
                << part[i].r(0) << " "
                << part[i].r(1) << " "
                << part[i].r(2) << "\n";
        }
    }

    int    get_N() const { return pars.Np; }
    double get_V() const { return V; }
    double get_L() const { return box.L; }
    double get_density() const { return box.density(pars.Np); }
    double get_energy() const { return U_current; }
};

#endif // SYSTEM_MUVT_HPP
