// ======================================================
//  EOS_from_SUS_muVT.cpp
//  Post-processing: da SUS (muVT) ricostruisce EOS P(rho)
//  usando reweighting in mu.
//
//  INPUT (in run_dir):
//    - logW_from_windows.dat   (# N logW(N) gauge: logW(Nmin)=0, tipicamente Nmin=0)
//
//  OUTPUT (in run_dir):
//    - EOS_mu_rho_P.dat        (# mu rho P Nmean lnXi)
//    - W_reweighted_muXXXX.dat (opzionale, per alcuni mu)
//
//  Uso:
//    ./eos_sus <run_dir> [mu_min mu_max dmu] [T V mu_ref] [dump_every]
//  Esempio:
//    ./eos_sus data/SUS/muVT_T1p40_V343_muRef_m3p50
// ======================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <filesystem>

namespace fs = std::filesystem;

struct Rec {
    int N = 0;
    double logW_ratio = 0.0; // log[ W(N|mu_ref) / W(Nmin|mu_ref) ], gauge
};

static double log_sum_exp(const std::vector<double>& a)
{
    double m = -std::numeric_limits<double>::infinity();
    for (double x : a) m = std::max(m, x);
    if (!std::isfinite(m)) return m;

    double s = 0.0;
    for (double x : a) s += std::exp(x - m);
    return m + std::log(s);
}

static bool read_logW_ratio(const fs::path& fname, std::vector<Rec>& out, int& Nmin, int& Nmax)
{
    std::ifstream fin(fname);
    if (!fin) return false;

    out.clear();
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        Rec r;
        if (!(iss >> r.N >> r.logW_ratio)) continue;
        if (!std::isfinite(r.logW_ratio)) continue;
        out.push_back(r);
    }
    if (out.empty()) return false;

    std::sort(out.begin(), out.end(), [](const Rec& a, const Rec& b){ return a.N < b.N; });
    Nmin = out.front().N;
    Nmax = out.back().N;
    return true;
}

static std::string mu_to_tag(double mu)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(2) << mu;
    std::string s = os.str();
    for (char& c : s) {
        if (c == '-') c = 'm';
        else if (c == '.') c = 'p';
    }
    return s;
}

int main(int argc, char* argv[])
{
    // -----------------------------
    // Default (coerenti col tuo caso)
    // -----------------------------
    double mu_min = -6.0;
    double mu_max = -2.0;
    double dmu    =  0.02;

    double T      =  1.4;
    double V      =  343.0;
    double mu_ref = -3.5;

    // dump_every = ogni quanti punti mu scrivere anche W(N;mu) per debug
    // 0 => non dumpare mai
    int dump_every = 0;

    if (argc < 2) {
        std::cerr << "Uso:\n"
                  << "  " << argv[0] << " <run_dir> [mu_min mu_max dmu] [T V mu_ref] [dump_every]\n"
                  << "Esempio:\n"
                  << "  " << argv[0] << " data/SUS/muVT_T1p40_V343_muRef_m3p50\n";
        return 1;
    }

    fs::path run_dir = fs::path(argv[1]);
    if (!fs::exists(run_dir) || !fs::is_directory(run_dir)) {
        std::cerr << "Errore: run_dir non valido: " << run_dir.string() << "\n";
        return 1;
    }

    // parsing argomenti opzionali
    int argi = 2;
    if (argc >= argi + 3) {
        mu_min = std::stod(argv[argi + 0]);
        mu_max = std::stod(argv[argi + 1]);
        dmu    = std::stod(argv[argi + 2]);
        argi  += 3;
    }
    if (argc >= argi + 3) {
        T      = std::stod(argv[argi + 0]);
        V      = std::stod(argv[argi + 1]);
        mu_ref = std::stod(argv[argi + 2]);
        argi  += 3;
    }
    if (argc >= argi + 1) {
        dump_every = std::stoi(argv[argi + 0]);
        argi += 1;
    }

    if (dmu <= 0.0) {
        std::cerr << "Errore: dmu deve essere > 0\n";
        return 1;
    }
    const double beta = 1.0 / T;

    // -----------------------------
    // Leggi logW(N)/W(Nmin)
    // -----------------------------
    const fs::path logW_file = run_dir / "logW_from_windows.dat";

    std::vector<Rec> data;
    int Nmin = 0, Nmax = 0;
    if (!read_logW_ratio(logW_file, data, Nmin, Nmax)) {
        std::cerr << "Errore: impossibile leggere " << logW_file.string() << "\n";
        return 1;
    }

    std::cout << "=== EOS from SUS (muVT) ===\n";
    std::cout << "run_dir: " << run_dir.string() << "\n";
    std::cout << "input:   " << logW_file.string() << "\n";
    std::cout << "N range: [" << Nmin << ", " << Nmax << "]  (#=" << data.size() << ")\n";
    std::cout << "T=" << T << "  beta=" << beta << "  V=" << V << "  mu_ref=" << mu_ref << "\n";
    std::cout << "mu grid: [" << mu_min << ", " << mu_max << "] step " << dmu << "\n";
    if (dump_every > 0) std::cout << "dump W(N;mu) every " << dump_every << " points\n";

    // -----------------------------
    // Output EOS
    // -----------------------------
    const fs::path out_eos = run_dir / "EOS_mu_rho_P.dat";
    std::ofstream fout(out_eos, std::ios::out | std::ios::trunc);
    if (!fout) {
        std::cerr << "Errore: impossibile scrivere " << out_eos.string() << "\n";
        return 1;
    }

    fout << "# EOS from SUS (muVT) via reweighting\n";
    fout << "# input_logW = " << logW_file.filename().string() << "\n";
    fout << "# gauge: logW(Nmin)=0 with Nmin=" << Nmin << "\n";
    fout << "# params: T=" << T << "  V=" << V << "  mu_ref=" << mu_ref << "\n";
    fout << "# columns: mu   rho(<N>/V)   P   <N>   lnXi\n";
    fout << std::scientific << std::setprecision(12);

    // workspace
    std::vector<double> a(data.size(), 0.0); // log-weights non normalizzati a_N(mu)
    std::vector<double> w(data.size(), 0.0); // W(N;mu) normalizzata

    int imu = 0;
    for (double mu = mu_min; mu <= mu_max + 1e-14; mu += dmu, ++imu)
    {
        // a_N(mu) = logW_ratio(N) + beta (mu-mu_ref) N
        // NB: logW_ratio è rispetto a Nmin; ok perché poi normalizziamo
        for (std::size_t i = 0; i < data.size(); ++i) {
            a[i] = data[i].logW_ratio + beta * (mu - mu_ref) * double(data[i].N);
        }

        // lnXi = ln sum_N exp(a_N)
        const double lnXi = log_sum_exp(a);
        if (!std::isfinite(lnXi)) continue;

        // W(N;mu) = exp(a_N - lnXi)
        double Nmean = 0.0;
        for (std::size_t i = 0; i < data.size(); ++i) {
            w[i] = std::exp(a[i] - lnXi);
            Nmean += double(data[i].N) * w[i];
        }

        const double rho = Nmean / V;
        const double P   = (1.0 / beta) * (lnXi / V); // kT/V * lnXi

        fout << mu << "  " << rho << "  " << P << "  " << Nmean << "  " << lnXi << "\n";

        // Dump distribuzione W(N;mu) per debug/plot (opzionale)
        if (dump_every > 0 && (imu % dump_every == 0)) {
            fs::path outW = run_dir / ("W_reweighted_mu_" + mu_to_tag(mu) + ".dat");
            std::ofstream fw(outW, std::ios::out | std::ios::trunc);
            fw << "# W(N;mu) reweighted from mu_ref\n";
            fw << "# mu=" << std::setprecision(12) << mu << "  mu_ref=" << mu_ref
               << "  T=" << T << "  V=" << V << "  lnXi=" << lnXi << "\n";
            fw << "# columns: N  W(N)\n";
            fw << std::scientific << std::setprecision(15);
            for (std::size_t i = 0; i < data.size(); ++i) {
                fw << data[i].N << " " << w[i] << "\n";
            }
        }
    }

    fout.close();
    std::cout << "OK. Scritto: " << out_eos.string() << "\n";
    return 0;
}
