// ======================================================
// analyze_SUS_muVT.cpp
// Legge i risultati SUS da cartelle windows/window_XXXX/window_summary.dat
// e ricostruisce logW(N) e W_norm(N).
//
// Uso:
//   ./analyze_SUS <run_dir>
// Esempio:
//   ./analyze_SUS data/SUS/muVT_T1p40_V343_muRef_m3p50
//
// Output (dentro run_dir):
//   - ratios_from_windows.dat
//   - logW_from_windows.dat
//   - W_from_windows_norm.dat
//   - bad_windows_from_windows.dat
// ======================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <limits>

namespace fs = std::filesystem;

struct WinRec {
    int  n   = -1;
    long Hn  = 0;
    long Hn1 = 0;
    double ratio = 0.0;
    bool ok = false;
    fs::path path;
};

static bool parse_window_summary(const fs::path& fname, WinRec& out)
{
    std::ifstream fin(fname);
    if (!fin) return false;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::istringstream iss(line);
        int n;
        long Hn, Hn1;
        double ratio;
        if (!(iss >> n >> Hn >> Hn1 >> ratio)) {
            continue;
        }

        out.n = n;
        out.Hn = Hn;
        out.Hn1 = Hn1;

        // ricalcolo ratio se possibile (più affidabile del campo salvato)
        if (Hn > 0 && Hn1 >= 0) {
            out.ratio = double(Hn1) / double(Hn);
        } else {
            out.ratio = ratio;
        }

        out.ok = (out.n >= 0 && out.Hn > 0 && out.Hn1 > 0 && out.ratio > 0.0
                  && std::isfinite(out.ratio));
        return true;
    }
    return false;
}

static fs::path find_latest_run_dir(const fs::path& root)
{
    // Cerca la sottocartella più “recente” (ultimo write time) che contenga windows/
    if (!fs::exists(root) || !fs::is_directory(root)) {
        throw std::runtime_error("Root non valida: " + root.string());
    }

    fs::path best;
    std::filesystem::file_time_type best_t;

    bool found = false;
    for (const auto& entry : fs::directory_iterator(root)) {
        if (!entry.is_directory()) continue;
        fs::path run_dir = entry.path();
        fs::path win_dir = run_dir / "windows";
        if (!fs::exists(win_dir) || !fs::is_directory(win_dir)) continue;

        auto t = fs::last_write_time(run_dir);
        if (!found || t > best_t) {
            best_t = t;
            best = run_dir;
            found = true;
        }
    }

    if (!found) {
        throw std::runtime_error("Nessun run_dir trovato dentro: " + root.string());
    }
    return best;
}

int main(int argc, char* argv[])
{
    try {
        fs::path run_dir;

        if (argc >= 2) {
            run_dir = fs::path(argv[1]);
        } else {
            // default: cerca dentro data/SUS il run più recente
            fs::path root = "data/SUS";
            run_dir = find_latest_run_dir(root);
            std::cout << "[info] Nessun argomento: uso run_dir più recente = "
                      << run_dir.string() << "\n";
        }

        fs::path windows_dir = run_dir / "windows";
        if (!fs::exists(windows_dir) || !fs::is_directory(windows_dir)) {
            throw std::runtime_error("Directory windows non trovata: " + windows_dir.string());
        }

        // --------------------------------------------------
        // 1) Scansiona window_XXXX e leggi window_summary.dat
        // --------------------------------------------------
        std::vector<WinRec> recs;

        for (const auto& entry : fs::directory_iterator(windows_dir)) {
            if (!entry.is_directory()) continue;
            fs::path wdir = entry.path();
            std::string name = wdir.filename().string();
            if (name.rfind("window_", 0) != 0) continue; // non inizia con window_

            fs::path summary = wdir / "window_summary.dat";
            if (!fs::exists(summary)) continue;

            WinRec r;
            r.path = summary;
            if (parse_window_summary(summary, r)) {
                recs.push_back(r);
            }
        }

        if (recs.empty()) {
            throw std::runtime_error("Nessun window_summary.dat trovato in: " + windows_dir.string());
        }

        // Ordina per n
        std::sort(recs.begin(), recs.end(),
                  [](const WinRec& a, const WinRec& b){ return a.n < b.n; });

        // Verifica n_min, n_max (atteso: n = 0..nmax-1)
        int n_min = recs.front().n;
        int n_max = recs.back().n; // ultimo n (finestre fino a n_max)
        if (n_min < 0) {
            throw std::runtime_error("n_min < 0: qualcosa non torna.");
        }

        // Serve W(N) da N=n_min fino a N=n_max+1
        int N_min = n_min;
        int N_max = n_max + 1;
        int sizeN = N_max - N_min + 1;

        std::vector<double> logW(sizeN, 0.0);
        std::vector<char>   has(sizeN, 0);

        // gauge: logW(N_min)=0
        logW[0] = 0.0;
        has[0] = 1;

        // Output file “ratios” aggregato
        {
            std::ofstream fr(run_dir / "ratios_from_windows.dat", std::ios::out | std::ios::trunc);
            fr << "# n  Hn  Hn1  ratio=Hn1/Hn  ok\n";
            for (const auto& r : recs) {
                fr << r.n << " " << r.Hn << " " << r.Hn1 << " "
                   << std::setprecision(15) << r.ratio << " "
                   << (r.ok ? 1 : 0) << "\n";
            }
        }

        // File bad windows
        std::ofstream fbad(run_dir / "bad_windows_from_windows.dat", std::ios::out | std::ios::trunc);
        fbad << "# n  Hn  Hn1  ratio  NOTE  path\n";

        // --------------------------------------------------
        // 2) Ricostruisci logW ricorsivo usando ratio
        // logW(n+1)=logW(n)+log(ratio(n))
        // --------------------------------------------------
        // Mappa rapida n -> record (assumiamo contigui ma gestiamo buchi)
        int expected_n = n_min;
        for (const auto& r : recs) {
            if (r.n < expected_n) continue;
            // Se ci sono buchi, segnala
            if (r.n != expected_n) {
                for (int miss = expected_n; miss < r.n; ++miss) {
                    fbad << miss << " 0 0 0  MISSING_WINDOW  (none)\n";
                    // lascia logW costante (placeholder)
                    int idx_miss = (miss + 1) - N_min; // logW(miss+1)
                    int idx_prev = miss - N_min;       // logW(miss)
                    if (0 <= idx_prev && idx_prev < sizeN && 0 <= idx_miss && idx_miss < sizeN && has[idx_prev]) {
                        logW[idx_miss] = logW[idx_prev];
                        has[idx_miss] = 1;
                    }
                }
            }

            int idx_n   = (r.n) - N_min;
            int idx_np1 = (r.n + 1) - N_min;

            if (idx_n < 0 || idx_np1 < 0 || idx_n >= sizeN || idx_np1 >= sizeN) {
                fbad << r.n << " " << r.Hn << " " << r.Hn1 << " "
                     << std::setprecision(15) << r.ratio
                     << "  OUT_OF_RANGE  " << r.path.string() << "\n";
                expected_n = r.n + 1;
                continue;
            }

            if (!has[idx_n]) {
                // non dovrebbe succedere se partiamo da N_min
                has[idx_n] = 1;
            }

            if (!r.ok) {
                fbad << r.n << " " << r.Hn << " " << r.Hn1 << " "
                     << std::setprecision(15) << r.ratio
                     << "  BAD_STATS  " << r.path.string() << "\n";
                // placeholder: non aggiornare (mantieni logW uguale)
                logW[idx_np1] = logW[idx_n];
                has[idx_np1] = 1;
            } else {
                logW[idx_np1] = logW[idx_n] + std::log(r.ratio);
                has[idx_np1] = 1;
            }

            expected_n = r.n + 1;
        }

        // --------------------------------------------------
        // 3) Scrivi logW
        // --------------------------------------------------
        {
            std::ofstream flog(run_dir / "logW_from_windows.dat", std::ios::out | std::ios::trunc);
            flog << "# N  logW(N)   (gauge: logW(" << N_min << ")=0)\n";
            flog << std::setprecision(15);
            for (int i = 0; i < sizeN; ++i) {
                int N = N_min + i;
                flog << N << " " << logW[i] << "\n";
            }
        }

        // --------------------------------------------------
        // 4) Normalizza in modo stabile: W ~ exp(logW - maxlog)
        // --------------------------------------------------
        double maxlog = -std::numeric_limits<double>::infinity();
        for (double v : logW) {
            if (std::isfinite(v) && v > maxlog) maxlog = v;
        }

        std::vector<double> W(sizeN, 0.0);
        double Z = 0.0;
        for (int i = 0; i < sizeN; ++i) {
            double x = logW[i] - maxlog;
            // clamp per evitare underflow estremo (non necessario, ma sicuro)
            if (x < -700.0) x = -700.0;
            W[i] = std::exp(x);
            Z += W[i];
        }
        if (Z <= 0.0) {
            throw std::runtime_error("Normalizzazione fallita: Z<=0");
        }
        for (double& v : W) v /= Z;

        // Scrivi W_norm
        {
            std::ofstream fw(run_dir / "W_from_windows_norm.dat", std::ios::out | std::ios::trunc);
            fw << "# N  W_norm(N)   (sum=1)\n";
            fw << std::setprecision(15);
            for (int i = 0; i < sizeN; ++i) {
                int N = N_min + i;
                fw << N << " " << W[i] << "\n";
            }
        }

        std::cout << "Analisi completata.\n";
        std::cout << "Run dir: " << run_dir.string() << "\n";
        std::cout << "Finestre lette: " << recs.size() << " (n_min=" << n_min << ", n_max=" << n_max << ")\n";
        std::cout << "Output:\n";
        std::cout << "  " << (run_dir / "ratios_from_windows.dat").string() << "\n";
        std::cout << "  " << (run_dir / "logW_from_windows.dat").string() << "\n";
        std::cout << "  " << (run_dir / "W_from_windows_norm.dat").string() << "\n";
        std::cout << "  " << (run_dir / "bad_windows_from_windows.dat").string() << "\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
