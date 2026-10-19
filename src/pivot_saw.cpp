// SAW planare con N passi e N+1 siti, reticolo di passo unitario.
// Una run evolve una sola catena mediante proposte Pivot: le configurazioni
// successive sono correlate. Le repliche e la loro aggregazione sono gestite
// dagli script esterni; qui non vengono eseguiti fit dell'esponente nu.
// Osservabili del confronto: Re2, rigetto della mossa, L e media(L/(N-k)).
// Il rigetto conserva la catena, mentre nel simple sampling si scarta una costruzione.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

struct Point {
    int x = 0;
    int y = 0;
};

// Due punti della griglia sono uguali se coincidono entrambe le coordinate.
static inline bool operator==(const Point& a, const Point& b) {
    return a.x == b.x && a.y == b.y;
}

struct PointHash {
    // Converte una coppia (x, y) in un hash a 64 bit: serve per cercare
    // rapidamente se un sito della griglia e' gia' occupato dalla camminata.
    std::size_t operator()(const Point& p) const noexcept {
        const std::uint64_t ux = static_cast<std::uint32_t>(p.x);
        const std::uint64_t uy = static_cast<std::uint32_t>(p.y);
        std::uint64_t z = (ux << 32) ^ uy;
        z ^= z >> 30;
        z *= 0xbf58476d1ce4e5b9ULL;
        z ^= z >> 27;
        z *= 0x94d049bb133111ebULL;
        z ^= z >> 31;
        return static_cast<std::size_t>(z);
    }
};

class SplitMix64 {
public:
    // Generatore usato solo per inizializzare bene lo stato del generatore
    // principale anche quando il seed scelto dall'utente e' semplice.
    explicit SplitMix64(std::uint64_t seed) : state_(seed) {}

    std::uint64_t next() {
        std::uint64_t z = (state_ += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

private:
    std::uint64_t state_;
};

class Xoshiro256StarStar {
public:
    // Generatore pseudo-casuale veloce; produce i pivot e le simmetrie casuali.
    explicit Xoshiro256StarStar(std::uint64_t seed) {
        SplitMix64 sm(seed);
        for (auto& v : s_) v = sm.next();
    }

    std::uint64_t next_u64() {
        const std::uint64_t result = rotl(s_[1] * 5ULL, 7) * 9ULL;
        const std::uint64_t t = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        return result;
    }

    std::uint64_t uniform_u64(std::uint64_t bound) {
        if (bound == 0) throw std::invalid_argument("uniform_u64 bound is zero");
        // Rejection sampling: evita il piccolo bias che si avrebbe usando solo r % bound.
        const std::uint64_t threshold = -bound % bound;
        while (true) {
            const std::uint64_t r = next_u64();
            if (r >= threshold) return r % bound;
        }
    }

private:
    static std::uint64_t rotl(const std::uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    std::array<std::uint64_t, 4> s_{};
};

enum class SymmetryClass {
    Rot90,
    Rot180,
    AxisReflection,
    DiagonalReflection
};

struct Symmetry {
    int a, b, c, d;
    const char* name;
    SymmetryClass cls;
};

static const std::array<Symmetry, 7> kD4NonIdentity = {{
    // Le 7 simmetrie non banali del quadrato: rotazioni e riflessioni.
    // L'identita' non viene proposta perche' lascerebbe la camminata invariata.
    {0, -1, 1, 0, "rot90", SymmetryClass::Rot90},
    {0, 1, -1, 0, "rot270", SymmetryClass::Rot90},
    {-1, 0, 0, -1, "rot180", SymmetryClass::Rot180},
    {1, 0, 0, -1, "reflect_x", SymmetryClass::AxisReflection},
    {-1, 0, 0, 1, "reflect_y", SymmetryClass::AxisReflection},
    {0, 1, 1, 0, "reflect_y_eq_x", SymmetryClass::DiagonalReflection},
    {0, -1, -1, 0, "reflect_y_eq_minus_x", SymmetryClass::DiagonalReflection},
}};

struct Options {
    // Parametri controllabili da riga di comando: lunghezze simulate,
    // numero di tentativi pivot, fase di termalizzazione e percorsi di output.
    std::vector<int> lengths{32, 64, 128};
    std::uint64_t attempts = 20000;
    std::uint64_t thermalization = 5000;
    std::uint64_t sample_every = 1;
    std::uint64_t seed = 12345;
    std::string out = "outputs/results.csv";
    std::string series_out;
    std::uint64_t validate_every = 0;
    bool autocorrelation = true;
    bool verbose = true;
};

struct RunningStats {
    std::uint64_t n = 0;
    long double mean = 0.0L;
    long double m2 = 0.0L;

    void push(long double x) {
        // Algoritmo online di Welford: aggiorna media e varianza senza
        // conservare tutti i campioni in memoria.
        ++n;
        const long double delta = x - mean;
        mean += delta / static_cast<long double>(n);
        const long double delta2 = x - mean;
        m2 += delta * delta2;
    }

    long double variance() const {
        return n > 1 ? m2 / static_cast<long double>(n - 1) : 0.0L;
    }
};

// Re2 e' la distanza quadratica tra estremi. Rg2 e Rm2 sono diagnostiche
// aggiuntive conservate nel programma, non le variabili del confronto in tesi.
struct Observables {
    long double re2 = 0.0L;
    long double rg2 = 0.0L;
    long double rm2 = 0.0L;
};

// Risultato di UNA proposta: L conta i siti controllati prima della prima
// collisione, nell'ordine dal pivot verso l'estremo. Non e' una distanza
// spaziale percorsa dal polimero ne' una lunghezza di persistenza.
struct PivotAttempt {
    bool accepted = false;
    int symmetry_index = 0;
    int tail_length = 0;
    int free_path = 0;
};

struct RunResult {
    int n = 0;
    std::uint64_t attempts = 0;
    std::uint64_t thermalization = 0;
    std::uint64_t samples = 0;
    std::uint64_t accepted = 0;
    std::array<std::uint64_t, 4> proposed_by_class{};
    std::array<std::uint64_t, 4> accepted_by_class{};
    RunningStats re2;
    RunningStats rg2;
    RunningStats rm2;
    RunningStats accept_indicator;
    // failed_*: medie condizionate al rigetto. attempt_*: medie su TUTTE
    // le proposte (quelle usate per il confronto con il simple sampling).
    RunningStats failed_free_path;
    RunningStats failed_free_path_fraction;
    RunningStats attempt_free_path;
    RunningStats attempt_free_path_fraction;
    double tau_re2 = std::numeric_limits<double>::quiet_NaN();
    double tau_rg2 = std::numeric_limits<double>::quiet_NaN();
    double tau_rm2 = std::numeric_limits<double>::quiet_NaN();
    double tau_accept = std::numeric_limits<double>::quiet_NaN();
    double seconds = 0.0;
};

static int class_index(SymmetryClass cls) {
    switch (cls) {
        case SymmetryClass::Rot90: return 0;
        case SymmetryClass::Rot180: return 1;
        case SymmetryClass::AxisReflection: return 2;
        case SymmetryClass::DiagonalReflection: return 3;
    }
    return 0;
}

static Point apply_symmetry(const Symmetry& g, const Point& rel) {
    // Applica la matrice 2x2 della simmetria a un punto espresso rispetto al pivot.
    return {g.a * rel.x + g.b * rel.y, g.c * rel.x + g.d * rel.y};
}

static std::vector<int> parse_lengths(const std::string& text) {
    // Converte una stringa tipo "200,400,800" nella lista di lunghezze da simulare.
    std::vector<int> result;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        const int value = std::stoi(item);
        if (value < 2) throw std::invalid_argument("all lengths must be >= 2");
        result.push_back(value);
    }
    if (result.empty()) throw std::invalid_argument("no lengths supplied");
    return result;
}

static Options parse_args(int argc, char** argv) {
    // Parser minimale degli argomenti: ogni opzione con valore consuma anche
    // argv[i+1], mentre --quiet e --help sono flag senza valore.
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto need_value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument(std::string("missing value for ") + name);
            return argv[++i];
        };

        if (key == "--lengths") opt.lengths = parse_lengths(need_value("--lengths"));
        else if (key == "--attempts") opt.attempts = std::stoull(need_value("--attempts"));
        else if (key == "--thermalization") opt.thermalization = std::stoull(need_value("--thermalization"));
        else if (key == "--sample-every") opt.sample_every = std::stoull(need_value("--sample-every"));
        else if (key == "--seed") opt.seed = std::stoull(need_value("--seed"));
        else if (key == "--out") opt.out = need_value("--out");
        else if (key == "--series-out") opt.series_out = need_value("--series-out");
        else if (key == "--validate-every") opt.validate_every = std::stoull(need_value("--validate-every"));
        else if (key == "--no-autocorrelation") opt.autocorrelation = false;
        else if (key == "--quiet") opt.verbose = false;
        else if (key == "--help") {
            std::cout
                << "Usage: pivot_saw.exe [options]\n"
                << "  --lengths 200,400,800\n"
                << "  --attempts 1000000\n"
                << "  --thermalization 200000\n"
                << "  --sample-every 1\n"
                << "  --seed 12345\n"
                << "  --out outputs/results.csv\n"
                << "  --series-out outputs/series_last_length.csv\n"
                << "  --validate-every 0  (0 disables invariant checks)\n"
                << "  --no-autocorrelation  (O(1) memory production mode)\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + key);
        }
    }
    if (opt.sample_every == 0) throw std::invalid_argument("--sample-every must be positive");
    return opt;
}

class PivotSaw2D {
public:
    explicit PivotSaw2D(int n) : walk_(static_cast<std::size_t>(n) + 1) {
        // Stato iniziale: camminata rettilinea lungo l'asse x, sicuramente self-avoiding.
        occupied_.reserve(walk_.size() * 2);
        for (int i = 0; i <= n; ++i) {
            const Point p{i, 0};
            walk_[static_cast<std::size_t>(i)] = p;
            occupied_.emplace(p, i);
            sum_x_ += p.x;
            sum_y_ += p.y;
            sum_r2_ += squared_norm(p);
        }
    }

    int length() const { return static_cast<int>(walk_.size()) - 1; }

    PivotAttempt attempt(Xoshiro256StarStar& rng) {
        const int n = length();
        // k uniforme in [1,N-1]: k=N non avrebbe coda, k=0 ruoterebbe
        // l'intera catena. Il primo legame resta quindi orientato lungo +x;
        // questo fissa l'orientazione globale senza cambiare Re2, che e'
        // invariante per rotazioni e riflessioni.
        const int k = 1 + static_cast<int>(rng.uniform_u64(static_cast<std::uint64_t>(n - 1)));
        PivotAttempt result;
        // Le sette matrici sono equiprobabili e l'insieme contiene le
        // inverse: con accettazione per sola autoesclusione, le proposte
        // soddisfano il bilancio dettagliato per pesi configurazionali uguali.
        result.symmetry_index = static_cast<int>(rng.uniform_u64(kD4NonIdentity.size()));
        result.tail_length = n - k;
        const Symmetry& g = kD4NonIdentity[static_cast<std::size_t>(result.symmetry_index)];
        const Point pivot = walk_[static_cast<std::size_t>(k)];

        // Buffer provvisorio: fino all'accettazione non modifichiamo walk_
        // o occupied_. In caso di collisione basta uscire, senza rollback.
        scratch_.clear();
        scratch_.reserve(static_cast<std::size_t>(n - k));

        // Trasforma la coda punto per punto. Se un candidato collide con la
        // parte fissa, il pivot viene rifiutato e la camminata originale resta valida.
        // La mappa persistente contiene tutta la SAW: gli indici > k appartengono
        // alla vecchia coda, che viene sostituita, e quindi non sono collisioni.
        // Una simmetria rigida e' biunivoca: la coda trasformata resta
        // internamente autoevitante. Basta cercare collisioni con il prefisso.
        for (int i = k + 1; i <= n; ++i) {
            const Point old = walk_[static_cast<std::size_t>(i)];
            const Point rel{old.x - pivot.x, old.y - pivot.y};
            const Point t = apply_symmetry(g, rel);
            const Point candidate{pivot.x + t.x, pivot.y + t.y};
            const auto hit = occupied_.find(candidate);
            if (hit != occupied_.end() && hit->second <= k) {
                result.accepted = false;
                // Esclude il sito che collide. Cambiare l'ordine dei
                // controlli cambierebbe questa diagnostica, non il rigetto.
                result.free_path = i - k - 1;
                return result;
            }
            scratch_.push_back(candidate);
        }

        // Prima rimuove tutte le vecchie coordinate. Questo ordine e' essenziale:
        // una nuova coordinata puo' coincidere con un sito della vecchia coda,
        // che non fa piu' parte della configurazione proposta.
        for (int i = k + 1; i <= n; ++i) {
            const std::size_t index = static_cast<std::size_t>(i);
            const std::size_t erased = occupied_.erase(walk_[index]);
            if (erased != 1) throw std::logic_error("persistent occupancy map is inconsistent");
        }

        // Aggiorna insieme coordinate, mappa di occupazione e somme geometriche.
        for (int i = k + 1; i <= n; ++i) {
            const std::size_t index = static_cast<std::size_t>(i);
            const Point old = walk_[index];
            const Point candidate = scratch_[static_cast<std::size_t>(i - k - 1)];
            sum_x_ += static_cast<long double>(candidate.x) - old.x;
            sum_y_ += static_cast<long double>(candidate.y) - old.y;
            sum_r2_ += squared_norm(candidate) - squared_norm(old);
            walk_[index] = candidate;
            const bool inserted = occupied_.emplace(candidate, i).second;
            if (!inserted) throw std::logic_error("accepted pivot produced a self-intersection");
        }
        result.accepted = true;
        // Tutta la coda e' stata verificata: L/(N-k)=1 per i successi.
        result.free_path = result.tail_length;
        return result;
    }

    Observables observables() const {
        // Le somme geometriche sono mantenute al momento del pivot: la misura e' O(1).
        const long double m = static_cast<long double>(walk_.size());
        const Point end = walk_.back();
        Observables obs;
        // r_0 resta nell'origine, quindi |r_N-r_0|^2 = |r_N|^2.
        obs.re2 = squared_norm(end);
        // Rg2 = media(|r|^2) - |baricentro|^2; Rm2 = media(|r-r_0|^2).
        obs.rg2 = sum_r2_ / m - (sum_x_ * sum_x_ + sum_y_ * sum_y_) / (m * m);
        obs.rm2 = sum_r2_ / m;
        return obs;
    }

    void validate_state() const {
        // Controllo costoso destinato ai test: self-avoidance, indice di ogni
        // monomero e somme incrementali devono coincidere con il calcolo diretto.
        if (occupied_.size() != walk_.size()) {
            throw std::logic_error("walk is not self-avoiding or occupancy size is wrong");
        }

        long double sum_x = 0.0L;
        long double sum_y = 0.0L;
        long double sum_r2 = 0.0L;
        for (std::size_t i = 0; i < walk_.size(); ++i) {
            const Point p = walk_[i];
            const auto hit = occupied_.find(p);
            if (hit == occupied_.end() || hit->second != static_cast<int>(i)) {
                throw std::logic_error("occupancy index does not match the walk");
            }
            sum_x += p.x;
            sum_y += p.y;
            sum_r2 += squared_norm(p);
        }
        if (sum_x != sum_x_ || sum_y != sum_y_ || sum_r2 != sum_r2_) {
            throw std::logic_error("incremental geometric sums are inconsistent");
        }
    }

private:
    static long double squared_norm(const Point& p) {
        return static_cast<long double>(p.x) * p.x + static_cast<long double>(p.y) * p.y;
    }

    std::vector<Point> walk_;
    std::vector<Point> scratch_;
    std::unordered_map<Point, int, PointHash> occupied_;
    long double sum_x_ = 0.0L;
    long double sum_y_ = 0.0L;
    long double sum_r2_ = 0.0L;
};

// FFT radix-2: funzione ausiliaria per l'autocorrelazione, non per
// generare le configurazioni. La dimensione passata e' una potenza di due.
static void fft(std::vector<std::complex<double>>& values, bool inverse) {
    const std::size_t n = values.size();
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(values[i], values[j]);
    }

    const double pi = std::acos(-1.0);
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double angle = (inverse ? 2.0 : -2.0) * pi / static_cast<double>(len);
        const std::complex<double> root(std::cos(angle), std::sin(angle));
        for (std::size_t start = 0; start < n; start += len) {
            std::complex<double> factor(1.0, 0.0);
            for (std::size_t j = 0; j < len / 2; ++j) {
                const std::complex<double> even = values[start + j];
                const std::complex<double> odd = values[start + j + len / 2] * factor;
                values[start + j] = even + odd;
                values[start + j + len / 2] = even - odd;
                factor *= root;
            }
        }
        if (len == n) break;  // Evita overflow di len <<= 1 per input estremi.
    }

    if (inverse) {
        const double scale = 1.0 / static_cast<double>(n);
        for (auto& value : values) value *= scale;
    }
}

static double integrated_autocorrelation_time(const std::vector<double>& x) {
    // Autocovarianze lineari via FFT: zero-padding a 2N evita il wrap-around.
    // La stessa finestra automatica della versione diretta preserva la
    // definizione statistica di tau_int, riducendo il costo a O(N log N).
    const std::size_t n = x.size();
    if (n < 8) return std::numeric_limits<double>::quiet_NaN();

    const double mean = std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(n);
    double c0_sum = 0.0;
    for (double v : x) {
        const double d = v - mean;
        c0_sum += d * d;
    }
    const double c0 = c0_sum / static_cast<double>(n);
    if (c0 <= 0.0) return 0.5;

    if (n > std::numeric_limits<std::size_t>::max() / 2) {
        throw std::length_error("time series is too long for FFT padding");
    }
    std::size_t fft_size = 1;
    while (fft_size < 2 * n) {
        if (fft_size > std::numeric_limits<std::size_t>::max() / 2) {
            throw std::length_error("FFT size overflow");
        }
        fft_size <<= 1;
    }

    std::vector<std::complex<double>> spectrum(fft_size);
    for (std::size_t i = 0; i < n; ++i) spectrum[i] = x[i] - mean;
    fft(spectrum, false);
    for (auto& value : spectrum) value = std::norm(value);
    fft(spectrum, true);

    // Convenzione tau_int = 1/2 + somma delle correlazioni normalizzate.
    // Per campioni indipendenti vale 1/2. Il taglio della somma e' una
    // stima pratica a finestra finita, non una misura esatta del mixing.
    // L'unita' e' l'intervallo tra campioni (sample_every tentativi).
    double tau = 0.5;
    const std::size_t max_lag = std::min<std::size_t>(n / 2, 10000);
    for (std::size_t lag = 1; lag <= max_lag; ++lag) {
        const double c = spectrum[lag].real() / static_cast<double>(n - lag);
        const double rho = c / c0;
        if (rho <= 0.0) break;
        tau += rho;
        if (static_cast<double>(lag) > 6.0 * tau) break;
    }
    return tau;
}

static RunResult run_one(int n, const Options& opt, std::uint64_t seed, bool save_series) {
    // Esegue una simulazione completa per una singola lunghezza N.
    PivotSaw2D saw(n);
    Xoshiro256StarStar rng(seed);
    RunResult result;
    result.n = n;
    result.attempts = opt.attempts;
    result.thermalization = opt.thermalization;

    std::vector<double> re2_series;
    std::vector<double> rg2_series;
    std::vector<double> rm2_series;
    std::vector<double> acc_series;
    // Disattivare l'autocorrelazione evita di conservare le serie temporali;
    // non rende indipendenti i campioni e non cambia le mosse della catena.
    if (opt.autocorrelation) {
        re2_series.reserve(static_cast<std::size_t>(opt.attempts / opt.sample_every + 1));
        rg2_series.reserve(re2_series.capacity());
        rm2_series.reserve(re2_series.capacity());
        acc_series.reserve(re2_series.capacity());
    }

    std::ofstream series_file;
    if (save_series && !opt.series_out.empty()) {
        // La serie temporale completa viene salvata solo quando richiesto,
        // normalmente per l'ultima lunghezza simulata.
        series_file.open(opt.series_out);
        if (!series_file) throw std::runtime_error("cannot open series output: " + opt.series_out);
        series_file << "attempt,accepted,free_path,free_path_fraction,Re2,Rg2,Rm2\n";
    }

    const auto t0 = std::chrono::steady_clock::now();
    for (std::uint64_t step = 1; step <= opt.thermalization + opt.attempts; ++step) {
        const PivotAttempt attempt = saw.attempt(rng);

        if (opt.validate_every > 0 && step % opt.validate_every == 0) {
            saw.validate_state();
        }

        // Si scarta il transiente della catena rettilinea iniziale.
        // La durata e' un parametro scelto dall'esterno: questo controllo
        // non dimostra automaticamente che l'equilibrio sia stato raggiunto.
        if (step <= opt.thermalization) continue;

        const int ci = class_index(kD4NonIdentity[static_cast<std::size_t>(attempt.symmetry_index)].cls);
        result.proposed_by_class[static_cast<std::size_t>(ci)]++;
        // Queste statistiche usano ogni tentativo di produzione,
        // indipendentemente da sample_every e dall'accettazione.
        result.attempt_free_path.push(attempt.free_path);
        if (attempt.tail_length > 0) {
            // Normalizzazione PRIMA della media: la coda N-k cambia
            // a ogni proposta. media(L/(N-k)) non e' media(L)/media(N-k).
            result.attempt_free_path_fraction.push(
                static_cast<long double>(attempt.free_path) / static_cast<long double>(attempt.tail_length));
        }
        if (attempt.accepted) {
            result.accepted++;
            result.accepted_by_class[static_cast<std::size_t>(ci)]++;
        } else {
            result.failed_free_path.push(attempt.free_path);
            if (attempt.tail_length > 0) {
                result.failed_free_path_fraction.push(
                    static_cast<long double>(attempt.free_path) / static_cast<long double>(attempt.tail_length));
            }
        }

        if ((step - opt.thermalization) % opt.sample_every == 0) {
            // Campionamento a intervalli fissi, anche dopo un rigetto:
            // la configurazione invariata deve mantenere il suo peso.
            // Il thinning riduce i campioni salvati, ma non garantisce
            // indipendenza statistica. Con sample_every=1 si misura sempre.
            const Observables obs = saw.observables();
            result.samples++;
            result.re2.push(obs.re2);
            result.rg2.push(obs.rg2);
            result.rm2.push(obs.rm2);
            result.accept_indicator.push(attempt.accepted ? 1.0L : 0.0L);

            if (opt.autocorrelation) {
                re2_series.push_back(static_cast<double>(obs.re2));
                rg2_series.push_back(static_cast<double>(obs.rg2));
                rm2_series.push_back(static_cast<double>(obs.rm2));
                acc_series.push_back(attempt.accepted ? 1.0 : 0.0);
            }

            if (series_file) {
                const double free_path_fraction = attempt.tail_length > 0
                    ? static_cast<double>(attempt.free_path) / attempt.tail_length
                    : 0.0;
                series_file << (step - opt.thermalization) << ',' << (attempt.accepted ? 1 : 0) << ','
                            << attempt.free_path << ',' << std::setprecision(17)
                            << free_path_fraction << ',' << static_cast<double>(obs.re2) << ','
                            << static_cast<double>(obs.rg2) << ','
                            << static_cast<double>(obs.rm2) << '\n';
            }
        }
    }
    if (opt.autocorrelation) {
        result.tau_re2 = integrated_autocorrelation_time(re2_series);
        result.tau_rg2 = integrated_autocorrelation_time(rg2_series);
        result.tau_rm2 = integrated_autocorrelation_time(rm2_series);
        result.tau_accept = integrated_autocorrelation_time(acc_series);
    }
    const auto t1 = std::chrono::steady_clock::now();
    result.seconds = std::chrono::duration<double>(t1 - t0).count();
    return result;
}

static void write_header(std::ofstream& out) {
    out << "N,attempts,thermalization,samples,accepted,acceptance,failed,failure,"
        << "Re2_mean,Re2_sd,Rg2_mean,Rg2_sd,Rm2_mean,Rm2_sd,"
        << "tau_Re2,tau_Rg2,tau_Rm2,tau_accept,"
        << "free_path_mean,free_path_fraction_mean,"
        << "failed_free_path_mean,failed_free_path_fraction_mean,"
        << "acc_rot90,acc_rot180,acc_axis_reflection,acc_diagonal_reflection,"
        << "fail_rot90,fail_rot180,fail_axis_reflection,fail_diagonal_reflection,"
        << "seconds\n";
}

static double ratio(std::uint64_t a, std::uint64_t b) {
    // Rapporto protetto: se il denominatore e' zero scrive NaN invece di dividere per zero.
    return b == 0 ? std::numeric_limits<double>::quiet_NaN()
                  : static_cast<double>(a) / static_cast<double>(b);
}

// Le colonne *_sd sono deviazioni standard delle configurazioni, NON
// errori sulle medie. Nel protocollo della tesi, i SEM del Pivot vengono
// stimati esternamente dalla dispersione delle medie fra repliche.
// Le accettazioni per classe sono condizionate alle proposte di quella
// classe: le quattro classi non hanno tutte la stessa molteplicita'.
static void write_result(std::ofstream& out, const RunResult& r) {
    out << r.n << ',' << r.attempts << ',' << r.thermalization << ',' << r.samples << ','
        << r.accepted << ',' << std::setprecision(17) << ratio(r.accepted, r.attempts) << ','
        << (r.attempts - r.accepted) << ',' << ratio(r.attempts - r.accepted, r.attempts) << ','
        << static_cast<double>(r.re2.mean) << ',' << std::sqrt(static_cast<double>(r.re2.variance())) << ','
        << static_cast<double>(r.rg2.mean) << ',' << std::sqrt(static_cast<double>(r.rg2.variance())) << ','
        << static_cast<double>(r.rm2.mean) << ',' << std::sqrt(static_cast<double>(r.rm2.variance())) << ','
        << r.tau_re2 << ',' << r.tau_rg2 << ',' << r.tau_rm2 << ',' << r.tau_accept << ','
        << static_cast<double>(r.attempt_free_path.mean) << ','
        << static_cast<double>(r.attempt_free_path_fraction.mean) << ','
        << static_cast<double>(r.failed_free_path.mean) << ','
        << static_cast<double>(r.failed_free_path_fraction.mean) << ','
        << ratio(r.accepted_by_class[0], r.proposed_by_class[0]) << ','
        << ratio(r.accepted_by_class[1], r.proposed_by_class[1]) << ','
        << ratio(r.accepted_by_class[2], r.proposed_by_class[2]) << ','
        << ratio(r.accepted_by_class[3], r.proposed_by_class[3]) << ','
        << ratio(r.proposed_by_class[0] - r.accepted_by_class[0], r.proposed_by_class[0]) << ','
        << ratio(r.proposed_by_class[1] - r.accepted_by_class[1], r.proposed_by_class[1]) << ','
        << ratio(r.proposed_by_class[2] - r.accepted_by_class[2], r.proposed_by_class[2]) << ','
        << ratio(r.proposed_by_class[3] - r.accepted_by_class[3], r.proposed_by_class[3]) << ','
        << r.seconds << '\n';
}

int main(int argc, char** argv) {
    try {
        const Options opt = parse_args(argc, argv);
        std::ofstream out(opt.out);
        if (!out) throw std::runtime_error("cannot open output: " + opt.out);
        write_header(out);

        for (std::size_t i = 0; i < opt.lengths.size(); ++i) {
            const int n = opt.lengths[i];
            // Ogni lunghezza riceve un seed derivato ma diverso, cosi' le run
            // non riusano esattamente la stessa sequenza casuale.
            const std::uint64_t seed = opt.seed + 0x9e3779b97f4a7c15ULL * (i + 1);
            if (opt.verbose) {
                std::cerr << "Running N=" << n << " attempts=" << opt.attempts
                          << " thermalization=" << opt.thermalization << "\n";
            }
            const bool save_series = (i + 1 == opt.lengths.size());
            const RunResult r = run_one(n, opt, seed, save_series);
            write_result(out, r);
            out.flush();
            if (opt.verbose) {
                std::cerr << "  acceptance=" << ratio(r.accepted, r.attempts)
                          << " Re2=" << static_cast<double>(r.re2.mean)
                          << " Rg2=" << static_cast<double>(r.rg2.mean)
                          << " seconds=" << r.seconds << "\n";
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
