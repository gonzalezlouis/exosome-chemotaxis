#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <fstream>

const double PI = 3.14159265358979323846;

struct Params {
    double v;                    // Source & follower speed (μm/min)
    double nu;                   // Molecules per minute (fixed)
    double alpha0;               // Signal strength
    double K;                    // Hill threshold
    double tau;                  // Memory window (min)
    double fixed_n_bar;          // Mean molecules per exosome
    double D;                    // Diffusion coefficient (μm^2/min)
    double dt_exo;               // Exosome update time step (min)
    double dt_follower;          // Follower update time step (min)
    double detection_radius;     // Detection radius (μm)
    double molecule_lifetime;    // Molecule lifetime (min)
    int num_trials;              // Repetitions per condition
    std::vector<double> D_values;
    std::vector<double> T_values;
};

std::vector<double> linspace(double start, double stop, int num) {
    std::vector<double> result(num);
    double step = (stop - start) / (num - 1);
    for (int i = 0; i < num; ++i)
        result[i] = start + i * step;
    return result;
}

double simulate_trial_2d(const Params& params, double n_bar, double T, int H) {
    int num_steps      = static_cast<int>(T / params.dt_exo);
    int follower_steps = static_cast<int>(params.dt_follower / params.dt_exo);
    double dx_mol      = std::sqrt(4 * params.D * params.dt_exo);

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> uni(0.0, 1.0);
    std::uniform_real_distribution<> ang(0.0, 2 * PI);
    std::poisson_distribution<> pd_nbar(n_bar);

    double source_x   = 10.0, source_y = 0.0;
    double follower_x =  0.0, follower_y = 0.0;

    std::vector<std::tuple<double,double,double,int>> particles;
    std::vector<std::tuple<double,double,int>> detections;

    for (int step = 0; step < num_steps; ++step) {
        double t = step * params.dt_exo;

        // Move source
        source_x += params.v * params.dt_exo;

        // Emit exosomes
        double mean_emit = (params.nu / n_bar) * params.dt_exo;
        std::poisson_distribution<> pd_emit(mean_emit);
        int n_emit = pd_emit(gen);
        for (int i = 0; i < n_emit; ++i) {
            int n = pd_nbar(gen);
            if (n > 0)
                particles.emplace_back(source_x, source_y, t, n);
        }

        // // Remove old particles
        // particles.erase(std::remove_if(particles.begin(), particles.end(),
        //     [t,&params](auto& p){ return t - std::get<2>(p) >= params.molecule_lifetime; }),
        //     particles.end());

        // Diffuse
        for (auto& p : particles) {
            double theta = ang(gen);
            std::get<0>(p) += dx_mol * std::cos(theta);
            std::get<1>(p) += dx_mol * std::sin(theta);
        }

        // Detect & remove
        for (auto it = particles.begin(); it != particles.end(); ) {
            double dx = std::get<0>(*it) - follower_x;
            double dy = std::get<1>(*it) - follower_y;
            if (std::hypot(dx, dy) < params.detection_radius) {
                double theta = std::atan2(dy, dx);
                detections.emplace_back(t, theta, std::get<3>(*it));
                it = particles.erase(it);
            } else {
                ++it;
            }
        }

        // NOTE: Do not hard-remove detections by an absolute cutoff time.
        // Memory is implemented via exponential decay (exp(-dt/params.tau)).
        // Keeping older detection events and letting the decay factor down-weight
        // their contributions provides a soft memory and avoids double-thresholding.

        // Update follower only every dt_follower
        if (step % follower_steps == 0) {
            double sum_alpha = 0.0, net_x = 0.0, net_y = 0.0;
            for (auto& [ti, theta, n] : detections) {
                double dt = t - ti;
                double hill = std::pow(n, H) / (std::pow(n, H) + std::pow(params.K, H));
                double alpha = params.alpha0 * hill * std::exp(-dt/params.tau);
                sum_alpha += alpha;
                net_x += alpha * std::cos(theta);
                net_y += alpha * std::sin(theta);
            }
            double p = sum_alpha / (1.0 + sum_alpha),
                   move_theta = (uni(gen) < p && (net_x||net_y)) 
                            ? std::atan2(net_y, net_x) : ang(gen);
            follower_x += params.v * params.dt_follower * std::cos(move_theta);
            follower_y += params.v * params.dt_follower * std::sin(move_theta);
        }
    }

    return (follower_x - 0.0) / T;
}

void run_simulation(const Params& base) {
    std::ofstream out("results_vs_D_T.txt");
    out << "n_bar = " << base.fixed_n_bar << "\n";

    const int H = 3;
    // progress helper
    auto print_progress = [](double progress, std::chrono::steady_clock::time_point start) {
        const int barWidth = 50;
        int pos = static_cast<int>(barWidth * progress);
        auto now = std::chrono::steady_clock::now();
        double elapsed_sec = std::chrono::duration<double>(now - start).count();
        double eta_sec = (progress > 0.0) ? elapsed_sec * (1.0 / progress - 1.0) : 0.0;
        int eta_min = static_cast<int>(eta_sec / 60);
        int eta_s = static_cast<int>(std::round(eta_sec)) % 60;
        std::cout << "[";
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0)
                  << "% (ETA: " << eta_min << "m " << eta_s << "s)\r";
        std::cout.flush();
    };

    size_t total_tasks = base.T_values.size() * base.D_values.size() * static_cast<size_t>(base.num_trials);
    size_t completed = 0;
    auto start_time = std::chrono::steady_clock::now();

    for (double T : base.T_values) {
        out << "T = " << T << ":\n";
        for (double D : base.D_values) {
            Params p = base;
            p.D = D;
            std::vector<double> vs;
            for (int i = 0; i < p.num_trials; ++i) {
                vs.push_back(simulate_trial_2d(p, p.fixed_n_bar, T, H));
                ++completed;
                print_progress(static_cast<double>(completed) / total_tasks, start_time);
            }

            double mean = std::accumulate(vs.begin(), vs.end(), 0.0) / vs.size();
            double sq = std::inner_product(vs.begin(), vs.end(), vs.begin(), 0.0);
            double sd = std::sqrt(sq/vs.size() - mean*mean);
            double se = sd / std::sqrt(vs.size());
            out << D << " " << mean << " " << se << "\n";
        }
        out << "\n";
    }
    std::cout << std::endl;
    std::cout << "Results saved to results_vs_D_T.txt\n";
}

int main() {
    Params params;
    params.v                = 3.4/60;
    params.nu               = 350.0;
    params.alpha0           = 1.0;
    params.K                = 25.0;
    params.tau              = 5.0;
    params.fixed_n_bar      = 50.0;
    params.dt_exo           = 0.1;
    params.dt_follower      = 50.0;
    params.molecule_lifetime= 7.0;
    params.detection_radius = 10.0;
    params.num_trials       = 35;
    params.D_values         = linspace(0.0, 1000.0, 50);
    params.T_values         = {180, 300, 600, 1000, 1440};

    auto t0 = std::chrono::high_resolution_clock::now();
    run_simulation(params);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Completed in "
              << std::chrono::duration<double>(t1-t0).count()
              << " s\n";
    return 0;
}