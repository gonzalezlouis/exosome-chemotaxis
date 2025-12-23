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
    double nu;                   // Molecules per minute
    double alpha0;               // Signal strength
    double K;                    // Hill threshold
    double tau;                  // Memory decay window (min)
    double dt_exo;               // Exosome update time step (min)
    double dt_follower;          // Follower update time step (min)
    double detection_radius;     // Detection radius (μm)
    int num_trials;              // Trials per condition
};

double simulate_trial_2d(const Params& params, double n_bar, double D, double T, int H) {
    int num_steps = static_cast<int>(T / params.dt_exo);
    int follower_steps = static_cast<int>(params.dt_follower / params.dt_exo);
    double dx_mol = std::sqrt(4 * D * params.dt_exo);

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> uni(0.0, 1.0);
    std::uniform_real_distribution<> ang(0.0, 2 * PI);
    std::poisson_distribution<> pd_nbar(n_bar);

    double source_x = 10.0, source_y = 0.0;
    double follower_x = 0.0, follower_y = 0.0;

    std::vector<std::tuple<double,double,double,int>> particles;
    std::vector<std::tuple<double,double,int>> detections;

    for (int step = 0; step < num_steps; ++step) {
        double t = step * params.dt_exo;

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

        // Diffuse all particles
        for (auto& p : particles) {
            double theta = ang(gen);
            std::get<0>(p) += dx_mol * std::cos(theta);
            std::get<1>(p) += dx_mol * std::sin(theta);
        }

        // Detect and remove particles
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

        // NOTE: Do not hard-prune detection events using an absolute cutoff time.
        // Memory is implemented via exponential decay (exp(-dt / params.tau)).
        // Keeping older detection events and letting the decay factor down-weight
        // their contribution provides a soft memory and avoids double-thresholding.

        // Move follower every dt_follower
        if (step % follower_steps == 0) {
            double sum_alpha = 0.0, net_x = 0.0, net_y = 0.0;
            for (auto& [ti, theta, n] : detections) {
                double dt = t - ti;
                double hill = std::pow(n, H) / (std::pow(n, H) + std::pow(params.K, H));
                double alpha = params.alpha0 * hill * std::exp(-dt / params.tau);
                sum_alpha += alpha;
                net_x += alpha * std::cos(theta);
                net_y += alpha * std::sin(theta);
            }
            double p = sum_alpha / (1.0 + sum_alpha);
            double move_theta = (uni(gen) < p && (net_x || net_y))
                           ? std::atan2(net_y, net_x) : ang(gen);
            follower_x += params.v * params.dt_follower * std::cos(move_theta);
            follower_y += params.v * params.dt_follower * std::sin(move_theta);
        }
    }

    return (follower_x - 0.0) / T;
}

void run_panel4c_comparison(const Params& base, const std::vector<double>& nbar_values, double T) {
    std::ofstream out("results_protocol_comparison.txt");
    const int H = 3;
    const double n0 = 1.0;
    const double D0 = 250.0;
    // Progress helpers
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

    size_t total_tasks = nbar_values.size() * 2 * static_cast<size_t>(base.num_trials);
    size_t completed = 0;
    auto start_time = std::chrono::steady_clock::now();

    for (double n_bar : nbar_values) {
        out << "n_bar = " << n_bar << ":\n";

        // Protocol A: Constant D
        {
            std::vector<double> resultsA;
            for (int i = 0; i < base.num_trials; ++i) {
                resultsA.push_back(simulate_trial_2d(base, n_bar, D0, T, H));
                ++completed;
                print_progress(static_cast<double>(completed) / total_tasks, start_time);
            }

            double mean = std::accumulate(resultsA.begin(), resultsA.end(), 0.0) / resultsA.size();
            double sq = std::inner_product(resultsA.begin(), resultsA.end(), resultsA.begin(), 0.0);
            double sd = std::sqrt(sq / resultsA.size() - mean * mean);
            double se = sd / std::sqrt(resultsA.size());
            out << "ConstantD " << mean << " " << se << "\n";
        }

        // Protocol B: Stokes-Einstein D
        {
            double D_scaled = D0 * std::pow(n0 / n_bar, 1.0 / 3.0);
            std::vector<double> resultsB;
            for (int i = 0; i < base.num_trials; ++i) {
                resultsB.push_back(simulate_trial_2d(base, n_bar, D_scaled, T, H));
                ++completed;
                print_progress(static_cast<double>(completed) / total_tasks, start_time);
            }

            double mean = std::accumulate(resultsB.begin(), resultsB.end(), 0.0) / resultsB.size();
            double sq = std::inner_product(resultsB.begin(), resultsB.end(), resultsB.begin(), 0.0);
            double sd = std::sqrt(sq / resultsB.size() - mean * mean);
            double se = sd / std::sqrt(resultsB.size());
            out << "StokesEinsteinD " << mean << " " << se << "\n";
        }

        out << "\n";
    }
    // Ensure progress bar finishes on its own line
    std::cout << std::endl;

    std::cout << "Results saved to results_protocol_comparison.txt\n";
}

int main() {
    Params params;
    params.v                = 3.4 / 60.0;
    params.nu               = 350.0;
    params.alpha0           = 1.0;
    params.K                = 25.0;
    params.tau              = 5.0;
    params.dt_exo           = 0.1;
    params.dt_follower      = 40.0;
    params.detection_radius = 10.0;
    params.num_trials       = 50;

    std::vector<double> nbar_values = {1, 10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250};
    double T = 600.0;

    auto t0 = std::chrono::high_resolution_clock::now();
    run_panel4c_comparison(params, nbar_values, T);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Completed in "
              << std::chrono::duration<double>(t1 - t0).count()
              << " s\n";
    return 0;
}
