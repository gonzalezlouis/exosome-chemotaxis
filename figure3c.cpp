#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>

const double PI = 3.14159265358979323846;

struct Params {
    double v;
    double nu;
    double alpha0;
    double K;
    double tau;
    double D;
    double dt_exo;
    double dt_follower;
    double T;
    double detection_radius;
    double molecule_lifetime;
    int num_trials;
    std::vector<double> n_bar_values;
};

std::vector<double> linspace(double start, double stop, int num) {
    std::vector<double> result(num);
    double step = (stop - start) / (num - 1);
    for (int i = 0; i < num; ++i) result[i] = start + i * step;
    return result;
}

void print_progress(double progress, std::chrono::steady_clock::time_point start) {
    const int barWidth = 50;
    int pos = static_cast<int>(barWidth * progress);

    auto now = std::chrono::steady_clock::now();
    double elapsed_sec = std::chrono::duration<double>(now - start).count();
    double eta_sec = (progress > 0.0) ? elapsed_sec * (1.0 / progress - 1.0) : 0.0;

    int eta_min = static_cast<int>(eta_sec / 60);
    int eta_s   = static_cast<int>(std::round(eta_sec)) % 60;

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0)
              << "% (ETA: " << eta_min << "m " << eta_s << "s)\r";
    std::cout.flush();
}

double simulate_trial_2d(const Params& params, double n_bar, int H, std::mt19937_64& gen) {
    int total_steps = static_cast<int>(params.T / params.dt_exo);
    int follower_interval = static_cast<int>(params.dt_follower / params.dt_exo);
    double dx_molecule = std::sqrt(4 * params.D * params.dt_exo);

    std::uniform_real_distribution<> uniform01(0.0, 1.0);
    std::uniform_real_distribution<> angle_dis(0.0, 2 * PI);
    std::poisson_distribution<> poisson_n_bar(n_bar);

    double source_x = 10.0, source_y = 0.0;
    double follower_x = 0.0, follower_y = 0.0;

    std::vector<std::tuple<double, double, double, int>> particles;
    std::vector<std::tuple<double, double, int>> detection_events;

    for (int step = 0; step < total_steps; ++step) {
        double t = step * params.dt_exo;
        source_x += params.v * params.dt_exo;

        // Exosome emission
        double mean_emit = (params.nu / n_bar) * params.dt_exo;
        std::poisson_distribution<> poisson_emit(mean_emit);
        int n_emit = poisson_emit(gen);
        for (int i = 0; i < n_emit; ++i) {
            int n = poisson_n_bar(gen);
            if (n > 0)
                particles.emplace_back(source_x, source_y, t, n);
        }

        // // Remove old particles
        // particles.erase(std::remove_if(particles.begin(), particles.end(),
        //     [t, &params](const auto& p) {
        //         return t - std::get<2>(p) >= params.molecule_lifetime;
        //     }), particles.end());

        // Diffusion
        for (auto& p : particles) {
            double theta = angle_dis(gen);
            std::get<0>(p) += dx_molecule * std::cos(theta);
            std::get<1>(p) += dx_molecule * std::sin(theta);
        }

        // Detection
        auto it = particles.begin();
        while (it != particles.end()) {
            double dx = std::get<0>(*it) - follower_x;
            double dy = std::get<1>(*it) - follower_y;
            if (std::hypot(dx, dy) < params.detection_radius) {
                double theta = std::atan2(dy, dx);
                detection_events.emplace_back(t, theta, std::get<3>(*it));
                it = particles.erase(it);
            } else {
                ++it;
            }
        }

        // Follower motion every dt_follower
        if (step % follower_interval == 0) {
            double sum_alpha = 0.0, net_x = 0.0, net_y = 0.0;

            // Do not hard-prune detection events by a fixed cutoff time.
            // Memory is implemented via exponential decay (exp(-dt / tau)).
            // Keeping older events and letting the decay factor down-weight them
            // provides a soft memory and avoids double-thresholding.

            for (const auto& [ti, theta, n] : detection_events) {
                double dt = t - ti;
                double hill = std::pow(n, H) / (std::pow(n, H) + std::pow(params.K, H));
                double alpha = params.alpha0 * hill * std::exp(-dt / params.tau);
                sum_alpha += alpha;
                net_x += alpha * std::cos(theta);
                net_y += alpha * std::sin(theta);
            }

            double move_theta;
            if ((net_x != 0.0 || net_y != 0.0) && uniform01(gen) < (sum_alpha / (1.0 + sum_alpha))) {
                move_theta = std::atan2(net_y, net_x);
            } else {
                move_theta = angle_dis(gen);
            }

            follower_x += params.v * params.dt_follower * std::cos(move_theta);
            follower_y += params.v * params.dt_follower * std::sin(move_theta);
        }
    }

    return follower_x / params.T;
}

void run_velocity_scan(const Params& params, int H, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << "\n";
        return;
    }

    out << "K = " << params.K << "\n";
    out << "H = " << H << "\n";
    out << "tau = " << params.tau << "\n";
    out << "nu = " << params.nu << "\n";
    out << std::fixed << std::setprecision(6);

    std::mt19937_64 gen(std::random_device{}());

    size_t total = params.n_bar_values.size();
    size_t current = 0;
    auto start_time = std::chrono::steady_clock::now();

    for (double n_bar : params.n_bar_values) {
        double sum_v = 0.0;
        for (int trial = 0; trial < params.num_trials; ++trial) {
            sum_v += simulate_trial_2d(params, n_bar, H, gen);
        }
        double avg_v = sum_v / params.num_trials;
        out << n_bar << " " << avg_v << "\n";

        current++;
        print_progress(static_cast<double>(current) / total, start_time);
    }

    std::cout << std::endl;
    out.close();
}

int main() {
    Params params;
    params.v = 3.4 / 60.0;
    params.alpha0 = 1.0;
    params.K = 25;
    params.tau = 5.0;
    params.D = 250.0;
    params.dt_exo = 0.1;
    params.dt_follower = 40.0;
    params.T = 450.0;
    params.detection_radius = 10.0;
    params.molecule_lifetime = 7.0;
    params.num_trials = 50;
    params.n_bar_values = linspace(1.0, 350.0, 75);

    std::vector<double> nu_values = {60.0, 100.0, 250.0, 600.0};

    for (double nu : nu_values) {
        params.nu = nu;
        std::string filename = "results_2d_nu" + std::to_string(static_cast<int>(nu)) + ".txt";
        std::cout << "Running nu = " << nu << "...\n";
        run_velocity_scan(params, 3, filename);
    }

    return 0;
}


