#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <string>

const double PI = 3.14159265358979323846;

struct Params {
    double v;                    // Source & follower speed
    double nu;                   // Molecules per minute
    double alpha0;               // Signal strength
    double K;                    // Hill threshold
    double tau;                  // Memory window
    double fixed_n_bar;          // Mean molecules per exosome
    double D;                    // Diffusion coefficient (varying)
    double dt_exo;               // Exosome update step
    double dt_follower;          // Follower update step
    double detection_radius;     // Detection radius
    double molecule_lifetime;    // Molecule lifetime
    int num_trials;              // Repetitions
    std::vector<double> D_values;
};

std::vector<double> linspace(double start, double stop, int num) {
    std::vector<double> result(num);
    double step = (stop - start) / (num - 1);
    for (int i = 0; i < num; ++i)
        result[i] = start + i * step;
    return result;
}

// Added wall_half_width argument. Pass -1.0 for no walls.
double simulate_trial_2d(const Params& params, double T, int H, double wall_half_width) {
    int num_steps      = static_cast<int>(T / params.dt_exo);
    int follower_steps = static_cast<int>(params.dt_follower / params.dt_exo);
    double dx_mol      = std::sqrt(4 * params.D * params.dt_exo);

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> uni(0.0, 1.0);
    std::uniform_real_distribution<> ang(0.0, 2 * PI);
    std::poisson_distribution<> pd_nbar(params.fixed_n_bar);

    double source_x   = 10.0, source_y = 0.0;
    double follower_x =  0.0, follower_y = 0.0;

    // x, y, emission_time, molecule_count
    std::vector<std::tuple<double,double,double,int>> particles;
    // time, angle, molecule_count
    std::vector<std::tuple<double,double,int>> detections;

    for (int step = 0; step < num_steps; ++step) {
        double t = step * params.dt_exo;

        // 1. Move Source
        source_x += params.v * params.dt_exo;

        // 2. Emit Exosomes
        double mean_emit = (params.nu / params.fixed_n_bar) * params.dt_exo;
        std::poisson_distribution<> pd_emit(mean_emit);
        int n_emit = pd_emit(gen);
        for (int i = 0; i < n_emit; ++i) {
            int n = pd_nbar(gen);
            if (n > 0)
                particles.emplace_back(source_x, source_y, t, n);
        }

        // 3. Diffuse & Reflect
        for (auto& p : particles) {
            double theta = ang(gen);
            std::get<0>(p) += dx_mol * std::cos(theta);
            std::get<1>(p) += dx_mol * std::sin(theta);

            // --- REFLECTIVE WALLS LOGIC ---
            if (wall_half_width > 0.0) {
                double& y = std::get<1>(p);
                // Reflect off Top Wall (+W)
                if (y > wall_half_width) {
                    y = 2 * wall_half_width - y;
                } 
                // Reflect off Bottom Wall (-W)
                else if (y < -wall_half_width) {
                    y = -2 * wall_half_width - y;
                }
            }
        }

        // 4. Detect & Remove
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

        // 5. Move Follower
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
            double p = sum_alpha / (1.0 + sum_alpha);
            double move_theta;
            
            if (uni(gen) < p && (std::abs(net_x) > 1e-9 || std::abs(net_y) > 1e-9)) {
                move_theta = std::atan2(net_y, net_x);
            } else {
                move_theta = ang(gen);
            }

            follower_x += params.v * params.dt_follower * std::cos(move_theta);
            follower_y += params.v * params.dt_follower * std::sin(move_theta);
        }
    }

    return follower_x / T; // Avg velocity
}

void run_simulation(const Params& base) {
    std::ofstream out("results_walls.txt");
    out << "n_bar = " << base.fixed_n_bar << "\n";

    const int H = 3;
    const double T = 1000.0; // Fixed Time for comparison
    const double cell_length = 20.0; // Diameter (2 * Radius)

    // Distances: Control (-1), 20L, 10L, 5L
    // Note: distance is from center, so wall_half_width = distance * cell_length
    std::vector<double> wall_distances_L = {-1.0, 100.0, 10.0, 5.0}; 

    // Progress bar helper
    auto print_progress = [](double progress) {
        int barWidth = 50;
        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
    };

    size_t total_tasks = wall_distances_L.size() * base.D_values.size() * base.num_trials;
    size_t completed = 0;

    for (double L : wall_distances_L) {
        double wall_half_width = (L < 0) ? -1.0 : (L * cell_length);
        
        if (L < 0) out << "Wall Dist = Infinity:\n";
        else       out << "Wall Dist = " << L << " cells:\n";

        for (double D : base.D_values) {
            Params p = base;
            p.D = D;
            std::vector<double> vs;
            for (int i = 0; i < p.num_trials; ++i) {
                vs.push_back(simulate_trial_2d(p, T, H, wall_half_width));
                ++completed;
                if (completed % 10 == 0) print_progress((double)completed / total_tasks);
            }

            double mean = std::accumulate(vs.begin(), vs.end(), 0.0) / vs.size();
            double sq = std::inner_product(vs.begin(), vs.end(), vs.begin(), 0.0);
            double sd = std::sqrt(sq/vs.size() - mean*mean);
            double se = sd / std::sqrt(vs.size());
            out << D << " " << mean << " " << se << "\n";
        }
        out << "\n";
    }
    std::cout << std::endl << "Done. Saved to results_walls.txt" << std::endl;
}

int main() {
    Params params;
    params.v                = 3.4/60; // 3.4 um/min
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

    run_simulation(params);
    return 0;
}