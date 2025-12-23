#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>

const double PI = 3.14159265358979323846;

struct Params {
    double v;
    double nu;
    double alpha0;
    double K;
    int H;
    double tau;
    double D;
    double n_bar;
    double dt_exo;
    double dt_follower;
    double T;
    double detection_radius;
    double molecule_lifetime;
    std::string file_suffix;
};

void run_snapshot_simulation(const Params& params) {
    double v_source = params.v;
    double nu = params.nu;
    double alpha0 = params.alpha0;
    double K = params.K;
    int H = params.H;
    double tau = params.tau;
    double D = params.D;
    double dt_exo = params.dt_exo;
    double dt_follower = params.dt_follower;
    double T = params.T;
    double detection_radius = params.detection_radius;
    double molecule_lifetime = params.molecule_lifetime;
    double n_bar = params.n_bar;

    int num_steps = static_cast<int>(T / dt_exo);
    int follower_interval = static_cast<int>(dt_follower / dt_exo);
    double dx_molecule = std::sqrt(4 * D * dt_exo);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_real_distribution<> angle_dis(0.0, 2 * PI);
    std::poisson_distribution<> poisson_n_bar(n_bar);

    double source_x = 10.0, source_y = 0.0;
    double follower_x = 0.0, follower_y = 0.0;

    std::vector<std::tuple<double, double, double, int>> particles;
    std::vector<std::tuple<double, double, int>> detection_events;
    std::vector<std::pair<double, double>> follower_trajectory;
    follower_trajectory.emplace_back(follower_x, follower_y);

    for (int step = 0; step < num_steps; ++step) {
        double current_time = step * dt_exo;

        source_x += v_source * dt_exo;

        double expected_exos = (nu / n_bar) * dt_exo;
        std::poisson_distribution<> poisson_nu_exo(expected_exos);
        int num_new = poisson_nu_exo(gen);

        for (int i = 0; i < num_new; ++i) {
            int n = poisson_n_bar(gen);
            if (n > 0) {
                particles.emplace_back(source_x, source_y, current_time, n);
            }
        }

        // particles.erase(std::remove_if(particles.begin(), particles.end(),
        //     [current_time, molecule_lifetime](const auto& p) {
        //         return current_time - std::get<2>(p) >= molecule_lifetime;
        //     }), particles.end());

        if (D > 0.0) {
            for (auto& p : particles) {
                double theta = angle_dis(gen);
                std::get<0>(p) += std::cos(theta) * dx_molecule;
                std::get<1>(p) += std::sin(theta) * dx_molecule;
            }
        }

        auto it = particles.begin();
        while (it != particles.end()) {
            double dx = std::get<0>(*it) - follower_x;
            double dy = std::get<1>(*it) - follower_y;
            if (std::sqrt(dx * dx + dy * dy) < detection_radius) {
                double theta = std::atan2(dy, dx);
                detection_events.emplace_back(current_time, theta, std::get<3>(*it));
                it = particles.erase(it);
            } else {
                ++it;
            }
        }

        // NOTE: Do not hard-prune detection events by an absolute cutoff time.
        // Memory is implemented via exponential decay (exp(-delta_t / tau)).
        // Removing the strict time-based erase avoids double-thresholding memory
        // and ensures old events still contribute (with small weight) via decay.

        if (step % follower_interval == 0) {
            double sum_alpha = 0.0, net_x = 0.0, net_y = 0.0;
            for (const auto& [t_i, theta_i, n_i] : detection_events) {
                double delta_t = current_time - t_i;
                double hill = std::pow(n_i, H) / (std::pow(n_i, H) + std::pow(K, H));
                double alpha_i = alpha0 * hill * std::exp(-delta_t / tau);
                sum_alpha += alpha_i;
                net_x += alpha_i * std::cos(theta_i);
                net_y += alpha_i * std::sin(theta_i);
            }

            double net_theta = std::atan2(net_y, net_x);
            double p_theta = sum_alpha / (1.0 + sum_alpha);
            double move_theta = (dis(gen) < p_theta && (net_x != 0.0 || net_y != 0.0)) ?
                net_theta : angle_dis(gen);

            follower_x += std::cos(move_theta) * v_source * dt_follower;
            follower_y += std::sin(move_theta) * v_source * dt_follower;
            follower_trajectory.emplace_back(follower_x, follower_y);
        }
    }

    std::ofstream traj("follower_trajectory_" + params.file_suffix + ".txt");
    for (const auto& [x, y] : follower_trajectory)
        traj << x << " " << y << "\n";
    traj.close();

    std::ofstream exo("exosomes_final_" + params.file_suffix + ".txt");
    for (const auto& [x, y, t, n] : particles)
        exo << x << " " << y << "\n";
    exo.close();

    std::ofstream leader("leader_final_" + params.file_suffix + ".txt");
    leader << source_x << " " << source_y << "\n";
    leader.close();

    std::cout << "Snapshot complete for D = " << D << "\n";
}

int main() {
    Params base;
    base.v = 3.4 / 60;
    base.nu = 300.0;
    base.alpha0 = 1.0;
    base.K = 25;
    base.H = 3;
    base.tau = 5.0;
    base.dt_exo = 0.1;
    base.dt_follower = 40.0;
    base.T = 1440.0;
    base.detection_radius = 10.0;
    // base.molecule_lifetime = 7.0;
    base.n_bar = 25.0;

    std::vector<double> D_values = {0.0, 50.0, 250.0, 500.0, 1000.0};
    for (double D_val : D_values) {
        base.D = D_val;
        base.file_suffix = std::to_string(static_cast<int>(D_val));
        run_snapshot_simulation(base);
    }

    return 0;
}


