#ifndef KALMAN_HPP
#define KALMAN_HPP

#include <iostream>
#include "ft_vec.hpp"
#include "ft_mat.hpp"
#include <iomanip>
#include <sys/time.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

typedef ft::vector<double> vec;
typedef ft::matrix<double> mat;

static const size_t kStateDim = 9;

static vec forward_from_euler(const vec &dir) {
    const double pitch = dir[1];
    const double yaw = dir[2];
    return vec({
        std::cos(pitch) * std::cos(yaw),
        std::cos(pitch) * std::sin(yaw),
        -std::sin(pitch)
    });
}

static vec project_on_dir(const vec &v, const vec &dir_unit) {
    return dir_unit * dot(v, dir_unit);
}

struct State {

    mat v;
    mat P;

    State() : v(kStateDim, 1), P(kStateDim) {}

    State(const vec &pos, double speed, const vec &acc, const vec &dir) : v(kStateDim, 1), P(kStateDim) {
        const vec dir_unit = forward_from_euler(dir).normalize();
        vec vel = dir_unit * (speed  * 1000 / 3600);//km/h en m/s
        vec acc_aligned = project_on_dir(acc, dir_unit);
        vec tmp({pos[0], pos[1], pos[2],
                 vel[0], vel[1], vel[2],
                 acc_aligned[0], acc_aligned[1], acc_aligned[2]});
        v = mat(tmp).transposed();
    }

    vec pos() {return vec({v[0][0], v[1][0], v[2][0]});}

    vec vel() {return vec({v[3][0], v[4][0], v[5][0]});}

    vec acc() {return vec({v[6][0], v[7][0], v[8][0]});}
};


class KalmanFilter {
public:
    mat F; // state transition matrix

    mat Q; // process noise covariance

    mat H_pos; // measurement matrix
    mat H_acc;

    mat R_pos; // measurement noise covariance
    mat R_acc;
    double acc_var;
    double gyr_var;
    double gps_var_base;

    State state;

    double thresh;
    double alpha;

    bool print;
    std::ofstream variances;
    size_t gps_updates;
    size_t gps_scaled;
    size_t gps_gated;
    double gps_last_scale;
    double nis_sum;
    double nis_max;
    double gps_res_sum;
    double gps_res_max;
    bool initialized;


    KalmanFilter(double process_noise = 1.0, double n = 1, bool print = false,
                    double accNoise = 1e-4, double gyrNoise = 1e-4, double gpsNoise = 0.01,
                    double thresh = 100, double alpha = 0.01
    ) : F(kStateDim), Q(kStateDim),
    H_pos(3, kStateDim, 0), H_acc(3, kStateDim, 0),
    R_pos(3), R_acc(3), acc_var(0.0), gyr_var(0.0), gps_var_base(0.0),
    thresh(thresh), alpha(alpha), print(print),
    gps_updates(0), gps_scaled(0), gps_gated(0), gps_last_scale(1.0),
    nis_sum(0.0), nis_max(0.0), gps_res_sum(0.0), gps_res_max(0.0), initialized(false) {
        H_pos = mat({
            {1, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 1, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 1, 0, 0, 0, 0, 0, 0}
        });
        H_acc = mat({
            {0, 0, 0, 0, 0, 0, 1, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 1, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 1}
        });

        const double tol = (thresh > 0.0) ? thresh : 1.0;
        gps_var_base = gpsNoise * n * n * (5.0 / tol) * (5.0 / tol);
        R_pos *= gps_var_base;

        acc_var = accNoise * n * n;
        R_acc *= acc_var;
        gyr_var = gyrNoise;

        const double dt = 0.01;
        const double q = acc_var * process_noise;
        const double dt2 = dt * dt;
        const double dt3 = dt2 * dt;
        const double dt4 = dt2 * dt2;
        const double dt5 = dt4 * dt;
        const double q11 = q * dt5 / 20.0;
        const double q12 = q * dt4 / 8.0;
        const double q13 = q * dt3 / 6.0;
        const double q22 = q * dt3 / 3.0;
        const double q23 = q * dt2 / 2.0;
        const double q33 = q * dt;
        Q = mat({
            {q11, 0,   0,   q12, 0,   0,   q13, 0,   0},
            {0,   q11, 0,   0,   q12, 0,   0,   q13, 0},
            {0,   0,   q11, 0,   0,   q12, 0,   0,   q13},
            {q12, 0,   0,   q22, 0,   0,   q23, 0,   0},
            {0,   q12, 0,   0,   q22, 0,   0,   q23, 0},
            {0,   0,   q12, 0,   0,   q22, 0,   0,   q23},
            {q13, 0,   0,   q23, 0,   0,   q33, 0,   0},
            {0,   q13, 0,   0,   q23, 0,   0,   q33, 0},
            {0,   0,   q13, 0,   0,   q23, 0,   0,   q33}
        });

        F = mat({
            {1, 0, 0, dt, 0, 0, 0.5 * dt * dt, 0, 0},
            {0, 1, 0, 0, dt, 0, 0, 0.5 * dt * dt, 0},
            {0, 0, 1, 0, 0, dt, 0, 0, 0.5 * dt * dt},
            {0, 0, 0, 1, 0, 0, dt, 0, 0},
            {0, 0, 0, 0, 1, 0, 0, dt, 0},
            {0, 0, 0, 0, 0, 1, 0, 0, dt},
            {0, 0, 0, 0, 0, 0, 1, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 1, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 1}
        });

        if (print) {
            variances.open("variances.txt");
            if (!variances.is_open())
                std::cerr << "Error opening variances.txt" << std::endl;
        }
    }

    ~KalmanFilter() {
        if (print)
            variances.close();
    }

    void init(State _state) {
        state = _state;
        initialized = true;

        const double pos_var = 1e-6;
        const double speed = norm(state.vel());
        const double sigma_dir = std::sqrt(gyr_var);
        double vel_var = speed * speed * sigma_dir * sigma_dir;
        if (vel_var < 1e-6)
            vel_var = 1e-6;

        state.P = mat(kStateDim, kStateDim, 0.0);
        state.P.set_diag(vec({pos_var, pos_var, pos_var}), 0);
        state.P.set_diag(vec({vel_var, vel_var, vel_var}), 3);
        state.P.set_diag(vec({acc_var, acc_var, acc_var}), 6);

        gps_last_scale = 1.0;
    }

    void predict() {
        state.v = F * state.v;
        state.P = F * state.P * F.transposed() + Q;
    }

    void align_longitudinal(const vec &dir_unit) {
        const double speed = norm(state.vel());
        const vec v_aligned = dir_unit * speed;
        state.v.set(3, 0, mat(v_aligned).transposed());

        const vec a_aligned = project_on_dir(state.acc(), dir_unit);
        state.v.set(6, 0, mat(a_aligned).transposed());
    }

    void update(const mat &measurement, const mat &H, const mat &R, bool is_gps = false) {
        const mat Ht = H.transposed();
        mat y = measurement;
        y -= H * state.v;
        const double innov = norm(y);

        mat R_eff = R;
        double nis = 0.0;
        if (is_gps) {
            mat S_nom = H * state.P * Ht + R;
            mat S_nom_inv = S_nom.inv();
            nis = (y.transposed() * S_nom_inv * y)[0][0];
            const double dof = 3.0;
            const double use_tol = thresh > 0.0;
            const double tol = use_tol * thresh + (1.0 - use_tol);
            const double tol_weight = use_tol * (5.0 / tol) + (1.0 - use_tol);
            const double res_scale2 = use_tol * (innov / tol) * (innov / tol);
            const double nis_scale = (nis / dof) * tol_weight;
            ++gps_updates;
            const double scale = std::max(1.0, std::max(res_scale2, nis_scale));
            R_eff *= std::max(scale, nis_scale * nis_scale * (gps_updates == 1));
#ifdef BONUS
            gps_last_scale = scale; 
            if (scale > 1.0)
                ++gps_scaled;
            nis_sum += nis;
            nis_max = std::max(nis_max, nis);
            gps_res_sum += innov;
            gps_res_max = std::max(gps_res_max, innov);
#endif
        }

        mat S = H * state.P * Ht + R_eff;
        mat K = state.P * Ht * S.inv();
        state.v = state.v + K * y;
        const mat I = mat(kStateDim);
        mat A = I;
        A -= K * H;
        state.P = A * state.P * A.transposed() + K * R_eff * K.transposed();
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    std::string pos_to_string() {
        const mat& v = state.v;
        if (print) {
            const mat& P = state.P;
            std::cout << std::fixed << "x : " << v[0][0] << ", y : " << v[1][0] << ", z : " << v[2][0] << std::endl;
            variances << std::fixed << "x : " << P[0][0] << ", y : " << P[1][1] << ", z : " << P[2][2] << std::endl;
        }

        return std::to_string(v[0][0]) + " " + std::to_string(v[1][0]) + " " + std::to_string(v[2][0]);
    }


    std::vector<std::string> split_set(std::string s, std::string delimiter) {
        std::vector<std::string> ret;
        size_t pos_start = 0, pos_end = 0;
        std::string token;

        while ((pos_start = s.find_first_not_of(delimiter, pos_end)) != std::string::npos) {
            pos_start = s.find_first_not_of(delimiter, pos_end);
            pos_end = s.find_first_of(delimiter, pos_start);
            ret.push_back(s.substr(pos_start, pos_end - pos_start));
        }
        return ret;
    }

    vec get_val(const std::string &msg, const std::string &name) {
        vec ret;
        std::string tmp = msg.substr(msg.find(name) + name.length());

        tmp = tmp.substr(tmp.find_first_of("0123456789.-\n "), tmp.find_first_not_of("0123456789.-\n "));

        std::vector<std::string> vals = split_set(tmp, "\n ");

        for (size_t i = 0; i < vals.size(); i++)
            ret.push_back(std::stod(vals[i]));
        return ret;
    }

    void parse(std::string msg) {
        while (msg.find("[") != std::string::npos)
            msg = msg.substr(0, msg.find("[")) + msg.substr(msg.find("]") + 1);

        const bool has_true_pos = msg.find("TRUE POSITION") != std::string::npos;
        if (has_true_pos && !initialized) {
            const vec acc = get_val(msg, "ACCELERATION");
            const vec dir = get_val(msg, "DIRECTION");
            init(State(get_val(msg, "TRUE POSITION"), get_val(msg, "SPEED")[0], acc, dir));
            return;
        }

        predict();

        const vec dir = get_val(msg, "DIRECTION");
        const vec dir_unit = forward_from_euler(dir).normalize();
        const vec acc = get_val(msg, "ACCELERATION");
        const vec acc_aligned = project_on_dir(acc, dir_unit);
        update(mat(acc_aligned).transposed(), H_acc, R_acc);

        if (msg.find("POSITION") != std::string::npos || has_true_pos) {
            const vec gps = msg.find("POSITION") != std::string::npos
                ? get_val(msg, "POSITION")
                : get_val(msg, "TRUE POSITION");
            update(mat(gps).transposed(), H_pos, R_pos, true);
        }
        align_longitudinal(dir_unit);
    }

    void log_gps_stats() const {
        if (gps_updates == 0) {
            std::cerr << "GPS updates: 0" << std::endl;
            return;
        }
        std::cerr << "GPS updates: " << gps_updates
                  << ", scaled: " << gps_scaled
                  << ", gated: " << gps_gated
                  << ", last_scale: " << gps_last_scale
                  << ", mean NIS: " << (nis_sum / gps_updates)
                  << ", max NIS: " << nis_max
                  << ", mean |gps-res|: " << (gps_res_sum / gps_updates)
                  << ", max |gps-res|: " << gps_res_max
                  << std::endl;
    }

};

#endif
