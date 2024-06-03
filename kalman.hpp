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

typedef ft::vector<double> vec;
typedef ft::matrix<double> mat;

double huberLoss(const mat& residual, double delta) {
    double sum = 0.0;
    for (size_t i = 0; i < residual.rows(); ++i) {
        for (size_t j = 0; j < residual.cols(); ++j) {
            double r = residual[i][j];
            if (abs(r) <= delta) {
                sum += 0.5 * r * r;
            } else {
                sum += delta * (abs(r) - 0.5 * delta);
            }
        }
    }
    return sum;
}

struct State {

    mat v;
    mat P;

    State() : v(12, 1), P(12) {}

    State(const vec &pos, double speed, const vec &acc, const vec &dir) : v(12, 1), P(12) {
        vec vel = direction(dir).normalize() * (speed  * 1000 / 3600);//km/h en m/s
        vec tmp({pos[0], pos[1], pos[2],
                 vel[0], vel[1], vel[2],
                 acc[0], acc[1], acc[2],
                 dir[0], dir[1], dir[2]});
        v = mat(tmp).transposed();
    }

    vec pos() {return vec({v[0][0], v[1][0], v[2][0]});}

    vec vel() {return vec({v[3][0], v[4][0], v[5][0]});}

    vec acc() {return vec({v[6][0], v[7][0], v[8][0]});}

    vec dir() {return vec({v[9][0], v[10][0], v[11][0]});}
};


class KalmanFilter {
public:
    mat F; // state transition matrix

    mat Q; // process noise covariance
    mat B; // control matrix

    mat H_pos; // measurement matrix
    mat H_acc;
    mat H_dir;

    mat R_pos; // measurement noise covariance
    mat R_acc;
    mat R_dir;

    State state;

    double thresh;
    double alpha;

    bool print;
    std::ofstream variances;


    KalmanFilter(double process_noise, double n = 1, bool print = false,
                    double accNoise = 10E-6, double gyrNoise = 10E-4, double gpsNoise = 1800,
                    double thresh = 100, double alpha = 0.001
    ) : F(12), Q(12), B(12, 12),
    H_pos(3, 12, 0), H_acc(3, 12, 0), H_dir(3, 12, 0), R_pos(3), R_acc(3), R_dir(3),
    thresh(thresh), alpha(alpha), print(print) {
        H_pos = mat({
            {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0},
            {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0},
            {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ,0}
        });

        H_acc = mat({
            {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 ,0},
            {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ,0},
            {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ,0}
        });

        H_dir = mat({
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ,0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ,0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,1}
        });

        R_pos *= gpsNoise * n * n * n * n * n * n * n;

        R_acc *= accNoise * n;

        R_dir *= gyrNoise * n;

        Q *= process_noise* n;

        double dt = 0.01;

        F = mat({
            {1,  0,  0,  dt, 0,   0,  0.5 * dt * dt,   0, 0, 0, 0, 0 },
            {0,  1,  0,  0,  dt,  0,  0,  0.5 * dt * dt,  0, 0, 0, 0 },
            {0,  0,  1,  0,  0,  dt,  0,  0,  0.5 * dt * dt, 0, 0, 0 },
            {0,  0,  0,  1,  0,   0,  dt, 0,   0,            0, 0, 0 },
            {0,  0,  0,  0,  1,   0,  0,  dt,  0,            0, 0, 0 },
            {0,  0,  0,  0,  0,   1,  0,  0,  dt,            0, 0, 0 },
            {0,  0,  0,  0,  0,   0,  1,  0,   0,            0, 0, 0 },
            {0,  0,  0,  0,  0,   0,  0,  1,   0,            0, 0, 0 },
            {0,  0,  0,  0,  0,   0,  0,  0,   1,            0, 0, 0 },
            {0,  0,  0,  0,  0,   0,  0,  0,   0,            1, 0, 0 },
            {0,  0,  0,  0,  0,   0,  0,  0,   0,            0, 1, 0 },
            {0,  0,  0,  0,  0,   0,  0,  0,   0,            0, 0, 1 }
        });

        B[4][0] = dt; B[5][1] = dt; B[6][2] = dt;

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

        state.P = mat({
            {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0},
            {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0},
            {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ,0},

            {0, 0, 0, R_dir[0][0], 0, 0, 0, 0, 0, 0, 0 ,0},
            {0, 0, 0, 0, R_dir[0][0], 0, 0, 0, 0, 0, 0 ,0},
            {0, 0, 0, 0, 0, R_dir[0][0], 0, 0, 0, 0, 0 ,0},

            {0, 0, 0, 0, 0, 0, R_acc[0][0], 0, 0, 0, 0 ,0},
            {0, 0, 0, 0, 0, 0, 0, R_acc[1][1], 0, 0, 0 ,0},
            {0, 0, 0, 0, 0, 0, 0, 0, R_acc[2][2], 0, 0 ,0},

            {0, 0, 0, 0, 0, 0, 0, 0, 0, R_dir[0][0], 0 ,0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, R_dir[1][1], 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, R_dir[2][2]}
         });
    }

    void predict() {
        vec tmp(direction(state.dir()));

        mat u = mat({0, 0, 0 , tmp[0], tmp[1], tmp[2], 0, 0, 0, 0, 0, 0}).transposed();

        state.v = F * state.v + B * u;
        state.P = F * state.P * F.transposed() + Q;
    }


    void update(const mat &measurement, const mat &H, mat &R) {
        mat S = H * state.P * H.transposed() + R;
        mat K = state.P * H.transposed() * S.inv();
        mat y = measurement.transposed() - H * state.v;


        double innov = 0.0;
        for (size_t i = 0; i < y.rows(); ++i)
            innov += y[i][0] * y[i][0];
        innov = sqrt(innov);
#ifdef BONUS
        double thresh = 3;
        if (innov < thresh) {
            R = R *  std::max(10E-2, R_pos[0][0] * (1.0 - (alpha / 1000) * (thresh - innov))) / R_pos[0][0];
        } else {
            R = R *  std::min(5000.0, R_pos[0][0] + alpha * 1000  * (innov - thresh)) / R_pos[0][0];
        }

        if (std::min(sqrt((y.transposed() * S.inv() * y)[0][0]), huberLoss(y, thresh)) > thresh * 10) {
            std::cerr << "Outlier detected, skipping update" << std::endl;
            return;
        }
#endif
        state.v = state.v + K * y;
        state.P = (mat(12) - K * H) * state.P;
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

        if (msg.find("TRUE POSITION") != std::string::npos) {
            init(State(get_val(msg, "TRUE POSITION"), get_val(msg, "SPEED")[0], get_val(msg, "ACCELERATION"), get_val(msg, "DIRECTION")));
        } else if (msg.find("POSITION")!= std::string::npos) {
            update(get_val(msg  , "POSITION"), H_pos, R_pos);
            update(get_val(msg, "ACCELERATION"), H_acc, R_acc);
            update(get_val(msg, "DIRECTION"), H_dir, R_dir);
        } else {
            update(get_val(msg, "ACCELERATION"), H_acc, R_acc);
            update(get_val(msg, "DIRECTION"), H_dir, R_dir);
        }
    }

};

#endif
