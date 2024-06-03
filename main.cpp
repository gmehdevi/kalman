#include <iostream>
#include <string>
#include "boost/asio.hpp"
#include "boost/array.hpp"
#include <sstream>
#include <iomanip>
#include <chrono>
#include "kalman.hpp"
#include <fstream>

#define PORT 4242
#define HANDSHAKE "READY"
#define TIMEOUT_SECONDS 5

using namespace std;
using namespace boost;
using namespace std::chrono;


tuple<double, double, bool, double, double, double, double, double> parse_arguments(int ac, const char *av[]) {
    double process_noise = 0.01;
    double n = 1;
    bool print = false;
    double accNoise = 10E-6;
    double gyrNoise = 10E-4;
    double gpsNoise = 1800;
    double thresh = 100;
    double alpha = 0.001;
    try {
        for (int i = 1; i < ac; i++) {
            if (strcmp(av[i], "--process_noise") == 0) {
                process_noise = stod(av[++i]);
            } else if (strcmp(av[i], "-n") == 0) {
                n = stod(av[++i]);
            } else if (strcmp(av[i], "-p") == 0) {
                print = true;
            } else if (strcmp(av[i], "--acc") == 0) {
                accNoise = stod(av[++i]);
            } else if (strcmp(av[i], "--gyr") == 0) {
                gyrNoise = stod(av[++i]);
            } else if (strcmp(av[i], "--gps") == 0) {
                gpsNoise = stod(av[++i]);
            } else if (strcmp(av[i], "-t") == 0) {
                thresh = stod(av[++i]);
            } else if (strcmp(av[i], "-a") == 0) {
                alpha = stod(av[++i]);
            } else {
                cerr << "Invalid flag: " << av[i] << endl;
                exit(1);
            }
        }
    } catch (std::exception &e) {
        cerr << "Invalid argument: " << e.what() << endl;
        exit(1);
    }
    return make_tuple(process_noise, n, print, accNoise, gyrNoise, gpsNoise, thresh, alpha);
}

int main(int ac, const char *av[]) {

    auto [process_noise, n, print, accNoise, gyrNoise, gpsNoise, thresh, alpha] = parse_arguments(ac, av);

    asio::io_service io_service;
    asio::ip::udp::socket socket(io_service, asio::ip::udp::endpoint(asio::ip::udp::v4(), 0));
    asio::ip::udp::endpoint remote_endpoint = asio::ip::udp::endpoint(asio::ip::address::from_string("127.0.0.1"), PORT);
    asio::ip::udp::endpoint local_endpoint;
    stringstream msg;
    boost::array<char, 4096> recv_buf;
    system::error_code error;
    duration<double> timeout(TIMEOUT_SECONDS);

    KalmanFilter kalman(process_noise, n, print, accNoise, gyrNoise, gpsNoise, thresh, alpha);

    std::cout << "Sending handshake..." << std::endl;

    auto startTime = high_resolution_clock::now();

    while (socket.available() == 0) {
        socket.send_to(asio::buffer(HANDSHAKE, strlen(HANDSHAKE)), remote_endpoint, 0, error);
        if (error)
            std::cerr << "Error sending handshake: " << error.message() << std::endl;

        auto currentTime = high_resolution_clock::now();
        if (currentTime - startTime >= timeout) {
            std::cout << "Handshake timeout reached." << std::endl;
            socket.close();
            return 1;
        }
        sleep(1);
    }

    ofstream file;
    file.open("log.txt", ofstream::out | ofstream::trunc);
    startTime = high_resolution_clock::now();
    size_t msgs = 0;

    while (true) {

        if (socket.available() > 0) {
            size_t len = socket.receive_from(asio::buffer(recv_buf), remote_endpoint, 0, error);
            if (error && error != asio::error::message_size)
                throw system::system_error(error);
            recv_buf[len] = '\0';
            msg << recv_buf.data();
            recv_buf.data()[0] = '\0';
        }
        if (msg.str().find("MSG_END") != string::npos) {
            startTime = high_resolution_clock::now();

            kalman.parse(msg.str());

            kalman.predict();
            socket.send_to(asio::buffer(kalman.pos_to_string()), remote_endpoint, 0, error);
            if (error)
                throw system::system_error(error);


            file << "\nMessage: " << msgs << "\n" << msg.str();
            msgs++;
            msg.str("");
        }

        auto currentTime = high_resolution_clock::now();
        if (currentTime - startTime >= timeout / 2) {
            std::cout << "Timeout reached. Exiting." << std::endl;
            break;
        }
    }
    while (socket.available() > 0)
       socket.receive_from(asio::buffer(recv_buf), remote_endpoint, 0, error), recv_buf.data()[0] = '\0';

    file << "\nReceived " << msgs << " messages." << std::endl;

    #ifdef BONUS
        std::cerr << "GPS variance: " << kalman.R_pos[0][0] << std::endl;
        std::cerr << "Acc variance: " << kalman.R_acc[0][0] << std::endl;
        std::cerr << "Gyr variance: " << kalman.R_dir[0][0] << std::endl;
    #endif
    
    file.close();
    socket.close();
    return 0;
}
