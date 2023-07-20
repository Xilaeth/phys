#include <complex>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <ctime>

#define EMAT Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>
#define EVEC Eigen::Vector<complex<double>, Eigen::Dynamic>
#define CLD complex<double>

using namespace std;

const int limit = 1000; 
const int time_steps = 1000;
const int x_steps = 100;
const CLD x_len = 1.0;
const CLD t_len = 0.0025;
const CLD kx = 0;
const CLD ky = 500;
const CLD dt = t_len / (CLD)time_steps;
const CLD dx = x_len / (CLD)x_steps;

const int k_steps = 100;
const CLD k_from = -750;
const CLD k_to = 750;
const CLD k_len = k_to - k_from;
const CLD dk = k_len / (CLD)k_steps;

const CLD sigma = 0.0025;
const CLD mean_x = 0.5;
const CLD mean_y = 0.1;

const CLD l = (CLD)1i * dt / ((CLD)4*dx*dx);
const CLD pot_fac = (CLD)(1i)*(dt/(CLD)2.0);

const double wall_y = 0.45;
const double wall_x1 = 0.470;
const double wall_x2 = 0.478;
const double wall_x3 = 0.522;
const double wall_x4 = 0.530;
const double wall_width = 0.02;

EMAT pot;

void gauss_seidel_tri(EMAT &vec) {
    CLD a;
    EMAT b;
    EMAT new_vec;
    vec.resize(x_steps, x_steps);
    new_vec.resize(x_steps, x_steps);
    b.resize(x_steps, x_steps);
    b.setZero();

    for (int ix = 1; ix < x_steps-1; ix++) {
        for (int iy = 1; iy < x_steps-1; iy++) {
            b(ix, iy) = (vec(ix-1, iy) + vec(ix+1, iy) + vec(ix, iy-1) + vec(ix, iy+1)) * l + vec(ix, iy) * ((CLD)1 - (CLD)4.0*l + pot(ix, iy));
       }
    }

    for (int rep = 1; rep <= limit; rep++) {
        new_vec.setZero();
        for (int ix = 1; ix < x_steps-1; ix++) {
            for (int iy = 1; iy < x_steps-1; iy++) {
                if (real((CLD)iy*dx) > wall_y && real((CLD)iy*dx) < (wall_y + wall_width)) {
                    if (real((CLD)ix*dx) < wall_x1 || real((CLD)ix*dx) > wall_x4 || (real((CLD)ix*dx) > wall_x2 && real((CLD)ix*dx) < wall_x3)) {
                        continue;
                    }
                } 
                a = -l * (new_vec(ix-1, iy) + vec(ix+1, iy) + new_vec(ix, iy-1) + vec(ix, iy+1));
                new_vec(ix, iy) = (b(ix, iy) - a) / ( (CLD)1 + (CLD)4*l - pot(ix, iy));
            }
        }

        if ((new_vec - vec).norm() < 1e-10) {
            cout << rep << " | Error: " << (new_vec - vec).norm() << endl;
            break;
        }

        if (rep == limit) {
            cout << "limit reached" << endl;
            exit(0);
        }
        vec = new_vec;
    }
}

void fourier_trafo(EMAT &k_vec, EMAT &vec) {
    CLD pre_fac = (dx*dx)/(CLD)(sqrt(2.0*M_PI));
    for (int ikx = 0; ikx < k_steps; ikx++) {
        for (int iky = 0; iky < k_steps; iky++) {
            for (int ix = 1; ix < x_steps-1; ix++) {
                for (int iy = 1; iy < x_steps-1; iy++) {
                    k_vec(ikx, iky) += pre_fac * vec(ix, iy) * exp( -(CLD)1i * ((((CLD)(k_from + (CLD)ikx*dk)) * (dx*(CLD)ix)) + (((CLD)(k_from + (CLD)iky*dk)) * (dx*(CLD)iy))));
                }
            }
        }
    }
}

int main() {
    EMAT vec;

    pot.resize(x_steps, x_steps);
    pot.setZero();
    vec.resize(x_steps, x_steps);
    vec.setZero();

    for (int i = 1; i < x_steps-1; i++) {
        for (int j = 1; j < x_steps-1; j++) {
            vec(i, j) = exp((-(pow(((CLD)i*dx - mean_x), 2) + pow((CLD)j*dx - mean_y, 2))/sigma) + ((CLD)1i * (kx * (CLD)i*dx + ky * (CLD)j*dx)));
        }
    }

    

    ofstream dat;
    ofstream dat2;
    dat.open("nicolson2d-k-space.txt");
    EMAT v0;
    v0 = vec;

    EMAT k_vec;
    k_vec.resize(k_steps, k_steps);
    k_vec.setZero();

    cout << k_vec << endl;

    for (int tim = 0; tim < time_steps; tim++) {
        if (tim % 100 == 0) {
            fourier_trafo(k_vec, vec);
            for (int i = 0; i < x_steps; i++) {
                for (int j = 0; j < x_steps; j++) {
                    dat << real( (k_vec.array().abs().square())(i, j) ) << "|";
                }
                dat << ";";
            }
            dat << endl;
            k_vec.setZero();
        }
        cout << "Step: " << tim << " | ";
        gauss_seidel_tri(vec);
        cout << tim << " | " << (vec.array().abs().square()).sum() << endl;
    }

    dat << ">";
    dat.close();

    return 0;
}