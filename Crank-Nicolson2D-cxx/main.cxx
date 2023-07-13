#include <numeric>
#include <complex>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>

#define EMAT Eigen::Matrix<complex<long double>, Eigen::Dynamic, Eigen::Dynamic>
#define EVEC Eigen::Vector<complex<long double>, Eigen::Dynamic>
#define CLD complex<long double>

using namespace std;

const int limit = 1000; 
const int time_steps = 1000;
const int x_steps = 400;
const CLD x_len = 1.0;
const CLD t_len = 0.004;
const CLD m = 1;
const CLD hbar = 1;
const CLD k = 500;
const CLD dt = t_len / (CLD)time_steps;
const CLD dx = x_len / (CLD)x_steps;

const CLD sigma = 0.001;
const CLD mean_x = 0.5;

const CLD lambda = (CLD)1i * (dt*hbar) / ((CLD)2*m*pow(dx, 2));
const CLD pot_fac = (CLD)(1i)/hbar;

EVEC gauss_seidel_tri(EVEC v, CLD l, int limit) {
    EVEC vec;
    EVEC a;
    EVEC b;
    EVEC new_vec;
    vec.resize(x_steps);
    new_vec.resize(x_steps);
    a.resize(x_steps);
    b.resize(x_steps);
    vec.setZero();
    a.setZero();
    b.setZero();
    for (int rep = 1; rep <= limit; rep++) {
        new_vec.setZero();
        for (int nth = 1; nth < x_steps-1; nth++) {
            if (rep == 1) {
                b(nth) = (v(nth-1) + v(nth+1)) * (l/(CLD)2.0) + (v(nth) * ((CLD)1 - (CLD)2*l));
            }
            a(nth) = (l/(CLD)2) * (new_vec(nth-1) + vec(nth+1));
        }
        new_vec = (a - b) / ( (CLD)1 + (CLD)2*l );
        if ((new_vec - vec).norm() < pow(10, -10)) {
            cout << rep << " | Error: " << (new_vec - vec).norm() << endl;
            break;
        }
        if (rep == limit) {
            cout << "limit reached" << endl;
            exit(0);
        }
        vec = new_vec;
    }
    return vec;
}

int main() {
    EVEC vec;
    vec.resize(x_steps);
    for (int i = 0; i < x_steps; i++) {
        vec(i) = (exp((-pow(((CLD)i*dx - mean_x), 2)/sigma) + ((CLD)1i * k * (CLD)i*dx)));
    }

    vector<EVEC> vs;
    for (int tim = 0; tim < time_steps; tim++) {
        vs.push_back(vec);
        EVEC v = gauss_seidel_tri(vec, lambda, limit);
        vec = v;
    }

    cout << vs[time_steps-1] << endl;

    ofstream dat;
    dat.open("data.txt");
    for (EVEC v2 : vs) {
        for (int i = 0; i < x_steps; i++) {
            dat << pow(abs(v2(i)), 2) << "|";
        }
        dat << endl;
    }
    dat << ">";
    dat.close();
    return 0;
}
