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
const int time_steps = 5000;
const int x_steps = 1000;
const CLD x_len = 1.0;
const CLD t_len = 0.004;
const CLD m = 1;
const CLD hbar = 1;
const CLD k = 500;
const CLD dt = t_len / (CLD)time_steps;
const CLD dx = x_len / (CLD)x_steps;

const CLD sigma = 0.001;
const CLD mean_x = 0.5;

const CLD lambda = (CLD)1i * (dt*hbar) / ((CLD)2*m*dx*dx);
const CLD lamb02 = (CLD)1i * (dt*hbar) / ((CLD)4*m*dx*dx);
const CLD pot_fac = (CLD)(1i)/hbar;

const CLD wall_width = 0.03;
const CLD wall_start = 0.7;

CLD potential(CLD x) {
    if (real(x) >= real(wall_start) && real(x) <= real(wall_start + wall_width)) {
        return -pot_fac*(CLD)70000;
    } else {
        return (CLD)0.0;
    }
}

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
                b(nth) = (v(nth-1) + v(nth+1)) * lamb02 + (v(nth) * (((CLD)1 - (CLD)2*l) + (dt/(CLD)2)*potential((CLD)nth*dx)));
            }
            a(nth) = lamb02 * (new_vec(nth-1) + vec(nth+1));
            new_vec(nth) = (a(nth) - b(nth)) / ( (CLD)1 + (CLD)2*l - (dt/(CLD)2)*potential((CLD)nth*dx) );
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
        cout << "Step: " << tim << " | ";
        EVEC v = gauss_seidel_tri(vec, lambda, limit);
        vec = v;
    }

    ofstream dat;
    dat.open("nicolson.txt");
    int count = 0;
    for (EVEC v2 : vs) {
        if (count % 5 == 0) {
            EVEC v_temp = v2.array().abs().square();
            //cout << v_temp << endl;
            for (int i = 0; i < x_steps; i++) {
                dat << real(v_temp(i)) << "|";
            }
            dat << endl;
        }
        count++;
    }
    dat << ">";
    dat.close();
    return 0;
}
