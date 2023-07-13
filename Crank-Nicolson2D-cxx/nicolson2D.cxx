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
const int time_steps = 100;
const int x_steps = 50;
const CLD x_len = 1.0;
const CLD t_len = 0.004;
const CLD m = 1;
const CLD hbar = 1;
const CLD kx = 500;
const CLD ky = 500;
const CLD dt = t_len / (CLD)time_steps;
const CLD dx = x_len / (CLD)x_steps;

const CLD sigma = 0.001;
const CLD mean_x = 0.5;
const CLD mean_y = 0.5;

const CLD lambda = (CLD)1i * (dt*hbar) / ((CLD)2*m*dx*dx);
const CLD lamb02 = (CLD)1i * (dt*hbar) / ((CLD)4*m*dx*dx);
const CLD pot_fac = (CLD)(1i)/hbar;

const CLD wall_width = 0.03;
const CLD wall_start = 0.7;

CLD potential(CLD x, CLD y) {
    return 0;
}

EMAT gauss_seidel_tri(EMAT v, CLD l, int limit) {
    EMAT vec;
    EMAT a;
    EMAT b;
    EMAT new_vec;
    vec.resize(x_steps, x_steps);
    new_vec.resize(x_steps, x_steps);
    a.resize(x_steps, x_steps);
    b.resize(x_steps, x_steps);
    vec.setZero();
    a.setZero();
    b.setZero();
    for (int rep = 1; rep <= limit; rep++) {
        new_vec.setZero();
        for (int nth = 1; nth < x_steps-1; nth++) {
            for (int nth2 = 1; nth2 < x_steps-1; nth2++) {
                if (rep == 1) {
                    b(nth, nth2) = (v(nth-1, nth2) + v(nth+1, nth2) + v(nth, nth2-1) + v(nth, nth2+1)) * lamb02 + (v(nth, nth2) * (((CLD)1 - (CLD)2*l) + (dt/(CLD)2)*potential((CLD)nth*dx, (CLD)nth2*dx)));
                }
                a(nth, nth2) = lamb02 * (new_vec(nth-1, nth2) + vec(nth+1, nth2) + new_vec(nth, nth2-1) + vec(nth, nth2+1));
                new_vec(nth, nth2) = (b(nth, nth2) - a(nth, nth2)) / ( (CLD)1 + (CLD)2*l - (dt/(CLD)2)*potential((CLD)nth*dx, (CLD)nth2*dx) );
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
    return vec;
}

int main() {
    EMAT vec;
    vec.resize(x_steps, x_steps);
    for (int i = 0; i < x_steps; i++) {
        for (int j = 0; j < x_steps; j++) {
            vec(i, j) = (exp((-(pow(((CLD)i*dx - mean_x), 2) + pow((CLD)j*dx, 2))/sigma) + ((CLD)1i * (kx * (CLD)i*dx + ky * (CLD)j*dx))));
        }
    }

    //cout << vec << endl;
    //cout << vec.norm() << endl;

    //EMAT v = gauss_seidel_tri(vec, lambda, limit);
    //cout << v << endl;

    vector<EMAT> vs;
    for (int tim = 0; tim < time_steps; tim++) {
        vs.push_back(vec);
        cout << "Step: " << tim << " | ";
        EMAT v = gauss_seidel_tri(vec, lambda, limit);
        cout << tim << endl;
        vec = v;
    }

    //EMAT v_temp = vs[0].array().abs().square();

    //cout << v_temp << endl;

    ofstream dat;
    dat.open("data.txt");
    int count = 0;
    for (EMAT v2 : vs) {
        if (count % 5 == 0) {
            EMAT v_temp = v2.array().abs().square();
            cout << v_temp << endl;
            for (int i = 0; i < x_steps; i++) {
                for (int j = 0; j < x_steps; j++) {
                    dat << real(v_temp(i, j)) << "|";
                }
                dat << ";";
            }
            dat << endl;
        }
        count++;
    }
    dat << ">";
    dat.close();
    return 0;
}
