#include <complex>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>

#define EVEC Eigen::Vector<complex<long double>, Eigen::Dynamic>
#define CLD complex<long double>

using namespace std;

const int limit = 1000; 
const int time_steps = 1000;
const int x_steps = 300;
const CLD x_len = 1.0;
const CLD t_len = 0.004;
const CLD m = 1;
const CLD hbar = 1;
const CLD k = 500;
const CLD dt = t_len / (CLD)time_steps;
const CLD dx = x_len / (CLD)x_steps;

const CLD sigma = 0.001;
const CLD mean_x = 0.5;

const CLD l = (CLD)1i * (dt*hbar) / ((CLD)2*m*dx*dx);
const CLD lamb02 = (CLD)1i * (dt*hbar) / ((CLD)4*m*dx*dx);
const CLD pot_fac = (CLD)(1i)/hbar;

const CLD wall_width = 0.03;
const CLD wall_start = 0.7;

EVEC pot;

CLD potential(CLD x) {
    if (real(x) >= real(wall_start) && real(x) <= real(wall_start + wall_width)) {
        return pot_fac*(CLD)70000;
    } else {
        return (CLD)0.0;
    }
}

void gauss_seidel_tri(EVEC &vec) {
    CLD a;
    EVEC b;
    EVEC new_vec;
    new_vec.resize(x_steps);
    new_vec.setZero();
    b.resize(x_steps);
    b.setZero();

    for (int rep = 1; rep <= limit; rep++) {
        for (int nth = 1; nth < x_steps-1; nth++) {
            b(nth) = (vec(nth-1) + vec(nth+1)) * lamb02 + (vec(nth) * (((CLD)1 - (CLD)2*l) + (dt/(CLD)2)*pot(nth)));
        }
    }

    for (int rep = 1; rep <= limit; rep++) {
        for (int nth = 1; nth < x_steps-1; nth++) {
            a = -lamb02 * (new_vec(nth-1) + vec(nth+1));
            new_vec(nth) = (b(nth) - a) / ( (CLD)1 + (CLD)2*l - (dt/(CLD)2)*pot(nth) );
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

//void integrate(double from, double to, EVEC func, :w


int main() {
    EVEC vec;
    vec.setZero();
    vec.resize(x_steps);
    for (int i = 1; i < x_steps-1; i++) {
        vec(i) = (exp((-pow(((CLD)i*dx - mean_x), 2)/sigma) + ((CLD)1i * k * (CLD)i*dx)));
    }

    cout << vec << endl;

    pot.resize(x_steps);
    pot.setZero();

    for (int nth = 1; nth < x_steps-1; nth++) {
        pot(nth) = potential((CLD)nth*dx);
    }

    ofstream dat;
    dat.open("nicolson-k-space.txt");
    int count = 0;

    for (int tim = 0; tim < time_steps; tim++) {
        cout << "Step: " << tim << " | ";
        if (count % 10 == 0) {
            for (int i = 0; i < x_steps; i++) {
                dat << real(vec.array().abs().square())(i) << "|";
            }
            dat << endl;
        }
        gauss_seidel_tri(vec);
        count++;
    }

    dat << ">";
    dat.close();

    return 0;
}
