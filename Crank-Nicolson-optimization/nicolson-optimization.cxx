#include <complex>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <ctime>

#define EMAT Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>
#define CLD complex<double>

using namespace std;

int limit = 1000; 
int time_steps = 50;
int x_steps = 400;
const CLD x_len = 1.0;
const CLD t_len = 0.0016;
const CLD kx = 0;
const CLD ky = 500;
CLD dt = t_len / (CLD)time_steps;
const CLD dx = x_len / (CLD)x_steps;

long iter_amount = 0;
//const int x_steps_steps = 50;
const int fraction = 20;
const int time_steps_steps = 50;
const int time_steps_max = 2000;

const CLD sigma = 0.0025;
const CLD mean_x = 0.5;
const CLD mean_y = 0.1;

CLD l = (CLD)1i * dt / ((CLD)4*dx*dx);
CLD pot_fac = (CLD)(1i)*(dt/(CLD)2.0);

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
            iter_amount += rep;
            break;
        }

        if (rep == limit) {
            cout << "limit reached" << endl;
            exit(0);
        }
        vec = new_vec;
    }
}

int main() {
    srand(time(0));

    pot.resize(x_steps, x_steps);
    pot.setZero();

    EMAT vec;
    EMAT vec_start;

    vec.resize(x_steps, x_steps);
    vec.setZero();
    vec_start.resize(x_steps, x_steps);
    vec_start.setZero();

    for (int i = 1; i < x_steps-1; i++) {
        for (int j = 1; j < x_steps-1; j++) {
            vec(i, j) = exp((-(pow(((CLD)i*dx - mean_x), 2) + pow((CLD)j*dx - mean_y, 2))/sigma) + ((CLD)1i * (kx * (CLD)i*dx + ky * (CLD)j*dx)));
        }
    }

    vec_start = vec;

    ofstream dat;
    dat.open("gauss-seidel-comp.txt");
    
    cout << "X-Steps: " << x_steps <<  " | X-len: " << real(x_len) << " | T-len: " << real(t_len) << endl;

    for (time_steps; time_steps <= time_steps_max; time_steps += time_steps_steps) {
        dt = t_len / (CLD)time_steps;
        l = (CLD)1i * dt / ((CLD)4*dx*dx);
        pot_fac = (CLD)(1i)*(dt/(CLD)2.0);
        for (int tim = 0; tim < (time_steps / fraction); tim++) {
            gauss_seidel_tri(vec);
        }
        long out = iter_amount * fraction * (x_steps * x_steps);
        cout << "T-Steps: " << time_steps << " | Iter-Amount: " <<  out << endl;
        dat << out;
        dat << "|";
        vec = vec_start;
        iter_amount = 0;
    }

    dat << ">";
    dat.close();

    return 0;
}
