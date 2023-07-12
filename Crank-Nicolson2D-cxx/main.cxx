#include <numeric>
#include <complex>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace std;

const int limit = 1000; 
const int time_steps = 10000;
const int x_steps = 400;
const complex<long double> x_len = 1.0;
const complex<long double> t_len = 0.004;
const complex<long double> m = 1;
const complex<long double> hbar = 1;
const complex<long double> k = 500;
const complex<long double> dt = t_len / (complex<long double>)time_steps;
const complex<long double> dx = x_len / (complex<long double>)x_steps;

const complex<long double> sigma = 0.001;
const complex<long double> mean_x = 0.5;

const complex<long double> lambda = (complex<long double>)1i * (dt*hbar) / ((complex<long double>)2*m*pow(dx, 2));
const complex<long double> pot_fac = (complex<long double>)(1i)/hbar;

int mat_to_ind(int i, int j) {
    return (int)(i * x_steps + j);
}

vector<complex<long double>> add_v(vector<complex<long double>> v1, vector<complex<long double>> v2) {
    vector<complex<long double>> v_out(0, v1.size());
    for (int i = 0; i < v1.size(); i++) {
        v_out.push_back(v1[i] + v2[i]);
    }
    return v_out;
}

vector<complex<long double>> sub_v(vector<complex<long double>> v1, vector<complex<long double>> v2) {
    vector<complex<long double>> v_out(0, v1.size());
    for (int i = 0; i < v1.size(); i++) {
        v_out.push_back(v1[i] - v2[i]);
    }
    return v_out;
}

vector<complex<long double>> div_v_s(vector<complex<long double>> v, complex<long double> scal) {
    vector<complex<long double>> v_out(0, v.size());
    for (int i = 0; i < v.size(); i++) {
        v_out.push_back(v[i] / scal);
    }
    return v_out;
}

long double norm_v(vector<complex<long double>> v1) {
    long double out = 0;
    for (int i = 0; i < v1.size(); i++) {
        out += pow(abs(v1[i]),2);
    }
    return sqrt(out);
}

vector<complex<long double>> gauss_seidel_tri(vector<complex<long double>> v, complex<long double> l, int limit) {
    vector<complex<long double>> vec(x_steps, 0);
    vector<complex<long double>> b(x_steps, 0);
    for (int rep = 1; rep <= limit; rep++) {
        vector<complex<long double>> new_vec(x_steps, 0);
        for (int nth = 1; nth < x_steps-1; nth++) {
            if (rep == 1) {
                b[nth] = (v[nth-1] + v[nth+1]) * (l/(complex<long double>)2.0) + (v[nth] * ((complex<long double>)1 - (complex<long double>)2*l));
            }
            complex<long double> a = (l/(complex<long double>)2) * (new_vec[nth-1] + vec[nth+1]);
            new_vec[nth] = (b[nth] - a) / ((complex<long double>)1 + (complex<long double>)2*l);
        }
        if (norm_v(sub_v(new_vec, vec)) < pow(10, -5)) {
            //cout << rep << " | Error: " << norm_v(sub_v(new_vec, vec)) << endl;
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
    vector<complex<long double>> vec;
    for (int i = 0; i < x_steps; i++) {
        vec.push_back(exp((-pow(((complex<long double>)i*dx - mean_x), 2)/sigma) + ((complex<long double>)1i * k * (complex<long double>)i*dx)));
    }
    vector<vector<complex<long double>>> vs;
    for (int tim = 0; tim < time_steps; tim++) {
        vs.push_back(vec);
        vector<complex<long double>> v = gauss_seidel_tri(vec, lambda, limit);
        vec = v;
    }

    ofstream dat;
    dat.open("data.txt");
    for (auto v2 : vs) {
        for (int i = 0; i < x_steps; i++) {
            dat << pow(abs(v2[i]), 2) << "|";
        }
        dat << endl;
    }
    dat << ">";
    dat.close();
    return 0;
}
