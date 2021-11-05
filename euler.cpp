#include<iostream>
#include<fstream>

int main() {
    const int nx = 250; // x steps
    const double h = 1./nx;
    double dt = 0.000004; // time steps in 2tau
    const int nt = 1250000; // time steps in 5tau
    const int nt2 = 2 * nt / 5;
    auto *u = new double[nx + 1]; // target vector - function u(x)
    auto *tmp = new double[nx + 1]; // for matrix multiplication
    // initial conditions
    for (int i = 0; i <= nx ; i++){
        u[i] = 0; tmp[i] = 0;
    }
    for (int i = 17 * nx / 30; i <= 23 * nx / 30; ++i) {
        u[i] = 1; tmp[i] = 1;
    }
    // consts for matrix multiplication
    const double b0 = -2 / h / (h + 2);
    const double a0 = -(2 + 3 * h) / (h + 2);
    const double K =(double) dt / h / h;
    std::ofstream fx("fx.txt");
    const int x = nx * 2 / 3; // calculation in x=2l/3
    std::ofstream f2("f2.txt");
    // t = 0...2tau
    for (int i = 0; i <= nt2; ++i) {
        for (int k = 0; k <= nx; ++k) {
            tmp[k] = u[k];    // tmp = u
        }
        u[0] += (a0 * tmp[0] + tmp[1]) * K + b0 * dt;
        for (int j = 1; j < nx; ++j) { // multiplying vector u(x) by sparse matrix
            u[j] += (tmp[j - 1] - 2 * tmp[j] + tmp[j + 1]) * K;  // multiplying vector u(x) by tridiagonal matrix
        }
        u[nx] += (tmp[nx - 1] - tmp[nx]) * K;
        fx << i * dt << ' ' << u[x] << '\n';  // x=2L/3
    }
    for (int i = 0; i <= nx; ++i)
    {
        f2 << i * h << ' ' << u[i] << '\n';   // result in t = 2tau
    }
    std::ofstream f5("f5.txt");
    // t = 2tau...5tau
    for (int i = nt2 + 1; i <= nt; ++i)
    {
        //tmp = u;
        for (int k = 0; k <= nx; k++) {
            tmp[k] = u[k];
        }
        u[0] += (a0 * tmp[0] + tmp[1]) * K + b0 * dt;
        for (int j = 1; j < nx; j++) {
            u[j] += (tmp[j - 1] - 2 * tmp[j] + tmp[j + 1]) * K; // multiplying vector u(x) by tridiagonal matrix
        }
        u[nx] += (tmp[nx - 1] - tmp[nx]) * K;
        fx << i * dt << ' ' << u[x] << '\n';  // x=2L/3
    }
    for (int i = 0; i <= nx; ++i) {
        f5 << i * h <<' ' << u[i] << '\n';  // result in t = 5tau
    }

    f2.close();
    f5.close();
    fx.close();
    delete[]u;
    delete[]tmp;
    return 0;
}