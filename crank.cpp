#include<iostream>
#include<fstream>

int main() {
    const int nx = 250; // x steps
    const double h = 1./nx;
    double dt = 0.000004; // time steps in 2tau
    const int nt = 1250000; // time steps in 5tau
    const int nt2 = 2 * nt / 5;
    auto *u = new double[nx + 1]; // target vector - function u(x)
    //vectors for tridiagonal matrix algorithm
    auto *g = new double[nx + 1];
    auto *q = new double[nx];
    auto *p = new double[nx];
    //initial conditions
    for (int i = 0; i <= nx; i++) {
        u[i] = 0;
    }
    for (int i = 17 * nx / 30; i <= 23 * nx / 30; i++) {
        u[i] = 1;
    }
    const double a0 = -(2 + 3 * h) / (h + 2);
    const double b0 = -2 / h / (h + 2);
    const double K = (double) dt / h / h;
    const int x = nx * 2 / 3; // calculation in x=2l/3
    // consts and calculations for tridiagonal matrix algorithm
    const double a1 = 1 - a0 * K / 2;
    const double an = 1 + K / 2;
    const double ai = 1 + K;
    const double bi = -K / 2;
    q[0] = -bi / a1;
    for (int i = 1; i < nx; i++) {
        q[i] = -bi / (ai + q[i - 1] * bi);
    }
    std::ofstream fx("fx.txt");
    std::ofstream f2("f2.txt");
    // t = 2tau
    for (int i = 0; i <= nt2; i++) {
        //g=u(t)+A*u(t)*dt/2+B
        g[0] = u[0] + (a0 * u[0] + u[1]) * K / 2 + b0 * dt;
        for (int k = 1; k < nx; k++) {
            g[k] = u[k] + (u[k - 1] - 2 * u[k] + u[k + 1]) * K / 2;
        }
        g[nx] = u[nx] + (u[nx - 1] - u[nx]) * K / 2;
        // Tridiagonal matrix algorithm
        // pk(forward steps)
        p[0] = g[0] / a1;
        for (int k = 1; k < nx; k++)
            p[k] = -q[k] * (g[k] - bi * p[k - 1]) / bi;
        // back substitution
        u[nx] = (g[nx] - p[nx - 1] * bi) / (an + bi * q[nx - 1]);
        for (int k = nx - 1; k >= 0; k--) {
            u[k] = q[k] * u[k + 1] + p[k];
        }
        fx << i * dt << ' ' << u[x] << '\n';  // x=2l/3
    }
    for (int i = 0; i <= nx; i++) {
        f2 << i * h << ' ' << u[i] << '\n';  // result in t = 2tau
    }
    std::ofstream f5("f5.txt");
    // t = 5tau
    for (int i = nt2 + 1; i <= nt; i++) {
        g[0] = u[0] + (a0*u[0] + u[1])*K / 2 + b0*dt;
        //g=u(t)+Bdt
        for (int k = 1; k < nx; k++) {
            g[k] = u[k] + (u[k - 1] - 2 * u[k] + u[k + 1]) * K / 2;
        }
        g[nx] = u[nx] + (u[nx - 1] - u[nx]) * K / 2;
        // Tridiagonal matrix algorithm
        // pk-forward steps
        p[0] = g[0] / a1;
        for (int k = 1; k < nx; k++) {
            p[k] = -q[k] * (g[k] - bi * p[k - 1]) / bi;
        }
        // back substitution
        u[nx] = (g[nx] - p[nx - 1] * bi) / (an + bi * q[nx - 1]);
        for (int k = nx - 1; k >= 0; --k) {
            u[k] = q[k] * u[k + 1] + p[k];
        }
        fx << i * dt << ' ' << u[x] << '\n';  // x=2l/3
    }
    for (int i = 0; i <= nx; i++)
        f5 << i * h << ' ' << u[i] << '\n'; // result in t = 5tau

    f2.close();
    f5.close();
    fx.close();
    delete[]u;
    delete[]g;
    delete[]p;
    delete[]q;
    return 0;
}