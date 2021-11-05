#include<iostream>
#include<cmath>
#include<fstream>

const double pi = 3.141592653589793238463;
double func(double x)
{
    double f = 1 / x - tan(x);
    return f;
}
// searching root of tg(y)=1/y in [a, b]
double binsearch(double(*f)(double), double a, double b, double eps)
{
    if (f(a) == 0) {
        return a;
    }
    if (f(b) == 0) {
        return b;
    }
    double mid;
    while (std::abs(b - a) > eps){
        mid = (a + b) / 2;
        if (f(mid) == 0) {
            return mid;
        }
        if (f(mid) * f(a) < 0) {
            b = mid;
        }
        else {
            a = mid;
        }
    }
    return mid;
}
// Fourier series
double sum(int N, const double* yk, const double* ck, double x, double t) {
    double s = 0;
    for (int i = 0; i <= N; ++i) {
        s += cos(yk[i] * (x - 1)) * ck[i] * exp(-yk[i] * yk[i] * t);
    }
    return s;
}

int main() {
    const int nx = 250; // x steps
    const int nt = 625000; // time steps in 5tau
    const int N = 50; // max sum index. Number of terms required
    const double eps = 2e-6; // binary search precision for eigenvalues
    auto *yk = new double[N + 1]; // eigenvalues
    auto *ck = new double[N + 1]; // fourier coeffs
    for (int i = 0; i <= N; i++) {
        yk[i] = binsearch(func, pi*i, pi*(i + 1), eps); // binary search of eigenvalues
        ck[i] = 2 * (sin(yk[i]) + 2 * sin(yk[i] / 10) * cos(yk[i] / 3)) / (yk[i] * (1 + sin(yk[i]) * sin(yk[i])));
    }
    std::ofstream f2("f2.txt");
    const double dx = 1./nx; //x step
    double x = 0; // current x
    // t = 2tau; 1 term in sum required
    for (int i = 0; i <= nx; ++i) {
        f2 << x << ' ' << sum(1, yk, ck, x, 2) - 1 << '\n';
        x += dx;
    }
    std::ofstream f5("f5.txt");
    x = 0;
    // t = 5tau; 1 term in sum required
    for (int i = 0; i <= nx; ++i) {
        f5 << x << ' ' << sum(1, yk, ck, x, 5) - 1 << '\n';
        x += dx;
    }
    const double dt = 5. / nt;  // time step
    std::ofstream fx("fx.txt");
    x = 2. / 3;
    double t = 0; //current time
    // x=2l/3; N=50 terms in sum required
    for (int i = 0; i <= nt; ++i) {
        fx << t << ' ' << sum(N, yk, ck, x, t) - 1 << '\n';
        t += dt;
    }
    f2.close();
    f5.close();
    fx.close();
    delete[]yk;
    delete[]ck;
    return 0;
}