#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>
#include <utility>
#include <functional>

#include "./constants.hpp"

using namespace std;
using namespace constants;

namespace EoS
{
    // Constant density equation of state
    std::function<long double(long double)> equation_of_state_const(long double rho_0)
    {
        return [rho_0](long double P) -> long double
        { return rho_0; };
    }

    // Polytropic equation state
    constexpr long double polytropic_equation_of_state(long double P, long double n, long double K)
    {
        return powl(P / K, n / (n + 1.0L)); // Basic polytropic equation of state
    }
    // Non-relativistic white dwarf polytropic equation state
    const long double mu_e_wd = 2.0L;
    const long double n_wd_nonR = 1.5L;
    const long double K_wd_nonR = 1.0L / 20.0L * powl((3.0L / __pi), (2.0L / 3.0L)) * ((__planck * __planck) / __mass_electron) * powl((1.0L / (mu_e_wd * __mass_proton)), (5.0L / 3.0L));
    constexpr long double equation_of_state_wd_non_r(long double P)
    {
        return polytropic_equation_of_state(P, n_wd_nonR, K_wd_nonR);
    }
    // Relativistic white dwarf polytropic equation state
    const long double n_wd_R = 3.0L;
    const long double K_wd_R = 1.0L / 8.0L * powl((3.0L / __pi), (2.0L / 3.0L)) * (__planck * __c) * powl((1.0L / (mu_e_wd * __mass_proton)), (4.0L / 3.0L));
    constexpr long double equation_of_state_wd_r(long double P)
    {
        return polytropic_equation_of_state(P, n_wd_R, K_wd_R);
    }
    // Interpolated numerical equations of state
    constexpr long double SLy4(long double P)
    {
        // csv data format:
        // a,b,d,c,x_min,discard
        // linear sorted x_min
        // goto minimum x_min < x
        // return a+b*x+c*x^2+d*x^3
        return 0;
    }
    constexpr long double FPS(long double P)
    {
        return 0;
    }

    // General white dwarf and neutron star equations of state

    // Neutron star

    long double chi(double x)
    {
        long double sqrt_term = sqrtl(1.0L + x * x);
        return (1.0L / (8.0L * __pi * __pi)) * ((x * ((2.0L * x * x) / 3.0L - 1.0L) * sqrt_term) - logl(x + sqrt_term));
    }

    // Neutron starn and White dwarf

    long double phi(double x)
    {
        long double sqrt_term = sqrtl(1.0L + x * x);
        return (1.0L / (8.0L * __pi * __pi)) * ((x * ((2.0L * x * x) / 3.0L - 1.0L) * sqrt_term) + logl(x + sqrt_term));
    }

    // Newton-Raphson solver

    double solve_for_x(std::function<double(double)> f,
                       double x0 = 1.0)
    {
        const double tol = 1e-12;
        const int max_iter = 100;
        double x = x0;

        for (int i = 0; i < max_iter; ++i)
        {

            double fx = f(x);

            double h = 1e-6;
            double dfx = (f(x + h) - f(x - h)) / (2.0 * h);

            double dx = -fx / dfx;
            x += dx;

            if (std::fabs(dx) < tol)
                return x;

            if (x <= 0.0)
                x = tol;
        }

        throw std::runtime_error("Newton solver failed to converge");
    }

    // NS

    double equation_of_state_neutron(double P)
    {
        if (P <= 0.0)
            return 0.0;

        const double A = 2.0e36;

        double target = P / A;

        double x = solve_for_x(
            [&](double xx)
            { return phi(xx) - target; });

        return A * chi(x) / (__c * __c);
    }

    // WD

    double equation_of_state_white_dwarf(double P)
    {
        if (P <= 0.0)
            return 0.0;

        const double A = 1.42180e25;
        const double B = 1.0088e-2;

        double target = P / A;

        double x = solve_for_x(
            [&](double xx)
            { return phi(xx) - target; });

        return std::pow(x / B, 3.0);
    }

}
