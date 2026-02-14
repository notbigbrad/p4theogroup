#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>
#include <utility>
#include <functional>
#include <fstream>
#include <string_view>

#include "./constants.hpp"

using namespace std;
using namespace constants;

namespace EoS
{
    // Constant density equation of state
    std::function<long double(long double)> equation_of_state_const(long double rho_0)
    {
        return [rho_0](long double P) -> long double { return rho_0; };
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

    // General white dwarf and neutron star equations of state

    // Neutron starn and White dwarf

    long double chi(double x)
    {
        long double sqrt_term = sqrtl(1.0L + x * x);
        return (1.0L / (8.0L * __pi * __pi)) * ((x * ((2.0L * x * x) / 3.0L - 1.0L) * sqrt_term) - logl(x + sqrt_term));
    }

    long double phi(double x)
    {
        long double sqrt_term = sqrtl(1.0L + x * x);
        return (1.0L / (8.0L * __pi * __pi)) * ((x * ((2.0L * x * x) / 3.0L - 1.0L) * sqrt_term) + logl(x + sqrt_term));
    }

    // Newton-Raphson solver

    double solve_for_x(std::function<long double(long double)> f, long double x0 = 1.0L)
    {
        const long double tol = 1e-12;
        const int max_iter = 100;
        long double x = x0;
        long double h = 1e-6;

        for (int i = 0; i < max_iter; ++i)
        {

            long double fx = f(x);

            long double dfx = (f(x + h) - f(x - h)) / (2.0L * h);

            long double dx = -fx / dfx;
            x += dx;

            if (std::fabsl(dx) < tol) return x;

            if (x <= 0.0L) x = tol;
        }

        throw std::runtime_error("Newton solver failed to converge");
    }

    // NS

    double equation_of_state_neutron(double P)
    {
        long double target = P / 2.0e36;

        long double x = solve_for_x( [&](long double xx) { return phi(xx) - target; } );

        return 2.0e36 * chi(x) / (__c * __c);
    }

    // WD

    double equation_of_state_white_dwarf(double P)
    {
        long double target = P / 1.42180e25;

        long double x = solve_for_x( [&](long double xx) { return phi(xx) - target; } );

        return std::pow(x / 1.0088e-2, 3.0);
    }

    // Read CSV
    vector<vector<long double>> read_csv(string f_name)
    {
        vector<vector<long double>> data;   // Data Vector
        ifstream file(f_name.c_str());      // Open File

        string line;

        while(getline(file, line))
        {
            vector<long double> line_data;
            size_t start = 0;

            for (size_t i = 0; i < line.length(); i++)
            {
                if (line[i] == ',')
                {
                    line_data.emplace_back(stold(line.substr(start, i - start)));
                }
            }
            data.emplace_back(line_data);
        }

        return data;
    }

    class interpolating_polynomial
    {
        public:
            interpolating_polynomial(string f_name)
            {
                data = read_csv(f_name);
            }
            interpolating_polynomial(char *f_char_arr)
            {
                string f_name(f_char_arr);
                data = read_csv(f_name);
            }
            long double density(long double P)
            {
                long double rho;

                // constant, linear, quadratic, cubic, (a, b)

                size_t pos;

                for (size_t i = 0; i < data.size(); i++)
                {
                    if (P < data[pos][4])
                    {
                        pos = i - 1; // If the current minimum pressure is lower than desired pressure we have gone too far so stop
                        i = data.size(); // TESTING CONDITION REMOVE IN FINAL VERSION
                        break;
                    }
                }

                rho = data[pos][0] + data[pos][1] * P + data[pos][2] * P * P + data[pos][3] * P * P * P;

                return rho;
            }
        private:
            vector<vector<long double>> data;
    };

    // Interpolated numerical equations of state
    interpolating_polynomial SLy4_poly("./Jupyter/SLy4_coefficients.csv");
    long double SLy4(long double P)
    {
        return SLy4_poly.density(P);
    }
    interpolating_polynomial FPS_poly("./Jupyter/FPS_coefficients.csv");
    long double FPS(long double P)
    {
        return FPS_poly.density(P);
    }

}
