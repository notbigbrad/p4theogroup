#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>
#include <utility>

#include "./constants.hpp"
#include "./butcher_tableau.hpp"

using namespace std;
using namespace project;

namespace project
{
    // POSSIBLY BREAK THESE OUT INTO SEPARATE FILES AND HAVE THEM BE PARAMETERS IN FUNCTION CALL
    long double equation_of_state(long double P, long double n, long double K)
    {
        return powl(P/K, n / (n+1.0L)); // Basic polytropic equation of state
    };
    long double mass_continuity(long double r, long double rho)
    {
        return __4pi * r * r * rho;
    };
    long double tolman_oppenheimer_volkoff(long double r, long double P, long double M, long double rho)
    {
        long double term1 = (rho + P / __c2);
        long double term2 = (M + 4.0L * __pi * r * r * r * P / __c2);
        long double term3 = r * (r - 2.0L * __gravitation * M / __c2);
        return -(__gravitation * term1 * term2) / term3;
    };

    // Apply the equation of state to create generic first order couple ODEs for the solver
    // Solver must take ODEs in this form

    const long double K_white_dwarf_non_rel = 1.0L / 20.0L * powl((3.0L / __pi), (2.0L / 3.0L)) * ((__planck * __planck) / __mass_electron) * powl((1.0L / (2.0L * __mass_proton)), (5.0L / 3.0L));
    const long double n_white_dwarf_non_rel = 1.5L;

    // F1 = dy1/dt(t, y1, y2)
    long double generic_first_order_ODE_1(long double r, long double M, long double P)
    {
        return mass_continuity(r, equation_of_state(P, n_white_dwarf_non_rel, K_white_dwarf_non_rel));
    };
    // F2 = dy2/dt(t, y1, y2)
    long double generic_first_order_ODE_2(long double r, long double M, long double P)
    {
        return tolman_oppenheimer_volkoff(r, P, M, equation_of_state(P, n_white_dwarf_non_rel, K_white_dwarf_non_rel));
    };

    // Use generic solver to solve TOV
    pair<vector<long double>, vector<long double>> solve_TOV(vector<long double> inital_conditions, long double h, long double n, long double P_c, long double K)
    {
        // Initial conditions to turn ODE into IVP
        long double P_c = P_c;
        long double M   = 0.0L;
    };

    // void MR_TOV_over_density_logspace(long double n, long double h, long double K, long double rho_min, long double rho_max, int density_runs)
    // {
    //     vector<long double> _xi;
    //     vector<long double> _theta;

    //     tie(_xi, _theta) = project::solve_TOV(h, n, 0);

    //     vector<long double> rho_c;

    //     for (int i = 0; i < density_runs; i++)
    //     {
    //         rho_c.push_back(rho_min * powl(rho_max / rho_min, static_cast<long double>(i) / (density_runs - 1)));
    //     }

    //     for (long double rho : rho_c)
    //     {
    //         long double alpha = sqrtl((n + 1.0L) * K * powl(rho, (1.0L / n - 1)) / (4.0L * __pi * __gravitation));
    //         long double R = alpha * _xi.back();
    //         long double M = 0.0L;
    //         for (int j = 1; j < _xi.size(); j++)
    //         {
    //             long double r = alpha * _xi[j];
    //             long double dr = alpha * (_xi[j] - _xi[j - 1]);
    //             long double rho_r = rho * powl(_theta[j], n);

    //             M += 4.0L * __pi * (r * r) * rho_r * dr;
    //         }

    //         cout << M << "," << R << endl;
    //     }
    // }
}