#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>
#include <utility>
#include <functional>
#include <fstream>
#include <string>

#include "./constants.hpp"
#include "../gen/Runge_Kutta.hpp"

using namespace std;
using namespace constants;

namespace project
{
    //// Wrapper to use the generic runge kutta scheme in the fourth order to solve LE



    // Re-wrap functions to be in the desired generic form
    long double theta(long double xi, array<long double, 2> y)
    {
        return y[0];
    }
    std::function<long double(long double, std::array<long double, 2>)> get_dtheta(long double n)
    {
        return [n](long double xi, std::array<long double, 2> y) -> long double {
            return -2.0L / xi * y[0] - (powl(y[1], n));
        };
    }



    // Generalised solver for LE
    pair< vector<long double>, array<vector<long double>, 2> >  LE_RK_wrapper(long double h, int n, function<bool( long double r, array<long double, 2> y )> solved_condition)
    {
        // Initial conditions to turn ODE into IVP
        long double theta_0 = {1.0L - (h * h) / 6.0L}; // Avoid singularity by using xi0=h
        long double dtheta_0 = -h / 3.0L;
        long double init[2] = { dtheta_0, theta_0 };

        // List of ODEs
        array<function<long double( long double r, array<long double, 2> y )>, 2> ODES = {get_dtheta(n), &theta};

        // Apply RK scheme to problem
        auto t = Runge_Kutta::RK4_2( h, init, ODES, solved_condition );

        return t;
    }



    // Demonstrative for arbitrary n plot
    bool LE_criterion( long double xi, array<long double, 2> y)
    {
        if (xi > 20.0L) { return true; }
        else { return false; }
    }

    // Run for arbitrary n values
    void LE_n(long double h, long double n_min, long double n_max, long double n_h, long double xi_max, string f_name)
    {
	ofstream file(f_name);

        vector<string> csv;

        vector<long double> r;
        array<vector<long double>, 2> y;

        // Unwrap answer to ODE
        tie(r, y) = LE_RK_wrapper(h, n_min, LE_criterion);

        // Resize csv string array
        csv.resize(y[1].size());

        // Grab desired values
        for (int j = 0; j < size(y[1]); j++)
        {
            csv[j] += to_string(r[j]) + "," + to_string(y[1][j]);
        }

        for (long double i = n_min + n_h; i <= n_max; i += n_h)
        {
            vector<long double> r;
            array<vector<long double>, 2> y;

            // Unwrap answer to ODE
            tie( r, y ) = LE_RK_wrapper(h, i, LE_criterion);

            // Grab desired values
            for (int j = 0; j < size(y[1]); j++)
            {
                csv[j] += "," + to_string(y[1][j]);
            }
        }
        for (auto& entry : csv)
        {
            file << entry << endl;
        }
    }



    // Solves for the radial extent of star
    bool LE_criterion2( long double xi, array<long double, 2> y)
    {
        if (y[1] <= 0 || isnan(y[1]) || xi > 200.0L) { return true; }
        else { return false; }
    }

    void LE_mass_radius(long double h, long double n, long double K, long double P_c_min, long double P_c_max, int runs, function<long double(long double)> equation_of_state, string f_name)
    {
	ofstream file(f_name);

        vector<long double> _xi;
        array<vector<long double>, 2> _thetas;
        vector<long double> P_c_list;

        // Use Lane-Emden solver
        tie(_xi, _thetas) = LE_RK_wrapper(h, n, LE_criterion2);

        // Get theta from solver
        vector<long double> _theta = _thetas[1];

        // Create log distributed central pressures to test
        for (int i = 0; i < runs; i++)
        {
            P_c_list.push_back(P_c_min * powl(P_c_max / P_c_min, static_cast<long double>(i) / (runs - 1)));
        }

        // Loop through all possible central pressures
        for (long double P_c : P_c_list)
        {
            long double rho = equation_of_state(P_c);
            // long double rho = P_c;
            long double alpha = sqrtl((n + 1.0L) * K * powl(rho, (1.0L / n - 1)) / (4.0L * __pi * __gravitation));
            long double R = alpha * _xi.back();
            long double M = 0.0L;
            for (int j = 1; j < _xi.size(); j++)
            {
                long double r = alpha * _xi[j];
                long double dr = alpha * (_xi[j] - _xi[j - 1]);
                long double rho_r = rho * powl(_theta[j], n);

                M += 4.0L * __pi * (r * r) * rho_r * dr;
            }

            file << M << "," << R << endl;
        }
    }
}
