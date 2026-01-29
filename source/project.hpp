#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>

using namespace std;

namespace project
{
    // Lane-Emden
    pair<vector<long double>, vector<long double>> solve_LE(long double h, long double n);
    void solve_LE_over_n(long double h, long double n_min, long double n_max, long double n_h, long double end_point);
    void MR_LE_over_density_logspace(long double n, long double h, long double K, long double rho_min, long double rho_max, int density_runs);
    
    // Tolman-Oppenheimer-Volkoff
    pair<vector<long double>, vector<long double>> solve_TOV(long double h, long double n);
    void MR_TOV_over_density_logspace(long double n, long double h, long double K, long double rho_min, long double rho_max, int density_runs);

    // Temporary
    long double generic_first_order_ODE_1(long double r, long double M, long double P);
    long double generic_first_order_ODE_2(long double r, long double M, long double P);
}