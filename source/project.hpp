#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>

using namespace std;

namespace project
{
    // Lane-Emden
    // pair<vector<long double>, vector<long double>> solve_LE(long double h, long double n); // OLD
    // void solve_LE_over_n(long double h, long double n_min, long double n_max, long double n_h, long double end_point); // OLD
    // void MR_LE_over_density_logspace(long double n, long double h, long double K, long double rho_min, long double rho_max, int density_runs); // OLD

    // Lane-Emden
    void LE_mass_radius(long double n, long double K, long double P_c_min, long double P_c_max, int runs);
    void LE_n(long double h, long double n_min, long double n_max, long double n_h, long double xi_max);
    
    // Tolman-Oppenheimer-Volkoff
    void TOV_mass_radius(long double P_c_min, long double P_c_max, int runs);
    void TOV_Pr(long double h, long double P_c);
}