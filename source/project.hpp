#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>
#include <functional>

using namespace std;

namespace project
{
    // Lane-Emden
    void LE_mass_radius(long double n, long double K, long double P_c_min, long double P_c_max, int runs, function<long double(long double)> equation_of_state);
    void LE_n(long double h, long double n_min, long double n_max, long double n_h, long double xi_max);
    
    // Tolman-Oppenheimer-Volkoff
    void TOV_mass_radius(long double P_c_min, long double P_c_max, int runs, function<long double(long double)> equation_of_state);
    void TOV_Pr(long double h, long double P_c, function<long double(long double)> equation_of_state);
    void TOV_Pr_const_rho(long double h, long double P_c);
}