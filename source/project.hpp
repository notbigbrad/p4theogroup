#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include <thread>
#include <functional>
#include <string>

using namespace std;

namespace project
{
    // Lane-Emden
    void LE_mass_radius(long double h, long double n, long double K, long double P_c_min, long double P_c_max, int runs, function<long double(long double)> equation_of_state, string f_name);
    void LE_n(long double h, long double n_min, long double n_max, long double n_h, long double xi_max, string f_name);

    // Tolman-Oppenheimer-Volkoff
    void TOV_mass_radius(long double h, long double P_c_min, long double P_c_max, int runs, function<long double(long double)> equation_of_state, string f_name);
    void TOV_Pr(long double h, long double P_c, function<long double(long double)> equation_of_state, string f_name);

    // Hydrostatic Equilibrium
    void HE_mass_radius(long double h, long double P_c_min, long double P_c_max, int runs, function<long double(long double)> equation_of_state, string f_name);
    void HE_Pr(long double h, long double P_c, function<long double(long double)> equation_of_state, string f_name);
}
