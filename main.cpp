#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <vector>
#include <random>
#include <map>
#include <set>
#include <mutex>
#include <atomic>
#include <thread>
#include <functional>
#include <filesystem>
#include <ctime>
#include <utility>

#include "./source/project.hpp"
#include "./source/equations_of_state.hpp"

using namespace std;
using namespace std::chrono;
using namespace project;
using namespace EoS;

int main(int argc, char **argv)
{   
    // Solve Lane-Emden for n in (-2,10)
    freopen("./output/laneEmdenSolutions.csv", "w", stdout);
    cout << setprecision(4);
    project::LE_n(0.01L, -2.0L, 10.0L, 0.5L, 20.0L);




    // Non-relativistic white dwarfs
    long double min_central_pressure = 1000000.0L;
    long double max_central_pressure = 10000000000000.0L;

    // Mass-Radius using LE
    freopen("./output/LEmassRadius_wd_non_r.csv", "w", stdout);
    cout << setprecision(std::numeric_limits<long double>::digits10 + 1);
    project::LE_mass_radius( n_wd_nonR, K_wd_nonR, min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_wd_non_r );

    // Mass-Radius using TOV
    freopen("./output/TOVmassRadius_wd_non_r.csv", "w", stdout);
    cout << setprecision(std::numeric_limits<long double>::digits10 + 1);
    project::TOV_mass_radius( min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_wd_non_r );




    // Relativistic white dwarfs

    // Mass-Radius using LE
    freopen("./output/LEmassRadius_wd_r.csv", "w", stdout);
    cout << setprecision(std::numeric_limits<long double>::digits10 + 1);
    project::LE_mass_radius( n_wd_R, K_wd_R, min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_wd_r );

    // Mass-Radius using TOV
    freopen("./output/TOVmassRadius_wd_r.csv", "w", stdout);
    cout << setprecision(std::numeric_limits<long double>::digits10 + 1);
    project::TOV_mass_radius( min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_wd_r );




    // P(r) for Constant Rho using TOV
    freopen("./output/TOVconstantdensity.csv", "w", stdout);
    cout << setprecision(std::numeric_limits<long double>::digits10 + 1);
    long double central_const_rho = 10000000.0L;
    project::TOV_Pr_const_rho( 1000000.0L, central_const_rho);




    // Pressure(Radius) TOV plots for white dwarf
    long double central_pressure = 10000000.0L;

    // Non-relativistic white dwarfs
    freopen("./output/Pr_TOV_wd_non_r.csv", "w", stdout);
    cout << setprecision(std::numeric_limits<long double>::digits10 + 1);
    project::TOV_Pr( 1000000.0L, central_pressure, equation_of_state_wd_non_r );

    // Relativistic white dwarfs
    freopen("./output/Pr_TOV_wd_r.csv", "w", stdout);
    cout << setprecision(std::numeric_limits<long double>::digits10 + 1);
    project::TOV_Pr( 1000000.0L, central_pressure, equation_of_state_wd_r );


    // Get maximum thread count
    // unsigned int max_threads = thread::hardware_concurrency();
    // cout << max_threads << " threads available." << endl << setprecision(std::numeric_limits<long double>::digits10 + 1);

    return 0;
}