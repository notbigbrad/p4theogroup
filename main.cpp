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
	// Solutions check-list

	// Equation of State
	// P(r)
	// M-R

	// Constant Density
	// LE
	//

	// Non-relativistic PeS
	// LE TOV
	// LE TOV

	// Relativistic PeS
	// LE TOV
	//

	//// NO MORE LANE-EMDEN

	// General White Dwarf
	//
	//

	// Basic Neutron Star
	//
	//

	// SLy4
	//
	//

	// FPS
	//
	//

	// ABS
	//
	//

	//// Main execution body

	long double TOV_h = 1000000.0L;
	long double LE_h = 0.001L;

	long double central_pressure = 10000000.0L;

	long double min_central_pressure = 100000.0L;
	long double max_central_pressure = 100000000000000.0L;



	//// Basic Lane-Emden

	cout << "Solving Lane-Emden for n in [-2,10]" << endl;
	project::LE_n(LE_h, -2.0L, 10.0L, 0.5L, 20.0L, "./output/laneEmdenSolutions.csv");



	//// Polytropic Equation of State

	// Non-relativistic white dwarfs

	cout << "Mass-Radius using LE non-relativistic PeS" << endl;
	project::LE_mass_radius(LE_h, n_wd_nonR, K_wd_nonR, min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_wd_non_r, "./output/LEmassRadius_wd_non_r.csv");

	cout << "Mass-Radius using TOV non-relativistic PeS" << endl;
	project::TOV_mass_radius(TOV_h, min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_wd_non_r, "./output/TOVmassRadius_wd_non_r.csv");
	cout << "TOV P(r) Non-relativistic white dwarfs PeS" << endl;
	project::TOV_Pr(TOV_h, central_pressure, equation_of_state_wd_non_r, "./output/Pr_TOV_wd_non_r.csv");

	// Relativistic white dwarfs

	cout << "Mass-Radius using LE relativistic PeS" << endl;
	project::LE_mass_radius(LE_h, n_wd_R, K_wd_R, min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_wd_r, "./output/LEmassRadius_wd_r.csv");

	cout << "Mass-Radius using TOV relativistic Pes" << endl;
	project::TOV_mass_radius(TOV_h, min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_wd_r, "./output/TOVmassRadius_wd_r.csv");
	cout << "TOV P(r) Relativistic white dwarfs PeS" << endl;
	project::TOV_Pr(TOV_h, central_pressure, equation_of_state_wd_r, "./output/Pr_TOV_wd_r.csv");

	// Constant density case for TOV

	cout << "P(r) for Constant Rho using TOV" << endl;
	long double central_pressure_c = 1000000.0L;
	project::TOV_Pr(TOV_h, central_pressure_c, equation_of_state_const(central_pressure_c), "./output/TOVconstantdensity.csv");



	//// General Equations of State, No More Lane-Emden


	// Neutron Star

	// TOV
	cout << "TOV neutron star" << endl;
	project::TOV_mass_radius(TOV_h, min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_neutron, "./output/TOVmassRadius_neutron.csv");
	project::TOV_Pr(TOV_h, central_pressure, equation_of_state_neutron, "./output/Pr_TOV_neutron.csv");

	// HE
	cout << "HE neutron star" << endl;
	project::HE_mass_radius(TOV_h, min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_neutron, "./output/HEmassRadius_neutron.csv");
	project::HE_Pr(TOV_h, central_pressure, equation_of_state_neutron, "./output/Pr_HE_neutron.csv");


	// White Dwarf

	// TOV
	cout << "TOV white dwarf" << endl;
	project::TOV_mass_radius(TOV_h, min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_white_dwarf, "./output/TOVmassRadius_white_dwarf.csv");
	project::TOV_Pr(TOV_h, central_pressure, equation_of_state_white_dwarf, "./output/Pr_TOV_white_dwarf.csv");

	// HE
	cout << "HE white dwarf" << endl;
	project::HE_mass_radius(TOV_h, min_central_pressure, max_central_pressure, (int)(1000), equation_of_state_white_dwarf, "./output/HEmassRadius_white_dwarf.csv");
	project::HE_Pr(TOV_h, central_pressure, equation_of_state_white_dwarf, "./output/Pr_HE_white_dwarf.csv");



	//// Interpolated Numerical Equations of State


	// SLy4

	// TOV
	cout << "TOV SLy4" << endl;
	project::TOV_mass_radius(TOV_h, min_central_pressure, max_central_pressure, (int)(1000), SLy4, "./output/TOVmassRadius_SLy4.csv");
	project::TOV_Pr(TOV_h, central_pressure, SLy4, "./output/Pr_TOV_SLy4.csv");

	// HE
	cout << "HE FPS" << endl;
	project::HE_mass_radius(TOV_h, min_central_pressure, max_central_pressure, (int)(1000), SLy4, "./output/HEmassRadius_SLy4.csv");
	project::HE_Pr(TOV_h, central_pressure, SLy4, "./output/Pr_HE_SLy4.csv");


	// FPS

	// TOV
	cout << "TOV FPS" << endl;
	project::TOV_mass_radius(TOV_h, min_central_pressure, max_central_pressure, (int)(1000), FPS, "./output/TOVmassRadius_FPS.csv");
	project::TOV_Pr(TOV_h, central_pressure, FPS, "./output/Pr_TOV_FPS.csv");

	// HE
	cout << "HE FPS" << endl;
	project::HE_mass_radius(TOV_h, min_central_pressure, max_central_pressure, (int)(1000), FPS, "./output/HEmassRadius_FPS.csv");
	project::HE_Pr(TOV_h, central_pressure, FPS, "./output/Pr_HE_FPS.csv");

	return 0;
}
