#ifndef _HYDRO_H_
#define _HYDRO_H_

#include <string>
#include <vector>

// UNITS
#define K_BOLTZMANN (1.38065e-16)   // boltzmann constant [erg K^-1]
#define M_PROTON (1.672622e-24)     // proton mass [g]
#define JOULE (1e7)                 // 1 Joule of energy [erg]
#define KILOMETER (1e5)             // 1 km of distance [cm]
#define GRAV_ACCEL (980)            // Earth gravitational accel [cm s^-2]
#define YEAR (3600*24*365.24)

#define CFL (.2)


class Cell
{
public:
    double R; // current radius

    // primitives
    double V;       // specific volume [cm^3 g^-1]
    double u;       // fluid velocity [cm s^-1]

    double P;       // pressure [dyne cm^-2]

    // useful
    double r;          // stretched reference coordinate
    double V_0;        // specific density; uniform in reference coords
    double V_initial;  // specific density initially =/= V_0, for eq. of state

    double dR;         // current zone spacing, ENTIRE ZONE [cm]
    double P_initial;        // initial pressure [dyne cm^-2], for eq. of state

    double gamma;

    double sound_speed() const;
    void apply_Equation_Of_State();

};

class Domain
{
public:
    Domain(const int N_cells, const double height, 
           const double gamma, const double mu);
    
    void print_to_file(const std::string filename) const;

    void add_blast(const double E_blast, const double altitude);

    void timestep();

    double get_current_time() const;

    void interpolate();

    void enforce_boundary_conditions();

    void calc_dR();

protected:
    double determine_dt() const;

private:
    const int N_cells;
    std::vector<Cell> cells;
    std::vector<Cell> cells_old;
    double time;

    const double gamma; // adiabatic index
    const double mu;    // mean molecular weight


};


#endif