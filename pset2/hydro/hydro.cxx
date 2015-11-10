
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <sstream>

#include "hydro.H"


int main(int argc, char const *argv[])
{


    // const double c_s = std::sqrt((7*K_BOLTZMANN*300)/(5*28.97*M_PROTON));
    // std::cout << c_s << std::endl;
    // const double total_sound_crossing_time = 30 * KILOMETER / c_s;
    // const double time_final = 100.5 * total_sound_crossing_time; // [s]
    // const double time_between_checkpoints = 10 * total_sound_crossing_time;
    // const std::string checkpoint_prefix("equilibrium_");
    
    const double time_final = 10.5; // [s]
    const double time_between_checkpoints = 1;
    const std::string checkpoint_prefix("shock_");

    
    const int N_cells = 100-1;
    const double height = 30 * KILOMETER;
    const double gamma = 7./5;
    const double mu = 28.97;

    Domain domain(N_cells, height, gamma, mu);


    unsigned long int previous_checkpoint = -1;   
    const unsigned long int k_max = 1000000;
    for (unsigned long int k=0; k<k_max; ++k)
    {
        const double time = domain.get_current_time();
        if (time > time_final) break;

        if (time >= (previous_checkpoint+1)*time_between_checkpoints)
        {
            const unsigned int long current_checkpoint = time / time_between_checkpoints; 
            std::stringstream ss;
            ss << std::setw(4) << std::setfill('0') << current_checkpoint;
            std::string checkpoint_str;
            ss >> checkpoint_str;

            std::cout << "Checkpoint: " << current_checkpoint << std::endl;
            domain.print_to_file(checkpoint_prefix + "checkpoint_"+ checkpoint_str + ".dat");
            previous_checkpoint = current_checkpoint;
        }

        domain.timestep();
    }

    return 0;
}


//******* DOMAIN ****//////
Domain::Domain(const int N_cells, const double height,
               const double gamma, const double mu)
    : N_cells(N_cells)
    , gamma(gamma)
    , mu(mu)
{
    if (N_cells%2 == 0)
    {
        throw std::invalid_argument("Need an odd N_cells");
    }

    cells.resize(N_cells);

    time = 0;

    const double T_0 = 300; // [K]
    const double V_0 = 1./0.001225; // from Wolfram Alpha

    const double scale_height = (K_BOLTZMANN * T_0) / (mu * M_PROTON * GRAV_ACCEL);

    std::cout << "Scale height: " << scale_height/1e5 << " km" << std::endl;

    const double r_min = 0; // by definition
    const double r_max = scale_height * (1-std::exp(-height/scale_height));
    const double dr = (r_max - r_min)/(N_cells-2);

    const double dR = height/(N_cells-2);
    for (int i = 0; i < N_cells; ++i)
    {

        Cell* c      = &(cells[i]);

        // c->r = (i-1) * dr;
        // c->R = -scale_height * std::log(1-(c->r/scale_height));
        // c->V_0 = V_0;
        // c->V   = V_0 / (1 - (c->r/scale_height));
        // c->V_initial = c->V;

        c->R         = (i-1) * dR;
        c->V_0       = V_0;
        c->V         = V_0 * std::exp( (c->R * mu * M_PROTON * GRAV_ACCEL)
                                     /(K_BOLTZMANN * T_0));
        c->V_initial = c->V;
        c->r         = (1-(V_0/c->V)) 
                        * ( K_BOLTZMANN * T_0 / ( mu * M_PROTON * GRAV_ACCEL));

        if(i%2==0)
        {
            // ZONES
            c->P = (K_BOLTZMANN * T_0) / (c->V * mu * M_PROTON);
            c->P_initial = c->P;

            c->gamma = gamma;
        }
        else
        {
            // EDGES
            c->u   = 0;
        }
    }

    this->calc_dR();

    this->add_blast(1e6 * JOULE, 10*KILOMETER);

    this->enforce_boundary_conditions();
    this->interpolate();
    cells_old = cells;

    this->print_to_file(std::string("initial.dat"));

}

void Domain::timestep()
{

    const double dt = this->determine_dt();


    // ACCELERATE
    for (int i=1; i<N_cells-1; i+=2)
    {
        cells[i].u += -dt*cells[i].V_0*(cells[i+1].P - cells[i-1].P)
        /(cells[i+1].r - cells[i-1].r) - dt*GRAV_ACCEL;
    }
    // // // boundary condition (absorbing)
    {
        const int i = N_cells-2;

        // const double dr_13 = cells[i-1].r - cells[i-3].r;
        // const double dr_35 = cells[i-3].r - cells[i-5].r;
        // const double tol=1e-3;
        // if ( std::abs(1-(dr_13/dr_35)) > tol) throw std::logic_error("need constant dr");
        // const double dr = dr_13/2;
        // const double dP_dr = (- (     cells[i-1].P  )
                              // + (3. * cells[i-3].P/2)
                              // - (     cells[i-5].P/2)) / dr;

        const double dP_dr = (cells[i-1].P - cells[i-3].P)/(cells[i-1].r - cells[i-3].r);

        cells[i].u = cells_old[i].u - dt*cells[i].V_0*dP_dr - dt*GRAV_ACCEL;
    }

    // MOVE
    // do this after accelerating to benefit from an implicit solve
    for (int i=1; i<N_cells-1; i+=2)
    {
        cells[i].R += dt*cells[i].u;
    }
    

    // DENSITIES
    for (int i=2; i<N_cells-2; i+=2)
    {
        cells[i].V = cells[i].V_0 * (cells[i+1].R - cells[i-1].R)
                                   /(cells[i+1].r - cells[i-1].r);
    }

    // EQUATION OF STATE
    for (int i=2; i<N_cells-2; i+=2)
    {
        cells[i].apply_Equation_Of_State();
    }

    this->calc_dR();

    this->enforce_boundary_conditions();
    this->interpolate(); // for convenience with plotting

    cells_old = cells;
    time += dt;

}

void Domain::calc_dR()
{
    // mainly for use in determine_dt
    for (int i=2; i<N_cells-2; i+=2)
    {
        cells[i].dR = cells[i+1].R - cells[i-1].R;
        if (cells[i].dR <= 0) 
        {
            std::cout << "Zone collision for zone " << i << std::endl;
            throw std::runtime_error("zones collided");
        }
    }
}

void Domain::interpolate()
{
    for (int i=2; i<N_cells-2; ++i)
    {
        if (i%2==0)
        {
            // ZONES
            cells[i].R = (cells[i+1].R + cells[i-1].R)/2;
            cells[i].u = (cells[i+1].u + cells[i-1].u)/2;
        }
        else
        {
            // EDGES
            cells[i].P = (cells[i+1].P + cells[i-1].P)/2;
            cells[i].V = (cells[i+1].V + cells[i-1].V)/2;
        }
    }
}

void Domain::enforce_boundary_conditions()
{
    cells[1].u = 0; // hard wall at [1]

    // not used if we're also using absorbing boundary conditions
    cells[N_cells-1].P = -cells[N_cells-3].P; // vacuum at [N_cells-2]
}

double Domain::determine_dt() const
{

    double c_s = cells[2].sound_speed();
    double u_total = std::abs(cells[1].u) + std::abs(cells[3].u) + c_s;
    double dt = cells[2].dR / u_total;

    for (int i=2; i<N_cells-2; i+=2)
    {
        u_total = std::abs(cells[i-1].u) + std::abs(cells[i+1].u) 
                + cells[i].sound_speed();

        double dt_tmp = cells[i].dR / u_total;
        if (dt_tmp < dt)  dt = dt_tmp;
    }

    std::cout << "dt: " << std::scientific << dt << std::endl;

    return CFL * dt;
}

double Domain::get_current_time() const
{
    return time;
}

void Domain::print_to_file(const std::string filename) const
{
    std::ofstream file;
    file.open(filename);

    file << std::scientific << std::left;
    file << "# t: " << time << " [s]" << std::endl;

    file  << std::left << std::setw(16) << "# R [cm]" 
                       << std::setw(16) << "r [cm]" 
                       << std::setw(16) << "rho [g cm^-3]" 
                       << std::setw(16) << "u [cm s^-1]"
                       << std::setw(16) << "P [dyne cm^-2]"  << std::endl;
    for (int i=2; i<N_cells-2; ++i)
    {
        file << cells[i].R << "\t" 
             << cells[i].r << "\t" 
             << 1/cells[i].V << "\t" 
             << cells[i].u << "\t" << cells[i].P << std::endl;
    }
    file.close();
}

void Domain::add_blast(const double E_blast, const double altitude)
{
    // find center of blast:
    int i_blast = 0;
    for (int i=1; i<N_cells-1; i+=2)
    {
        if (cells[i].R > altitude)
        {
            i_blast = i-1; //  previous zone
            break;
        }
    }
    if (i_blast == 0) throw std::invalid_argument("no valid zone at given altitude");
    std::cout << "Adding blast to zone: " << i_blast << std::endl;
    

    // ADD E_BLAST
    const double area=1;

    const double mass = area*cells[i_blast].dR/cells[i_blast].V;
    const double prev_energy = mass * cells[i_blast].P * cells[i_blast].V 
                                    / (gamma-1);
    const double new_energy = E_blast + prev_energy;
    const double new_pressure =   (gamma-1) * new_energy 
                                / (mass * cells[i_blast].V);
    cells[i_blast].P = new_pressure;
    cells[i_blast].P_initial = cells[i_blast].P;
    cells[i_blast].V_initial = cells[i_blast].V;
}


//******* Cell ****//////

void Cell::apply_Equation_Of_State()
{
    P = P_initial * std::pow(V_initial / V, gamma );
}

double Cell::sound_speed() const
{
    const double sound_speed = std::sqrt(gamma * P * V);
    return sound_speed;
}
