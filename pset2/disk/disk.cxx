#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>

void print(const std::vector<double>& Rs, const std::vector<double>& Sigmas, 
    const double time, const std::string filename)
{
    std::ofstream file;
    file.open(filename);

    file << "# time = " << time << " [s]" << std::endl;
    file << "# R Sigma" << std::endl;

    for(int i=0; i<Rs.size(); ++i)
    {
        file << Rs[i] << " " << Sigmas[i] << std::endl;
    }

    file.close();
}

void apply_boundary_conditions( std::vector<double>& Sigmas)
{
    Sigmas[0] = 0; // dirichlet
    Sigmas[Sigmas.size()-1] = 0; // dirichlet

}

int main(int argc, char const *argv[])
{
    
    const double G = 6.67e-8; // [dyne cm^2 g^-2]
    const double alpha = .001;
    const double T_0 = 150; // [K] MMSN Temperature at 1 AU (Chiang & Goldreich 1997)
    const double M_sun = 1.99e33; // [g]
    const double M_earth = 5.67e27; // [g]

    const double k_boltzmann = 1.38e-16; // [erg K^-1]
    const double mu = .67; // mean molecular weight
    const double m_proton = 1.67e-24; // [g]
    const double c_s_0 = std::sqrt(k_boltzmann * T_0 / (mu * m_proton));

    const double C_diffusive = 2 * alpha * std::pow(c_s_0,2) / std::sqrt(G * M_sun);

    const double AU = 1.496e13; // [cm]

    const double R_min = .1 * AU;
    const double R_max = 10 * AU;
    const double R_earth = 1 * AU;

    const double N_zones = 100;


    std::vector<double> Rs (N_zones);
    std::vector<double> Sigmas (N_zones);

    const double dR = (R_max - R_min) / (N_zones-3);

    for (int i=0; i<N_zones; ++i)
    {
        Rs[i] = dR * (i-.5) + R_min;
        if (Rs[i] <= 0) throw std::logic_error("Can't have negative radii");
    }

    for (int i=0; i<N_zones; ++i)
    {
        Sigmas[i] = (M_earth / std::pow(AU,2)) 
                    * std::exp( -std::pow(Rs[i] - R_earth, 2) / (.05*AU*AU) ) ;
    }
    apply_boundary_conditions(Sigmas);
    std::vector<double> Sigmas_old = Sigmas;

    print(Rs, Sigmas, 0, std::string("initial"));

    const int k_max = 10000;
    const int checkpoints = 20;
    const int print_every = k_max / checkpoints;
    int previous_checkpoint = -1;
    double time = 0;
    for (int k=0; k<k_max; ++k)
    {

        if (k >= ((previous_checkpoint+1)*print_every))
        {
            const int current_checkpoint = k / print_every; 
            std::stringstream ss;
            ss << std::setw(4) << std::setfill('0') << current_checkpoint;
            std::string checkpoint_str;
            ss >> checkpoint_str;

            std::cout << "Checkpoint: " << current_checkpoint << std::endl;
            print(Rs, Sigmas, time, "checkpoint_"+ checkpoint_str + ".dat");
            previous_checkpoint = current_checkpoint;
        }

        // DETERMINE TIME STEP
        const double CFL = .0002;
        const double dt = CFL * std::sqrt(dR) / C_diffusive;

        // EVOLVE
        for (int i=1; i<N_zones-1; ++i)
        {

            Sigmas[i] = Sigmas_old[i] + 
                (C_diffusive*dt/(Rs[i]*dR*dR))
                * ( std::sqrt(Rs[i]+(.5*dR))*( std::pow(Rs[i+1],2)*Sigmas_old[i+1] 
                                              -std::pow(Rs[i  ],2)*Sigmas_old[i  ])
                  - std::sqrt(Rs[i]-(.5*dR))*( std::pow(Rs[i  ],2)*Sigmas_old[i  ]
                                              -std::pow(Rs[i-1],2)*Sigmas_old[i-1]) );

            if (Sigmas[i] < 0) throw std::runtime_error("negative density");
        }

        apply_boundary_conditions(Sigmas);

        Sigmas_old = Sigmas;
        time += dt;


        // if (time >= (previous_checkpoint+1)*time_between_checkpoints)
        // {
        //     const unsigned int long current_checkpoint = time / time_between_checkpoints; 
        //     std::stringstream ss;
        //     ss << std::setw(4) << std::setfill('0') << current_checkpoint;
        //     std::string checkpoint_str;
        //     ss >> checkpoint_str;

        //     std::cout << "Checkpoint: " << current_checkpoint << std::endl;
        //     domain.print_to_file(checkpoint_prefix + "checkpoint_"+ checkpoint_str + ".dat");
        //     previous_checkpoint = current_checkpoint;
        // }

    }



    return 0;
}