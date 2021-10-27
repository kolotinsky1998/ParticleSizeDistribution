/////////////////////////////////////
// This program is written by ///////
// Daniil Kolotinskii and     ///////
// and Vladislav Nikolaev     ///////
// to reproduce results from  ///////
// the article "Nucleation and///////
// growth of lead oxide       ///////
// particles in liquid        ///////
// lead-bismuth eutectic"     ///////
// published by               ///////
// Kristof Gladinez           ///////
// //////////////////////////////////

#include <vector>
#include "math.h"
#include "iostream"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "read.hpp" 

using std::cout;
using std::vector;
using std::endl;

int main(int argc, char *argv[]){

        ///////////////////////////////////////
        //////// set calculation flags ////////
        ///////////////////////////////////////

	int self_consistent_flag = 1; //0 - oxygen concentration is read from file, 1 - oxygen concentration is calculated self-consistently
        int distribution_out_flag = 1; // 0 - print time dependent parameters, 1 - print particle size distribution function

        ///////////////////////////////////////
        ////// load oxygen Gladinez data //////
        ///////////////////////////////////////

        ReadC_O readC_O;
        gsl_interp_accel *acc;
        gsl_spline *spline;
        acc = readC_O.acc;
	spline = readC_O.spline;	

	///////////////////////////////////////
	////// used parameters ////////////////
	///////////////////////////////////////
	
<<<<<<< HEAD
	double J_0 = 1e35; //How many nuclei are generated in unit volume in unit time (pre-exponential factor) [nuclei/m**3/s]
	const double rho_PbO = 9450; //Density of solid PbO [kg/m**3]
	const double M_PbO = 0.22320; //Molecular mass of PbO [kg/mole]
	const double V_PbO = M_PbO/rho_PbO; //Molar volume of PbO [m**3/mole]
	const double M_O = 0.016; // Molar mass of O [kg/mole]
	const double pi = 3.1415; // Pi constant
	double r_c = 1e-9; //Critical radius of nuclei [m]
	const double dt = 0.25; // Integration time step [s]
	const double tmax = 13.33*60*60; //Hole time of the system evolution [s]
	const int N = (int) (tmax/dt); // Number of timesteps
	const int n = 200; // Number of grid nodes of discrete particle size distribution function 
	const double T_0 = 673; //Initial temperature [K]
	const double T_end = 564; //Final temperature [K]
	const double R = 8.31; //Universal gas constant [J/mole/K]
	const double k_B = 1.38e-23; //Boltzmann's constant [J/K/mole]
	const double C_O_0 = 8.46e-5; //Initial oxygen concentration [wt %]
	const double M_LBE = 0.208; //Molar mass of LBE [kg/mole]
	const double N_A = 6.02e23; // Avogadro's constant [1/mole]
	double Z; //Zeldovich constant 
	double k = 1.2; //
=======
	double J_0 = 1e35; //nuclei/m**3/s
	double rho_PbO = 9450; //kg/m**3
	double M_PbO = 0.22320; //kg/mol
	double V_PbO = M_PbO/rho_PbO; //m**3/mole
	double M_O = 0.016; // kg/mol
	double pi = 3.1415;
	double r_c = 1e-9; //m //requires special attention
	double dt = 1.0; //s
	double tmax = 13.33*60*60; //s
	int N = (int) (tmax/dt); // number of timesteps
	int n = 200; // number of grid nodes
	double T_0 = 723; //K
	double T_end = 563; //K
	double R = 8.31; //J/mole/K
	double k_B = 1.38e-23; // J/K/mol
	double C_O_0 = 1e-4; //wt %
	double M_LBE = 0.208; // kg/mole
	double N_A = 6.02e23;
	double Z;
	double k = 1.2;
>>>>>>> 6a35c52ed0c30cba8d3a2171009b71f736ceaecf

	///////////////////////////////////////
	////// used time dependent parameters /
	///////////////////////////////////////
	
	vector<double> t(N); //time [s]
	vector<double> T(N); //temperature [K]
	vector<double> D_O(N); //Oxygen diffusion in PbO [m**2/s]
	vector<double> rho_LBE(N); //density of LBE [kg/m**3]
	vector<double> C_Os(N); //Oxygen concentration of saturation [wt %]
	vector<double> C_O(N); //Oxygen concentration [wt %]
	vector<double> sigma(N); //surface energe for nuclei formation [J/m**2]
	vector<double> V_LBE(N); //LBE molar volume [m**3/mole]
	vector<double> mu(N); //Dynamics viscocity of LBE [N*/m**2]

	
	for(int p = 0; p<N; p++){
		t[p] = p*dt;
	}

	for (int p = 0; p<N; p++){
		T[p] = T_0  - (T_0-T_end)*t[p]/tmax;
	}

	for (int p = 0; p<N; p++){
		C_Os[p] = pow(10, 2.64-4426.0/T[p]);
	}

	for (int p = 0; p<N; p++){
		D_O[p] = 2.391e-8*exp(-43073/R/T[p]);
	}

	for (int p = 0; p<N; p++){
		C_O[p] = C_O_0;
	}

	for (int p = 0; p<N; p++){
		rho_LBE[p] = 11096.0 - 1.3236*T[p];
		V_LBE[p] = M_LBE/rho_LBE[p];
	}

	for (int p = 0; p<N; p++){
		sigma[p] = 0.1976 -2.4075e-4*T[p];
	}

	for (int p = 0; p<N; p++){
		mu[p] = 4.91e-4*exp(754.1/T[p]);
	}
	
	///////////////////////////////////////
	////// arrays to find PSD /////////////
	///////////////////////////////////////
	
	vector<double> particles_amount; //Amount of particles in one population
	vector<double> particles_radius; //Radius of particles in one population [m]
	vector<double> particles_rate; 
	vector<double> grid(n); //Grid for particle size distribution function bining
	vector<double> distribution(n-1); //Particle size distribution function
	double v;
	grid[0] = r_c;
	for (int j = 1; j<n; j++){
			grid[j] = grid[j-1]*pow(k,1.0/3.0);
	}
	for (int j = 0; j<n-1; j++){
			distribution[j] = 0;
	}
	
	////////////////////////////////////////
	//////// main solution block ///////////
	////////////////////////////////////////
	double I; //PbO formation rate [J/K/mole]
	double dG_v; //Volume energy of nuclei formation [J/m**3]
	double dG_c; //Surface energy of nucles formation [J/m**2]
	double J; //Rate of nuclei formation [nuclei/m**3/s]
	double S; //Store oxygen concentration in PbO nuclei form [wt %]
	int size; //Number of particle populations
	for (int p = 0; p<N-1; p++){
		//If calculation is launched with given C_O time function, get C_O from data interpolation
		if (self_consistent_flag == 0){
		    C_O[p] = gsl_spline_eval(spline, T[p], acc)*C_Os[p];
		}
	        
		//If C_O greater that C_Os (saturation) than add population with given calculated radius and amount of particles
		if (C_O[p] > C_Os[p]){
		    I = R*T[p]*log(C_O[p]/C_Os[p])/2.0;
		    dG_v = -2*I/(V_PbO);
		    r_c = -2*sigma[p]/dG_v;
		    dG_c = 16*pi*pow(sigma[p],3)/(3*dG_v*dG_v);
		    J_0 = N_A/V_LBE[p];
		    J_0 *= k_B*N_A*T[p]/(3*pi*V_LBE[p]*mu[p]);
		    Z = 1/(2*pi*r_c*r_c)*V_PbO/N_A*sqrt(sigma[p]/(k_B*T[p]));
		    J_0 *= Z;
		    J = J_0*exp(-dG_c/(k_B*T[p]));
		    particles_amount.push_back(J*dt);
		    particles_radius.push_back(r_c);
                    size = particles_amount.size()-1;
		//Do not create new population
		} else {
		    J = 0;
		    size = particles_amount.size();
                }
                //Update particle position
		for (int i=0; i<size; i++) {
		    particles_radius[i] = sqrt(pow(particles_radius[i],2) + 2*D_O[p]*(C_O[p]-C_Os[p])*rho_LBE[p]*V_PbO/M_O*dt);
		}
                //If output flag set NOT to print particle size distribution
                if (distribution_out_flag == 0){
		    cout << t[p] << "\t"<< T[p] << "\t"<< C_O[p] << "\t"<< C_Os[p] << "\t" << C_O[p]/C_Os[p] << "\t" << J << "\t" << endl;
		}

		//Update concentration of dissolved oxygen
		S = 0;
		size = particles_amount.size();
		for (int i=0; i< size; i++){
		    S += 4.0/3.0*pi*pow(particles_radius[i],3)*particles_amount[i]*rho_PbO*M_O/rho_LBE[p]/M_PbO*100;
		}	
                //If oxygen concentration is calculated self-consistenly
                if (self_consistent_flag == 1){
		    C_O[p+1] = C_O_0 - S;
		}
		//If output flag set to print particle size distribution
		if (distribution_out_flag == 1){
		    if (p % int(N/360) == 0) {
			for (int j=0; j<n-1; j++) {
			    distribution[j] = 0;
			}
			for (int i=0; i<size; i++) {
                	    for (int j=0; j<n-1; j++) {
                                if ((particles_radius[i]>grid[j])&&(particles_radius[i]<grid[j+1])){
                                distribution[j] += particles_amount[i];
                        	}
                	    }
			}

			for (int j=0; j<n-1; j++) {
                	    cout << grid[j] << "\t" << distribution[j] << endl;
        		}
			cout << "\n\n";
		    }
		}	
	}
	
	return 0;
}
