#include <vector>
#include "math.h"
#include "iostream"

using namespace std;

int main(int argc, char *argv[]){
	

	///////////////////////////////////////
	////// used constants /////////////////
	///////////////////////////////////////
	
	double J_0 = 1e35; //nuclei/m**3/s
	double rho_PbO = 9450; //kg/m**3
	double M_PbO = 0.22320; //kg/mol
	double V_PbO = M_PbO/rho_PbO; //m**3/mole
	double M_O = 0.016; // kg/mol
	double pi = 3.1415;
	double r_c = 1e-9; //m
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

	///////////////////////////////////////
	////// time dependent arrays //////////
	///////////////////////////////////////
	
	vector<double> t(N); //s
	vector<double> T(N); //K
	vector<double> D_O(N); //m**2/s
	vector<double> rho_LBE(N); //kg/m**3
	vector<double> C_Os(N); //wt %
	vector<double> C_O(N); //wt %
	vector<double> sigma(N); //J/m**2
	vector<double> V_LBE(N); //m**3/mole
	vector<double> mu(N);

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
		D_O[p] = 2.391e-10*exp(-43073/R/T[p]);
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
	
	//for (int p = 0; p<N; p++){
	//	cout << t[p] << "\t" << T[p] << "\t" << C_Os[p] << "\t" << D_O[p] << "\t" << C_O[p] << "\t" << rho_LBE[p] << endl;
	//}

	
	///////////////////////////////////////
	////// arrays to find PSD /////////////
	///////////////////////////////////////
	
	vector<double> particles_amount;
	vector<double> particles_radius;
	vector<double> particles_rate;
	vector<double> grid(n);
	vector<double> distribution(n-1);
	double v;
	//v = D_O[N-1]*(C_O_0-C_Os[N-1])*rho_LBE[N-1]*V_PbO/M_O/2;
	grid[0] = r_c;
	for (int j = 1; j<n; j++){
			//grid[j] = (double)j/(double)n*sqrt(pow(r_c,2) + 2*tmax*v*100);
			grid[j] = grid[j-1]*pow(k,1.0/3.0);
	}
	for (int j = 0; j<n-1; j++){
			distribution[j] = 0;
	}
	
	////////////////////////////////////////
	//////// main solution block ///////////
	////////////////////////////////////////
	double I;
	double dG_v;
	double dG_c;
	double J;
	int counter;
	double S;
	int size;
	for (int p = 0; p<N-1; p++){
		if (C_O[p] > C_Os[p]){
		    I = R*T[p]*log(C_O[p]/C_Os[p])/2.0;
		    dG_v = -2*I/(2*V_PbO);
			r_c = -2*sigma[p]/dG_v;
		    dG_c = 16*pi*pow(sigma[p],3)/(3*dG_v*dG_v);
			J_0 = N_A/V_LBE[p];
			J_0 *= k_B*N_A*T[p]/(3*pi*V_LBE[p]*mu[p]);
			Z = 1/(2*pi*r_c*r_c)*V_PbO/N_A*sqrt(sigma[p]/(k_B*T[p]));
			J_0 *= Z;
		    J = J_0*exp(-dG_c/(k_B*T[p]));
			//cout << t[p] << "\t" << J << endl;
		    	
		    particles_amount.push_back(J*dt);
		    particles_radius.push_back(r_c);
		    v = D_O[p]*(C_O[p]-C_Os[p])*rho_LBE[p]*V_PbO/M_O/r_c;
		    particles_rate.push_back(v);
            size = particles_amount.size()-1;
		} else {
			J = 0;
			size = particles_amount.size();
        }
		for (int i=0; i<size; i++) {
		    particles_radius[i] += particles_rate[i]*dt;
            v = D_O[p]*(C_O[p]-C_Os[p])*rho_LBE[p]*V_PbO/M_O/particles_radius[i];
		    particles_rate[i] = v;
		}
        //counter = 0;
		/*
                for (int i=0; i<size; i++) {
                    if (particles_rate[i] == 0){
                        counter ++;
                    }
                }
                for (int i=0; i<counter; i++){
                    particles_amount.pop_back();
                    particles_radius.pop_back();
                    particles_rate.pop_back();
                }
		*/
		S = 0;
		size = particles_amount.size();
		for (int i=0; i< size; i++){
			S += 4.0/3.0*pi*pow(particles_radius[i],3)*particles_amount[i]*rho_PbO*M_O/rho_LBE[p]/M_PbO/100.0;
		}	
		C_O[p+1] = C_O_0 - S;
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
	
	size = particles_amount.size();

	for (int i=0; i<size; i++) {
		for (int j=0; j<n-1; j++) {
			if ((particles_radius[i]>grid[j])&&(particles_radius[i]<grid[j+1])){
				distribution[j] += particles_amount[i];
			}
		}
	}
	
	for (int p=0; p<N; p++) {
		cout << t[p] << "\t" << C_O[p] << "\t" << C_O[p] / C_Os[p] << endl;
	}

	cout << "\n\n";
	
	for (int j=0; j<n-1; j++) {
		cout << grid[j] << "\t" << distribution[j] << endl;
	}
	
	return 0;
}
