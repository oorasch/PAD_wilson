//Massless Schwinger model with Path Activation Determinants, gauge action = Wilson

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <string>
#include <math.h>
#include <limits>
#include <gsl/gsl_statistics.h>
#include "initialize.h"
#include "determinant.h"
#include "constants.h"
#include "updates.h"
#include "debug.h"
#include "observables.h"


using namespace std;

int main(){

	//Auxiliary arrays
 	int neib[constants::V][4];		//Neighbour field
	vector<vector<int>> k_link (constants::V, vector<int> (2, 0)); //link variables for fermion loops
	vector<vector<int>> sig_link (constants::V, vector<int> (2, 0));	//auxiliary sigma link variables, from Hubbard-Stratonovich-Trafo
	vector<int> s_site (constants::V, 1); //loop occupation number of sites
	vector<int> plaq_occ (constants::V, 0); //plaquette occupation number, x = lower left corner
	vector<double> M (constants::V*constants::V, 0.0); // fermion matrix
	vector<double> I_bessel; // storage for ratios of Bessel functions
	//vector<double> M (constants::V*constants::V, 0.0);
	double detM = 0.0; //log of fermion matrix
	
	auto start = std::chrono::high_resolution_clock::now();
	
	neib_init(constants::Nt, constants::Ns, neib);//initialization of the neighbour field
		
	do
	{
		sig_init(sig_link);
		calc_M_array(neib, sig_link, s_site, M); //fill the fermion matrix with the initial configuration
		detM = det(constants::V, M);//calculate the log of the determinant for the initial configuration
	}
	while( detM == 0.0);
	
	//cout << detM << endl;
	
	ofstream outfile1, outfile2;
	outfile1.open("plaq_occ_beta_full3.dat");
	outfile2.open("topo_chrg_beta_full3.dat");
		
	//for(double theta_bar = -1.5; theta_bar <= 1.6; theta_bar += 0.05)
	for(double beta = 0.1; beta <= 3.1; beta += 0.1)
	{
		//constants::theta_bar = theta_bar;
		constants::beta = beta;
		constants::eta = 0.5*constants::beta - 0.5*constants::theta_bar;
		constants::eta_bar = 0.5*constants::beta + 0.5*constants::theta_bar;
		
		bessel_init(I_bessel);
		
//		for(int i=0; i<constants::max_bessel; i++)
//		{
//			cout << beta << " " << I_bessel.at(3) << endl;
//		}
//		cout << "+++++++++++++++++++++++++++++++++++++" << endl << endl;
		
		do_updates(constants::nequi, neib, I_bessel, detM, s_site, plaq_occ, k_link, sig_link);
		for(int i = 0; i < constants::nmeas; i++)
		{
			do_updates(constants::nskip, neib, I_bessel, detM, s_site, plaq_occ, k_link, sig_link);
			//cout << constants::beta << " " << constants::eta << " " << meas_plaq(plaq_occ, I_bessel) << endl;
			outfile1 << meas_plaq(plaq_occ, I_bessel) << endl;
			outfile2 << meas_topo_chrg(plaq_occ, I_bessel) << endl;
		}
	}
	outfile1.close(); outfile2.close();

//	for(int i=0; i<constants::Ns; i++)
//	{
//		for(int j=0; j<constants::Nt; j++)
//		{
//			int k = j + i*constants::Nt;
//			
//			cout << k_link.at(k).at(0) << " " << k_link.at(k).at(1) << " --- " << plaq_occ.at(k) << endl; 
//		}
//		cout << endl;
//	}
	
	pauli_check(neib, k_link);
	flux_check(neib, k_link);
	plaquette_check(neib, k_link, plaq_occ);
	
	auto finish = std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "# Elapsed time: " << elapsed.count() << " s\n";

	return 0;
}
