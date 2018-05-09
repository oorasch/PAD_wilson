#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include "updates.h"
#include "initialize.h"
#include "determinant.h"
#include "constants.h"

const unsigned seed(std::chrono::system_clock::now().time_since_epoch().count());
////std::default_random_engine generator(seed);
std::ranlux48 generator(seed);

using namespace std;

/*######################################################################################################*/
/*############################################ Loop updates ############################################*/
/*######################################################################################################*/

void do_updates(const int n, const int neib[][4], const vector<double>& I_bessel, double& detM, vector<int>& s_site, vector<int>& plaq_occ, vector<vector<int>>& k_link, vector<vector<int>>& sig_link)
{
	for(int i = 0; i < n; i++)
	{		
		do_loop_updates(neib, I_bessel, detM, s_site, plaq_occ, k_link, sig_link);
		do_blanket_update(I_bessel, plaq_occ);
		do_sigma_update(neib, s_site, detM, sig_link);
		//cout << "-----------------------------------" << endl;
	}
}

void do_loop_updates(const int neib[][4], const vector<double>& I_bessel, double& detM, vector<int>& s_site, vector<int>& plaq_occ, vector<vector<int>>& k_link, vector<vector<int>>& sig_link)
{
	vector<int> coords (4, 0);
	int site_check = 0;
	int link_check = 0;
	
	for(int i = 0; i < constants::V; i++)
	{
		coords.at(0) = i;
		coords.at(1) = neib[i][1];
		coords.at(2) = neib[coords.at(1)][0];
		coords.at(3) = neib[i][0];
				
		//cout << site_check << " " << link_check << endl;
		
		site_check = s_site.at(coords.at(0)) + s_site.at(coords.at(1)) + s_site.at(coords.at(2)) + s_site.at(coords.at(3));
		link_check = abs(k_link[coords.at(0)][1]) + abs(k_link[coords.at(1)][0]);
		link_check += abs(k_link[coords.at(3)][1]) + abs(k_link[coords.at(0)][0]);
		
		if(site_check == 4 && link_check == 0)//insert a plaquette
		{
			plaquette_ins_update(neib, coords, sig_link, I_bessel, detM, s_site, plaq_occ, k_link);
		}
		else if(site_check == 0 && link_check == 4)//delete a plaquette
		{
			plaquette_del_update(neib, coords, sig_link, I_bessel, detM, s_site, plaq_occ, k_link);
		}
		else if(site_check == 2 && link_check == 1)//expand a loop fragment
		{
			loop_exp_update(neib, coords, sig_link, I_bessel, detM, s_site, plaq_occ, k_link);
		}
		else if(site_check == 0 && link_check == 3)//collapse a loop fragment
		{
			loop_col_update(neib, coords, sig_link, I_bessel, detM, s_site, plaq_occ, k_link);
		}
		//else if(site_check == 0 && link_check == 2)//rotate loop
		else if(k_link[coords.at(0)][1]*k_link[coords.at(3)][1] < 0 || k_link[coords.at(0)][0]*k_link[coords.at(1)][0] < 0)//rotate loop 
		{
			loop_rot_update(coords, I_bessel, plaq_occ, k_link);
		}
	}
}

//insertion of a plaquette
void plaquette_ins_update(const int neib[][4], const vector<int>& coords, const vector<vector<int>>& sig_link, const vector<double>& I_bessel, double& detM, vector<int>& s_site, vector<int>& plaq_occ, vector<vector<int>>& k_link)
{
	std::uniform_int_distribution<int> delta_distribution(0,1);
	std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);
	int delta = 2*delta_distribution(generator) - 1; //+-1, the plaquette occupation number changes by -+1
	vector<int> s_site_prime = s_site;
	vector<double> M_prime (constants::V*constants::V, 0.0);
	double detM_prime = 0.0;
	double rho = 1.0;
	//int delta_L = 4;
	
	//change the loop occupation variables
	s_site_prime.at(coords.at(0)) = 0;
	s_site_prime.at(coords.at(1)) = 0;
	s_site_prime.at(coords.at(2)) = 0;
	s_site_prime.at(coords.at(3)) = 0;
	
	//fill the fermion matrix
	calc_M_array(neib, sig_link, s_site_prime, M_prime);
	
	//calculate the log of the determinant
	detM_prime = det(constants::V, M_prime);
	
	//cout << plaq_occ.at(coords.at(0)) << " ";

	//cout << detM_prime/detM << endl;
	
	if(detM_prime != 0.0) //catch the case detM = 0
	{
	
		//calculate the Metropolis weight
		rho /= 2.0*2.0*2.0*2.0;
		rho *= detM_prime/detM;
		
		if(delta == 1)
		{
			rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
			rho *= sqrt(constants::eta_bar/constants::eta); //when delta = 1, delta_p = - 1
		}
		else
		{
			rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
			rho *= sqrt(constants::eta/constants::eta_bar);
		}
	
		//cout << rho << endl;
	
		if(rho >= 1.0)
		{
			s_site.at(coords.at(0)) = 0;
			s_site.at(coords.at(1)) = 0;
			s_site.at(coords.at(2)) = 0;
			s_site.at(coords.at(3)) = 0;
		
			plaq_occ.at(coords.at(0)) -= delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa
		
			detM = detM_prime; // store the fermion determinant
		
			//change the loop configuration
			k_link[coords.at(0)][1] += delta;
			k_link[coords.at(1)][0] += delta;
			k_link[coords.at(3)][1] -= delta;
			k_link[coords.at(0)][0] -= delta;
		}
		else if(uniform_distribution(generator) < rho)
		{
			s_site.at(coords.at(0)) = 0;
			s_site.at(coords.at(1)) = 0;
			s_site.at(coords.at(2)) = 0;
			s_site.at(coords.at(3)) = 0;
		
			plaq_occ.at(coords.at(0)) -= delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa
		
			detM = detM_prime; // store the fermion determinant
		
			//change the loop configuration
			k_link[coords.at(0)][1] += delta;
			k_link[coords.at(1)][0] += delta;
			k_link[coords.at(3)][1] -= delta;
			k_link[coords.at(0)][0] -= delta;
		}
	}
	
}

//deletion of a plaquette
void plaquette_del_update (const int neib[][4], const vector<int>& coords, const vector<vector<int>>& sig_link, const vector<double>& I_bessel, double& detM, vector<int>& s_site, vector<int>& plaq_occ, vector<vector<int>>& k_link)
{
	std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);
	vector<int> s_site_prime = s_site;
	int delta = k_link[coords.at(0)][1];
	vector<double> M_prime (constants::V*constants::V, 0.0);
	double detM_prime = 0.0;
	double rho = 1.0;
	//int delta_L = -4;
	
	//change the loop occupation variables
	s_site_prime.at(coords.at(0)) = 1;
	s_site_prime.at(coords.at(1)) = 1;
	s_site_prime.at(coords.at(2)) = 1;
	s_site_prime.at(coords.at(3)) = 1;
	
	//fill the fermion matrix
	calc_M_array(neib, sig_link, s_site_prime, M_prime);
	
	//calculate the log of the determinant
	detM_prime = det(constants::V, M_prime);
	
	if(detM_prime != 0.0 ) //catch the case detM = 0
	{
		//calculate the Metropolis weight
		rho *= 2.0*2.0*2.0*2.0;
		rho *= detM_prime/detM;
		
		if(delta == 1)
		{
			rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
			rho *= sqrt(constants::eta/constants::eta_bar);
		}
		else
		{
			rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
			rho *= sqrt(constants::eta_bar/constants::eta);
		}
	
		if(rho  >= 1.0)
		{
			s_site.at(coords.at(0)) = 1;
			s_site.at(coords.at(1)) = 1;
			s_site.at(coords.at(2)) = 1;
			s_site.at(coords.at(3)) = 1;
		
			//is this true?!?
			plaq_occ.at(coords.at(0)) += delta; //
		
			detM = detM_prime; // store the fermion determinant
		
			//change the loop configuration
			k_link[coords.at(0)][1] = 0;
			k_link[coords.at(1)][0] = 0;
			k_link[coords.at(3)][1] = 0;
			k_link[coords.at(0)][0] = 0;
		}
		else if(uniform_distribution(generator) < rho)
		{
			s_site.at(coords.at(0)) = 1;
			s_site.at(coords.at(1)) = 1;
			s_site.at(coords.at(2)) = 1;
			s_site.at(coords.at(3)) = 1;
		
			plaq_occ.at(coords.at(0)) += delta; //
		
			detM = detM_prime; // store the fermion determinant
		
			//change the loop configuration
			k_link[coords.at(0)][1] = 0;
			k_link[coords.at(1)][0] = 0;
			k_link[coords.at(3)][1] = 0;
			k_link[coords.at(0)][0] = 0;
		}
	}
	
}

//expand a loop
void loop_exp_update(const int neib[][4], const vector<int>& coords, const vector<vector<int>>& sig_link, const vector<double>& I_bessel, double& detM, vector<int>& s_site, vector<int>& plaq_occ, vector<vector<int>>& k_link)
{
	std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);
	vector<int> s_site_prime = s_site;
	int delta;
	vector<double> M_prime (constants::V*constants::V, 0.0);
	double detM_prime = 0.0;
	double rho = 1.0;
	//int delta_L = 2;
	
	if(k_link[coords.at(0)][0] != 0)//loop in time direction
	{
		delta = k_link[coords.at(0)][0];
		s_site_prime.at(coords.at(1)) = 0;
		s_site_prime.at(coords.at(2)) = 0;
		
		//fill the fermion matrix
		calc_M_array(neib, sig_link, s_site_prime, M_prime);

		//calculate the log of the determinant
		detM_prime = det(constants::V, M_prime);

		if(detM_prime != 0.0 ) //catch the case detM = 0
		{
		
			rho /= 2.0*2.0;
			rho *= detM_prime/detM;
		
			if(delta == 1)
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta_bar/constants::eta);
			}
			else
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta/constants::eta_bar);
			}
	
			if(rho >= 1.0)
			{
				s_site.at(coords.at(1)) = 0;
				s_site.at(coords.at(2)) = 0;
	
				plaq_occ.at(coords.at(0)) -= delta; //
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = delta;
				k_link[coords.at(1)][0] = delta;
				k_link[coords.at(3)][1] = -delta;
				k_link[coords.at(0)][0] = 0; //we have to delete this link
			}
			else if(uniform_distribution(generator) < rho)
			{
				s_site.at(coords.at(1)) = 0;
				s_site.at(coords.at(2)) = 0;
	
				plaq_occ.at(coords.at(0)) -= delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = delta;
				k_link[coords.at(1)][0] = delta;
				k_link[coords.at(3)][1] = -delta;
				k_link[coords.at(0)][0] = 0; //we have to delete this link
			}
		}
	}
	else if(k_link[coords.at(0)][1] != 0)//loop in spatial dierction
	{
		delta = k_link[coords.at(0)][1];
		s_site_prime.at(coords.at(2)) = 0;
		s_site_prime.at(coords.at(3)) = 0;
		
		//fill the fermion matrix
		calc_M_array(neib, sig_link, s_site_prime, M_prime);

		//calculate the log of the determinant
		detM_prime = det(constants::V, M_prime);

		if(detM_prime != 0.0 ) //catch the case detM = 0
		{
			rho /= 2.0*2.0;
			rho *= detM_prime/detM;
		
			if(delta == 1)
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta/constants::eta_bar);
			}
			else
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta_bar/constants::eta);
			}
			
			if(rho >= 1.0)
			{
				s_site.at(coords.at(2)) = 0;
				s_site.at(coords.at(3)) = 0;
	
				//is this true?!?
				plaq_occ.at(coords.at(0)) += delta; //
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = 0;
				k_link[coords.at(1)][0] = -delta;
				k_link[coords.at(3)][1] = delta;
				k_link[coords.at(0)][0] = delta; //we have to delete this link
			}
			else if(uniform_distribution(generator) < rho)
			{
				s_site.at(coords.at(2)) = 0;
				s_site.at(coords.at(3)) = 0;
	
				plaq_occ.at(coords.at(0)) += delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = 0;
				k_link[coords.at(1)][0] = -delta;
				k_link[coords.at(3)][1] = delta;
				k_link[coords.at(0)][0] = delta; //we have to delete this link
			}
			
		}
	}
	else if(k_link[coords.at(3)][1] != 0) //loop in spatial direction
	{
		delta = k_link[coords.at(3)][1];
		s_site_prime.at(coords.at(0)) = 0;
		s_site_prime.at(coords.at(1)) = 0;
		
		//fill the fermion matrix
		calc_M_array(neib, sig_link, s_site_prime, M_prime);

		//calculate the log of the determinant
		detM_prime = det(constants::V, M_prime);

		if(detM_prime != 0.0 ) //catch the case detM = 0
		{
			rho /= 2.0*2.0;
			rho *= detM_prime/detM;
		
			if(delta == 1)
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta_bar/constants::eta);
			}
			else
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta/constants::eta_bar);
			}
			
			if(rho >= 0.0)
			{
				s_site.at(coords.at(0)) = 0;
				s_site.at(coords.at(1)) = 0;
	
				//is this true?!?
				plaq_occ.at(coords.at(0)) -= delta; //
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = delta;
				k_link[coords.at(1)][0] = delta;
				k_link[coords.at(3)][1] = 0;
				k_link[coords.at(0)][0] = -delta; //we have to delete this link
			}
			else if(uniform_distribution(generator) < rho)
			{
				s_site.at(coords.at(0)) = 0;
				s_site.at(coords.at(1)) = 0;
	
				plaq_occ.at(coords.at(0)) -= delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = delta;
				k_link[coords.at(1)][0] = delta;
				k_link[coords.at(3)][1] = 0;
				k_link[coords.at(0)][0] = -delta; //we have to delete this link
			}
			
		}
	}
	else if(k_link[coords.at(1)][0] != 0)//loop in time direction//loop in temporal direction
	{
		delta = k_link[coords.at(1)][0];
		s_site_prime.at(coords.at(0)) = 0;
		s_site_prime.at(coords.at(3)) = 0;
		
		//fill the fermion matrix
		calc_M_array(neib, sig_link, s_site_prime, M_prime);

		//calculate the log of the determinant
		detM_prime = det(constants::V, M_prime);

		if(detM_prime != 0.0) //catch the case detM = 0
		{
			rho /= 2.0*2.0;
			rho *= detM_prime/detM;
		
			if(delta == 1)
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta/constants::eta_bar);
			}
			else
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta_bar/constants::eta);
			}
			
			if(rho >= 0.0)
			{
				s_site.at(coords.at(0)) = 0;
				s_site.at(coords.at(3)) = 0;
	
				//is this true?!?
				plaq_occ.at(coords.at(0)) += delta; //
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = -delta;
				k_link[coords.at(1)][0] = 0;
				k_link[coords.at(3)][1] = delta;
				k_link[coords.at(0)][0] = delta; //we have to delete this link
			}
			else if(uniform_distribution(generator) < rho)
			{
				s_site.at(coords.at(0)) = 0;
				s_site.at(coords.at(3)) = 0;
	
				plaq_occ.at(coords.at(0)) += delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = -delta;
				k_link[coords.at(1)][0] = 0;
				k_link[coords.at(3)][1] = delta;
				k_link[coords.at(0)][0] = delta; //we have to delete this link
			}
			
		}
	}
}

//collapse a loop
void loop_col_update(const int neib[][4], const vector<int>& coords, const vector<vector<int>>& sig_link, const vector<double>& I_bessel, double& detM, vector<int>& s_site, vector<int>& plaq_occ, vector<vector<int>>& k_link)
{
	std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);
	vector<int> s_site_prime = s_site;
	int delta;
	vector<double> M_prime (constants::V*constants::V, 0.0);
	double detM_prime = 0.0;
	double rho = 1.0;
	//int delta_L = -2;
	
	if(k_link[coords.at(0)][0] == 0)//no loop in time direction 
	{
		delta = k_link[coords.at(1)][0];
		s_site_prime.at(coords.at(1)) = 1;
		s_site_prime.at(coords.at(2)) = 1;
		
		//fill the fermion matrix
		calc_M_array(neib, sig_link, s_site_prime, M_prime);

		//calculate the log of the determinant
		detM_prime = det(constants::V, M_prime);

		if(detM_prime != 0.0) //catch the case detM = 0
		{
			rho *= 2.0*2.0;
			rho *= detM_prime/detM;
		
			if(delta == 1)
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta/constants::eta_bar);
			}
			else
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta_bar/constants::eta);
			}
			
			if(rho >= 1.0)
			{
				s_site.at(coords.at(1)) = 1;
				s_site.at(coords.at(2)) = 1;
	
				//is this true?!?
				plaq_occ.at(coords.at(0)) += delta; //
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = 0;
				k_link[coords.at(1)][0] = 0;
				k_link[coords.at(3)][1] = 0;
				k_link[coords.at(0)][0] = delta; //we have to delete this link
			}
			else if(uniform_distribution(generator) < rho)
			{
				s_site.at(coords.at(1)) = 1;
				s_site.at(coords.at(2)) = 1;
	
				plaq_occ.at(coords.at(0)) += delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = 0;
				k_link[coords.at(1)][0] = 0;
				k_link[coords.at(3)][1] = 0;
				k_link[coords.at(0)][0] = delta; //activate this link
			}
		}
	}
	else if(k_link[coords.at(0)][1] == 0)//no loop in spatial direction
	{
		delta = k_link[coords.at(3)][1];
		s_site_prime.at(coords.at(2)) = 1;
		s_site_prime.at(coords.at(3)) = 1;
		
		//fill the fermion matrix
		calc_M_array(neib, sig_link, s_site_prime, M_prime);

		//calculate the log of the determinant
		detM_prime = det(constants::V, M_prime);

		if(detM_prime != 0.0) //catch the case detM = 0
		{
			rho *= 2.0*2.0;
			rho *= detM_prime/detM;
		
			if(delta == 1)
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta_bar/constants::eta);
			}
			else
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta/constants::eta_bar);
			}
			
			if(rho >= 1.0)
			{
				s_site.at(coords.at(2)) = 1;
				s_site.at(coords.at(3)) = 1;
	
				//is this true?!?
				plaq_occ.at(coords.at(0)) -= delta; //
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = delta;
				k_link[coords.at(1)][0] = 0;
				k_link[coords.at(3)][1] = 0;
				k_link[coords.at(0)][0] = 0; //we have to delete this link
			}
			else if(uniform_distribution(generator) < rho)
			{
				s_site.at(coords.at(2)) = 1;
				s_site.at(coords.at(3)) = 1;
	
				plaq_occ.at(coords.at(0)) -= delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = delta;
				k_link[coords.at(1)][0] = 0;
				k_link[coords.at(3)][1] = 0;
				k_link[coords.at(0)][0] = 0; //we have to delete this link
			}
			
		}
	}
	else if(k_link[coords.at(3)][1] == 0)//no loop in spatial direction
	{
		delta = k_link[coords.at(0)][1];
		s_site_prime.at(coords.at(0)) = 1;
		s_site_prime.at(coords.at(1)) = 1;
		
		//fill the fermion matrix
		calc_M_array(neib, sig_link, s_site_prime, M_prime);

		//calculate the log of the determinant
		detM_prime = det(constants::V, M_prime);

		if(detM_prime != 0.0) //catch the case detM = 0
		{
			rho *= 2.0*2.0;
			rho *= detM_prime/detM;
		
			if(delta == 1)
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta/constants::eta_bar);
			}
			else
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta_bar/constants::eta);
			}
			
			if(rho >= 1.0)
			{
				s_site.at(coords.at(0)) = 1;
				s_site.at(coords.at(1)) = 1;
	
				//is this true?!?
				plaq_occ.at(coords.at(0)) += delta; //
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = 0;
				k_link[coords.at(1)][0] = 0;
				k_link[coords.at(3)][1] = delta;
				k_link[coords.at(0)][0] = 0; //we have to delete this link
			}
			else if(uniform_distribution(generator) < rho)
			{
				s_site.at(coords.at(0)) = 1;
				s_site.at(coords.at(1)) = 1;
	
				plaq_occ.at(coords.at(0)) += delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = 0;
				k_link[coords.at(1)][0] = 0;
				k_link[coords.at(3)][1] = delta;
				k_link[coords.at(0)][0] = 0; //we have to delete this link
			}
			
		}
	}
	else if(k_link[coords.at(1)][0] == 0)//no loop in temporal direction
	{
		delta = k_link[coords.at(0)][0];
		s_site_prime.at(coords.at(0)) = 1;
		s_site_prime.at(coords.at(3)) = 1;
		
		//fill the fermion matrix
		calc_M_array(neib, sig_link, s_site_prime, M_prime);

		//calculate the log of the determinant
		detM_prime = det(constants::V, M_prime);

		if( detM_prime != 0.0 ) //catch the case detM = 0
		{
			rho *= 2.0*2.0;
			rho *= detM_prime/detM;
		
			if(delta == 1)
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta_bar/constants::eta);
			}
			else
			{
				rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
				rho *= sqrt(constants::eta/constants::eta_bar);
			}
			
			
			if(rho >= 0.0)
			{
				s_site.at(coords.at(0)) = 1;
				s_site.at(coords.at(3)) = 1;
	
				//is this true?!?
				plaq_occ.at(coords.at(0)) -= delta; //
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = 0;
				k_link[coords.at(1)][0] = delta;
				k_link[coords.at(3)][1] = 0;
				k_link[coords.at(0)][0] = 0; //we have to delete this link
			}
			else if(uniform_distribution(generator) < rho)
			{
				s_site.at(coords.at(0)) = 1;
				s_site.at(coords.at(3)) = 1;
	
				plaq_occ.at(coords.at(0)) -= delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa
	
				detM = detM_prime; // store the fermion determinant
	
				//change the loop configuration
				k_link[coords.at(0)][1] = 0;
				k_link[coords.at(1)][0] = delta;
				k_link[coords.at(3)][1] = 0;
				k_link[coords.at(0)][0] = 0; //we have to delete this link
			}
			
		}
	}
	
}


//join two loops/ cut a loop by rotation of loop fragments
void loop_rot_update(const vector<int>& coords, const vector<double>& I_bessel, vector<int>& plaq_occ, vector<vector<int>>& k_link)
{
	std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);
	int delta;
	double rho = 1.0;
	
	if(k_link[coords.at(0)][1] != 0)
	{
		delta = k_link[coords.at(0)][1];
		if(delta == 1)
		{
			rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
			rho *= sqrt(constants::eta_bar/constants::eta);
		}
		else
		{
			rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
			rho *= sqrt(constants::eta/constants::eta_bar);
		}
		
		if(rho >= 1.0)
			{
				plaq_occ.at(coords.at(0)) += delta; //
	
				//change the loop configuration
				k_link[coords.at(0)][1] = 0;
				k_link[coords.at(1)][0] = -delta;
				k_link[coords.at(3)][1] = 0;
				k_link[coords.at(0)][0] = delta; //we have to delete this link
			}
			else if(uniform_distribution(generator) < rho)
			{
				plaq_occ.at(coords.at(0)) += delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa
	
				//change the loop configuration
				k_link[coords.at(0)][1] = 0;
				k_link[coords.at(1)][0] = -delta;
				k_link[coords.at(3)][1] = 0;
				k_link[coords.at(0)][0] = delta; //we have to delete this link
			}
	}
	else
	{
		delta = k_link[coords.at(0)][0];
		if(delta == 1)
		{
			rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) - 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
			rho *= sqrt(constants::eta_bar/constants::eta);
		}
		else
		{
			rho *= I_bessel.at(abs(plaq_occ.at(coords.at(0)) + 1))/I_bessel.at(abs(plaq_occ.at(coords.at(0))));
			rho *= sqrt(constants::eta/constants::eta_bar);
		}
		
		if(rho >= 0.0)
		{
			plaq_occ.at(coords.at(0)) -= delta; //

			//change the loop configuration
			k_link[coords.at(0)][1] = delta;
			k_link[coords.at(1)][0] = 0;
			k_link[coords.at(3)][1] = -delta;
			k_link[coords.at(0)][0] = -0; //we have to delete this link
		}
		else if(uniform_distribution(generator) < rho)
		{
			plaq_occ.at(coords.at(0)) -= delta; //if delta introduces clock-wise flux, plaquette occupation must compensate that => anti-clock-wise and vice versa

			//change the loop configuration
			k_link[coords.at(0)][1] = delta;
			k_link[coords.at(1)][0] = 0;
			k_link[coords.at(3)][1] = -delta;
			k_link[coords.at(0)][0] = 0; //we have to delete this link
		}
	}
	
}

/*######################################################################################################*/
/*########################################## Blanket update ############################################*/
/*######################################################################################################*/

void do_blanket_update(const vector<double>& I_bessel, vector<int>& plaq_occ)
{
	std::uniform_int_distribution<int> delta_distribution(0,1);
	std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);
	int delta = 2*delta_distribution(generator) - 1; //+-1, the plaquette occupation number changes by -+1
	double rho = 1.0, aux;
	
	aux = sqrt(constants::eta/constants::eta_bar);
	for(auto it: plaq_occ) 
	{
		rho *= I_bessel.at(abs(it + delta))/I_bessel.at(abs(it));
	}
	
	if(delta == 1) rho *= pow(aux, 1.0*constants::V);
	else					 rho *= pow(aux, (-1.0)*constants::V);

	//cout << rho << endl;
	
	if(rho >= 1.0)
	{
		for(auto& it: plaq_occ) it += delta;
	}
	else if(uniform_distribution(generator) < rho)
	{
		for(auto& it: plaq_occ) it += delta;
	}
	
	//cout << "->" << plaq_occ.at(0) << endl;
}

/*######################################################################################################*/
/*########################################### Sigma update #############################################*/
/*######################################################################################################*/

void do_sigma_update(const int neib [][4], const vector<int> s_site, double detM, vector<vector<int>>& sig_link)
{
	std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);
	//vector<vector<int>> sig_link_prime = sig_link;
	vector<double> M_prime (constants::V*constants::V, 0.0);
	double detM_prime = 0.0;
	double rho = 1.0;

	
	for(int i = 0; i < constants::V; i++)
	{
		for(int mu = 0; mu <= 1; mu++)
		{
			vector<vector<int>> sig_link_prime = sig_link;
			sig_link_prime.at(i).at(mu) *= -1;
			
			//fill the fermion matrix
			calc_M_array(neib, sig_link_prime, s_site, M_prime);
	
			//calculate the log of the determinant
			detM_prime = det(constants::V, M_prime);
			
			if(detM_prime != 0.0) //catch the case detM = 0
			{
				rho = detM_prime/detM;
	
				if(rho >= 1.0)
				{
					sig_link.at(i).at(mu) *= -1;
					detM = detM_prime;
				}
				else if(uniform_distribution(generator) < rho)
				{
					sig_link.at(i).at(mu) *= -1;
					detM = detM_prime;
				}
			}
		}
	}
}

