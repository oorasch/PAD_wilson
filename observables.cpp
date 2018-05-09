#include <iostream>
#include <vector>
#include <cmath>
#include "constants.h"
#include "observables.h"

using namespace std;

double meas_plaq(const vector<int>& plaq_occ, const vector<double>& I_bessel)
{
	double tmp = 0.0;
	double aux1 = 0.25*(constants::eta + constants::eta_bar)/sqrt(constants::eta*constants::eta_bar),
				 aux2 = (0.25/constants::eta - 0.25/constants::eta_bar);

	for(auto it: plaq_occ)
	{
		if(it == 0)
		{
			tmp += 2.0*aux1*I_bessel.at(1)/I_bessel.at(0);
		}
		else
		{
			tmp += aux1*(I_bessel.at(abs(it)+1) + I_bessel.at(abs(it)-1))/I_bessel.at(abs(it));
			tmp += it*aux2;
		} 				
		
	}
	
	return tmp/constants::V;	
}

double meas_topo_chrg(const vector<int>& plaq_occ, const vector<double>& I_bessel)
{	
	double tmp = 0.0;
	double aux1 =  (0.125/constants::pi)*(constants::eta - constants::eta_bar)/sqrt(constants::eta*constants::eta_bar);
	double aux2 = (0.125/constants::eta + 0.125/constants::eta_bar)/constants::pi;
	
	for(auto it: plaq_occ)
	{
		if(it == 0) 
		{
			tmp += 2.0*aux1*I_bessel.at(1)/I_bessel.at(0);
		}
		else
		{
			tmp += aux1*(I_bessel.at(abs(it)+1) + I_bessel.at(abs(it)-1))/I_bessel.at(abs(it));
			tmp -= it*aux2;
		}
		
	}
	
	return tmp/constants::V;
}


